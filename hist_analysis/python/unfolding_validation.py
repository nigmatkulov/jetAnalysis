"""Matched-only geometric unfolding and conditional covariance diagnostics."""

from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
import ROOT

from hist_analysis.python.forward_backward import fold_eta_histogram, forward_backward_ratio
from hist_analysis.python.histogram_io import load_object, open_root_file
from hist_analysis.python.unfolding import (
    as_pt_intervals, build_roounfold_response, flatten_pt_eta_projections,
    flatten_sparse_response, histogram_fingerprint, project_eta_by_pt,
    validate_common_binning, write_unfolding_output,
)


def edges(axis):
    return np.array([axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)])


def values(h):
    return np.array([h.GetBinContent(i) for i in range(1, h.GetNbinsX() + 1)])


def variances(h):
    return np.array([h.GetBinError(i)**2 for i in range(1, h.GetNbinsX() + 1)])


def matrix_values(h, *, errors=False):
    get = h.GetBinError if errors else h.GetBinContent
    return np.array([[get(i, j) for j in range(1, h.GetNbinsY()+1)]
                     for i in range(1, h.GetNbinsX()+1)])


def root_matrix(a):
    result = ROOT.TMatrixD(*a.shape)
    for i, row in enumerate(a):
        for j, value in enumerate(row):
            result[i][j] = value
    return result


def numpy_matrix(m):
    return np.array([[m[i][j] for j in range(m.GetNcols())] for i in range(m.GetNrows())])


def divide(a, b):
    return np.divide(a, b, out=np.full(np.broadcast_shapes(np.shape(a), np.shape(b)), np.nan),
                     where=b > 0)


def histogram(template, name, y, variance):
    h = template.Clone(name)
    h.SetDirectory(0)
    h.Reset("ICES")
    if not h.GetSumw2N():
        h.Sumw2()
    for i, (v, s) in enumerate(zip(y, variance), 1):
        h.SetBinContent(i, v)
        h.SetBinError(i, np.sqrt(max(s, 0)))
    return h


def _same_edges(a, b):
    return len(a) == len(b) and np.allclose(a, b, atol=1e-9, rtol=0)


def _check_pt(axis, intervals):
    source = edges(axis)
    for boundary in np.ravel(intervals):
        if not np.any(np.isclose(source, boundary, rtol=0, atol=1e-9)):
            raise ValueError(f"pT boundary {boundary:g} splits a source bin or includes flow")


def flatten_spectrum(h, intervals, eta_range, *, name, eta_axis=None):
    """Validate source edges and independently check direct weighted TH2 sums."""
    if h.GetDimension() != 2 or not h.GetSumw2N():
        raise ValueError(f"{h.GetName()}: need a weighted TH2 with Sumw2")
    pt, eta = h.GetXaxis(), h.GetYaxis()
    if any(axis.TestBit(ROOT.TAxis.kAxisRange) for axis in (pt, eta)):
        raise ValueError(f"{h.GetName()}: remove active source axis ranges before projection")
    transposed = pt.GetXmax() <= 20
    if transposed:
        pt, eta = eta, pt
    _check_pt(pt, intervals)
    if eta_axis is not None and not _same_edges(edges(eta), edges(eta_axis)):
        raise ValueError(f"{h.GetName()}: eta axis differs from response")
    template = project_eta_by_pt(h, intervals, name_prefix=name, eta_range=eta_range)
    flat, layout = flatten_pt_eta_projections(template, name=name, pt_bins=intervals)
    eta_ids = [eta.FindBin(x) for x in (edges(template[0].GetXaxis())[:-1]
                                        + np.diff(edges(template[0].GetXaxis()))/2)]
    direct = []
    for low, high in intervals:
        for y in eta_ids:
            cells = [h.GetBin(y, x) if transposed else h.GetBin(x, y)
                     for x in range(pt.FindBin(low+0.001), pt.FindBin(high-0.001)+1)]
            direct.append([sum(h.GetBinContent(c) for c in cells),
                           sum(h.GetBinError(c)**2 for c in cells)])
    np.testing.assert_allclose(np.array(direct).T, [values(flat), variances(flat)],
                               rtol=1e-10, atol=0)
    return flat, template, layout


def _audit_sparse(sparse, layout, eta_edges, matrix):
    """Cross-check ROOT projections against a single scan of stored coordinates."""
    from array import array
    maps = []
    for pt_dim, eta_dim in ((0, 1), (2, 3)):
        pt, eta = sparse.GetAxis(pt_dim), sparse.GetAxis(eta_dim)
        pt_ids = {b: i for i, (lo, hi) in enumerate(layout.pt_intervals)
                  for b in range(pt.FindBin(lo+0.001), pt.FindBin(hi-0.001)+1)}
        eta_ids = {eta.FindBin(x): i for i, x in enumerate(
            (eta_edges[1:]+eta_edges[:-1])/2)}
        maps.append((pt_ids, eta_ids))
    coordinates = array("i", [0]*4)
    expected = np.zeros((2, layout.n_global_bins, layout.n_global_bins))
    dropped = dict(gen_only=0., reco_only=0., neither=0.)
    for i in range(sparse.GetNbins()):
        weight = sparse.GetBinContent(i, coordinates)
        error2 = sparse.GetBinError2(i)
        if not np.isfinite(weight) or not np.isfinite(error2) or weight < 0 or error2 < 0:
            raise ValueError("Invalid or negative-weight sparse response cell")
        selected = []
        for dim, (pt_ids, eta_ids) in zip((0, 2), maps):
            a, b = pt_ids.get(coordinates[dim]), eta_ids.get(coordinates[dim+1])
            selected.append(None if a is None or b is None else a*layout.n_eta_bins+b)
        j, k = selected
        if j is not None and k is not None:
            expected[:, k, j] += (weight, error2)
        else:
            dropped["gen_only" if j is not None else "reco_only" if k is not None else "neither"] += weight
    np.testing.assert_allclose(matrix_values(matrix), expected[0], rtol=1e-11, atol=0)
    np.testing.assert_allclose(matrix_values(matrix, errors=True)**2, expected[1],
                               rtol=1e-11, atol=0)
    return dropped


def _fraction(matched, total):
    """Delta method on the disjoint matched and complement weighted yields."""
    a, t = values(matched), values(total)
    va, vt = variances(matched), variances(total)
    for part, whole in ((a, t), (va, vt)):
        if (not np.isfinite([part, whole]).all() or np.any(part < 0)
                or np.any(whole < 0) or np.any(part > whole*(1+1e-5))):
            raise ValueError("Matched and inclusive content or Sumw2 violate the weighted subset relation")
    f = np.minimum(divide(a, t), 1.)
    # Work with relative errors: a common normalization does not change this variance.
    error = np.hypot((1-f)*divide(np.sqrt(va), t),
                     f*divide(np.sqrt(np.maximum(vt-va, 0)), t))
    return f, error**2


@dataclass
class Training:
    matrix: object
    matched_reco: object
    matched_gen: object
    reco: object
    gen: object
    purity: object
    inverse_efficiency: object
    layout: object
    templates: list
    metadata: dict

    @property
    def probabilities(self):
        m = matrix_values(self.matrix)
        total = values(self.matched_gen)
        return np.divide(m, total, out=np.zeros_like(m), where=total > 0)


def make_training(matrix, reco, gen, layout, templates, metadata, *, cached=None):
    """Derive matched marginals only from the nonflow ROOT matrix projections."""
    nx, ny = matrix.GetNbinsX(), matrix.GetNbinsY()
    if (not matrix.GetSumw2N() or nx != ny or nx != layout.n_global_bins
            or any(matrix.GetBinContent(i, j) != 0 or matrix.GetBinError(i, j) != 0
                   for i, j in ([(i, j) for i in (0, nx+1) for j in range(ny+2)]
                                + [(i, j) for j in (0, ny+1) for i in range(1, nx+1)]))):
        raise ValueError("Unexpected flattened response size, flow, or missing Sumw2")
    m = matrix_values(matrix)
    w2 = matrix_values(matrix, errors=True)**2
    if (not np.isfinite(m).all() or not np.isfinite(w2).all()
            or np.any(m < 0) or np.any((m == 0) & (w2 != 0))):
        raise ValueError("Bayesian migration probabilities require nonnegative finite bins")
    mr = matrix.ProjectionX("matched_reco", 1, ny, "e")
    mg = matrix.ProjectionY("matched_gen", 1, nx, "e")
    for side, matched, total, axis in ((0, mg, gen, matrix.GetYaxis()),
                                      (1, mr, reco, matrix.GetXaxis())):
        matched.SetDirectory(0)
        validate_common_binning([matched, total], label="training projections")
        np.testing.assert_allclose(values(matched), m.sum(axis=side), rtol=1e-11, atol=0)
        np.testing.assert_allclose(variances(matched), w2.sum(axis=side), rtol=1e-11, atol=0)
        if axis.TestBit(ROOT.TAxis.kAxisRange):
            raise ValueError("Flattened response axis must not have an active range")
    if cached is None:
        p, vp = _fraction(mr, reco)
        e, ve = _fraction(mg, gen)
        inverse = divide(1., e)
        vinverse = divide(ve, e**4)
        purity = histogram(reco, "purity", p, vp)
        inv = histogram(gen, "inverse_efficiency", inverse, vinverse)
    else:
        purity, inv = cached
        validate_common_binning([purity, reco], label="cached purity")
        validate_common_binning([inv, gen], label="cached inverse efficiency")
    return Training(matrix, mr, mg, reco, gen, purity, inv, layout, templates, metadata)


def load_training(filename, keys, pt_bins, eta_range, *, selection, cache=None, reuse=False):
    """Read coherent geometric inputs or reuse factors tied to the source ROOT revision."""
    filename = Path(filename).resolve()
    if cache is not None and Path(cache).resolve() == filename:
        raise ValueError("Correction cache must be separate from source MC")
    intervals = as_pt_intervals(pt_bins)
    root_file = open_root_file(str(filename))
    try:
        key_cycles = {key: int(root_file.GetKey(key).GetCycle()) for key in keys
                      if root_file.GetKey(key)}
        if len(key_cycles) != 3:
            raise KeyError(f"Training keys missing: {set(keys)-key_cycles.keys()}")
        uuid = str(root_file.GetUUID().AsString())
    finally:
        root_file.Close()
    reco_source, gen_source = (load_object(str(filename), k) for k in keys[:2])
    reco, templates, layout = flatten_spectrum(reco_source, intervals, eta_range,
                                               name="training_reco")
    gen, gen_templates, _ = flatten_spectrum(gen_source, intervals, eta_range,
                                              name="training_gen")
    validate_common_binning([templates[0], gen_templates[0]], label="Reco/Gen eta")
    eta_edges = edges(templates[0].GetXaxis())
    source_stat = filename.stat()
    identity = dict(schema=1, source=str(filename), uuid=uuid, keys=key_cycles,
                    roles=list(keys), source_size=source_stat.st_size,
                    source_mtime_ns=source_stat.st_mtime_ns, selection=selection,
                    pt_intervals=[list(p) for p in intervals], eta_range=list(eta_range),
                    eta_edges=eta_edges.tolist(), normalization="original producer weights",
                    reco=histogram_fingerprint(reco), gen=histogram_fingerprint(gen))
    if reuse:
        if cache is None:
            raise ValueError("Select a correction cache to reuse")
        f = open_root_file(str(cache))
        try:
            info = f.Get("unfoldingConfiguration")
            metadata = json.loads(str(info.GetTitle())) if info else {}
        finally:
            f.Close()
        if metadata.get("identity") != identity:
            raise ValueError("Incompatible cached response or corrections")
        matrix, p, inv = (load_object(str(cache), key) for key in
                          ("matched_response", "purity", "inverse_efficiency"))
        if [histogram_fingerprint(h) for h in (matrix, p, inv)] != metadata["digests"]:
            raise ValueError("Cached response or corrections failed checksum verification")
        return make_training(matrix, reco, gen, layout, templates, metadata, cached=(p, inv))

    sparse = load_object(str(filename), keys[2])
    if (not sparse.InheritsFrom("THnSparse") or sparse.GetNdimensions() != 4
            or not sparse.GetCalculateErrors()):
        raise ValueError("Expected a four-coordinate geometric response with Sumw2")
    for d, source in ((0, gen_source), (2, reco_source)):
        pt, eta = source.GetXaxis(), source.GetYaxis()
        if pt.GetXmax() <= 20:
            pt, eta = eta, pt
        for offset, axis in enumerate((pt, eta)):
            stored = sparse.GetAxis(d+offset)
            if (not _same_edges(edges(stored), edges(axis))
                    or stored.TestBit(ROOT.TAxis.kAxisRange)):
                raise ValueError("Response axes do not match the inclusive spectra")
    matrix, _ = flatten_sparse_response(sparse, intervals, name="matched_response",
                                        layout=layout, eta_range=eta_range)
    dropped = _audit_sparse(sparse, layout, eta_edges, matrix)
    result = make_training(matrix, reco, gen, layout, templates,
                           dict(identity=identity, matched_weight_excluded=dropped))
    if cache is not None:
        objects = (matrix, result.purity, result.inverse_efficiency)
        result.metadata["digests"] = [histogram_fingerprint(h) for h in objects]
        write_unfolding_output(cache, histograms=(*objects, result.matched_reco,
                                                 result.matched_gen), metadata=result.metadata)
    return result


def _bayes_derivative(a, measured, prior, iterations):
    """Independent iterative-Bayes recursion and derivative at fixed training."""
    prior = prior / prior.sum()
    derivative = None
    for _ in range(iterations):
        prediction = a @ prior
        reciprocal = np.divide(1., prediction, out=np.zeros_like(prediction),
                               where=prediction > 0)
        kernel = prior[:, None] * a.T * reciprocal
        result = kernel @ measured
        if derivative is None:
            derivative = kernel
        else:
            correction = (np.diag(a.T @ (measured * reciprocal))
                          - (kernel * (measured * reciprocal)) @ a)
            derivative = kernel + correction @ derivative
        # An unnormalized next prior is equivalent and avoids differentiating
        # the arbitrary common normalization of the prior.
        prior = result
    return result, derivative


@dataclass
class Spectrum:
    histogram: object
    covariance: np.ndarray
    supported: np.ndarray


@dataclass
class ValidationResult:
    reco: Spectrum
    signal: Spectrum
    unfolded: Spectrum
    refolded_signal: Spectrum
    refolded: Spectrum
    refolded_signal_cross: np.ndarray
    refolded_cross: np.ndarray
    bundle: object
    scale: float


def run_unfolding(training, measured, *, iterations, scale=None, flat_prior=False,
                  measured_covariance=None):
    """Propagate measurement covariance at fixed migration, purity and efficiency."""
    if not isinstance(iterations, int) or iterations < 1:
        raise ValueError("Use a positive integer iteration count")
    t = training
    validate_common_binning([measured, t.reco], label="measurement and response")
    m = matrix_values(t.matrix)
    a = t.probabilities
    positive = m[m > 0]
    if not positive.size:
        raise ValueError("The matched response is empty")
    # The local RooUnfold sanitizes isolated entries below an absolute 1e-9.
    # Use a small safety margin only if needed; always verify actual matrix cells.
    scale = max(1., 1e-8/positive.min()) if scale is None else scale
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Use a finite positive response scale")
    p = values(t.purity)
    inv = values(t.inverse_efficiency)
    valid_reco = np.isfinite(p) & (p > 0)
    valid_gen = np.isfinite(inv)
    p, inv = np.where(np.isfinite(p), p, 0.), np.where(valid_gen, inv, 0.)
    reco_values = values(measured)
    vr = np.diag(variances(measured)) if measured_covariance is None else measured_covariance
    if (vr.shape != (len(p), len(p)) or not np.isfinite(vr).all()
            or not np.allclose(vr, vr.T, atol=0, rtol=1e-10) or np.any(vr.diagonal() < 0)
            or not np.isfinite(reco_values).all() or np.any(reco_values < 0)):
        raise ValueError("Invalid measurement contents or covariance")
    signal = p * reco_values
    if signal.sum() <= 0:
        raise ValueError("No positive signal on trained response support")
    vs = p[:, None] * vr * p
    bundle = build_roounfold_response(
        ROOT, t.matched_gen, t.matched_reco, t.matrix, scale=scale,
        matched_only=True, require_fakes=False, name="matched_only",
    )
    if bundle.response.UseOverflowStatus():
        raise RuntimeError("Matched response must exclude flow")
    # RooUnfold f9fd812's empty-template ProjectionX misses the last regular
    # Gen column. Restore the actual projection without inferring fake entries.
    roo_reco = bundle.response.Hmeasured()
    for i in range(roo_reco.GetNbinsX()+2):
        roo_reco.SetBinContent(i, bundle.scaled_measured.GetBinContent(i))
        roo_reco.SetBinError(i, bundle.scaled_measured.GetBinError(i))
    for actual, expected in ((roo_reco, bundle.scaled_measured),
                             (bundle.response.Htruth(), bundle.scaled_truth)):
        np.testing.assert_allclose(values(actual), values(expected), rtol=1e-10, atol=0)
        np.testing.assert_allclose(variances(actual), variances(expected), rtol=1e-10, atol=0)
    raw = numpy_matrix(bundle.response.Mresponse(False))/scale
    # Isolated zero-prior cells can acquire RooUnfold dummy diagonal entries;
    # its Bayes loop skips all zero-prior columns.
    checked = raw[:, valid_gen]
    np.testing.assert_allclose(checked, m[:, valid_gen], rtol=1e-10, atol=0,
                               err_msg="RooUnfold changed supported migration cells")
    prior = valid_gen.astype(float) if flat_prior else values(t.matched_gen)
    prior /= prior.sum()
    scaled = histogram(measured, "roo_measured", signal*scale, vs.diagonal()*scale**2)
    algo = ROOT.RooUnfoldBayes(bundle.response, scaled, iterations, False, False)
    algo.SetVerbose(-1)
    reference, jac = _bayes_derivative(a, signal, prior, iterations)
    if flat_prior:
        # Diagnostic only. The local SetPriors checks an uninitialized size
        # before the first unfold. Keep the nominal result on RooUnfold.
        x, vx = reference, jac @ vs @ jac.T
    else:
        if measured_covariance is not None:
            algo.SetMeasuredCov(root_matrix(vs*scale**2))
        roo = algo.Hunfold()
        x = values(roo)/scale
        vx = numpy_matrix(algo.Eunfold(ROOT.RooUnfolding.kCovariance))/scale**2
        # The independent recursion also provides cross covariance with reco.
        np.testing.assert_allclose(x, reference, rtol=1e-9, atol=1e-14*signal.max())
        np.testing.assert_allclose(vx, jac @ vs @ jac.T, rtol=1e-7,
                                   atol=1e-12*np.max(vx.diagonal()))
    corrected = inv*x
    vc = inv[:, None]*vx*inv
    # Efficiency cancels the inverse efficiency on supported truth bins.
    matched_refold = a @ x
    fold_jac = a @ jac
    vf = a @ vx @ a.T
    reciprocal = np.divide(1., p, out=np.zeros_like(p), where=valid_reco)
    inclusive_refold = reciprocal * matched_refold
    vi = reciprocal[:, None]*vf*reciprocal
    cross_signal = fold_jac @ vs
    cross_inclusive = reciprocal[:, None]*((fold_jac * p) @ vr)

    def spectrum(template, name, content, covariance, mask):
        return Spectrum(histogram(template, name, content, covariance.diagonal()),
                        covariance, mask)

    return ValidationResult(
        Spectrum(measured, vr, np.ones(len(p), dtype=bool)),
        spectrum(measured, "purity_corrected_reco", signal, vs, np.isfinite(values(t.purity))),
        spectrum(t.gen, "corrected_unfolded", corrected, vc, valid_gen),
        spectrum(measured, "matched_refolded", matched_refold, vf, valid_reco),
        spectrum(measured, "refolded", inclusive_refold, vi, valid_reco),
        cross_signal, cross_inclusive, bundle, float(scale),
    )


def project(spectrum, transform, template, *, name):
    """Apply an aggregation or selection: y'=Ay, V'=AVA^T."""
    transform = np.asarray(transform)
    content = transform @ values(spectrum.histogram)
    covariance = transform @ spectrum.covariance @ transform.T
    supported = ~np.any((transform != 0) & ~spectrum.supported, axis=1)
    return Spectrum(histogram(template, name, content, covariance.diagonal()),
                    covariance, supported)


def ratio(numerator, denominator, *, name, cross=None):
    a, b = values(numerator.histogram), values(denominator.histogram)
    valid = numerator.supported & denominator.supported & (b > 0)
    y = divide(a, b)
    g = np.divide(1., b, out=np.zeros_like(b), where=valid)
    h = -g*g*a
    c = g[:, None]*numerator.covariance*g + h[:, None]*denominator.covariance*h
    if cross is not None:
        term = g[:, None]*cross*h
        c += term + term.T
    y[~valid] = np.nan
    return Spectrum(histogram(numerator.histogram, name, y, c.diagonal()), c, valid)


def forward_backward(spectrum, *, name):
    """Use the project's mirrored ROOT division and add the known cross term."""
    f, b = fold_eta_histogram(spectrum.histogram, name=name)
    central = forward_backward_ratio(f, b, name=name)
    half = f.GetNbinsX()
    pos, neg = np.arange(half, 2*half), np.arange(half-1, -1, -1)
    c = spectrum.covariance
    fs = Spectrum(f, c[np.ix_(pos, pos)], spectrum.supported[pos])
    bs = Spectrum(b, c[np.ix_(neg, neg)], spectrum.supported[neg])
    result = ratio(fs, bs, name=name, cross=c[np.ix_(pos, neg)])
    for i in range(half):
        if result.supported[i]:
            # Keep ROOT's central division; use full covariance for the errors.
            result.histogram.SetBinContent(i+1, central.GetBinContent(i+1))
    return result


def agreement(prediction, observed, cross, *, support=None):
    """Return residuals, pulls and a pseudoinverse quadratic form on support."""
    mask = prediction.supported & observed.supported
    if support is not None:
        mask &= support
    delta = values(prediction.histogram) - values(observed.histogram)
    c = prediction.covariance + observed.covariance - cross - cross.T
    c = (c + c.T)/2
    scale = np.max(observed.covariance.diagonal(), initial=0.)
    error = np.sqrt(np.maximum(c.diagonal(), 0.))
    pulls = divide(delta, error)
    delta[~mask] = np.nan
    pulls[~mask] = np.nan
    eig, vectors = np.linalg.eigh(c[np.ix_(mask, mask)])
    cutoff = max(scale*1e-12, np.max(eig, initial=0.)*1e-10)
    keep = eig > cutoff
    modes = vectors.T @ delta[mask]
    denominator = np.abs(values(observed.histogram)[mask]).sum()
    return dict(residual=delta, pulls=pulls, statistic=float(np.sum(modes[keep]**2/eig[keep])),
                rank=int(keep.sum()), null_residual=float(np.linalg.norm(modes[~keep])),
                eigenvalue_cutoff=float(cutoff), bins=int(mask.sum()),
                relative_l1=float(np.nansum(np.abs(delta))/denominator) if denominator else None)


def validation_views(training, result, truth=None, *, reco_truth_cross=None,
                     unfolded_truth_cross=None):
    """Build flattened and pT-sliced diagnostics from the same covariance."""
    t, r = training, result
    n = t.layout.n_global_bins
    zero = np.zeros((n, n))
    closure_cross = zero if unfolded_truth_cross is None else unfolded_truth_cross

    def fixed(h):
        return Spectrum(h, np.diag(np.nan_to_num(variances(h))), np.isfinite(values(h)))

    truth = fixed(truth) if truth is not None else None
    quantities = {
        "original": {"Reco": r.reco},
        "eta": {"Reco": r.reco, "Unfolded": r.unfolded},
        "refolded": {"Reco": r.reco, "Refolded": r.refolded},
        "refolded_ratio": {"Refolded / Reco": ratio(r.refolded, r.reco,
                                                  name="refolded_over_reco", cross=r.refolded_cross)},
        "matched_refolded": {"Corrected Reco": r.signal, "Matched refolded": r.refolded_signal},
        "matched_refolded_ratio": {"Matched refolded / Corrected Reco": ratio(
            r.refolded_signal, r.signal, name="matched_refolded_over_corrected_reco",
            cross=r.refolded_signal_cross)},
        "purity": {"Purity": fixed(t.purity)},
        "inverse_efficiency": {"Inverse efficiency": fixed(t.inverse_efficiency)},
    }
    if truth is not None:
        for kind in ("original", "eta"):
            quantities[kind]["Gen"] = truth
        quantities["to_gen"] = {
            "Reco / Gen": ratio(r.reco, truth, name="reco_over_gen", cross=reco_truth_cross),
            "Unfolded / Gen": ratio(r.unfolded, truth, name="unfolded_over_gen",
                                    cross=unfolded_truth_cross),
        }
    slices = {key: {} for key in (*quantities, "fb")}
    metrics = {
        "matched": agreement(r.refolded_signal, r.signal, r.refolded_signal_cross),
        "inclusive": agreement(r.refolded, r.reco, r.refolded_cross),
    }
    if truth is not None:
        metrics["closure"] = agreement(r.unfolded, truth, closure_cross)
    for i, interval in enumerate(t.layout.pt_intervals):
        local = np.arange(i*t.layout.n_eta_bins, (i+1)*t.layout.n_eta_bins)
        selector = np.eye(n)[local]
        template = t.templates[i]

        def block(s, label):
            # Direct indexing is the exact AVA^T for a one-block selection.
            v = s.covariance[np.ix_(local, local)]
            y = values(s.histogram)[local]
            return Spectrum(histogram(template, f"{label}_pt{i}", y, v.diagonal()),
                            v, s.supported[local])

        for kind, curves in quantities.items():
            if kind != "original":
                slices[kind][interval] = {label: block(s, s.histogram.GetName())
                                          for label, s in curves.items()}
        slices["fb"][interval] = {label: forward_backward(s, name=f"fb_{label}_pt{i}")
                                  for label, s in slices["eta"][interval].items()}
        if i == 0:
            projected = project(r.unfolded, selector, template, name="projection_check")
            direct = slices["eta"][interval]["Unfolded"]
            np.testing.assert_array_equal(values(projected.histogram), values(direct.histogram))
            np.testing.assert_array_equal(projected.covariance, direct.covariance)
        for label, kind, pair, cross in (
            ("inclusive", "refolded", ("Refolded", "Reco"), r.refolded_cross),
            ("matched", "matched_refolded", ("Matched refolded", "Corrected Reco"),
             r.refolded_signal_cross),
        ):
            group = slices[kind][interval]
            metrics[f"{label}_pt{i}"] = agreement(group[pair[0]], group[pair[1]],
                                                  cross[np.ix_(local, local)])
        if truth is not None:
            group = slices["eta"][interval]
            local_cross = closure_cross[np.ix_(local, local)]
            metrics[f"closure_pt{i}"] = agreement(group["Unfolded"], group["Gen"],
                                                  local_cross)
            if unfolded_truth_cross is not None:
                def fb_derivative(label):
                    y = values(group[label].histogram)
                    half = y.size // 2
                    denominator = y[:half][::-1]
                    gradient = np.zeros((half, y.size))
                    pos, neg = np.arange(half, y.size), np.arange(half-1, -1, -1)
                    reciprocal = np.divide(1., denominator, out=np.zeros(half),
                                           where=denominator > 0)
                    gradient[np.arange(half), pos] = reciprocal
                    gradient[np.arange(half), neg] = -y[half:] * reciprocal**2
                    return gradient

                fb = slices["fb"][interval]
                fb_cross = fb_derivative("Unfolded") @ local_cross @ fb_derivative("Gen").T
                metrics[f"fb_closure_pt{i}"] = agreement(fb["Unfolded"], fb["Gen"], fb_cross)
    return quantities, slices, metrics


def truth_refold_check(training, truth, measured):
    """Compare a known truth forward prediction before any fit to measured Reco."""
    validate_common_binning([truth, measured, training.gen], label="truth refold check")
    p = values(training.purity)
    inv = values(training.inverse_efficiency)
    good_r = np.isfinite(p) & (p > 0)
    good_g = np.isfinite(inv)
    efficiency = np.divide(1., inv, out=np.zeros_like(inv), where=good_g)
    matched = training.probabilities @ (efficiency * values(truth))
    prediction = divide(matched, p)
    reco = values(measured)
    total = reco[good_r].sum()
    return dict(relative_l1=float(np.abs(prediction[good_r]-reco[good_r]).sum()/total),
                integral_ratio=float(prediction[good_r].sum()/total),
                unsupported_gen_weight=float(values(truth)[~good_g].sum()))


def neighboring_factors(training):
    """Tabulate unsmoothed eta steps with local MC effective entries and migration."""
    t = training
    m = matrix_values(t.matrix)
    for side, h, matched in (("purity", t.purity, t.matched_reco),
                             ("inverse_efficiency", t.inverse_efficiency, t.matched_gen)):
        f, v = values(h), variances(h)
        w, w2 = values(matched), variances(matched)
        neff = divide(w**2, w2)
        same = divide(m.diagonal(), w)
        for k, (low, high) in enumerate(t.layout.pt_intervals):
            eta = edges(t.templates[k].GetXaxis())
            for j in range(1, t.layout.n_eta_bins):
                i = k*t.layout.n_eta_bins+j
                if np.isfinite(f[[i-1, i]]).all():
                    yield dict(factor=side, pt_low=low, pt_high=high, eta_boundary=eta[j],
                               left=f[i-1], right=f[i], difference=f[i]-f[i-1],
                               local_stat_scale=np.sqrt(v[i]+v[i-1]),
                               left_neff=neff[i-1], right_neff=neff[i],
                               left_same_bin_fraction=same[i-1], right_same_bin_fraction=same[i])


def save_validation(filename, training, result, flat, slices, metrics, *, provenance):
    """Write curves, full core covariances, refolding cross terms and bin diagnostics."""
    objects = {s.histogram.GetName(): s.histogram
               for group in (*flat.values(), *(x for kind in slices.values() for x in kind.values()))
               for s in group.values()}
    matrices = {
        "reco_covariance": result.reco.covariance,
        "signal_covariance": result.signal.covariance,
        "unfolded_covariance": result.unfolded.covariance,
        "matched_refolded_covariance": result.refolded_signal.covariance,
        "refolded_covariance": result.refolded.covariance,
        "matched_refolded_signal_cross": result.refolded_signal_cross,
        "refolded_reco_cross": result.refolded_cross,
    }
    for kind in ("eta", "fb", "refolded", "refolded_ratio", "matched_refolded",
                 "matched_refolded_ratio", "to_gen"):
        for group in slices.get(kind, {}).values():
            for s in group.values():
                matrices[f"{s.histogram.GetName()}_covariance"] = s.covariance
    metadata = dict(training=training.metadata, provenance=provenance, scale=result.scale,
                    covariance="measurement only, conditional on trained response and factors",
                    unsupported={s.histogram.GetName(): (np.flatnonzero(~s.supported)+1).tolist()
                                 for group in flat.values() for s in group.values()},
                    metrics={key: {k: v for k, v in info.items() if k not in {"pulls", "residual"}}
                             for key, info in metrics.items()})
    write_unfolding_output(filename, histograms=objects.values(),
                           response=result.bundle.response, metadata=metadata)
    f = open_root_file(str(filename), "UPDATE")
    try:
        f.cd()
        for name, a in matrices.items():
            if root_matrix(a).Write(name) <= 0:
                raise OSError(f"Could not save {name}")
        for name, info in metrics.items():
            for kind in ("residual", "pulls"):
                a = info[kind]
                h = ROOT.TH1D(f"{name}_{kind}", "", len(a), .5, len(a)+.5)
                h.SetDirectory(0)
                for i, y in enumerate(a, 1):
                    h.SetBinContent(i, y)
                if h.Write() <= 0:
                    raise OSError(f"Could not save {h.GetName()}")
    finally:
        f.Close()
