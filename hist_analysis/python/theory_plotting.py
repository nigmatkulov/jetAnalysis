"""Validation of stored theory envelopes; no uncertainty calculation lives here."""
import math


def validate_total_band(central, band, *, context):
    """Require a one-to-one match to all finite central TH1 bins, including edges."""
    if not central.InheritsFrom('TH1') or central.GetDimension() != 1:
        raise TypeError(f'{context}: expected a one-dimensional central TH1')
    if not band.InheritsFrom('TGraphAsymmErrors'):
        raise TypeError(f'{context}: expected TGraphAsymmErrors total band')
    bins = [i for i in range(1, central.GetNbinsX() + 1)
            if math.isfinite(central.GetBinContent(i))]
    if not bins or len(bins) != band.GetN():
        raise ValueError(f'{context}: central/band valid-point mismatch')
    records = []
    for point, bin_index in enumerate(bins):
        x, y = band.GetPointX(point), band.GetPointY(point)
        errors = (band.GetErrorXlow(point), band.GetErrorXhigh(point),
                  band.GetErrorYlow(point), band.GetErrorYhigh(point))
        if not all(math.isfinite(v) for v in (x, y, *errors)) or any(v < 0 for v in errors):
            raise ValueError(f'{context}: invalid point {point}')
        axis = central.GetXaxis()
        pairs = ((x, axis.GetBinCenter(bin_index)),
                 (y, central.GetBinContent(bin_index)),
                 (x-errors[0], axis.GetBinLowEdge(bin_index)),
                 (x+errors[1], axis.GetBinUpEdge(bin_index)))
        if not all(math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9) for a, b in pairs):
            raise ValueError(f'{context}: central/band alignment mismatch at point {point}')
        records.append((x, y, errors[2], errors[3]))
    print(f'{context}: {band.ClassName()}, {len(records)} valid points; '
          f'(x, central, total down, total up) = {records}')
    return records
