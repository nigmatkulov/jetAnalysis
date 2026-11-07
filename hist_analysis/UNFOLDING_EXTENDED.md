# Extended unfolding producer contract

Implemented 2026-10-12 in `processForestSimple.C`, in addition to the existing
geometrically matched fiducial response. This is a producer artifact family;
existing unfolding notebooks do not yet consume it.

## Observable and event correspondence

Truth is the global leading/subleading gen pair among jets supplied by the
forest. Reco is the independently ranked leading/subleading good-reco pair,
after configured jet ID / track-max pruning and the requested JER variation.
Both observables use the pair's pT average and arithmetic mean etaCM.

By explicit analysis decision, response entries connect the two pairs from the
SAME EVENT, without geometric matching or a positive-reference-pT requirement.
Thus rank/composition changes, and even selected reco jets without references,
are modeled as migrations whenever at least two gen jets exist. Truth is NOT
the reference pair attached to the reco-leading jets: that would change the
truth observable when a third jet enters the reco-leading pair.

The extended population remains conditional on the producer's existing event
selection and gen/reco overweight veto. It is not inclusive of rejected MC
events, and it cannot recover jets omitted upstream from the forest.

## ROOT keys and axes

For each key base below, append `_<i>` for nominal, `JerDefExtra_<i>` for
eta-dependent default, or `ControlledSample_<i>` / `StudiedSample_<i>` for the
nominal complementary split. Here `i` is the existing eta-cut index. There are
224 new histogram keys (7 roles × 8 cuts × 4 families), plus the MC-only
`unfoldingExtendedContract` TNamed metadata. Data output is unchanged.

| Key base | Type | Axes |
|---|---|---|
| `hUnfoldingExtendedResponse` | THnSparseD, 5D | gen pTave, gen etaCM, reco pTave, reco etaCM, truth category |
| `hUnfoldingExtendedTruth` | THnSparseD, 3D | gen pTave, gen etaCM, truth category |
| `hUnfoldingExtendedMatchedTruth` | THnSparseD, 3D | same as truth |
| `hUnfoldingExtendedMiss` | THnSparseD, 3D | same as truth |
| `hUnfoldingExtendedMeasured` | TH2D | reco pTave, reco etaCM |
| `hUnfoldingExtendedMatchedMeasured` | TH2D | reco pTave, reco etaCM |
| `hUnfoldingExtendedFake` | TH2D | reco pTave, reco etaCM |

Truth category is explicit: **1 = outside fiducial, 2 = fiducial**. Fiducial means
that the global gen pair passes the existing leading/subleading pT, delta-phi
and BOTH individual-jet etaCM cuts for this index. The category is not inferred
from dijet etaCM or pTave. Changing the eta-cut index can change the same truth
event's category.

Coordinate binning is the existing pTave [0,1000) GeV in 100 bins and etaCM
[-3.6,3.6) in 72 bins. Events outside these storage ranges are retained in ROOT
flow bins, not dropped. Flow bins aggregate coordinates beyond the range;
a consumer must explicitly retain them as guard components or account for their
feed-in/inefficiency when restricting the model. They cannot resolve shapes
beyond the storage range. The category itself always lies in a visible bin.

All fills retain event weights and Sumw2. Each family is filled once per event,
using only the first reco smearing replica. The zero/one-good-reco path still
fills extended truth and misses, including outside-fiducial truth. Nominal
controlled/studied objects remain complementary, using the existing event split.

## Exclusive accounting

| At least two gen jets | Selected reco leading pair | Extended contribution |
|---|---|---|
| yes | yes | truth, response, matched truth, measured, matched measured |
| yes | no | truth, miss |
| no | yes | measured, fake |
| no | no | none |

These identities hold per truth category/bin, including coordinate flow bins:

- truth = response truth marginal + miss;
- measured = response reco marginal + fake;
- directly filled matched marginals = the respective response projections.

Here “matched” denotes the same-event correspondence, not a geometric jet match.
The names distinguish marginals from inclusive distributions and are specific
to this extended family. The existing fiducial misses/fakes and pair diagnostics
keep their earlier geometric matching definitions. In particular, outside-
fiducial feed-in is an extended RESPONSE component, not an extended fake.

## How to unfold and report fiducial truth

A new consumer must flatten **(gen pTave, gen etaCM, category)** into truth bins,
and **(reco pTave, reco etaCM)** into measured bins. This generally makes a
rectangular response with twice as many visible truth bins as measured bins.
Do not pass the 5D response through the existing four-coordinate flattening
helper: that would integrate out the category and lose the fiducial distinction.
Do not combine it with the old fiducial fake/purity correction, which would
subtract outside-fiducial events already modeled by the extended response.

Build the inclusive detector response using the EXTENDED truth, measured and
response objects, accounting for misses and extended fakes once. Alternatively,
use the matched-only migration kernel with extended purity and category-resolved
efficiency, applying each correction exactly once. The truth denominator must
include pairs with no selected reco, rather than only response-producing truth.

After unfolding, select category 2, preserving its full covariance submatrix,
and report only this fiducial component. Category 1 remains a guard/nuisance
truth population, to be constrained by data where possible and otherwise by a
stated model. If coordinates are restricted, derive effective misses/fakes for
that restricted domain; stored full-domain components alone are insufficient.

An explicit category makes the observable well-defined but does not guarantee
that data distinguish fiducial and nonfiducial contributions. Sparse or nearly
degenerate response columns require model variations, guard-population constraints
or coarser binning. Zero-response truth bins remain unmeasured and require masks
or a separate extrapolation. Assess independent closure, prior/guard sensitivity,
finite-MC uncertainty and covariance before interpreting this as a measurement.

## Validation

Run from the repository root:

```bash
root -l -b -q -e 'gROOT->ProcessLine(".L JetCorrector.cc"); gROOT->ProcessLine(".L JetUncertainty.cc"); gROOT->ProcessLine(".L processForestSimple.C"); gROOT->ProcessLine(".L tests/test_unfolding_extended.C"); test_unfolding_extended();'
```

The fixture checks category membership, threshold/angular/eta failures, rank
changes, absent pairs, missing references, positive/negative weights, orientations,
splits, coordinate flows, independent response marginals, miss/fake accounting
and Sumw2. It also exercises production callers in all smearing modes and verifies
persisted keys/metadata. Generated test output is written to `/tmp`.

Validation results: Release configure/build and CTest (1/1) passed, with only
the two existing unused-parameter warnings. The final extended fixture passed
288 scenario/frame/split/weight combinations across all eta cuts, plus production
caller checks for all smearing modes and 0/1/2 good-reco jets. Existing pair-
population, seven-category, missing-reco and JER-cache ROOT regressions passed;
all 57 Python tests passed. A deterministic 96-event synthetic comparison against
the immediately preceding producer found all 756 existing ROOT objects exactly
identical, with 224 added histograms and one metadata object. `git diff --check`
passed. No full production was regenerated and no extended unfolding notebook
was executed: this change supplies and validates producer ingredients only.
