# Physics and implementation review

Reviewed 2026-10-10 against the working version following `4a8ab4c`.
The original review added comments and documentation only. A subsequent fix
restores missing-reco truth accounting (P02) and changes pPb pT-hat event and
weight intervals to `[low, high)` (boundary portion of P04). Calibration files,
random numbers, histogram keys and the job-normalization denominator are unchanged.
Findings below distinguish demonstrated implementation defects from conventions
that need an analysis decision. Passing regression tests does not certify the
calibration or statistical model for a measurement.

## Scope and reading guide

The detailed review covers `processForestSimple.C`, `AnalysisContract.h`, the
JEC/JEU helpers, histogram I/O/projection/normalization, hemisphere folding,
closure preparation, trigger stitching, systematic extraction, and the Python
and compiled unfolding helpers. Notebook configuration and ratio usage, local
runners, and processing documentation were cross-checked. External calibration
payloads were treated as data; generated files and notebook outputs were not
rewritten. The standalone PYTHIA source tree and every historical plotting
routine were not audited line by line.

Read the producer in this order:

1. `setupChains` / `setupBranches`: validate inputs and bind aligned event trees.
2. `processEvents`: orient vertex z, apply event selection and weights, load jets,
   veto MC tails, and dispatch physics fills.
3. `processGenDijets` / `processRefDijets`: truth and matched-reference populations.
4. `processRecoJets`: fill inclusive ID diagnostics, remove ID-failing jets, and
   prepare data energy-scale shifts.
5. `processRecoDijets`: reselect the two leading jets for each momentum variation.
6. `prepareUnfolding`: classify selected pairs and fill matched/miss/fake objects.
7. `writeOutput`: persist weighted yields and JER-mode metadata.

Downstream, `histogram_production.py` caches raw yields, `projections.py` selects
complete bins, `forward_backward.py` folds signed eta, `histogram_ops.py` creates
display shapes and ratios, and `unfolding.py` handles response accounting and
correction covariance. `regularization_scan.cpp` accelerates toy ensemble moments.

## Confirmed implementation issues

### P01 — High: all MC selects the embedding JEC branch

In the entry point, `isMc = (mcType != 0)` is Boolean, but
`isEmbedding = (isMc == 1)`. Thus both `mcType=1` and `mcType=2` enter the
Embedded calibration branch; the Unembedded branch is unreachable for MC.
The nominal PYTHIA and pp CP5 runners therefore do not select the unembedded
payload claimed by the earlier README. The file names themselves distinguish
embedded and unembedded calibrations.

This is a control-flow defect. Its yield/shape impact depends on the actual
payloads and input sample; the present review does not assign a numerical bias.
A proposed fix is to test `mcType == 1`, followed by explicit decisions about
which payload pp should use. Validate printed calibration paths, corrected jet
pT, threshold migrations, and response accounting before regenerating outputs.
No calibration was switched. As a concrete payload comparison, a 100 GeV
jet at eta 0 in the Pb-going branch is corrected to 105.000 GeV by Embedded
and 105.504 GeV by Unembedded; at eta +1.5 the values are 111.738 and
113.252 GeV. These are evaluations of the current files, not production yield
bias estimates. The difference can affect threshold acceptance and pair ordering.

### P02 — Fixed: truth misses with fewer than two selected reco jets

`processRecoJets` removes ID-failing jets. Before the fix, `processRecoDijets`
returned immediately when fewer than two remained, before either call to
`prepareUnfolding`. A selected generator dijet with zero or one selected reco
jet was therefore missing from the explicit miss and classification objects.
`prepareUnfolding` itself correctly supports this case when called directly.
The existing seven-category unit check calls that helper directly and does not
exercise the early return in its caller.

A synthetic pp event with gen pT 100/60 GeV, eta +/-0.2, delta phi 3.14, and no
reco jets gave stored miss integral **0** through `processRecoDijets`, but **1**
through a direct `prepareUnfolding` call with unit weight at eta-cut index 4.

The factorized Python workflow derives effective misses from inclusive truth
minus the matched response marginal, so it can recover the missing yield.
However, its residual over stored misses then includes this caller omission,
not just migrations outside the flattened acceptance. Direct consumers of
older stored misses can underestimate inefficiency. The caller now invokes
`prepareUnfolding` before the early return for both nominal and eta-dependent
response families. Selected truth contributes one event weight to misses and
gen-pass/reco-fail classification, with Sumw2 receiving the squared weight.
The nominal controlled/studied split is preserved. No smearing replica, measured
yield, fake, or matched-response entry is added. Truth-failing events add nothing.
Events with at least two selected reco jets follow the existing path.
`tests/test_unfolding_missing_reco.C` exercises the actual caller across frames,
directions, smearing modes, eta cuts, splits and positive/negative weights.

### P03 — High for stochastic systematics: sorting changes the random-draw identity

The JER hash excludes the variation index, but its `round` includes `iJet`.
`prepareUnfolding`, called at scale indices 4 and 7, sorts the shared reco
collection by varied momentum. If a pair changes rank, the same physical jet
gets a different index in subsequent variations and therefore a different
normal deviate. The claimed common draw across all Down/Default/Up templates
is not guaranteed. Hybrid jets with good matches do not use the deviate, but
unmatched/stochastic jets and pseudo-reco diagnostic replicas can be affected.

A synthetic hybrid pair with reco/reference pT 100/120 and 99/81 GeV swaps rank
under the nominal central SF. The collection remained reordered after the
systematic loop. For the original 100 GeV jet in event `(run,lumi,evt)=(1,1,17)`,
changing the hash round from 0 to 1 changed the deviate from **-0.726256** to
**-0.802291**. This demonstrates the two ingredients of the defect, without
claiming a measured production-spectrum difference.

Storing `refPtSmearedDef`, `refPtSmearedUp`, and `refPtSmearedDown` on each
`RecoJet` can solve this if all are computed from one cached normal deviate
before any sorting. The values then travel with the jet. Reinitializing the same
seed alone is insufficient if its hash includes a changing rank. Each variation
still needs its own pair selection; eta-dependent branches and pseudo-reco
replicas need the same stable-identity treatment. No JER change is applied here.
A fix should preserve a stable original jet identifier or cache each replica's
draw before any sorting. Validate the identity of each jet's draw, including
rank swaps, unmatched jets, and all systematic branches.

### P04 — Medium/high for absolute MC yields: pPb normalization depends on job partitioning

`eventWeight` divides legacy pT-hat coefficients times reference generated
counts by the total entries in the current invocation, before cuts. Processing
one sample as a single chain and processing its disjoint parts separately then
adding their outputs are not equivalent. For identical event distributions,
two equal partitions each receive twice the per-event weight of the full chain;
the merged yield is twice the full-chain yield. Different-sized partitions can
also alter the merged shape. One invocation per complete pT-hat sample is the
local runner convention; arbitrary MC sublists need a shared denominator or
an explicitly documented merge rescaling.

The boundary inconsistency is now fixed: both pPb event selection and the
weight-table lookup use `[low, high)`. For example, pT-hat exactly 50 GeV is
rejected by the 30–50 sample and accepted by the 50–80 sample with its own
coefficient. This is independent of the invocation-size denominator: changing
interval endpoints cannot make separately normalized jobs additive.

For N accepted entries with coefficient C and unit vertex correction, one full
job gives each event C/N and total C. Two equal jobs each give C/(N/2), and their
sum is 2C. More generally, the denominator is input entries before selection;
each job contributes C times its own selected fraction. Unequal job sizes weight
those fractions (and bin shapes) differently from a full-chain invocation.
A shared full-sample denominator, or multiplying each output by N_job/N_sample
before merging, restores additivity under the same coefficient/selection scheme.
Neither normalization prescription is implemented without campaign metadata.
Sample cross sections, reference counts, and units are not fully traceable from
code alone. pp uses its stored weight directly and has no such entry-count
normalization. Boundary tests cover every configured pPb interval in both
directions and confirm that the pp range bypass remains unchanged.

### P07 — Medium: a legacy reco/gen shape ratio requests binomial errors

In `macro/plotMcClosures.C`, `hInclusiveJetReco2GenEtaStdBinsPtBins` divides
independently normalized reconstructed and generator eta shapes with `"B"`.
Reco bins are not a pass subset of the same truth-coordinate denominator:
reconstruction introduces migrations, and shape normalization changes both
samples. The binomial statistical model is therefore inappropriate, even
though ROOT accepts the call. This affects displayed uncertainty, not the
central ratio. Other notebook comparison options expose `B`; a configurable
option alone is not a validated efficiency model.

Use standard propagation for a descriptive ratio, or construct a correlated
covariance from shared-event samples for quantitative fits. A true efficiency
requires matched numerator and inclusive denominator defined in the SAME truth
bins. Ordinary independent errors also omit shared-event covariance; switching
away from `B` does not by itself solve that problem.

### P08 — Medium: generic sparse projections retain selection state

`projections.project_histogram` calls `SetRange` on a supplied THnSparse and does
not restore its ranges. A later projection with no `axis_ranges` can silently
reuse the previous selection. A two-entry sparse example gave integrals **1**
for the restricted projection, **1** when the next call omitted ranges, and
**2** after explicitly resetting the source axis. This is a confirmed API state hazard when an
object is reused; current impact depends on the calling workflow. The dedicated
unfolding flattener resets its pT ranges in a `finally` block, so this finding
is about the generic projection helper. A future fix should preserve/restore
all original axis ranges or operate on an independent object.

The generic bin-range helper also swaps indices when the selected range is
empty or reversed. Validate configured intervals and avoid interpreting that
fallback as a physically meaningful empty selection.

## Conventions and statistical limitations requiring analysis judgment

### P05 — pp calibrations and signed vertex diagnostics

`collisionSystem='pp'` disables the eta boost/flip and pPb event weighting/skim
cuts, but does not select a separate pp JEC or resolution fit. The same
pPb/direction-based calibration selection remains, including P01. Reusing it
might be intentional for a controlled comparison, but needs validation for an
absolute pp measurement. The legacy vertex-z sign flip also lacks the pp guard:
pp eta is unflipped while signed vz depends on the direction flag. This does
not change its even |vz| event cut or stored pp weight, but changes diagnostics.

### P06 — pseudorapidity, x estimates, and tail cleaning

`etaCM = etaLab - 0.465` is the conventional shifted-pseudorapidity observable.
An additive longitudinal boost is exact for rapidity, not for pseudorapidity
of a massive jet; see the [PDG kinematics review, section 49.5.2](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-kinematics.pdf). The dijet eta is the arithmetic mean of the two jet etas,
not eta of the summed four-vector. These definitions can be correct for a
measurement of that observable, but should not be labelled an exact massive
four-vector CM transformation.

The stored x estimators use the massless leading-order two-jet relation
`x+/- = sum(pT_i exp(+/-eta_i))/8160`. They are not generator-record Bjorken x:
jet mass, eta versus rapidity, additional radiation, and missing out-of-acceptance
energy are approximations. In pp, the historical xPb key describes the negative-z
proton estimator, not a nuclear momentum fraction.

The MC overweight veto removes entire events before jet ID and fiducial cuts,
including truth events based on reconstructed-jet tails. It is an analysis
sample-cleaning choice with possible detector-dependent truth bias; it is not
a generic QCD kinematic bound. Its applicability to every production/tune
requires a dedicated stability study. Trigger ID 0 has no explicit MB-HLT check;
any MB selection must already be supplied by the chosen data sample.

### P09 — covariance and replica accounting

- Unit-area `TH1::Scale` treats the normalization factor as fixed when scaling
  errors. The fitted normalization is derived from the same bins, so normalized
  shapes have correlations that a TH1 error array does not represent.
- Trigger integral ratios use independent errors; overlapping trigger events
  and chained shared references can introduce covariance. Scaling spectra uses
  central stitching factors only; their uncertainties are common nuisances.
- Generator smearing fills 20 dijet replicas and 10 inclusive replicas at full
  event weight. Mode-2 reco diagnostics use 10 replicas; official hybrid and
  stochastic production use one. Replica normalization must be accounted for
  in yield comparisons, and replicas are not independent collisions.
- Effective miss/fake subtraction uses ROOT's independent `Add` errors although
  the matched population is a subset. Its bin contents are accounting inputs;
  diagnostic errors are not covariance-aware inefficiency estimates.
- Factorized purity and efficiency are fixed MC inputs. Returned unfolded
  covariance propagates measured uncertainty and efficiency scaling; finite-MC
  response uncertainty requires separate toys/variations.
- The compiled toy scan accumulates diagonal ensemble moments, not full
  cross-bin covariance. Nonnegative probability-based response helpers do not
  generally support arbitrary signed-weight unfolding, even though producer
  histograms and hemisphere folding can retain signed weights.
- F/B uses raw yields and independent ROOT division. For a zero denominator,
  ROOT returns content/error zero; that bin is an undefined ratio, not a
  physical measurement of zero. Binned folding retains the previously accepted
  exact-edge differences and cannot reconstruct event-level entries/moments.

These are scope limitations, not evidence that every displayed result is wrong.
They matter when moving from diagnostic shapes to significance, chi-square,
or a final covariance-bearing measurement.

## ROOT documentation checks

Context7 resolved CERN ROOT as `/root-project/root`; queries checked weighted
histogram errors, division compatibility, and Python/C++ ownership. Its returned
snippets were incomplete for `TH1::Divide`, so the official ROOT implementation
was also consulted rather than treating snippets as authoritative formulas.

- [ROOT TH1 implementation and API comments](https://github.com/root-project/root/blob/master/hist/hist/src/TH1.cxx): default division assumes uncorrelated inputs, `B` uses binomial errors, and denominator-zero cells are set to zero.
- [ROOT 6.40 ownership notes](https://github.com/root-project/root/blob/master/README/ReleaseNotes/v640/index.md): retain Python objects referenced by C++, and handle ownership explicitly to avoid deletion errors.

ROOT sources describe library behavior; the physics findings above follow from
this repository's definitions and control flow, not from Context7 certifying an
analysis convention.

## Original review validation and remaining work

Checks run successfully:

- Release configure/build and CTest: 1/1 passed. Two pre-existing unused-parameter
  warnings remain (`isPbGoing`, `eta`).
- Python compileall and the full regression suite: 49/49 passed, with no skips.
- ROOT macro checks: pp/pPb branch/weight/event selection, all seven unfolding
  pair categories in both orientations, and generator x-correlation checks.
- Synthetic review diagnostics reproduced P02, the rank/hash mechanism in P03,
  sparse range persistence in P08, and ROOT zero-denominator ratio behavior.
- Python AST and C++ token comparisons against HEAD confirmed that all 19
  modified source files have unchanged executable syntax.
- All 31 notebook JSON files and ordinary Python code-cell syntax were valid;
  notebook files were not edited.
- `git diff --check` passed.

No full production or data-dependent notebook was executed in the original
annotation-only review; numerical physics fixes required a separate
calibration/input convention decision. The synthetic diagnostics expose existing
behavior rather than claiming it has been repaired. Full production regeneration, external payload certification,
and execution of all data-dependent notebooks are outside this review. Fixing
P01–P04 changes physics outputs and needs its own explicitly validated change;
the follow-up below addresses P02 and only the boundary portion of P04.

## Follow-up validation

The missing-reco and half-open-boundary fixes are executable changes, unlike
the original annotation-only review. Release build, CTest, Python regressions
and ROOT caller/boundary/category/x-correlation checks were rerun. The caller
regression covers 960 combinations, checking both miss families and absence of
measured/fake/response fills. On 1000-entry real-input subsets, all 564 pPb objects were identical. In pp,
559/564 were identical; only the expected five miss/classification objects
changed, each by a weighted integral of 0.00333732014. Selected-event counts
remained 831 (pPb) and 974 (pp), and histogram keys, measured spectra, fake
spectra and matched responses were unchanged. These subsets do not establish
the frequency of exact pT-hat boundaries; synthetic tests exercise those edges.
CTest passed 1/1, Python regressions passed 49/49, and syntax checks passed.
Full production and notebook outputs are not regenerated.
