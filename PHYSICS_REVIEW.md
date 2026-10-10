# Physics and implementation review

Reviewed 2026-10-10 against the working version following `4a8ab4c`.
This review adds comments and documentation only. Selections, calibration files,
weights, random numbers, histogram keys, and numerical algorithms are unchanged.
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
No calibration was switched in this review.

### P02 — High: the caller drops truth misses with fewer than two selected reco jets

`processRecoJets` removes ID-failing jets. `processRecoDijets` then returns
immediately when fewer than two remain, before either call to
`prepareUnfolding`. A selected generator dijet with zero or one selected reco
jet is therefore missing from the explicit miss and classification objects.
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
stored misses can underestimate inefficiency. A fix must invoke response
accounting on every MC event without indexing a missing reco pair, then check
`truth = matched truth + misses` and `measured = matched reco + fakes` over the
same acceptance.

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

Additionally, pPb event selection includes `[low,high]`, whereas weight-table
intervals use `(low,high]`. An entry exactly at a lower boundary uses the
preceding weight interval (pT-hat exactly 15 uses the fallback coefficient).
Whether such events occur, or adjacent generated samples overlap there, was
not established. Sample cross sections, reference counts, and units are not
fully traceable from code alone. Validate campaign metadata and boundary
conventions before using absolute yields. pp uses its stored weight directly
and does not have this entry-count normalization.

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

## Validation and remaining work

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

No full production or data-dependent notebook was executed: this change has no
executable edits, and confirming numerical physics fixes requires a separate
calibration/input convention decision. The synthetic diagnostics expose existing
behavior rather than claiming it has been repaired. Full production regeneration, external payload certification,
and execution of all data-dependent notebooks are outside this review. Fixing
P01–P04 changes physics outputs and needs its own explicitly validated change.
