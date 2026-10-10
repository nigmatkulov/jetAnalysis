# Physics and implementation review

Reviewed 2026-10-10 against the working version following `4a8ab4c`.
The original review added comments and documentation only. A subsequent fix
restores missing-reco truth accounting (P02) and changes pPb pT-hat event and
weight intervals to `[low, high)` (boundary portion of P04). A further fix
caches JER momenta with stable draw identity across variations and replicas
(P03). Calibration files, histogram keys and the job-normalization denominator
are unchanged; stochastic outputs after rank swaps can change.
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

### P03 — Fixed: sorting changed the random-draw identity

Before the fix, the JER hash excluded the variation index, but its `round`
included `iJet`. `prepareUnfolding`, called at scale indices 4 and 7, sorted the
shared reco collection by varied momentum. If a pair changed rank, the same
physical jet got a different index in subsequent variations and therefore a
different normal deviate. The claimed common draw across all Down/Default/Up
templates was not guaranteed. Hybrid jets with good matches do not use the deviate, but
unmatched/stochastic jets and pseudo-reco diagnostic replicas can be affected.

A synthetic hybrid pair with reco/reference pT 100/120 and 99/81 GeV swaps rank
under the nominal central SF. The collection remained reordered after the
systematic loop. For the original 100 GeV jet in event `(run,lumi,evt)=(1,1,17)`,
changing the hash round from 0 to 1 changed the deviate from **-0.726256** to
**-0.802291**. This demonstrates the two ingredients of the defect, without
claiming a measured production-spectrum difference.

The implementation now freezes each jet's initial corrected-pT index and
computes one normal deviate per jet/replica before sorting. `RecoJet` stores
`refPtSmeared` (base unit SF), `refPtSmearedEta` (eta-dependent unit SF), and
`refPtSmearedJerDef/Up/Down` (configured-resolution SF variations). An internal
nine-value cache retains every historical template, including x2, x0.9 and
eta-dependent Down/Default/Up. All values travel with their jet. Each variation
still independently selects its leading pair; unfolding explicitly selects its
nominal or eta-dependent cache entry. The default direct-helper interface still
uses `refPtSmeared`. Mode 2 refreshes the cache per replica, retaining stable jet
identity across all ten replicas. No resolution constants or match cuts change.

`tests/test_jer_cached_momenta.C` exercises the actual caller, rank swaps,
unmatched jets, all modes, both resolution configurations, and the final
pseudo-reco replica. It compares every cached branch with the stable-identity
draw and verifies that response selection preserves the cache. Histogram keys,
binning and output paths are unchanged. Stochastic branches after a rank swap
can change numerically, as intended; these are common-draw corrections, not new
random seeds or resolution inputs.

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

## JER implementation and full pT-hat-80 audit (2026-10-10)

### Input and calibration routing (P01)

Ran the compiled executable on the complete, existing Pb-going pPb pT-hat-80
files, with trigger 0, jet-ID method 2, vertex filter 0, hybrid JER, and base
resolution. Outputs and logs were written to `/tmp`; production inputs and
histogram outputs were not overwritten.

| Mode | Input below `$HOME/cernbox/ana/pPb8160/` | Entries | Selected |
| --- | --- | ---: | ---: |
| Embedding (`mcType=1`) | `embedding/Pbgoing/forest/HiForestSkim_pPb_MC_pthat80_Pbgoing_embedded.root` | 783686 | 639706 |
| PYTHIA (`mcType=2`) | `pythia/Pbgoing/forest/HiForestSkim_pPb_MC_pthat80_Pbgoing_unembedded.root` | 984016 | 802116 |

Each input was accepted into the chains (accepted=1, skipped=0), and both runs
finished. The `JEC files` log selected
`aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt`
in BOTH runs. This is correct for embedding and inconsistent with the existing
Unembedded branch for standalone PYTHIA. The input routing is correct; P01 is
confirmed as a calibration-selector defect. The logged vector is passed to
`JetCorrector`, whose initializer constructs `SingleJetCorrector` for each file.
`tests/test_jec_payloads.C` independently loaded/evaluated both real Pb-going
payloads and reproduced the distinct correction values in P01. This is not a
validation of their physical suitability. No JEC selector or payload changed.
Both forests lacked a supported rho branch and used the documented rho=0 hash.

Reproduce the new-version runs from the repository root, choosing mcType 1 or 2:

```bash
build/processForestSimple "$HOME/cernbox/ana/pPb8160/embedding/Pbgoing/forest/HiForestSkim_pPb_MC_pthat80_Pbgoing_embedded.root" /tmp/jer-embedding-full.root 1 1 80 0 2 0 pPb
build/processForestSimple "$HOME/cernbox/ana/pPb8160/pythia/Pbgoing/forest/HiForestSkim_pPb_MC_pthat80_Pbgoing_unembedded.root" /tmp/jer-pythia-full.root 2 1 80 0 2 0 pPb
```

### Size of the previously omitted misses (P02)

Compared full-sample outputs with the pre-fix executable (`7602717` behavior).
The nominal truth, nominal unfolding-measured spectrum, and nominal matched
response were identical, including errors, at every eta-cut index. Thus the
nominal factorized workflow's efficiency and unfolding inputs do not change:
it already uses matched/inclusive truth, recovering missing events indirectly.
Stored misses and gen-pass/reco-fail diagnostics were incomplete instead.

The table integrates all pT/eta bins, including flow bins, of the nominal
objects. Percentages are weighted; the omitted counts are unweighted events.

| Sample | Jet `|etaCM|` cut | Omitted events | Omitted / selected truth | Omitted / corrected misses |
| --- | ---: | ---: | ---: | ---: |
| Embedding | 1.9 | 0 | 0% | 0% |
| Embedding | 3.0 | 11 | 0.002215% | 0.013745% |
| PYTHIA | 1.9 | 3 | 0.0007212% | 0.0088849% |
| PYTHIA | 3.0 | 231 | 0.036109% | 0.207226% |

For matched yield M, correct misses L, and omitted misses D, the correct
single-bin efficiency is M/(M+L), whereas an explicit-miss consumer previously
obtained M/(M+L-D). Dividing a measured yield by the latter undercorrects it by
D/(M+L). This gives a per-bin or global-efficiency estimate, not a measured
shift from a full multidimensional unfolding. For broad-acceptance PYTHIA, the
global efficiency changes from 0.826049 to 0.825750 under that explicit-miss
prescription. The current matched/inclusive-truth workflow already has the
latter value, subject to histogram floating-point precision.

Small integrated effects can hide larger shape effects. After integrating all
pT and grouping dijet eta into width-0.2 bins, the largest omitted/truth fraction
among bins with effective truth entries >=100 was at eta `[2.8,3.0)` for the
3.0 jet acceptance: 0.3559% for embedding and 3.4622% for PYTHIA. At acceptance
1.9 the maxima were 0% and 0.004202%, respectively. These are generator dijet
CM eta bins. Sparse pT/eta tails can behave differently. This one pT-hat sample
and the existing Embedded JEC choice do not establish an all-sample measurement
bias. No full notebook unfolding was rerun or inferred to have shifted.

### JER regression and output differences (P03)

Release configure/build passed (only the two existing unused-parameter
warnings). CTest passed 1/1; Python regressions passed 49/49; compileall and
`git diff --check` passed. ROOT tests passed for cached JER momenta, missing-reco
accounting, all seven response categories, pp/pPb weights and boundaries,
x-correlations, and the two real JEC payloads. Context7 ROOT documentation was
consulted for sparse-projection behavior; numerical conclusions use the local
histograms and implementation.

Compared the full embedding sample with the pre-JER-fix executable (`f89754e`):
483/564 objects were exactly identical. The 81 changed objects were confined to
JER Up and eta-dependent JER templates/response/miss/fake/classification objects,
which previously could use draws reassigned after nominal sorting. Nominal
response, nominal measured/truth spectra, baseline spectra, event selection,
and all keys were unchanged. On the 1000-entry pp subset, unchanged object
counts were 557/564 (hybrid), 483/564 (stochastic), and 422/564 (pseudo-reco).
Pseudo-reco diagnostic templates can change in later replicas too, since their
jet identity is now stable across replicas; its first-replica nominal unfolding
response remains unchanged. No data-dependent notebook was executed, and no
production histogram files were replaced.
