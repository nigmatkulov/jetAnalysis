# Step 00: Analysis contract and current-code audit

**Status:** In progress. Contract tests, demonstrated correctness fixes, and one
Pb-going embedded-MC execution are complete. Three analysis channels, full
histogram-schema review, and response-comparison closure remain outstanding.

## Goal and problem

Create one authoritative specification for coordinate conventions, inputs,
selections, matching, weights, and observables. Direction mistakes cannot be
repaired safely after the two beam orientations are combined, especially because
some MC and data productions were stored with flipped orientations.

## Facts from the current repository

- `processForestSimple.C` reads `ak4PFJetAnalyzer/t` and applies AK4PF JEC files
  selected by data/MC type and direction.
- `etaCM(...)` has separate data and MC sign branches and uses 0.465 in active
  dijet processing; a scan from 0.463 to 0.495 also exists.
- The agreed nominal cuts are present in the dijet functions.
- Multiple jet CM acceptances, including the nominal 1.9, are filled.
- Reference-dijet duplicate fills found during the audit have been removed, and
  the missing `refphi` branch enable was found and corrected in the working tree.
- Both response construction from generator jets and an alternative reference-jet
  response appear in the working copy. Their intended roles are not documented.

## Recommended solution

1. Write a convention table for each sample: physical direction, stored detector
   orientation, sign operation, CM boost, JEC set, and expected proton-going sign.
2. Express the observable and every cut mathematically and in code-level terms,
   including interval-boundary conventions.
3. Trace one synthetic or hand-selected event through every direction/data-MC
   branch and compare the computed signs with the convention table.
4. Audit nominal generator, reference, and reconstructed dijet definitions,
   including ordering, angular cuts, matching, misses, fakes, and double fills.
5. Record a histogram schema: name, axes, units, selection, weight, and purpose.
6. Treat the existing `etaLab(...)` and `etaCM(...)` implementations as the
   authoritative starting convention. Derive their four data/MC × beam-direction
   transformations explicitly and test that the rest of the processing follows
   those transformations consistently.
7. Implement and test both response interpretations during the audit:
   independently clustered generator dijets and reference-matched reconstructed
   dijets. Do not select the nominal construction until their efficiencies,
   matching assumptions, and independent closure behavior have been compared.
8. Define separate configurations and output namespaces for Pb-going MC, p-going
   MC, Pb-going data, and p-going data. Explicitly prohibit direction combination
   until the relevant direction-specific corrections and validations are complete.

This step should correct only demonstrated defects needed to establish the
contract. Broader refactoring is excluded.

## Validation and exit criteria

- Unit-level convention tests cover data/MC × p-going/Pb-going and pass.
- Configuration and output-schema checks distinguish all four analysis channels.
- A histogram audit finds no unintended duplicate fills.
- Generator, reconstructed, and response selections are explicitly compared.
- Nominal acceptance is index/value 1.9 without deleting diagnostic acceptances.
- The analysis contract and histogram schema are approved by the analyst.

## Required figures

- Direction-overlay sanity plots before and after orientation transformation,
  saved as PDF and PNG.
- A small signed-eta event example or diagram, saved as PDF and PNG.

## Decisions recorded

1. The current `etaCM(...)` and `etaLab(...)` functions are the reference for the
   analysis orientation. Their behavior must be documented and tested rather than
   replaced using an inferred filename convention.
2. Both independently clustered generator dijets and reconstructed-jet reference
   matches will be tested. The nominal response will be selected only after the
   scientific comparison documented in this step.

## Remaining inputs

No additional design decision is required to start this step. If the audit finds
that actual files contradict the reference functions, the discrepancy will be
documented with event-level evidence and returned to the analyst for a decision;
it will not be resolved by guessing from filenames.

## Results

Plan approved on 2026-07-14. Execution began on 2026-07-14.

### Authoritative orientation contract

The stored jet eta is passed directly to the transformations. Positive final
`etaCM` is defined as proton-going.

| Channel | `etaLab(eta)` | `etaCM(eta, shift)` | JEC selected by current code |
|---|---:|---:|---|
| Pb-going MC | `-eta` | `-(eta + shift)` | Pb-going embedded/unembedded MC |
| p-going MC | `eta` | `eta - shift` | p-going embedded/unembedded MC |
| Pb-going data | `eta` | `eta - shift` | p-going embedded MC L2Relative + data residual |
| p-going data | `-eta` | `-(eta + shift)` | Pb-going embedded MC L2Relative + data residual |

The nominal shift is 0.465. The data JEC direction reversal is retained as an
explicit property of the current implementation and is not inferred from file
names.

### Observable and nominal boundaries

The leading and subleading jets are ordered by the pT used by the applicable
variation. A pair passes when

`pT1 >= 50 GeV`, `pT2 >= 40 GeV`, `abs(deltaPhi) >= 2*pi/3`, and both
`abs(etaCM) <= etaMax`. The nominal `etaMax` is 1.9; the other values in
`etaCuts` remain diagnostics. Then `pTave=(pT1+pT2)/2` and
`etaDijet=(etaCM1+etaCM2)/2`.

### Audit findings and corrections

- `processRefDijets` filled full and forward/backward histograms twice for every
  qualifying pair. The duplicate fills were removed.
- `prepareUnfolding` did not compile, used a 50/50 threshold in one helper,
  compared subleading pT with an eta limit, referred to undefined variables,
  and never set the booleans used by its response fill. It now uses the shared
  nominal boundary function and explicit direct/swapped matching.
- Independently clustered generator and reconstructed-reference response
  interpretations remain separate. No nominal interpretation is selected.
- The response schema is recorded in
  [histogram_schema.csv](histogram_schema.csv). A complete inventory of legacy
  diagnostic histograms remains future work within this step.

### Validation record

- Base Git commit tested: `1288468` (`Code is cleaned up. Test is added`), plus
  the working-tree correction that enables the `refphi` branch.
- ROOT: 6.40.02.
- Baseline cached build: failed because it referenced removed Homebrew GCC 15
  executables.
- AppleClang Release configuration and build in `build/`: passed on 2026-07-14.
- `ctest --test-dir build --output-on-failure`: one of one tests passed. The
  focused test covers nominal pT, inclusive eta, and inclusive delta-phi
  boundaries. The authoritative `etaLab` and `etaCM` definitions remain in
  `processForestSimple.C`.
- Scoped Release execution through the configuration and runner in
  `pPb8160_analyzeMc_Pbgoing.py`: passed after running all 536,200 events in the
  Pb-going embedded ptHat-120 sample; 425,448 passed event selection and the
  configured ROOT output was written. Runtime was 2 minutes 39 seconds. The
  committed script currently lists all pt-hat samples, so the validation invoked
  its runner with `[120]` to avoid starting an unapproved full production.

Input provenance:

- Path: `/Users/gnigmat/cernbox/ana/pPb8160/embedding/Pbgoing/forest/HiForestSkim_pPb_MC_pthat120_Pbgoing_embedded.root`
- ROOT UUID: `13cc371e-6e47-11f0-9520-d1338e80beef`
- Size: 1,833,795,295 bytes; entries: 536,200.
- Configuration: embedding MC, Pb-going, ptHat 120, JEU off, JER systematics on,
  no trigger requirement, `jetIdTightLeptVeto`.

Output provenance:

- Path: `/Users/gnigmat/cernbox/ana/pPb8160/embedding/Pbgoing/embedding_Pbgoing_ptHat120.root`
- ROOT UUID: `ea714fef-f29b-47df-b915-949cf3ef0562`
- Size: 15,599,740 bytes.

Nominal eta-index 5 (`etaCut=1.9`) object audit:

| Object | Entries |
|---|---:|
| Generator dijets | 244,206 |
| Reconstructed dijets | 248,013 |
| Reconstructed-associated reference dijets | 248,000 |
| Independently clustered generator response | 233,001 |
| Generator misses | 11,205 |
| Reference-selected dijets | 248,066 |
| Reference response | 245,342 |

The initial zero reference-selected yield was explained by branch setup:
`refphi` had a `SetBranchAddress` call but no corresponding
`SetBranchStatus("refphi", 1)`. Its value therefore remained zero and all
reference pairs failed the delta-phi selection. Enabling the branch and rerunning
Release populated both reference histograms. This validates the defect diagnosis
but does not select either response construction as nominal.

### Exit-criteria assessment

| Criterion | Result |
|---|---|
| Four data/MC x direction convention tests | Not met: transformations are documented, but the unit test currently covers selection boundaries only |
| Four distinct channel configurations and outputs | Not met: only Pb-going embedded MC has been executed |
| No unintended duplicate fills | Partially met: known reference duplicates were removed; full histogram audit remains incomplete |
| Gen/reco/reference/response selections compared | Partially met: entry counts are recorded and both responses are populated, but efficiencies and independent closure are not yet compared |
| Nominal acceptance is index/value 1.9 | Met for the inspected output: index 5 is 1.9 and diagnostic acceptances remain |
| Contract and histogram schema approved | Contract/code changes approved by the analyst; full legacy histogram schema remains incomplete |

Synthetic contract figures:

- [Direction overlay](figures/direction_overlay.pdf) ([PNG](figures/direction_overlay.png))
- [Signed-eta example](figures/signed_eta_event.pdf) ([PNG](figures/signed_eta_event.png))

Step 00 is not complete until the remaining three direction/data channels are
run, the full histogram schema is reviewed, and both response interpretations
receive event-level efficiency and closure comparisons.
