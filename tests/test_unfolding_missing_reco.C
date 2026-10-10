// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>

// Exercise the real caller AFTER reco-ID pruning, not just prepareUnfolding.
// Momentum shifts cannot create a second jet from a 0/1-jet collection, so each
// selected truth pair must be recorded exactly once as a miss in both families.
void test_unfolding_missing_reco() {
    TH1::AddDirectory(false);
    Histograms hs;
    createHistograms(hs, true);
    JERSmearingHelper jer;
    auto check = [](bool pass, const char *message) {
        if (!pass) throw std::runtime_error(message);
    };
    auto reset = [&]() {
        for (int i = 0; i < nEtaCuts; ++i) {
            for (auto *h : {hs.hGenDijetPtEtaCMMissArr[i].get(),
                            hs.hGenDijetPtEtaCMMissJerDefExtraArr[i].get(),
                            hs.hGenDijetPtEtaCMMissControlledSampleArr[i].get(),
                            hs.hGenDijetPtEtaCMMissStudiedSampleArr[i].get(),
                            hs.hRecoDijetPtEtaCMFakeArr[i].get(),
                            hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[i].get(),
                            hs.hRecoDijetPtEtaCMFakeControlledSampleArr[i].get(),
                            hs.hRecoDijetPtEtaCMFakeStudiedSampleArr[i].get(),
                            hs.hRecoDijetPtEtaCMArr[i].get(),
                            hs.hRecoDijetPtEtaCMArrJerDefUnfold[i].get(),
                            hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold[i].get()}) h->Reset();
            for (auto *h : {hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[i].get(),
                            hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[i].get(),
                            hs.hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr[i].get(),
                            hs.hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr[i].get()}) {
                h->Reset(); h->Sumw2();
            }
            hs.hUnfoldingPairClassificationArr[i]->Reset();
            hs.hUnfoldingPairClassificationJerDefExtraArr[i]->Reset();
        }
    };
    const bool originalSplit = isFillControlledSample;
    for (bool pp : {false, true}) for (bool pbGoing : {false, true})
    for (auto mode : {JERSmearingMode::Hybrid, JERSmearingMode::Stochastic,
                      JERSmearingMode::TruthToPseudoReco})
    for (bool etaDependent : {false, true}) for (bool controlled : {false, true})
    for (int nReco : {0, 1}) for (int truthCase : {0, 1, 2, 3, 4})
    for (float weight : {2.5f, -1.25f}) {
        reset();
        isFillControlledSample = controlled;
        // case 0: accepted everywhere; 1: eta acceptance varies with iCut;
        // 2: pT fail; 3: delta-phi fail; 4: no truth pair.
        std::vector<GenJet> gen{{100.f, .2f, 0.f, 100.f},
                                {60.f, -.2f, 3.14f, 60.f}};
        if (truthCase == 1) gen[0].eta = 1.8f;
        if (truthCase == 2) gen[1].pt = 39.f;
        if (truthCase == 3) gen[1].phi = .2f;
        if (truthCase == 4) gen.clear();
        std::vector<RecoJet> reco;
        if (nReco) {
            RecoJet j{};
            j.recoPt = 80.f; j.refPt = 75.f;
            j.recoEta = j.refEta = .2f;
            j.refPtSmeared = -1.f; // no smearing is needed on the missing-pair path
            reco.push_back(j);
        }
        processRecoDijets(pp, pbGoing, true, weight, reco, gen, hs, &jer,
                          mode, etaDependent);
        for (int i = 0; i < nEtaCuts; ++i) {
            const bool pass = gen.size() >= 2 && analysis_contract::passesDijetSelection(
                gen[0].pt, gen[1].pt,
                etaCM(gen[0].eta, analysis_contract::kNominalEtaShift, pbGoing, true, pp),
                etaCM(gen[1].eta, analysis_contract::kNominalEtaShift, pbGoing, true, pp),
                deltaPhi(gen[0].phi, gen[1].phi), etaCuts[i]);
            const double expected = pass ? weight : 0.;
            auto *nominal = hs.hGenDijetPtEtaCMMissArr[i].get();
            auto *extra = hs.hGenDijetPtEtaCMMissJerDefExtraArr[i].get();
            check(nominal->Integral() == expected && extra->Integral() == expected,
                  "caller must fill both missing-truth families exactly once");
            double variance = 0.;
            for (int b = 0; b < nominal->GetNcells(); ++b)
                variance += nominal->GetBinError(b) * nominal->GetBinError(b);
            check(std::abs(variance - (pass ? weight * weight : 0.)) < 1.e-12,
                  "missing-truth Sumw2 must retain the event weight exactly once");
            auto *chosen = controlled ? hs.hGenDijetPtEtaCMMissControlledSampleArr[i].get()
                                      : hs.hGenDijetPtEtaCMMissStudiedSampleArr[i].get();
            auto *other = controlled ? hs.hGenDijetPtEtaCMMissStudiedSampleArr[i].get()
                                     : hs.hGenDijetPtEtaCMMissControlledSampleArr[i].get();
            check(chosen->Integral() == expected && other->Integral() == 0.,
                  "missing truth must enter only the selected event subsample");
            for (auto *classification : {hs.hUnfoldingPairClassificationArr[i].get(),
                                         hs.hUnfoldingPairClassificationJerDefExtraArr[i].get()}) {
                check(classification->Integral() == expected &&
                      classification->GetBinContent(unfolding_diagnostics::kGenPassRecoFail) == expected,
                      "missing reco pair must classify as gen-pass/reco-fail");
            }
            check(hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[i]->GetEntries() == 0. &&
                  hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[i]->GetEntries() == 0. &&
                  hs.hRecoDijetPtEtaCMFakeArr[i]->Integral() == 0. &&
                  hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[i]->Integral() == 0. &&
                  hs.hRecoDijetPtEtaCMArr[i]->Integral() == 0. &&
                  hs.hRecoDijetPtEtaCMArrJerDefUnfold[i]->Integral() == 0. &&
                  hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold[i]->Integral() == 0.,
                  "missing reco pair must not add response, measured or fake yields");
        }
    }
    // Data must not touch unallocated MC objects on the early-return path.
    reset();
    Histograms dataHs; // MC pointers intentionally remain null
    for (int nReco : {0, 1}) {
        std::vector<GenJet> gen;
        std::vector<RecoJet> reco(nReco);
        processRecoDijets(false, false, false, 1.f, reco, gen, dataHs, nullptr,
                          JERSmearingMode::Hybrid, false);
        check(hs.hGenDijetPtEtaCMMissArr[4]->Integral() == 0., "data acquired truth misses");
    }
    isFillControlledSample = originalSplit;
    std::cout << "unfolding caller: 0/1 reco jets, both response families, cuts, modes, splits and weights passed\n";
}
