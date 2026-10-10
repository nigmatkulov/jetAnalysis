// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>
void test_unfolding_processing() {
    TH1::AddDirectory(false);
    Histograms hs;
    createHistograms(hs, true);
    auto reco = [](float pt, float eta, float phi, bool matched = true) {
        RecoJet jet{};
        jet.recoPt = jet.refPtSmeared = pt;
        jet.recoEta = jet.refEta = eta;
        jet.recoPhi = jet.refPhi = phi;
        jet.refPt = matched ? pt : -1.f;
        return jet;
    };
    auto check = [](bool pass, const char *message) {
        if (!pass) throw std::runtime_error(message);
    };
    for (bool pp : {false, true}) for (bool pbGoing : {false, true}) {
        for (int category = 1; category <= unfolding_diagnostics::nPairCategories; ++category) {
            std::vector<GenJet> gen{{100.f, .2f, 0.f, 100.f}, {60.f, -.2f, 3.14f, 60.f}};
            std::vector<RecoJet> measured{reco(100.f,.2f,0.f), reco(60.f,-.2f,3.14f)};
            if (category == unfolding_diagnostics::kMatchedSwapped) {
                measured[0].refPtSmeared = 60.f;
                measured[1].refPtSmeared = 100.f;
            }
            if (category == unfolding_diagnostics::kGenPassRecoFail) measured.clear();
            if (category == unfolding_diagnostics::kRecoPassGenFail) gen.clear();
            if (category == unfolding_diagnostics::kSelectedPairMismatch) measured[0].recoEta = 1.f;
            if (category == unfolding_diagnostics::kMissingOneRecoReference) measured[0].refPt = -1.f;
            if (category == unfolding_diagnostics::kMissingBothRecoReferences)
                measured[0].refPt = measured[1].refPt = -1.f;
            for (int i = 0; i < nEtaCuts; ++i) {
                hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[i]->Reset();
                hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[i]->Sumw2();
                hs.hGenDijetPtEtaCMMissArr[i]->Reset();
                hs.hRecoDijetPtEtaCMFakeArr[i]->Reset();
                hs.hUnfoldingPairClassificationArr[i]->Reset();
            }
            prepareUnfolding(pp, pbGoing, true, 2.5f, gen, measured,
                hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr, hs.hGenDijetPtEtaCMMissArr,
                hs.hRecoDijetPtEtaCMFakeArr, hs.hUnfoldingPairClassificationArr);
            const bool matched = category == unfolding_diagnostics::kMatchedDirect ||
                                 category == unfolding_diagnostics::kMatchedSwapped;
            auto *response = hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[5].get();
            check(response->GetEntries() == (matched ? 1. : 0.), "response membership incorrect");
            check(hs.hGenDijetPtEtaCMMissArr[5]->Integral() ==
                  (!matched && category != unfolding_diagnostics::kRecoPassGenFail ? 2.5 : 0.),
                  "truth miss membership incorrect");
            check(hs.hRecoDijetPtEtaCMFakeArr[5]->Integral() ==
                  (!matched && category != unfolding_diagnostics::kGenPassRecoFail ? 2.5 : 0.),
                  "measured fake membership incorrect");
            check(hs.hUnfoldingPairClassificationArr[5]->GetBinContent(category) == 2.5,
                  "pair category or weight incorrect");
            if (matched) {
                Int_t coordinates[4];
                check(response->GetBinContent(Long64_t(0), coordinates) == 2.5,
                      "response weight incorrect");
                check(response->GetBinError2(Long64_t(0)) == 6.25, "response error incorrect");
            }
        }
    }
    std::cout << "unfolding: all seven pair categories and both CM orientations passed\n";
}
