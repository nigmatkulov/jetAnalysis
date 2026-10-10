// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>

void test_jer_cached_momenta() {
    TH1::AddDirectory(false);
    Histograms hs;
    createHistograms(hs, true);
    JERSmearingHelper jer;
    auto etaScale = std::make_unique<TF1>("testCachedEtaScale", "1.7", -5., 5.);
    jer.addJERSmearingFuncVsEta(0., 10000., std::move(etaScale));
    run = 1; lumi = 1; evt = 17; eventRho = 12.f;
    auto check = [](bool pass, const char *message) {
        if (!pass) throw std::runtime_error(message);
    };
    auto makeJet = [](float pt, float ref, float eta, float phi) {
        RecoJet j{};
        j.recoPt = pt; j.refPt = ref;
        j.recoEta = j.refEta = eta;
        j.recoPhi = j.refPhi = phi;
        return j;
    };
    // The matched central pair swaps rank under nominal hybrid scaling.
    // Unmatched and forward jets additionally exercise stochastic/eta branches.
    const std::vector<RecoJet> initial{
        makeJet(100.f, 120.f, .2f, 0.f),
        makeJet(99.f, 81.f, -.2f, 3.14f),
        makeJet(70.f, -1.f, 1.2f, 1.f),
        makeJet(65.f, 60.f, -1.4f, 2.f),
        makeJet(60.f, 5.f, .3f, 1.5f)};
    for (auto mode : {JERSmearingMode::Hybrid, JERSmearingMode::Stochastic,
                      JERSmearingMode::TruthToPseudoReco})
    for (bool etaDependent : {false, true}) {
        auto reco = initial;
        std::vector<GenJet> gen{{120.f, .2f, 0.f, 120.f}, {81.f, -.2f, 3.14f, 81.f}};
        processRecoDijets(true, true, true, 1.f, reco, gen, hs, &jer, mode, etaDependent);
        const int replica = mode == JERSmearingMode::TruthToPseudoReco ? 9 : 0;
        for (const auto &jet : reco) {
            check(jet.smearingIndex >= 0 && jet.smearingIndex < static_cast<int>(initial.size()), "lost jet draw identity");
            const auto &original = initial[jet.smearingIndex];
            check(original.recoPt == jet.recoPt, "sorting changed draw identity");
            const double normal = jer.normalDeviate(original.recoPt, original.recoEta,
                eventRho, run, lumi, evt, replica * static_cast<int>(initial.size()) + jet.smearingIndex);
            const bool hasReference = original.refPt > 0.;
            auto expected = [&](double sf, bool eta) {
                if (!hasReference && mode == JERSmearingMode::TruthToPseudoReco) return -1.f;
                const bool goodMatch = hasReference &&
                    std::abs(original.recoPt-original.refPt) <
                    3. * jer.resolution(original.recoPt, original.recoEta, eta) * original.recoPt;
                return jer.smearRecoMomentum(original.recoPt,
                    hasReference ? original.refPt : original.recoPt, original.recoEta,
                    sf, eta, mode, goodMatch, normal);
            };
            check(jet.refPtSmeared == expected(1., false), "base unit-SF field changed meaning");
            check(jet.refPtSmearedEta == expected(1., true), "eta unit-SF field incorrect");
            for (int scale = 0; scale < 9; ++scale) {
                const int variation = (scale == 3 || scale == 6) ? -1 :
                                      (scale == 4 || scale == 7) ? 0 : 1;
                const double sf = scale < 3 ? (scale == 1 ? 2. : scale == 2 ? 0.9f : 1.) :
                    jer.JERScaleFactorSystematics(original.recoEta, variation);
                check(jet.jerMomenta[scale] == expected(sf, etaDependent || scale >= 6),
                      "variation or replica used a different physical jet's draw");
            }
            check(jet.refPtSmearedJerDown == jet.jerMomenta[3] &&
                  jet.refPtSmearedJerDef == jet.jerMomenta[4] &&
                  jet.refPtSmearedJerUp == jet.jerMomenta[5], "named SF fields disagree with templates");
        }
        if (mode == JERSmearingMode::Hybrid) {
            check(reco[0].recoPt == 99.f, "test did not exercise a rank swap");
            // A response helper call must not overwrite any of the cached fields.
            const auto saved = reco;
            prepareUnfolding(true, true, true, 1.f, gen, reco,
                hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr, hs.hGenDijetPtEtaCMMissArr,
                hs.hRecoDijetPtEtaCMFakeArr, hs.hUnfoldingPairClassificationArr,
                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 4);
            for (const auto &jet : reco) {
                const auto found = std::find_if(saved.begin(), saved.end(), [&](const RecoJet &j) {
                    return j.smearingIndex == jet.smearingIndex;
                });
                check(found != saved.end() && found->jerMomenta == jet.jerMomenta,
                      "response selection mutated cached momenta");
            }
        }
    }
    std::cout << "JER cache: shared draws, rank swaps, unmatched jets, eta branches, all modes and replicas passed\n";
}
