// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>

// Exercise data nominal/JEU ordering, cuts, both directions and eta boundaries.
// The output can also be compared against the pre-migration producer.
void test_dijet_full_eta_processing(const char *output = "/tmp/fb-data-synthetic.root") {
    TH1::AddDirectory(false);
    Histograms hs;
    createHistograms(hs, false);
    std::vector<GenJet> gen;
    auto jet = [](float pt, float eta, float phi) {
        RecoJet j{};
        j.recoPt = pt;
        j.recoPtJeuUp = pt * 1.1f;
        j.recoPtJeuDown = pt * .9f;
        j.recoEta = eta;
        j.recoPhi = phi;
        j.refPt = -1.f;
        return j;
    };
    for (bool pp : {false, true}) for (bool pbGoing : {false, true}) {
        for (float eta : {-3.5f, -2.5f, -1.5f, -.15f, 0.f, .15f, 1.5f, 2.5f, 3.5f}) {
            for (float pt : {49.f, 50.f, 80.f}) {
                std::vector<RecoJet> reco{jet(pt, eta, 0.f), jet(45.f, eta, 3.14f), jet(30.f, .1f, 1.f)};
                processRecoDijets(pp, pbGoing, false, 1.f, reco, gen, hs,
                                  nullptr, JERSmearingMode::Hybrid, false);
            }
        }
    }
    if (!(hs.hRecoDijetPtEtaCMArr[4]->Integral() > 0.) ||
        !(hs.hRecoDijetPtEtaCMArrJeuUp[4]->Integral() > 0.) ||
        !(hs.hRecoDijetPtEtaCMArrJeuDown[4]->Integral() > 0.))
        throw std::runtime_error("Data nominal/JEU test has no selected dijets");
    TString path(output);
    writeOutput(path, hs, false, JERSmearingMode::Hybrid, false);
    std::cout << "data full eta: nominal and JEU variations filled in both systems/directions\n";
}
