// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>
void test_x_correlations() {
    TH1::AddDirectory(false);
    Histograms hs;
    createHistograms(hs, true);
    JERSmearingHelper jer;
    auto check = [](bool pass, const char *message) {
        if (!pass) throw std::runtime_error(message);
    };
    for (bool pp : {false, true}) {
        for (bool pbGoing : {false, true}) {
            for (int i = 0; i < nEtaCuts; ++i) {
                hs.hGenDijetPtEtaCMArr[i]->Reset();
                hs.hGenDijetPtEtaCMXPbArr[i]->Reset();
                hs.hGenDijetPtEtaCMXPbArr[i]->Sumw2();
                hs.hGenDijetPtEtaCMXpArr[i]->Reset();
                hs.hGenDijetPtEtaCMXpArr[i]->Sumw2();
                hs.hGenDijetPtEtaCMXPbOverXpArr[i]->Reset();
                hs.hGenDijetPtEtaCMXPbOverXpArr[i]->Sumw2();
            }
            // Unequal momenta catch using pTave in place of individual momenta.
            std::vector<GenJet> jets{{100.f, .8f, 0.f, 100.f},
                                     {60.f, -.2f, 3.14f, 60.f}};
            processGenDijets(pp, pbGoing, true, 2.5, jets, hs, &jer);
            const double eta1 = etaCM(.8f, analysis_contract::kNominalEtaShift, pbGoing, true, pp);
            const double eta2 = etaCM(-.2f, analysis_contract::kNominalEtaShift, pbGoing, true, pp);
            const double xp = (100. * std::exp(eta1) + 60. * std::exp(eta2)) / 8160.;
            const double xPb = (100. * std::exp(-eta1) + 60. * std::exp(-eta2)) / 8160.;
            const double coordinates[3] = {80., (eta1 + eta2) / 2., xPb};
            auto *h = hs.hGenDijetPtEtaCMXPbArr[5].get();
            check(h->GetBinContent(h->GetBin(coordinates)) == 2.5, "xPb coordinate or weight incorrect");
            check(std::abs(h->GetBinError(h->GetBin(coordinates)) - 2.5) < 1.e-12, "xPb Sumw2 incorrect");
            double values[3] = {80., (eta1 + eta2) / 2., xp};
            h = hs.hGenDijetPtEtaCMXpArr[5].get();
            check(h->GetBinContent(h->GetBin(values)) == 2.5, "xp coordinate incorrect");
            values[2] = xPb / xp;
            h = hs.hGenDijetPtEtaCMXPbOverXpArr[5].get();
            check(h->GetBinContent(h->GetBin(values)) == 2.5, "ratio coordinate incorrect");
            check(h->GetEntries() == 1., "x must be filled once, without smearing replicas");
            check(hs.hGenDijetPtEtaCMArr[5]->Integral() == 2.5, "nominal yield differs from x yield");
            jets = {{49.f, .8f, 0.f, 49.f}, {40.f, -.2f, 3.14f, 40.f}};
            processGenDijets(pp, pbGoing, true, 2.5, jets, hs, &jer);
            check(h->GetEntries() == 1., "failed dijet selection entered x histogram");
        }
    }
    TString output = Form("%s/jetAnalysis-x-validation-%d.root",
                          gSystem->TempDirectory(), gSystem->GetPid());
    writeOutput(output, hs, true, JERSmearingMode::Hybrid, false);
    TFile file(output, "READ");
    for (int i = 0; i < nEtaCuts; ++i) {
        for (const char *base : {"hGenDijetPtEtaCMXPb", "hGenDijetPtEtaCMXp",
                                 "hGenDijetPtEtaCMXPbOverXp"}) {
            auto *saved = dynamic_cast<THnSparseD *>(file.Get(Form("%s_%d", base, i)));
            check(saved && saved->GetNdimensions() == 3, "serialized x key or dimensions incorrect");
        }
    }
    file.Close();
    gSystem->Unlink(output);
    std::cout << "x correlations: coordinates, weights, errors, frames and selection passed\n";
}
