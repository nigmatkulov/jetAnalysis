#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <vector>

//________________
void pPb8160_jerScalingFactors() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    double fJerEtaBinEdges[] = {-5.191, -3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139, 5.191};
    const int fJerNEtaBins = sizeof(fJerEtaBinEdges) / sizeof(double) - 1;
    for (int i = 0; i < fJerNEtaBins; ++i) {
        if ( fJerEtaBinEdges[i] >= fJerEtaBinEdges[i + 1] ) {
            std::cerr << Form("Error: Eta bin edges are not in increasing order! Bin %d: %f, Bin %d: %f", i, fJerEtaBinEdges[i], i + 1, fJerEtaBinEdges[i + 1]) << std::endl;
            return;
        }
    }
    std::vector<float> fJerDef{1.1922, 1.1869, 1.7788, 1.3418, 1.2963, 1.1512, 1.1426, 1.1000, 1.1278, 1.1609, 1.1464, 1.1948, 1.15958, 1.15958, 1.1948, 1.1464, 1.1609, 1.1278, 1.1000, 1.1426, 1.1512, 1.2963, 1.3418, 1.7788, 1.1869, 1.1922};
    std::vector<float> fJerLow{1.0434, 1.0626, 1.578, 1.1327, 1.0592, 1.0372, 1.0212, 0.9921, 1.0292, 1.0584, 1.0832, 1.1296, 1.095, 1.095, 1.1296, 1.0832, 1.0584, 1.0292, 0.9921, 1.0212, 1.0372, 1.0592, 1.1327, 1.578, 1.0626, 1.0434};
    std::vector<float> fJerHi{1.341, 1.3112, 1.9796, 1.5509, 1.5334, 1.2652, 1.264, 1.2079, 1.2264, 1.2634, 1.2096, 1.26, 1.224, 1.224, 1.26, 1.2096, 1.2634, 1.2264, 1.2079, 1.264, 1.2652, 1.5334, 1.5509, 1.9796, 1.3112, 1.341};

    // std::cout << Form("Number of eta bins: %d", fJerNEtaBins) << std::endl;
    // std::cout << Form("Number of def SFs: %d", (int)fJerDef.size()) << std::endl;
    // std::cout << Form("Number of low SFs: %d", (int)fJerLow.size()) << std::endl;
    // std::cout << Form("Number of hi SFs: %d", (int)fJerHi.size()) << std::endl;

    auto hJerSFDef = std::make_unique<TH1D>("hJerSFDdef", "hJerSFDdef;#eta;Scaling Factor", fJerNEtaBins, fJerEtaBinEdges);
    hJerSFDef->Sumw2();
    auto hJerSFLow = std::make_unique<TH1D>("hJerSFLow", "hJerSFLow;#eta;Scaling Factor", fJerNEtaBins, fJerEtaBinEdges);
    hJerSFLow->Sumw2();
    auto hJerSFHi = std::make_unique<TH1D>("hJerSFHi", "hJerSFHi;#eta;Scaling Factor", fJerNEtaBins, fJerEtaBinEdges);
    hJerSFHi->Sumw2();

    hJerSFDef->SetMarkerStyle(20); hJerSFDef->SetMarkerColor(kBlack); hJerSFDef->SetLineColor(kBlack); hJerSFDef->SetLineWidth(2); hJerSFDef->SetMarkerSize(1.4);
    hJerSFLow->SetMarkerStyle(21); hJerSFLow->SetMarkerColor(kRed); hJerSFLow->SetLineColor(kRed); hJerSFLow->SetLineWidth(2); hJerSFLow->SetMarkerSize(1.4);
    hJerSFHi->SetMarkerStyle(22); hJerSFHi->SetMarkerColor(kBlue); hJerSFHi->SetLineColor(kBlue); hJerSFHi->SetLineWidth(2); hJerSFHi->SetMarkerSize(1.4);

    for (int i = 0; i < fJerNEtaBins; ++i) {
        hJerSFDef->SetBinContent(i + 1, fJerDef[i]);
        hJerSFLow->SetBinContent(i + 1, fJerLow[i]);
        hJerSFHi->SetBinContent(i + 1, fJerHi[i]);
    }
    
    auto c = std::make_unique<TCanvas>("c", "c", 1000, 1000);
    c->cd();
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetLeftMargin(0.15);

    hJerSFDef->Draw();
    hJerSFLow->Draw("same");
    hJerSFHi->Draw("same");
    hJerSFDef->GetYaxis()->SetRangeUser(0.8, 2.0);

    auto leg = std::make_unique<TLegend>(0.45, 0.7, 0.6, 0.8);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hJerSFDef.get(), "Default", "p");
    leg->AddEntry(hJerSFLow.get(), "Down", "p");
    leg->AddEntry(hJerSFHi.get(), "Up", "p");
    leg->Draw("same");
    c->SaveAs("pPb8160_jerScalingFactors.pdf");

    auto oFile = std::make_unique<TFile>("pPb8160_jerScalingFactors.root", "recreate", "", 208);
    hJerSFDef->Write();
    hJerSFLow->Write();
    hJerSFHi->Write();
    oFile->Write();
    oFile->Close();
}
