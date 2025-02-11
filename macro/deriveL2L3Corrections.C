// ROOT headers
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TRatioPlot.h"
#include "TPad.h"

// C++ headers
#include <iostream>
#include <vector>

//________________
void plotReco2GenPtRatios(TFile *f, int collSystem = 0, double energy = 5.02) {
    TH3D *h3D = (TH3D*)f->Get("hRecoInclusiveJetPtRawOverPtRefVsPtEta");
    if ( !h3D ) {
        std::cerr << "Histogram hRecoInclusiveJetPtRawOverPtRefVsPtEta not found" << std::endl;
        return;
    }

    // Standard binning

    // pt  : 150 bins from 5. to 1505. GeV
    // eta : 52 bins from -5.2 to 5.2

    // Eta binning start: etaStart + (jetEtaBinsLow[j] - 1) * etaStep
    // Eta binning end: etaStart + jetEtaBinsHi[j] * etaStep
    // -5.2, -3.6, -3.0, -2.8, -2.6, -2.4, -2.0, -1.6, 0., 1.6, 2.0, 2.4, 2.6, 2.8, 3.0, 3.6
    std::vector<int> jetEtaBinsLow {1,  9, 12, 13, 14, 15, 16, 17, 19, 27, 35, 37, 39, 40, 41, 42, 45};
    // -3.6, -3.0, -2.8, -2.6, -2.4, -2.0, -1.6, 0., 1.6, 2.0, 2.4, 2.6, 2.8, 3.0, 3.6, 5.2
    std::vector<int> jetEtaBinsHi  {8, 11, 12, 13, 14, 15, 16, 18, 26, 34, 36, 38, 39, 40, 41, 44, 52};

    std::vector<double> jetEtaBins { -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.879, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.314, 3.489, 3.664, 3.839, 4.013, 4.191, 4.363, 4.538, 4.716, 4.889, 5.191 };

}

//________________
void deriveL2Corrections() {

    // Macro derives L2L3 corrections from the ratio of reco raw pt / ref pt.
    // Raw pt is the reconstructed jet pt, ref pt is the generated jet pt.
    // The ratio is the L2L3 correction factor.
    // The correction factors are plotted for each eta bin as a function of pt,
    // then total integrals is splitted into several pT intervals with 5-10% width of the total area.
    // For each eta bin and pT interval, the mean value of the correction factor is calculated.
    // The mean values is plotted as a function of pT, and fitted with a function.
    // The correction values are stored as a 2D array C++ array, and written to a file.

    // Base style
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(kBird);

    // Username of the machine
    TString uname = gSystem->GetFromPipe("whoami");

    // Collision system: 0 = pp, 1 = pPb, 2 = PbPb
    int collSystem = 0;
    double energy = 5.02;

    // pPb8160 embedding
    // TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", username.Data()) );
    TFile *inputFile = TFile::Open( Form("/Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_Pbgoing_eta2.root", username.Data()) );
    if ( !inputFile || inputFile->IsZombie() ) {
        std::cerr << Form("File not found: /Users/%s/cernbox/ana/pPb8160/embedding/Pbgoing/oEmbedding_pPb8160_Pbgoing_ak4.root", username.Data()) << std::endl;
        return;
    }
    collSystem = 1;
    energy = 8.16;
}