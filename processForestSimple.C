#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"

#include "JetCorrector.h"

#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>

const char* RED    = "\033[1;31m";
const char* GREEN  = "\033[1;32m";
const char* RESET  = "\033[0m";

// Eta shifts for pPb 8.16 TeV collisions, used for etaCM calculation
const int nEtaShifts = 13;
static constexpr std::array<float, nEtaShifts> etaShift{0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495 };

// Range for the selection
float ptHatRange[2] { 0.f, 8160.f};

// Event info variables
float ptHat;
float vz;

// Skim/event filter variables
int pBeamScrapingFilter;
int pPAprimaryVertexFilter;
int HBHENoiseFilterResultRun2Loose;
int phfCoincFilter;
int pVertexFilterCutdz1p0;
int pVertexFilterCutGplus;
int pVertexFilterCutVtx1;

// Maximum number of jets to store info for (adjust as needed)
const int maxJets = 1500;

// Gen jet info
int   nGenJets;           // number of generated jets
float genJetPt[maxJets];  // generated jet pT
float genJetEta[maxJets]; // generated jet eta
float genJetPhi[maxJets]; // generated jet phi

// Reco and ref jet info
int   nRecoJets;                  // number of reconstructed jets

float recoJetPtRaw[maxJets];      // reconstructed jet pT without JEC
float recoJetTrackMaxPt[maxJets]; // reconstructed track pT max in the jet
float recoJetEta[maxJets];        // reconstructed jet eta
float recoJetPhi[maxJets];        // reconstructed jet phi
float recoJetPfNHF[maxJets];      // reconstructed jet neutral hadron fraction
float recoJetPfNEF[maxJets];      // reconstructed jet neutral electromagnetic fraction
float recoJetPfCHF[maxJets];      // reconstructed jet charged hadron fraction
float recoJetPfMUF[maxJets];      // reconstructed jet muon fraction
float recoJetPfCEF[maxJets];      // reconstructed jet charged electromagnetic fraction
float recoJetPfCHM[maxJets];      // reconstructed jet charged muon fraction
float recoJetPfCEM[maxJets];      // reconstructed jet charged electromagnetic fraction
float recoJetPfNHM[maxJets];      // reconstructed jet neutral hadron fraction
float recoJetPfNEM[maxJets];      // reconstructed jet neutral electromagnetic fraction
float recoJetPfMUM[maxJets];      // reconstructed jet muon fraction

float refJetPt[maxJets];   // reference jet pT
float refJetEta[maxJets];  // reference jet eta
float refJetPhi[maxJets];  // reference jet phi

//________________
struct GenJet {
    float pt;
    float eta;
    float phi;
};

//________________
struct RecoJet {
    float recoPtRaw;
    float recoPt;
    float recoPtExtra;
    float recoEta;
    float recoPhi;
    float recoPfNHF;
    float recoPfNEF;
    float recoPfCHF;
    float recoPfMUF;
    float recoPfCEF;
    float recoPfCHM;
    float recoPfCEM;
    float recoPfNHM;
    float recoPfNEM;
    float recoPfMUM;
    float refPt;
    float refEta;
};

//________________
struct Histograms {
    //
    // Event level histograms
    //
    std::unique_ptr<TH1D> hPtHatUnweighted;
    std::unique_ptr<TH1D> hPtHat;
    std::unique_ptr<TH1D> hVzUnweighted;
    std::unique_ptr<TH1D> hVz;

    //
    // Gen level histograms
    //

    // Gen jet histograms
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaLabUnflipped;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaLab;
    std::unique_ptr<TH1D> hGenInclusiveJetPt;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hGenInclusiveJetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaStdBins;

    // Gen dijet histograms
    std::unique_ptr<TH1D> hGenDijetPtAve;
    std::unique_ptr<TH1D> hGenDijetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMForwardShifted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMBackwardShifted[nEtaShifts];

    //
    // Reco level histograms
    //

    // Reco jet histograms
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaLabUnflipped;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaLab;
    std::unique_ptr<TH1D> hRecoInclusiveJetPt;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaStdBins;
    std::unique_ptr<TH3D> hRecoInclusiveJetJESPtEtaStdBins;
    std::unique_ptr<TH3D> hRecoInclusiveJetJESExtraPtEtaStdBins;
    std::unique_ptr<TH3D> hRecoInclusiveJetJECPtEtaStdBins;

    //
    // Ref level histograms
    //


};

//________________
// Event weight calculation for pPb8160 MC samples
double eventWeight(const float &ptHat, const float &vz, TF1 &fVzWeight, Long64_t eventsInSample) {

    double genWeight = 1.0;
    double vzWeight = 1./ fVzWeight.Eval( vz );
    // vzWeight = 1.;
    // Magic numbers are (cross section x Nevents generated)
    if (ptHat > 15.0 && ptHat <= 30.)       { genWeight = 1.0404701e-06 * 961104 ; }
    else if (ptHat > 30. && ptHat <= 50.)   { genWeight = 7.7966624e-08 * 952110 ; }
    else if (ptHat > 50. && ptHat <= 80.)   { genWeight = 1.0016052e-08 * 952554 ; }
    else if (ptHat > 80. && ptHat <= 120.)  { genWeight = 1.3018269e-09 * 996844 ; }
    else if (ptHat > 120.&& ptHat <= 170.)  { genWeight = 2.2648493e-10 * 964681 ; }
    else if (ptHat > 170. && ptHat <= 220.) { genWeight = 4.0879112e-11 * 999260 ; }
    else if (ptHat > 220. && ptHat <= 280.) { genWeight = 1.1898939e-11 * 964336 ; }
    else if (ptHat > 280. && ptHat <= 370.) { genWeight = 3.3364433e-12 * 995036 ; }
    else if (ptHat > 370. && ptHat <= 460.) { genWeight = 7.6612402e-13 * 958160 ; }
    else if (ptHat > 460. && ptHat <= 540.) { genWeight = 2.1341026e-13 * 981427 ; }
    else if (ptHat > 540.)                  { genWeight = 7.9191586e-14 * 1000000; }
    genWeight /= eventsInSample;

    return genWeight * vzWeight;
}

//________________
void createHistograms(Histograms &hs, const bool &isMc = false) {


    std::cout << "Creating histograms...";
    // Standard eta binning for L2L3 corrections
    double jetEtaL2L3StdVals[] = { -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, 
                                 -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, 
                                 -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
                                 -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, 
                                 -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, 
                                 -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
                                 -0.435, -0.348, -0.261, -0.174, -0.087,  0.000, 
                                  0.087,  0.174,  0.261,  0.348,  0.435,  0.522, 
                                  0.609,  0.696,  0.783,  0.879,  0.957,  1.044, 
                                  1.131,  1.218,  1.305,  1.392,  1.479,  1.566, 
                                  1.653,  1.740,  1.830,  1.930,  2.043,  2.172, 
                                  2.322,  2.500,  2.650,  2.853,  2.964,  3.139, 
                                  3.314,  3.489,  3.664,  3.839,  4.013,  4.191, 
                                  4.363,  4.538,  4.716,  4.889,  5.191 };
    int jetEtaL2L3StdBins = sizeof(jetEtaL2L3StdVals)/sizeof(double) - 1;

    // Single jet pt binning
    const int nJetPtBins = 100;
    double jetPtBins[] = { 0., 1000. };
    const int nJetEtaBins = 100;
    double jetEtaBins[] = { -3.6, 3.6 };
    const int nJetJESBins = 100;
    double jetJESBins[] = { 0., 2. };

    // Dijet binning
    const int nDijetPtBins = 50;
    double dijetPtBins[] = { 30., 530.};
    const int nDijetEtaBins = 60;
    double dijetEtaBins[] = { -3., 3. };
    const int nDijetEtaFBBins = 30;
    double dijetEtaFBBins[] = { 0., 3. };

    //
    // Event level histograms
    //
    hs.hPtHatUnweighted = std::make_unique<TH1D>("hPtHatUnweighted", 
                                              "#hat{p}_{T} unweighted;#hat{p}_{T} (GeV);dN/d#hat{p}_{T}", 
                                              100, 0., 1000.);
    hs.hPtHatUnweighted->Sumw2();
    hs.hPtHat = std::make_unique<TH1D>("hPtHat", 
                                              "#hat{p}_{T};#hat{p}_{T} (GeV);dN/d#hat{p}_{T}", 
                                              100, 0., 1000.);
    hs.hPtHat->Sumw2();

    hs.hVzUnweighted = std::make_unique<TH1D>("hVzUnweighted", 
                                           "vz unweighted;vz (cm);dN/dvz", 
                                           60, -15., 15.);
    hs.hVzUnweighted->Sumw2();
    hs.hVz = std::make_unique<TH1D>("hVz", 
                                 "vz;vz (cm);dN/dvz", 
                                 60, -15., 15.);
    hs.hVz->Sumw2();

    if (isMc) {

        //
        // Gen level histograms
        //

        // Gen jet histograms
        hs.hGenInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLabUnflipped", 
                                                "Gen jet pT vs eta (lab frame, unflipped);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hGenInclusiveJetEtaLabUnflipped = std::make_unique<TH1D>("hGenInclusiveJetEtaLabUnflipped", 
                                                "Gen jet eta (lab frame, unflipped);#eta;dN/d#eta", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaLabUnflipped->Sumw2();
        hs.hGenInclusiveJetEtaLab = std::make_unique<TH1D>("hGenInclusiveJetEtaLab", 
                                                "Gen jet eta (lab frame);#eta;dN/d#eta", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaLab->Sumw2();
        hs.hGenInclusiveJetPt = std::make_unique<TH1D>("hGenInclusiveJetPt", 
                                                "Gen jet p_{T};p_{T} (GeV);dN/dp_{T}", 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPt->Sumw2();
        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
            hs.hGenInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%d", iShift), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Sumw2();
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift) 


        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%.3f", etaShift[0]), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0]->Sumw2();
        hs.hGenInclusiveJetEtaCMShifted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%.3f", etaShift[0]), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShifted[0]->Sumw2();
        hs.hGenInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hGenInclusiveJetPtEtaStdBins", 
                                                "Gen jet pT vs eta (standard bins);#eta;p_{T} (GeV)",
                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals,
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtEtaStdBins->Sumw2();

        // Gen dijet histograms
        hs.hGenDijetPtAve = std::make_unique<TH1D>("hGenDijetPtAve", 
                                                "Gen dijet pT average;p_{T}^{ave} (GeV);dN/dp_{T}^{ave}", 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1]);
        hs.hGenDijetPtAve->Sumw2();
        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenDijetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMShiftedUnweighted_%d", iShift), 
                                                Form("Gen dijet eta (CM frame, shifted by %.3f, unweighted);#eta;dN/d#eta", etaShift[iShift]), 
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Sumw2();
            hs.hGenDijetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMShifted_%d", iShift), 
                                                Form("Gen dijet eta (CM frame, shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetEtaCMShifted[iShift]->Sumw2();
            hs.hGenDijetEtaCMForwardShifted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMForwardShifted_%d", iShift), 
                                                Form("Gen dijet eta (CM frame, forward shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetEtaCMForwardShifted[iShift]->Sumw2();
            hs.hGenDijetEtaCMBackwardShifted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMBackwardShifted_%d", iShift), 
                                                Form("Gen dijet eta (CM frame, backward shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetEtaCMBackwardShifted[iShift]->Sumw2();
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

    }

    //
    // Reco level histograms
    //

    hs.hRecoInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnflipped", 
                                              "Reco jet pT vs eta (lab frame, unflipped);#eta;p_{T} (GeV)", 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                              nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetEtaLabUnflipped = std::make_unique<TH1D>("hRecoInclusiveJetEtaLabUnflipped", 
                                              "Reco jet eta (lab frame, unflipped);#eta;dN/d#eta", 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetEtaLab = std::make_unique<TH1D>("hRecoInclusiveJetEtaLab", 
                                              "Reco jet eta (lab frame);#eta;dN/d#eta", 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaLab->Sumw2();
    hs.hRecoInclusiveJetPt = std::make_unique<TH1D>("hRecoInclusiveJetPt", 
                                              "Reco jet p_{T};p_{T} (GeV);dN/dp_{T}", 
                                              nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPt->Sumw2();
    for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
        hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
                                              Form("Reco jet eta (CM frame, shifted by %.3f, unweighted);#eta;dN/d#eta", etaShift[iShift]), 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
        hs.hRecoInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShifted_%d", iShift), 
                                              Form("Reco jet eta (CM frame, shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetEtaCMShifted[iShift]->Sumw2();
    } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

    hs.hRecoInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaStdBins", 
                                              "Reco jet pT vs eta (standard bins);#eta;p_{T} (GeV)",
                                              jetEtaL2L3StdBins, jetEtaL2L3StdVals,
                                              nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtEtaStdBins->Sumw2();

    if (isMc) {
        hs.hRecoInclusiveJetJESPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESPtEtaStdBins", 
                                                "Reco jet JES (reco/ref) vs pT vs eta (standard bins);JES;p_{T} (GeV);#eta",
                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESPtEtaStdBins->Sumw2();
        hs.hRecoInclusiveJetJESExtraPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESExtraPtEtaStdBins", 
                                                "Reco jet JES (reco/ref) vs pT vs eta (standard bins, extra);JES;p_{T} (GeV);#eta",
                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESExtraPtEtaStdBins->Sumw2();
        hs.hRecoInclusiveJetJECPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJECPtEtaStdBins",
                                                "Reco jet JEC (recoCorr/recoRaw) vs pT vs eta (standard bins, corrected);JEC;p_{T} (GeV);#eta",
                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJECPtEtaStdBins->Sumw2();
    }

    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
void setupChains(const TString &input, TChain &hltTree, TChain &eventTree, TChain &skimTree, 
                 TChain &jetTree, const bool &isMc) {

    // Check input exists
    if ( input.Length() <= 0 ) {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        exit(0);
    }
    // Regular input
    else {
        // If input is a single ROOT file
        if ( input.Index(".root") > 0 ) {
            std::cout << Form( "Adding %s file to chains\n", input.Data() );
            hltTree.Add( input.Data() );
            eventTree.Add( input.Data() );
            skimTree.Add( input.Data() );
            jetTree.Add( input.Data() );
        }
        // Assuming that list of files is provided instead of a single file
        else {
            std::ifstream inputStream( input.Data() );

            if ( !inputStream ) std::cout << Form( "ERROR: Cannot open file list: %s\n", input.Data() );
            Int_t nFiles = 0;
            std::string file;
            size_t pos;
            while ( getline( inputStream, file ) ) {
                // NOTE: our external formatters may pass "file NumEvents"
                //       Take only the first part
                //cout << "DEBUG found " <<  file << endl;
                pos = file.find_first_of(" ");
                if ( pos != std::string::npos ) file.erase( pos, file.length() - pos );
                //cout << "DEBUG found [" <<  file << "]" << endl;

                // Check that file is of a correct name
                if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos &&
                     file.find("AOD") != std::string::npos ) {
                    
                    // Open file
                    TFile* ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if ( ftmp && !ftmp->IsZombie() && ftmp->GetNkeys() ) {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str() );
                        // Adding file to chains
                        hltTree.Add( file.c_str() );
                        eventTree.Add( file.c_str() );
                        skimTree.Add( file.c_str() );
                        jetTree.Add( file.c_str() );
                        ++nFiles;
                    } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp) {
                        ftmp->Close();
                    } //if (ftmp)
                } //if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos && file.find("AOD") != std::string::npos )
            } //while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
        } // else {   if file list
    } // else {   if normal input

    // Connect chains between each other
    hltTree.AddFriend(&eventTree);
    hltTree.AddFriend(&skimTree);
    hltTree.AddFriend(&jetTree);
}

//________________
void setupBranches(TChain &mainTree, const bool &isMc) {

    std::cout << "Setting up branches...";

    //
    // Disable all branches first
    //
    mainTree.SetBranchStatus("*", 0);

    //
    // Enable only the branches needed for this analysis
    //

    // Event level branches
    mainTree.SetBranchStatus("vz", 1);
    mainTree.SetBranchStatus("pthat", 1);

    // Enable skim/event filter branches
    mainTree.SetBranchStatus("pBeamScrapingFilter", 1);
    mainTree.SetBranchStatus("pPAprimaryVertexFilter", 1);
    mainTree.SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    mainTree.SetBranchStatus("phfCoincFilter", 1);
    mainTree.SetBranchStatus("pVertexFilterCutdz1p0", 1);
    mainTree.SetBranchStatus("pVertexFilterCutGplus", 1);
    mainTree.SetBranchStatus("pVertexFilterCutVtx1", 1);

    // Jet level branches
    if (isMc) {
        mainTree.SetBranchStatus("ngen", 1);
        mainTree.SetBranchStatus("genpt", 1);
        mainTree.SetBranchStatus("geneta", 1);
        mainTree.SetBranchStatus("genphi", 1);
    }

    mainTree.SetBranchStatus("nref", 1);
    mainTree.SetBranchStatus("rawpt", 1);
    mainTree.SetBranchStatus("trackMax", 1);
    mainTree.SetBranchStatus("jteta", 1);
    mainTree.SetBranchStatus("jtPfNHF", 1);
    mainTree.SetBranchStatus("jtPfNEF", 1);
    mainTree.SetBranchStatus("jtPfCHF", 1);
    mainTree.SetBranchStatus("jtPfMUF", 1);
    mainTree.SetBranchStatus("jtPfCEF", 1);
    mainTree.SetBranchStatus("jtPfCHM", 1);
    mainTree.SetBranchStatus("jtPfCEM", 1);
    mainTree.SetBranchStatus("jtPfNHM", 1);
    mainTree.SetBranchStatus("jtPfNEM", 1);
    mainTree.SetBranchStatus("jtPfMUM", 1);

    mainTree.SetBranchStatus("refpt", 1);
    mainTree.SetBranchStatus("refeta", 1);

    // Set branch addresses
    mainTree.SetBranchAddress("vz", &vz);
    mainTree.SetBranchAddress("pthat", &ptHat);

    mainTree.SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);
    mainTree.SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter);
    mainTree.SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose);
    mainTree.SetBranchAddress("phfCoincFilter", &phfCoincFilter);
    mainTree.SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0);
    mainTree.SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus);
    mainTree.SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1);

    if (isMc) {
        mainTree.SetBranchAddress("ngen",   &nGenJets);
        mainTree.SetBranchAddress("genpt",  &genJetPt);
        mainTree.SetBranchAddress("geneta", &genJetEta);
        mainTree.SetBranchAddress("genphi", &genJetPhi);
    }

    mainTree.SetBranchAddress("nref",     &nRecoJets);
    mainTree.SetBranchAddress("rawpt",    &recoJetPtRaw);
    mainTree.SetBranchAddress("trackMax", &recoJetTrackMaxPt);
    mainTree.SetBranchAddress("jteta",    &recoJetEta);
    mainTree.SetBranchAddress("jtphi",    &recoJetPhi);
    mainTree.SetBranchAddress("jtPfNHF",  &recoJetPfNHF);
    mainTree.SetBranchAddress("jtPfNEF",  &recoJetPfNEF);
    mainTree.SetBranchAddress("jtPfCHF",  &recoJetPfCHF);
    mainTree.SetBranchAddress("jtPfMUF",  &recoJetPfMUF);
    mainTree.SetBranchAddress("jtPfCEF",  &recoJetPfCEF);
    mainTree.SetBranchAddress("jtPfCHM",  &recoJetPfCHM);
    mainTree.SetBranchAddress("jtPfCEM",  &recoJetPfCEM);
    mainTree.SetBranchAddress("jtPfNHM",  &recoJetPfNHM);
    mainTree.SetBranchAddress("jtPfNEM",  &recoJetPfNEM);
    mainTree.SetBranchAddress("jtPfMUM",  &recoJetPfMUM);

    mainTree.SetBranchAddress("refpt", &refJetPt);
    mainTree.SetBranchAddress("refeta", &refJetEta);
    mainTree.SetBranchAddress("refphi", &refJetPhi);

    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
void setupInput(const TString &input, TChain &hltTree, TChain &eventTree, TChain &skimTree, 
                TChain &jetTree, const bool &isMc) {
    setupChains(input, hltTree, eventTree, skimTree, jetTree, isMc);
    setupBranches(hltTree, isMc);
}

//________________
float etaCM(const float &etaLab, const float &etaShift, const bool &isPbGoing, const bool &isMc) {
    float etaCM;
    if (isPbGoing) {
        etaCM = -1.0 * (etaLab + etaShift);
    }
    else {
        etaCM = etaLab - etaShift;
    }
    return etaCM;
}

//________________
float etaLab(const float &eta, const bool &isPbGoing, const bool &isMc) {
    float etaL = eta;
    if ( isMc ) { // For embedding: Pb goes to negative, p goes to positive
        if (isPbGoing) {
            etaL = -etaL;
        }
    }
    else { // For data: p goes to negative, Pb goes to positive
        if (isPbGoing) {
        }
        else {
            etaL = -etaL;
        }
    }
    return etaL;
}

//________________
void setPtHatRange(const int &ptHatSample) {
    std::cout << "Setting pT hat range for sample " << ptHatSample;
    if (ptHatSample == 15)       { ptHatRange[0] = 15.; ptHatRange[1] = 30.; }
    else if (ptHatSample == 30)  { ptHatRange[0] = 30.; ptHatRange[1] = 50.; }
    else if (ptHatSample == 50)  { ptHatRange[0] = 50.; ptHatRange[1] = 80.; }
    else if (ptHatSample == 80)  { ptHatRange[0] = 80.; ptHatRange[1] = 120.; }
    else if (ptHatSample == 120) { ptHatRange[0] = 120.; ptHatRange[1] = 170.; }
    else if (ptHatSample == 170) { ptHatRange[0] = 170.; ptHatRange[1] = 220.; }
    else if (ptHatSample == 220) { ptHatRange[0] = 220.; ptHatRange[1] = 280.; }
    else if (ptHatSample == 280) { ptHatRange[0] = 280.; ptHatRange[1] = 370.; }
    else if (ptHatSample == 370) { ptHatRange[0] = 370.; ptHatRange[1] = 460.; }
    else if (ptHatSample == 460) { ptHatRange[0] = 460.; ptHatRange[1] = 540.; }
    else if (ptHatSample == 540) { ptHatRange[0] = 540.; ptHatRange[1] = std::numeric_limits<float>::max(); }
    std::cout << Form(": [%.1f, %.1f]", ptHatRange[0], ptHatRange[1]) << std::endl;
}

//________________
bool isGoodEvent(const bool &isPbGoing, const bool &isMc, const int &ptHatSample) {

    if (isMc) {
        if (ptHat < ptHatRange[0] || ptHat > ptHatRange[1]) {
            return false;
        }
    }

    return ( ( std::abs(vz) <= 15. ) && 
             ( pBeamScrapingFilter == 1 ) && 
             ( pPAprimaryVertexFilter == 1 ) && 
             ( HBHENoiseFilterResultRun2Loose == 1 ) && 
             ( phfCoincFilter == 1 ) && 
             ( pVertexFilterCutdz1p0 == 1 ) );
}

//________________
void processGenJets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                    std::vector<GenJet> &genJets, Histograms &hs) {

    for (int iGenJet{0}; iGenJet < nGenJets; ++iGenJet) {
        float genJetPtL = genJetPt[iGenJet];
        float genJetEtaL = etaLab(genJetEta[iGenJet], isPbGoing, isMc);
        float genJetPhiL = genJetPhi[iGenJet];

        hs.hGenInclusiveJetPtEtaLabUnflipped->Fill(genJetEtaL, genJetPtL);
        hs.hGenInclusiveJetEtaLabUnflipped->Fill(genJetEtaL);
        hs.hGenInclusiveJetEtaLab->Fill(etaLab(genJetEta[iGenJet], isPbGoing, isMc));
        hs.hGenInclusiveJetPt->Fill(genJetPtL);

        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Fill(etaCM(genJetEta[iGenJet], etaShift[iShift], isPbGoing, isMc));
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Fill(etaCM(genJetEta[iGenJet], etaShift[iShift], isPbGoing, isMc));
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        hs.hGenInclusiveJetPtEtaStdBins->Fill(etaLab(genJetEta[iGenJet], isPbGoing, isMc), genJetPtL);
    } // for (int iGenJet{0}; iGenJet < nGenJets; ++iGenJet)
}

//________________
void processGenDijets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                      std::vector<GenJet> &genJets, Histograms &hs) {

    // Must be at least 2 jets to form a dijet system
    if (genJets.size() < 2) return;

    // Sort jets by pT to identify leading and subleading jets
    std::sort(genJets.begin(), genJets.end(), [](const GenJet &a, const GenJet &b) { return a.pt > b.pt; });

    GenJet leadingJet = genJets[0];
    GenJet subleadingJet = genJets[1];

    // Make a dijet cut
    if (leadingJet.pt < 50. || subleadingJet.pt < 40.) return;
    float dphi = leadingJet.phi - subleadingJet.phi;
    if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
    if (std::abs(dphi) < TMath::TwoPi() / 3.) return;

    for (int iShift = 0; iShift < nEtaShifts; ++iShift) {

        float leadingEtaCMShifted = etaCM(leadingJet.eta, etaShift[iShift], isPbGoing, isMc);
        float subleadingEtaCMShifted = etaCM(subleadingJet.eta, etaShift[iShift], isPbGoing, isMc);
        float dijetEtaCMShifted = 0.5 * (leadingEtaCMShifted + subleadingEtaCMShifted);
        float dijetPtAve = 0.5 * (leadingJet.pt + subleadingJet.pt);

        if (std::abs(leadingEtaCMShifted) > 1.9 || std::abs(subleadingEtaCMShifted) > 1.9) continue;
        
        hs.hGenDijetPtAve->Fill(dijetPtAve, weight);

        // if (dijetPtAve < 60. || dijetPtAve > 80.) continue;
        hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Fill(dijetEtaCMShifted);
        hs.hGenDijetEtaCMShifted[iShift]->Fill(dijetEtaCMShifted, weight);
        (dijetEtaCMShifted > 0) ? hs.hGenDijetEtaCMForwardShifted[iShift]->Fill(dijetEtaCMShifted, weight) : 
                                hs.hGenDijetEtaCMBackwardShifted[iShift]->Fill( std::abs(dijetEtaCMShifted), weight);
    } // for (int iShift = 0; iShift < nEtaShifts; ++iShift)
}

//________________
void processRecoJets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                     std::vector<RecoJet> &recoJets, Histograms &hs, JetCorrector &jec) {
    
    float recoJetPt{0.};
    float recoJetPtExtraCorr{0.};
    float recoJetEtaLabFlipped{0.};
    float recoJetEtaCM{0.};
    float recoJetJECFactor{0.};
    float recoJetJES{0.};

    // Loop over reco jets
    for (int iRecoJet = 0; iRecoJet < nRecoJets; iRecoJet++) {

        recoJetPt = {0.};
        recoJetPtExtraCorr = {0.};
        recoJetEtaLabFlipped = {0.};
        recoJetEtaCM = {0.};
        recoJetJECFactor = {0.};
        recoJetJES = {0.};

        jec.SetJetPT( recoJetPtRaw[iRecoJet] );
        jec.SetJetEta( recoJetEta[iRecoJet] );
        jec.SetJetPhi( recoJetPhi[iRecoJet] );
        recoJetPt = jec.GetCorrectedPT();
        recoJetJECFactor = recoJetPt / recoJetPtRaw[iRecoJet];

        // Should be very careful with dropping low-pT jets
        if (recoJetPt < 15.) continue; // Skip low-pT jets

        

        RecoJet recoJet{recoJetPtRaw[iRecoJet], recoJetPt, recoJetPtExtraCorr,
                        recoJetEta[iRecoJet], recoJetPhi[iRecoJet],
                        recoJetPfNHF[iRecoJet], recoJetPfNEF[iRecoJet], recoJetPfCHF[iRecoJet], 
                        recoJetPfMUF[iRecoJet], recoJetPfCEF[iRecoJet], recoJetPfCHM[iRecoJet], 
                        recoJetPfCEM[iRecoJet], recoJetPfNHM[iRecoJet], recoJetPfNEM[iRecoJet], 
                        recoJetPfMUM[iRecoJet], refJetPt[iRecoJet], refJetEta[iRecoJet]};
        recoJets.push_back(recoJet);

        recoJetEtaLabFlipped = etaLab(recoJetEta[iRecoJet], isPbGoing, isMc);
        // recoJetEtaCM = etaCM(recoJetEta[iRecoJet], 0.465f, isPbGoing, isMc);

        // std::cout << Form("Reco jet %d: raw pT = %.1f GeV, corrected pT = %.1f GeV, ref pT: %.1f, eta = %.2f", 
        //                   iRecoJet, recoJetPtRaw[iRecoJet], recoJetPt, refJetPt[iRecoJet], recoJetEta[iRecoJet]) << std::endl;

        hs.hRecoInclusiveJetPtEtaLabUnflipped->Fill(recoJetEta[iRecoJet], recoJetPt, weight);
        hs.hRecoInclusiveJetEtaLabUnflipped->Fill(recoJetEta[iRecoJet], weight);
        hs.hRecoInclusiveJetEtaLab->Fill(recoJetEtaLabFlipped, weight);
        hs.hRecoInclusiveJetPt->Fill(recoJetPt, weight);
        // hs.hRecoInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
        // hs.hRecoInclusiveJetEtaCMShifted[nEtaShifts];

        hs.hRecoInclusiveJetPtEtaStdBins->Fill(recoJetEta[iRecoJet], recoJetPt, weight);
        hs.hRecoInclusiveJetJECPtEtaStdBins->Fill(recoJetJECFactor, recoJetPt, recoJetEta[iRecoJet], weight);

        if (isMc) {
            if ( refJetPt[iRecoJet] < 0.) continue;
            recoJetJECFactor = recoJetPt / recoJetPtRaw[iRecoJet];
            recoJetJES = recoJetPt / refJetPt[iRecoJet];
            hs.hRecoInclusiveJetJESPtEtaStdBins->Fill(recoJetJES, recoJetPt, recoJetEta[iRecoJet], weight);
            // hs.hRecoInclusiveJetJESExtraPtEtaStdBins;
        } // if (isMc)

    } // for (int iRecoJet = 0; iRecoJet < nRecoJets; iRecoJet++)
}

//________________
void processRecoDijets(const bool &isPbGoing, const bool &isMc, const float &weight, 
                       std::vector<RecoJet> &recoJets, Histograms &hs) {
    // Empty for now
}

//________________
void writeOutput(TString &oFileName, Histograms &hs, const bool &isMc) {

    std::cout << "Writing output to file: " << oFileName.Data();
    // Output file
    int compressionSetting = 208; // LZMA compression
    // TFile *fOut = TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting);
    auto fOut = std::unique_ptr<TFile>( TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting) );

    //
    // Event level histograms
    //
    hs.hPtHatUnweighted->Write();
    hs.hPtHat->Write();
    hs.hVzUnweighted->Write();
    hs.hVz->Write();

    //
    // Gen level histograms
    //
    if (isMc) {
        // Gen jets
        hs.hGenInclusiveJetPtEtaLabUnflipped->Write();
        hs.hGenInclusiveJetEtaLabUnflipped->Write();
        hs.hGenInclusiveJetEtaLab->Write();
        hs.hGenInclusiveJetPt->Write();
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Write();
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Write();
        }
        hs.hGenInclusiveJetPtEtaStdBins->Write();

        // Gen dijets
        hs.hGenDijetPtAve->Write();
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Write();
            hs.hGenDijetEtaCMShifted[iShift]->Write();
            hs.hGenDijetEtaCMForwardShifted[iShift]->Write();
            hs.hGenDijetEtaCMBackwardShifted[iShift]->Write();
        }
    }

    //
    // Reco level histograms
    //
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetEtaLab->Write();
    hs.hRecoInclusiveJetPt->Write();
    for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
        hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Write();
        hs.hRecoInclusiveJetEtaCMShifted[iShift]->Write();
    }
    hs.hRecoInclusiveJetPtEtaStdBins->Write();
    hs.hRecoInclusiveJetJECPtEtaStdBins->Write();
    if (isMc) {
        hs.hRecoInclusiveJetJESPtEtaStdBins->Write();
        hs.hRecoInclusiveJetJESExtraPtEtaStdBins->Write();
    }

    fOut->Close();
    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
/**
    List of input parameters for the analysis
    
    input: input file name (or file list)
    output: output file name
    mcType (default 2): 0 - data, 1 - embedding, 2 - pythia
    isPbGoingDir (default 1): 0 - p-going, 1 - Pb-going
    ptHatSample (default 30, for MC only): 15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540
    jeuSyst (default 0): 0 - default, 1 - JEU+, -1 - JEU-
    jerSyst (default other -99): 0 - std.JER, 1 - JER+, -1 - JER-, other - only JEC applied (no smearing)
    triggerId (default 0): 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100
    recoJetSelMethod (default 2): 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetIdTightLeptVeto, 3 - jetIdLoose
*/
#if defined(__CINT__) || defined(__CLING__)
void processForestSimple( const char* input = "", const char* output = "", 
                          int mcType = 2, int isPbGoingDir = 1, int ptHatSample = 30, int jeuSyst = 0,
                          int jerSyst = -99, int triggerId = 0, int recoJetSelMethod = 2 ) {
#else
int main(int argc, char* argv[]) {

    const char* input = "";
    const char* output = "";
    int mcType = 2;
    int isPbGoingDir = 1;
    int ptHatSample = 30;
    int jeuSyst = 0;
    int jerSyst = -99;
    int triggerId = 0;
    int recoJetSelMethod = 2;

    if (argc <= 1) {
        std::cerr << RED << "No input parameters provided. Terminating." << RESET << std::endl;
        std::cout << GREEN << "Usage: " << argv[0] << " [options]\n"
                  << "Options:\n"
                  << "  <input_file_or_list>  Input ROOT file or file list (default: empty)\n"
                  << "  <output_file>         Output ROOT file name (default: empty)\n"
                  << "  <0|1|2>               MC type (default 2): 0 - data, 1 - embedding, 2 - pythia\n"
                  << "  <0|1>                 Direction (default 1): 0 - p-going, 1 - Pb-going\n"
                  << "  <value>               pT hat sample for MC (default: 30)\n"
                  << "  <-1|0|1>              JEU systematic variation (default 0): -1 - JEU-, 0 - not applied, 1 - JEU+\n"
                  << "  <-99|0|1|-1>          JER systematic variation (default -99): -99 - only JEC applied, -1 - JER-, 0 - std.JER, 1 - JER+ \n"
                  << "  <value>               Trigger ID (default 0): 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100\n"
                  << "  <value>               Reco jet selection method (default 2): 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetIdTightLeptVeto, 3 - jetIdLoose\n"
                  << RESET;
        return 0;
    }

    input = argv[1];
    output = argv[2];
    mcType = std::atoi(argv[3]);
    isPbGoingDir = std::atoi(argv[4]);
    ptHatSample = std::atoi(argv[5]);
    jeuSyst = std::atoi(argv[6]);
    jerSyst = std::atoi(argv[7]);
    triggerId = std::atoi(argv[8]);
    recoJetSelMethod = std::atoi(argv[9]);
#endif

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    bool isMc = (mcType != 0);
    bool isPythia = ((mcType != 0) && (mcType == 2));
    bool isEmbedding = (isMc == 1);
    bool isPbGoing = (isPbGoingDir == 1);

    // Define the standard eta bins for L2L3 corrections
    TString generator = isPythia ? "pythia" : "embedding";
    TString tag = isPythia ? "unembedded" : "embedded";
    TString directionStr = isPbGoing ? "Pbgoing" : "pgoing";

    std::vector< std::string > fJECFileNames;

    if (isMc) {
        if (isPbGoingDir) {
            if (isEmbedding) {
                // PYTHIA+EPOS
                fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt");
            }
            else {
                // PYTHIA
                fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_Unembedded_MC_L2Relative_AK4PF.txt");
            }
        }
        else {
            if (isEmbedding) {
                // PYTHIA+EPOS
                fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt");
            }
            else {
                // PYTHIA
                fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF.txt");
            }
        }
    }
    else {
        if (isPbGoingDir) { // Remember to flip to p-going for data
            fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt");
        }
        else {
            fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt");
        }
        // Old correction
        // JECFileDataName = "Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt";
        // New correction (with updated residuals)
        fJECFileNames.push_back("../aux_files/pPb_8160/JEC/Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PF.txt");
        
        // JEUFileName = "../aux_files/pPb_8160/JEC/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF.txt";
    } // else


    // Set pT hat range for MC samples
    if (isMc) setPtHatRange(ptHatSample);

    // Setup input filename
    TString fname;
    if (input == nullptr || input[0] == '\0') {
        fname = Form("~/cernbox/ana/pPb8160/%s/%s/forest/HiForestSkim_pPb_MC_pthat%d_%s_%s.root", 
                     generator.Data(), directionStr.Data(), ptHatSample, directionStr.Data(), tag.Data());
    }
    else {
        fname = input;
    }
    std::cout << "Input file: " << GREEN << fname.Data() << RESET << std::endl;
 
    // Setup output filename
    TString oFileName;
    if (output == nullptr || output[0] == '\0') {
        oFileName = Form("macro/eta_shift/%s_%s_ptHat%d.root", generator.Data(), directionStr.Data(), ptHatSample);
    }
    else {
        oFileName = output;
    }
    std::cout << "Output file name: " << GREEN << oFileName.Data() << RESET << std::endl;


    // TFile *f = TFile::Open( fname );
    auto mainTree = std::make_unique<TChain>("hltanalysis/HltTree");
    auto eventTree = std::make_unique<TChain>("hiEvtAnalyzer/HiTree");
    auto skimTree = std::make_unique<TChain>("skimanalysis/HltTree");
    auto jetTree = std::make_unique<TChain>("ak4PFJetAnalyzer/t");

    // TChain *mainTree = new TChain("hltanalysis/HltTree");
    // TChain *eventTree = new TChain("hiEvtAnalyzer/HiTree");
    // TChain *skimTree = new TChain("skimanalysis/HltTree");
    // TChain *jetTree = new TChain("ak4PFJetAnalyzer/t");
    setupInput(fname, *mainTree, *eventTree, *skimTree, *jetTree, isMc);
    // setupChain(fname.Data(), mainTree, eventTree, skimTree, jetTree);

    auto fVzWeight = std::make_unique<TF1>("fVzWeight", "pol8", -15.1, 15.1);
    fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);

    ///////////////////////
    // Create histograms //
    ///////////////////////
    Histograms hs;
    createHistograms(hs, isMc);

    //////////////////////////////////
    // Setup Jet Energy Corrections //
    //////////////////////////////////
    auto jec = std::make_unique<JetCorrector>(fJECFileNames);

    double weight{0.};

    //////////////////////////////
    //        Processing        //
    //////////////////////////////

    Long64_t nEntries = mainTree->GetEntries();
    std::cout << Form("%s%s%s %s%s%s direction, ptHat%d sample, number of entries:%s %lld %s", 
                       GREEN, (mcType ? ( (isPythia) ? "Pythia" : "Embedding") : "Data" ), RESET,
                       GREEN, (isPbGoing ? "Pb-going" : "p-going"), RESET,
                       ptHatSample, 
                       GREEN, nEntries, RESET) << std::endl;

    int nEventsProcessed{0};
    int nGoodEvents{0};

    std::vector<GenJet> genJets;
    std::vector<RecoJet> recoJets;

    int printEvery{50000};
    Long64_t nextPrintAt{printEvery};
    auto startTime = std::chrono::steady_clock::now();

    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {

        mainTree->GetEntry(iEntry);
        nEventsProcessed++;
        Long64_t entriesPassed = iEntry + 1;
        if ((entriesPassed >= nextPrintAt) || (entriesPassed == nEntries)) {
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - startTime;
            double elapsedSeconds = elapsed.count();
            double etaSeconds = 0.0;
            if (entriesPassed > 0 && entriesPassed < nEntries) {
                etaSeconds = elapsedSeconds * static_cast<double>(nEntries - entriesPassed) / static_cast<double>(entriesPassed);
            }

            Long64_t elapsedTotalSeconds = static_cast<Long64_t>(std::llround(elapsedSeconds));
            Long64_t elapsedMinutes = elapsedTotalSeconds / 60;
            Long64_t elapsedRemainingSeconds = elapsedTotalSeconds % 60;

            Long64_t etaTotalSeconds = static_cast<Long64_t>(std::llround(etaSeconds));
            Long64_t etaMinutes = etaTotalSeconds / 60;
            Long64_t etaRemainingSeconds = etaTotalSeconds % 60;

            std::cout << Form("Processed %lld / %lld entries | elapsed %lld m %lld s | ETA %lld m %lld s",
                              entriesPassed, nEntries,
                              elapsedMinutes, elapsedRemainingSeconds,
                              etaMinutes, etaRemainingSeconds) << std::endl;

            nextPrintAt += printEvery;
        }

        // std::cout << "\n========================================\n";
        // std::cout << Form("Event: %d ptHat = %.1f GeV, vz = %.2f cm", iEntry, ptHat, vz) << std::endl;
        // std::cout << Form("%.1f <= ptHat <= %.1f ? %s", ptHatRange[0], ptHatRange[1], 
        //                  (ptHat >= ptHatRange[0] && ptHat <= ptHatRange[1]) ? "true" : "false") << std::endl;
        // std::cout << Form("|vz| <= 15 cm ? %s", (std::abs(vz) <= 15.) ? "true" : "false") << std::endl;

        // Check the event is satisfies basic event selection

        if (!isGoodEvent(isPbGoing, isMc, ptHatSample)) continue;

        nGoodEvents++;

        if (isPbGoing) vz = -vz; 
        weight = eventWeight(ptHat, vz, *fVzWeight, nEntries);

        hs.hPtHatUnweighted->Fill(ptHat);
        hs.hPtHat->Fill(ptHat, weight);
        hs.hVzUnweighted->Fill(vz);
        hs.hVz->Fill(vz, weight);

        //
        // Gen level processing
        //
        if (isMc) {
            if ( !genJets.empty() ) genJets.clear();
            processGenJets(isPbGoing, isMc, weight, genJets, hs);
            processGenDijets(isPbGoing, isMc, weight, genJets, hs);
            if ( !genJets.empty() ) genJets.clear();
        }

        //
        // Reco level processing
        //
        if ( !recoJets.empty() ) recoJets.clear();
        processRecoJets(isPbGoing, isMc, weight, recoJets, hs, *jec);
        processRecoDijets(isPbGoing, isMc, weight, recoJets, hs);
        if ( !recoJets.empty() ) recoJets.clear();
    } // for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)

    std::cout << Form("Processed %d events, of which %d passed selection.", nEventsProcessed, nGoodEvents) << std::endl;
    
    // Write output
    writeOutput(oFileName, hs, isMc);

    // if (eventTree) { delete eventTree; eventTree = nullptr;}
    // if (skimTree)  { delete skimTree;  skimTree = nullptr; }
    // if (jetTree)   { delete jetTree;   jetTree = nullptr;  }
    // if (mainTree)  { delete mainTree;  mainTree = nullptr; }
}
