#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TString.h"
#include "THnSparse.h"

#include "JetCorrector.h"
#include "JetUncertainty.h"
#include "AnalysisContract.h"

#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <cstring>


#include "Histograms.h"
namespace jetanalysis {
void createHistograms(Histograms &hs, const bool &isMc) {


    // Preserve Sumw2 before weighted filling, including sparse response cells.
    // TH2 conventions are usually x=pT, y=eta, but some inclusive response
    // diagnostics use other axes; downstream semantic projections check ranges.
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
    double jetPtBins[] = { 10., 1010. };
    const int nJetEtaBins = 160;
    double jetEtaBins[] = { -4.0, 4.0 };
    const int nJetJESBins = 50;
    double jetJESBins[] = { 0., 2. };

    // Dijet binning
    const int nDijetPtBins = 100;
    double dijetPtBins[] = { 0., 1000.};
    const int nDijetEtaBins = 72;
    double dijetEtaBins[] = { -3.6, 3.6 };

    const int nDimResponse = 4;
    int responseBins[nDimResponse] = { nDijetPtBins, nDijetEtaBins, nDijetPtBins, nDijetEtaBins };
    double responseBinsMin[nDimResponse] = { dijetPtBins[0], dijetEtaBins[0], dijetPtBins[0], dijetEtaBins[0] };
    double responseBinsMax[nDimResponse] = { dijetPtBins[1], dijetEtaBins[1], dijetPtBins[1], dijetEtaBins[1] };

    //
    // Event level histograms
    //
    hs.hVzRaw = std::make_unique<TH1D>("hVzRaw", 
                                 "vz raw;vz (cm);dN/dvz", 
                                 60, -15., 15.);
    hs.hVzRaw->Sumw2();
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
        // Event level histograms
        //
        hs.hGenDijetPtAveOverPtHatVsPtHat = std::make_unique<TH2D>("hGenDijetPtAveOverPtHatVsPtHat", 
                                            "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenDijetPtAveOverPtHatVsPtHat->Sumw2();
        hs.hGenDijetPtAveOverPtHatVsPtHatPass = std::make_unique<TH2D>("hGenDijetPtAveOverPtHatVsPtHatPass", 
                                            "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenDijetPtAveOverPtHatVsPtHatPass->Sumw2();
        hs.hGenLeadJetPtOverPtHatVsPtHat = std::make_unique<TH2D>("hGenLeadJetPtOverPtHatVsPtHat", 
                                            "Gen leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Gen p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenLeadJetPtOverPtHatVsPtHat->Sumw2();
        hs.hGenLeadJetPtOverPtHatVsPtHatPass = std::make_unique<TH2D>("hGenLeadJetPtOverPtHatVsPtHatPass", 
                                            "Gen leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Gen p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenLeadJetPtOverPtHatVsPtHatPass->Sumw2();
        hs.hRecoDijetPtAveOverPtHatVsPtHat = std::make_unique<TH2D>("hRecoDijetPtAveOverPtHatVsPtHat", 
                                            "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoDijetPtAveOverPtHatVsPtHat->Sumw2();
        hs.hRecoDijetPtAveOverPtHatVsPtHatPass = std::make_unique<TH2D>("hRecoDijetPtAveOverPtHatVsPtHatPass", 
                                            "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoDijetPtAveOverPtHatVsPtHatPass->Sumw2();
        hs.hRecoLeadJetPtOverPtHatVsPtHat = std::make_unique<TH2D>("hRecoLeadJetPtOverPtHatVsPtHat", 
                                            "Reco leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Reco p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoLeadJetPtOverPtHatVsPtHat->Sumw2();
        hs.hRecoLeadJetPtOverPtHatVsPtHatPass = std::make_unique<TH2D>("hRecoLeadJetPtOverPtHatVsPtHatPass", 
                                            "Reco leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Reco p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoLeadJetPtOverPtHatVsPtHatPass->Sumw2();

        //
        // Gen level histograms
        //

        // Gen jet histograms
        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLabUnflippedUnweighted", 
                                                "Gen jet #eta (lab frame, unflipped, unweighted) vs p_{T};p_{T} (GeV);#eta",
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted->Sumw2();
        hs.hGenInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLabUnflipped", 
                                                "Gen jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta",
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hGenInclusiveJetPtEtaLab = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLab", 
                                                "Gen jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta", 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaLab->Sumw2();
        hs.hGenInclusiveJetPtEtaCM = std::make_unique<TH2D>("hGenInclusiveJetPtEtaCM", 
                                                "Gen jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM}", 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaCM->Sumw2();
        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
                                                Form("Gen jet #eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
            hs.hGenInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%d", iShift), 
                                                Form("Gen jet #eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Sumw2();
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift) 


        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%.3f", etaShift[0]), 
                                                Form("Gen jet #eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0]->Sumw2();
        hs.hGenInclusiveJetEtaCMShifted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%.3f", etaShift[0]), 
                                                                    Form("Gen jet #eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShifted[0]->Sumw2();
        hs.hGenInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hGenInclusiveJetPtEtaStdBins", 
                                                                "Gen jet #eta (standard bins) vs p_{T};p_{T} (GeV);#eta",
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hs.hGenInclusiveJetPtEtaStdBins->Sumw2();
        hs.hGenInclusiveJetPtPtHat = std::make_unique<TH2D>("hGenInclusiveJetPtPtHat", 
                                                            "Gen jet #hat{p}_{T} vs p_{T};p_{T} (GeV);#hat{p}_{T} (GeV)",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtPtHat->Sumw2();

        hs.hGenInclusiveJetJESDefPtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefPtEta", 
                                                                "Gen jet JES default variation pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefPtEta->Sumw2();
        hs.hGenInclusiveJetJESDefExtraPtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefExtraPtEta", 
                                                                "Gen jet JES default variation with extra SF pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefExtraPtEta->Sumw2();
        hs.hGenInclusiveJetJESDefDoublePtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefDoublePtEta", 
                                                                "Gen jet JES default variation with double SF pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefDoublePtEta->Sumw2();
        hs.hGenInclusiveJetJESDefNinetyPtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefNinetyPtEta", 
                                                                "Gen jet JES default variation with 90%% SF pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefNinetyPtEta->Sumw2();

        hs.hGenInclusiveJetDefPtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefPtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefPtEtaFlipped->Sumw2();
        hs.hGenInclusiveJetDefExtraPtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefExtraPtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped, extra SF);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefExtraPtEtaFlipped->Sumw2();
        hs.hGenInclusiveJetDefDoublePtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefDoublePtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped, double SF);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefDoublePtEtaFlipped->Sumw2();
        hs.hGenInclusiveJetDefNinetyPtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefNinetyPtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped, 90%% SF);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefNinetyPtEtaFlipped->Sumw2();
        

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
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        hs.hGenDijetPtEtaLabUnflipped = std::make_unique<TH2D>("hGenDijetPtEtaLabUnflipped", 
                                                               "Gen dijet #eta (lab frame, unflipped) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta",
                                                               nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                               nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hGenDijetPtEtaLabUnflipped->Sumw2();
        hs.hGenDijetPtEtaLab = std::make_unique<TH2D>("hGenDijetPtEtaLab", 
                                                      "Gen dijet #eta (lab frame) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta",
                                                      nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                      nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hGenDijetPtEtaLab->Sumw2();
        hs.hGenDijetPtEtaCM = std::make_unique<TH2D>("hGenDijetPtEtaCM", 
                                                     "Gen dijet #eta (CM frame) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", 
                                                     nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                     nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hGenDijetPtEtaCM->Sumw2();

        // Loop over eta cuts for dijet histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hGenDijetPtEtaLabArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaLab_%d", iCut),
                                                Form("Gen dijet #eta (lab frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaLabArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMControlledSampleArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMControlledSample_%d", iCut),
                                                Form("Gen dijet #eta_{CM}, controlled sample (|#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMControlledSampleArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMStudiedSampleArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMStudiedSample_%d", iCut),
                                                Form("Gen dijet #eta_{CM}, studied sample (|#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMStudiedSampleArr[iCut]->Sumw2();
#ifdef DO_MIXED_PRIOR
            hs.hGenDijetPtEtaCMMixedPriorArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMMixedPrior_%d", iCut),
                                                Form("Gen dijet mixed prior [#mu_{L}=-0.5, #sigma_{L}=1.2#sigma_{0}; #mu_{R}=+0.7, #sigma_{R}=0.8#sigma_{0}], |#eta| < %.1f;p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMixedPriorArr[iCut]->Sumw2();
#endif

            // Smearing studies
            hs.hGenDijetDefPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (def smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefPtEtaCMArr[iCut]->Sumw2();

            hs.hGenDijetDefExtraPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefExtraPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (extra smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefExtraPtEtaCMArr[iCut]->Sumw2();

            hs.hGenDijetDefDoublePtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefDoublePtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (def 2x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefDoublePtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetDefNinetyPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefNinetyPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (def 0.9x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefNinetyPtEtaCMArr[iCut]->Sumw2();

            // Jet-based massless light-cone x estimates, not incoming parton x.
            auto makeXCorrelation = [&](const char *name, const char *label,
                                        double lowExponent, double highExponent) {
                const int bins[3] = {nDijetPtBins, nDijetEtaBins, 120};
                const double low[3] = {dijetPtBins[0], dijetEtaBins[0], std::pow(10., lowExponent)};
                const double high[3] = {dijetPtBins[1], dijetEtaBins[1], std::pow(10., highExponent)};
                auto histogram = std::make_unique<THnSparseD>(
                    Form("%s_%d", name, iCut),
                    Form("Gen jet-based estimate;p_{T}^{ave} (GeV);#eta_{CM}^{dijet};%s", label),
                    3, bins, low, high);
                double edges[121];
                for (int bin = 0; bin <= 120; ++bin)
                    edges[bin] = std::pow(10., lowExponent + (highExponent - lowExponent) * bin / 120.);
                histogram->SetBinEdges(2, edges);
                histogram->Sumw2();
                return histogram;
            };
            hs.hGenDijetPtEtaCMXPbArr[iCut] = makeXCorrelation("hGenDijetPtEtaCMXPb", "x_{Pb}^{jets}", -5., 1.);
            hs.hGenDijetPtEtaCMXpArr[iCut] = makeXCorrelation("hGenDijetPtEtaCMXp", "x_{p}^{jets}", -5., 1.);
            hs.hGenDijetPtEtaCMXPbOverXpArr[iCut] = makeXCorrelation("hGenDijetPtEtaCMXPbOverXp", "x_{Pb}^{jets}/x_{p}^{jets}", -4., 4.);

            auto createPairDiagnostics = [&](UnfoldingPairHistograms &pairs, const char *suffix) {
                pairs.hRecoDijetPtEtaCMMatchedGenSelected[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaCMMatchedGenSelected%s_%d", suffix, iCut),
                    "hRecoDijetPtEtaCMMatchedGenSelected;p_{T}^{ave} (GeV);#eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                pairs.hRecoDijetPtEtaCMMatchedGenSelected[iCut]->Sumw2();
                pairs.hRecoDijetPtEtaCMMatchedGenFound[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaCMMatchedGenFound%s_%d", suffix, iCut),
                    "hRecoDijetPtEtaCMMatchedGenFound;p_{T}^{ave} (GeV);#eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                pairs.hRecoDijetPtEtaCMMatchedGenFound[iCut]->Sumw2();
                pairs.hGenDijetPtEtaCMSelectedRecoFound[iCut] = std::make_unique<TH2D>(
                    Form("hGenDijetPtEtaCMSelectedRecoFound%s_%d", suffix, iCut),
                    "hGenDijetPtEtaCMSelectedRecoFound;p_{T}^{ave} (GeV);#eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                pairs.hGenDijetPtEtaCMSelectedRecoFound[iCut]->Sumw2();
                pairs.hGenDijetPtEtaCMSelectedRecoNotFound[iCut] = std::make_unique<TH2D>(
                    Form("hGenDijetPtEtaCMSelectedRecoNotFound%s_%d", suffix, iCut),
                    "hGenDijetPtEtaCMSelectedRecoNotFound;p_{T}^{ave} (GeV);#eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                pairs.hGenDijetPtEtaCMSelectedRecoNotFound[iCut]->Sumw2();
                pairs.hGenDijetPtEtaCMMatchedRecoSelected[iCut] = std::make_unique<TH2D>(
                    Form("hGenDijetPtEtaCMMatchedRecoSelected%s_%d", suffix, iCut),
                    "hGenDijetPtEtaCMMatchedRecoSelected;p_{T}^{ave} (GeV);#eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                pairs.hGenDijetPtEtaCMMatchedRecoSelected[iCut]->Sumw2();
                pairs.hRecoDijetPtEtaCMSelectedGenNotFound[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaCMSelectedGenNotFound%s_%d", suffix, iCut),
                    "hRecoDijetPtEtaCMSelectedGenNotFound;p_{T}^{ave} (GeV);#eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                pairs.hRecoDijetPtEtaCMSelectedGenNotFound[iCut]->Sumw2();
            };
            createPairDiagnostics(hs.unfoldingPairs, "");
            createPairDiagnostics(hs.unfoldingPairsJerDefExtra, "JerDefExtra");
            createPairDiagnostics(hs.unfoldingPairsControlled, "ControlledSample");
            createPairDiagnostics(hs.unfoldingPairsStudied, "StudiedSample");

            auto createExtended = [&](ExtendedUnfoldingHistograms &extended, const char *suffix) {
                const int bins5[5] = {nDijetPtBins, nDijetEtaBins, nDijetPtBins, nDijetEtaBins, 2};
                const double low5[5] = {dijetPtBins[0], dijetEtaBins[0], dijetPtBins[0], dijetEtaBins[0], .5};
                const double high5[5] = {dijetPtBins[1], dijetEtaBins[1], dijetPtBins[1], dijetEtaBins[1], 2.5};
                const int bins3[3] = {nDijetPtBins, nDijetEtaBins, 2};
                const double low3[3] = {dijetPtBins[0], dijetEtaBins[0], .5};
                const double high3[3] = {dijetPtBins[1], dijetEtaBins[1], 2.5};
                extended.response[iCut] = std::make_unique<THnSparseD>(
                    Form("hUnfoldingExtendedResponse%s_%d", suffix, iCut),
                    "Same-event leading-pair response;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM};Truth fiducial category",
                    5, bins5, low5, high5);
                extended.response[iCut]->Sumw2();
                extended.response[iCut]->GetAxis(4)->SetBinLabel(1, "Outside fiducial");
                extended.response[iCut]->GetAxis(4)->SetBinLabel(2, "Fiducial");
                extended.truth[iCut] = std::make_unique<THnSparseD>(
                    Form("hUnfoldingExtendedTruth%s_%d", suffix, iCut),
                    "Extended Truth;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Truth fiducial category",
                    3, bins3, low3, high3);
                extended.truth[iCut]->Sumw2();
                extended.truth[iCut]->GetAxis(2)->SetBinLabel(1, "Outside fiducial");
                extended.truth[iCut]->GetAxis(2)->SetBinLabel(2, "Fiducial");
                extended.matchedTruth[iCut] = std::make_unique<THnSparseD>(
                    Form("hUnfoldingExtendedMatchedTruth%s_%d", suffix, iCut),
                    "Extended MatchedTruth;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Truth fiducial category",
                    3, bins3, low3, high3);
                extended.matchedTruth[iCut]->Sumw2();
                extended.matchedTruth[iCut]->GetAxis(2)->SetBinLabel(1, "Outside fiducial");
                extended.matchedTruth[iCut]->GetAxis(2)->SetBinLabel(2, "Fiducial");
                extended.miss[iCut] = std::make_unique<THnSparseD>(
                    Form("hUnfoldingExtendedMiss%s_%d", suffix, iCut),
                    "Extended Miss;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Truth fiducial category",
                    3, bins3, low3, high3);
                extended.miss[iCut]->Sumw2();
                extended.miss[iCut]->GetAxis(2)->SetBinLabel(1, "Outside fiducial");
                extended.miss[iCut]->GetAxis(2)->SetBinLabel(2, "Fiducial");
                extended.measured[iCut] = std::make_unique<TH2D>(
                    Form("hUnfoldingExtendedMeasured%s_%d", suffix, iCut),
                    "Extended Measured;Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                extended.measured[iCut]->Sumw2();
                extended.matchedMeasured[iCut] = std::make_unique<TH2D>(
                    Form("hUnfoldingExtendedMatchedMeasured%s_%d", suffix, iCut),
                    "Extended MatchedMeasured;Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                extended.matchedMeasured[iCut]->Sumw2();
                extended.fake[iCut] = std::make_unique<TH2D>(
                    Form("hUnfoldingExtendedFake%s_%d", suffix, iCut),
                    "Extended Fake;Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                extended.fake[iCut]->Sumw2();
            };
            createExtended(hs.unfoldingExtended, "");
            createExtended(hs.unfoldingExtendedJerDefExtra, "JerDefExtra");
            createExtended(hs.unfoldingExtendedControlled, "ControlledSample");
            createExtended(hs.unfoldingExtendedStudied, "StudiedSample");

            // Response matrix and missing dijets histograms
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[iCut] = std::make_unique<THnSparseF>(Form("hGenDijetPtEtaCMVsRecoPtEtaCM_%d", iCut), 
                                                            Form("Gen dijet #eta_{CM} vs Reco dijet #eta_{CM} vs Gen dijet p_{T}^{ave} vs Reco dijet p_{T}^{ave};Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}"), 
                                                            nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr[iCut] = std::make_unique<THnSparseF>(
                                                            Form("hGenDijetPtEtaCMVsRecoPtEtaCMControlledSample_%d", iCut),
                                                            "Controlled-sample gen dijet vs reco dijet response;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                                                            nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr[iCut] = std::make_unique<THnSparseF>(
                                                            Form("hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSample_%d", iCut),
                                                            "Studied-sample gen dijet vs reco dijet response;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                                                            nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr[iCut]->Sumw2();
            // hs.hRefDijetPtEtaCMVsRecoPtEtaCMArr[iCut] = std::make_unique<THnSparseF>(Form("hRefDijetPtEtaCMVsRecoPtEtaCM_%d", iCut),
            //                                                 "Reference-matched dijet response;Ref p_{T}^{ave} (GeV);Ref #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
            //                                                 nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            // hs.hRefDijetPtEtaCMVsRecoPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMMissArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMMiss_%d", iCut), 
                                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMissArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMMissControlledSampleArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMMissControlledSample_%d", iCut),
                                                                Form("Controlled-sample gen dijet miss (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMissControlledSampleArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMMissStudiedSampleArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMMissStudiedSample_%d", iCut),
                                                                Form("Studied-sample gen dijet miss (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMissStudiedSampleArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMFakeArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMFake_%d", iCut),
                                                                Form("Unmatched reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMFakeArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMFakeControlledSampleArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMFakeControlledSample_%d", iCut),
                                                                Form("Controlled-sample unmatched reco dijet (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMFakeControlledSampleArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMFakeStudiedSampleArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMFakeStudiedSample_%d", iCut),
                                                                Form("Studied-sample unmatched reco dijet (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMFakeStudiedSampleArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[iCut] = std::make_unique<THnSparseF>(
                Form("hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCM_%d", iCut),
                "Gen dijet vs reco dijet with eta-dependent default JER;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMMissJerDefExtraArr[iCut] = std::make_unique<TH2D>(
                Form("hGenDijetPtEtaCMMissJerDefExtra_%d", iCut),
                Form("Gen dijet miss for eta-dependent default JER (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMissJerDefExtraArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[iCut] = std::make_unique<TH2D>(
                Form("hRecoDijetPtEtaCMFakeJerDefExtra_%d", iCut),
                Form("Reco dijet fake for eta-dependent default JER (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[iCut]->Sumw2();
            hs.hUnfoldingPairClassificationArr[iCut] = std::make_unique<TH1D>(
                Form("hUnfoldingPairClassification_%d", iCut),
                Form("Unfolding pair classification (|#eta| < %.1f);Pair category;Weighted events", etaCuts[iCut]),
                unfolding_diagnostics::nPairCategories, 0.5,
                unfolding_diagnostics::nPairCategories + 0.5);
            hs.hUnfoldingPairClassificationArr[iCut]->Sumw2();

            auto *classificationAxis = hs.hUnfoldingPairClassificationArr[iCut]->GetXaxis();
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMatchedDirect, "Matched direct");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMatchedSwapped, "Matched swapped");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kGenPassRecoFail, "Gen pass, reco fail");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kRecoPassGenFail, "Reco pass, gen fail");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kSelectedPairMismatch, "Selected-pair mismatch");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMissingOneRecoReference, "Missing one reco reference");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMissingBothRecoReferences, "Missing both reco references");
            hs.hUnfoldingPairClassificationJerDefExtraArr[iCut] = std::make_unique<TH1D>(
                Form("hUnfoldingPairClassificationJerDefExtra_%d", iCut),
                Form("Unfolding pair classification for eta-dependent default JER (|#eta| < %.1f);Pair category;Weighted events", etaCuts[iCut]),
                unfolding_diagnostics::nPairCategories, 0.5,
                unfolding_diagnostics::nPairCategories + 0.5);
            hs.hUnfoldingPairClassificationJerDefExtraArr[iCut]->Sumw2();
            auto *extraClassificationAxis = hs.hUnfoldingPairClassificationJerDefExtraArr[iCut]->GetXaxis();
            for (int category = 1; category <= unfolding_diagnostics::nPairCategories; ++category) {
                extraClassificationAxis->SetBinLabel(category, classificationAxis->GetBinLabel(category));
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)



        //
        // Ref level histograms
        //

        // Ref jet histograms
        hs.hRefInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetPtEtaLabUnflipped", 
                                                                      "Ref jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta;", 
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetPtEtaLab = std::make_unique<TH2D>("hRefInclusiveJetPtEtaLab", 
                                                              "Ref jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta;", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetPtEtaLab->Sumw2();
        hs.hRefInclusiveJetPtEtaCM = std::make_unique<TH2D>("hRefInclusiveJetPtEtaCM", 
                                                              "Ref jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM};", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetPtEtaCM->Sumw2();
        hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetTrkMaxPtEtaLabUnflipped", 
                                                                      "Ref jet track max p_{T} (lab frame, unflipped) vs jet p_{T} (TrkMax selection);p_{T}^{jet} (GeV);#eta;", 
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetTrkMaxPtEtaLab = std::make_unique<TH2D>("hRefInclusiveJetTrkMaxPtEtaLab", 
                                                                    "Ref jet track max p_{T} (lab frame) vs jet p_{T} (TrkMax selection);p_{T}^{jet} (GeV);#eta;", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetTrkMaxPtEtaLab->Sumw2();
        hs.hRefInclusiveJetTrkMaxPtEtaCM = std::make_unique<TH2D>("hRefInclusiveJetTrkMaxPtEtaCM", 
                                                                    "Ref jet track max p_{T} (CM frame) vs jet p_{T} (TrkMax selection);p_{T}^{jet} (GeV);#eta_{CM};", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetTrkMaxPtEtaCM->Sumw2();
        hs.hRefInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRefInclusiveJetPtEtaStdBins", 
                                                                "Ref jet #eta (standard bins) vs p_{T};p_{T} (GeV);#eta;",
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hs.hRefInclusiveJetPtEtaStdBins->Sumw2();
        hs.hRefInclusiveJetNoSelPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetNoSelPtEtaLabUnflipped", 
                                                                      "Ref jet #eta (lab frame, unflipped) vs p_{T} (no selection);p_{T} (GeV);#eta;", 
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetNoSelPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetNoSelPtEtaLab = std::make_unique<TH2D>("hRefInclusiveJetNoSelPtEtaLab", 
                                                              "Ref jet #eta (lab frame) vs p_{T} (no selection);p_{T} (GeV);#eta;", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetNoSelPtEtaLab->Sumw2();
        hs.hRefInclusiveJetNoSelPtEtaCM = std::make_unique<TH2D>("hRefInclusiveJetNoSelPtEtaCM", 
                                                              "Ref jet #eta (CM frame) vs p_{T} (no selection);p_{T} (GeV);#eta_{CM};", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetNoSelPtEtaCM->Sumw2();

        auto makeRefSelInclusiveJetHistogram = [&](const char *name, const char *title) {
            auto histogram = std::make_unique<TH2D>(
                name, title, nJetPtBins, jetPtBins[0], jetPtBins[1],
                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            histogram->Sumw2();
            return histogram;
        };
        hs.hRefSelInclusiveJetPtEtaLabUnflipped = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetPtEtaLabUnflipped", "Ref-selected jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta");
        hs.hRefSelInclusiveJetPtEtaLab = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetPtEtaLab", "Ref-selected jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta");
        hs.hRefSelInclusiveJetPtEtaCM = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetPtEtaCM", "Ref-selected jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM}");
        hs.hRefSelInclusiveJetTrkMaxPtEtaLabUnflipped = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetTrkMaxPtEtaLabUnflipped", "Ref-selected jet #eta (lab frame, unflipped) vs p_{T} (TrkMax);p_{T} (GeV);#eta");
        hs.hRefSelInclusiveJetTrkMaxPtEtaLab = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetTrkMaxPtEtaLab", "Ref-selected jet #eta (lab frame) vs p_{T} (TrkMax);p_{T} (GeV);#eta");
        hs.hRefSelInclusiveJetTrkMaxPtEtaCM = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetTrkMaxPtEtaCM", "Ref-selected jet #eta (CM frame) vs p_{T} (TrkMax);p_{T} (GeV);#eta_{CM}");
        hs.hRefSelInclusiveJetNoSelPtEtaLabUnflipped = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetNoSelPtEtaLabUnflipped", "Ref-selected jet #eta (lab frame, unflipped) vs p_{T} (NoSel);p_{T} (GeV);#eta");
        hs.hRefSelInclusiveJetNoSelPtEtaLab = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetNoSelPtEtaLab", "Ref-selected jet #eta (lab frame) vs p_{T} (NoSel);p_{T} (GeV);#eta");
        hs.hRefSelInclusiveJetNoSelPtEtaCM = makeRefSelInclusiveJetHistogram(
            "hRefSelInclusiveJetNoSelPtEtaCM", "Ref-selected jet #eta (CM frame) vs p_{T} (NoSel);p_{T} (GeV);#eta_{CM}");

        // Ref dijet histograms

        hs.hRefDijetPtEtaLabUnflipped = std::make_unique<TH2D>("hRefDijetPtEtaLabUnflipped", 
                                                               "Ref dijet #eta (lab frame, unflipped) vs p_{T}^{ave} (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta", 
                                                               nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                               nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRefDijetPtEtaLabUnflipped->Sumw2();
        hs.hRefDijetPtEtaLab = std::make_unique<TH2D>("hRefDijetPtEtaLab", 
                                                      "Ref dijet #eta (lab frame) vs p_{T}^{ave} (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta", 
                                                      nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                      nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRefDijetPtEtaLab->Sumw2();
        hs.hRefDijetPtEtaCM = std::make_unique<TH2D>("hRefDijetPtEtaCM", 
                                                     "Ref dijet #eta (CM frame) vs p_{T}^{ave} (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta_{CM}", 
                                                     nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                     nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRefDijetPtEtaCM->Sumw2();

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefDijetPtEtaLabArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaLab_%d", iCut),
                                            Form("Ref dijet #eta (lab frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaLabArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaLabPointingResArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaLabPointingRes_%d", iCut),
                                            Form("Reco dijet #eta (lab frame, matched pair, |#eta^{reco}| < %.1f) vs reco p_{T}^{ave} (pointing resolution);p_{T}^{ave} (GeV);#eta", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaLabPointingResArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaLabPointingResArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaLabPointingRes_%d", iCut),
                                            Form("Ref dijet #eta (lab frame, |#eta| < %.1f) vs p_{T}^{ave} (pointing resolution);p_{T}^{ave} (GeV);#eta", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaLabPointingResArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaLabPointingResRefPtAveArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaLabPointingResRefPtAve_%d", iCut),
                                            Form("Ref dijet #eta (lab frame, |#eta| < %.1f) vs ref p_{T}^{ave} (pointing resolution);p_{T}^{ave,ref} (GeV);#eta", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaLabPointingResRefPtAveArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCM_%d", iCut), 
                                                Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMControlledSampleArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMControlledSample_%d", iCut),
                                                Form("Ref dijet #eta_{CM}, controlled sample (|#eta| < %.1f) vs reco p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMControlledSampleArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMStudiedSampleArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMStudiedSample_%d", iCut),
                                                Form("Ref dijet #eta_{CM}, studied sample (|#eta| < %.1f) vs reco p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMStudiedSampleArr[iCut]->Sumw2();

            hs.hRefDijetPtEtaCMArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMJerUp_%d", iCut),
                                            Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMArrJerUp[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMJerDown_%d", iCut),
                                            Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMArrJerDown[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMJerDef_%d", iCut),
                                            Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMArrJerDef[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMPointingResArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMPointingRes_%d", iCut),
                                            Form("Reco dijet #eta_{CM} (matched pair, |#eta^{reco}| < %.1f) vs reco p_{T}^{ave} (pointing resolution);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMPointingResArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMPointingResArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMPointingRes_%d", iCut),
                                            Form("Ref dijet #eta_{CM} (CM frame, |#eta^{ref}| < %.1f) vs reco p_{T}^{ave} (pointing resolution);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMPointingResArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaCMPointingResRefPtAveArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMPointingResRefPtAve_%d", iCut),
                                            Form("Ref dijet #eta_{CM} (CM frame, |#eta^{ref}| < %.1f) vs ref p_{T}^{ave} (pointing resolution);p_{T}^{ave,ref} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMPointingResRefPtAveArr[iCut]->Sumw2();
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)

        //
        // Ref-selected level histograms
        //

        // Loop over eta cuts for dijet histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefSelDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hRefSelDijetPtEtaCM_%d", iCut), 
                                                Form("Ref-selected dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefSelDijetPtEtaCMArr[iCut]->Sumw2();
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    } // if (fIsMC)


    //
    // Reco level histograms
    //

    // Reco jet histograms
    hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnflippedUnweighted", 
                                                                    "Reco jet #eta (lab frame, unflipped, unweighted) vs p_{T};p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted->Sumw2();
    hs.hRecoInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnflipped", 
                                                                    "Reco jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetPtEtaLab = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLab", 
                                                          "Reco jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta", 
                                                          nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                          nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaLab->Sumw2();
    hs.hRecoInclusiveJetPtEtaCM = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaCM", 
                                                          "Reco jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM}", 
                                                          nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                          nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaCM->Sumw2();
    hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetTrkMaxPtEtaLabUnflipped", 
                                                                        "Reco jet #eta (lab frame, unflipped) vs p_{T} (TrkMax selection);p_{T} (GeV);#eta", 
                                                                        nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetTrkMaxPtEtaLab = std::make_unique<TH2D>("hRecoInclusiveJetTrkMaxPtEtaLab", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (TrkMax selection);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetTrkMaxPtEtaLab->Sumw2();
    hs.hRecoInclusiveJetTrkMaxPtEtaCM = std::make_unique<TH2D>("hRecoInclusiveJetTrkMaxPtEtaCM", 
                                                                    "Reco jet #eta (CM frame) vs p_{T} (TrkMax selection);p_{T} (GeV);#eta_{CM}", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetTrkMaxPtEtaCM->Sumw2();
    hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetNoSelPtEtaLabUnflipped", 
                                                                    "Reco jet #eta (lab frame, unflipped) vs p_{T} (no selection);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetNoSelPtEtaLab = std::make_unique<TH2D>("hRecoInclusiveJetNoSelPtEtaLab", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (no selection);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetNoSelPtEtaLab->Sumw2();
    hs.hRecoInclusiveJetNoSelPtEtaCM = std::make_unique<TH2D>("hRecoInclusiveJetNoSelPtEtaCM", 
                                                                    "Reco jet #eta (CM frame) vs p_{T} (no selection);p_{T} (GeV);#eta_{CM}", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetNoSelPtEtaCM->Sumw2();
    // for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
    //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
    //                                           Form("Reco jet #eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
    //                                           nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
    //     hs.hRecoInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShifted_%d", iShift), 
    //                                           Form("Reco jet #eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
    //                                           nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    //     hs.hRecoInclusiveJetEtaCMShifted[iShift]->Sumw2();
    // } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

    hs.hRecoInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaStdBins", 
                                                            "Reco jet #eta (standard bins) vs p_{T};p_{T} (GeV);#eta",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hs.hRecoInclusiveJetPtEtaStdBins->Sumw2();
    hs.hRecoInclusiveJetPtPtHat = std::make_unique<TH2D>("hRecoInclusiveJetPtPtHat", 
                                                            "Reco jet #hat{p}_{T} vs p_{T};p_{T} (GeV);#hat{p}_{T} (GeV)",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtPtHat->Sumw2();

    if (isMc) {
        hs.hRecoInclusiveJetPtEtaLabMatched = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabMatched", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (matched to ref);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetPtEtaLabMatched->Sumw2();
        hs.hRecoInclusiveJetPtEtaLabUnmatched = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnmatched", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (unmatched to ref);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetPtEtaLabUnmatched->Sumw2();
        hs.hRecoInclusiveJetJESPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESPtEta", 
                                                                    "Reco jet JES (reco/ref) vs pT vs #eta (standard bins);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESPtEta->Sumw2();

        hs.hRecoInclusiveJetJESDefNoSFPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefNoSFPtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins, no scaling);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefNoSFPtEta->Sumw2();
        hs.hRecoInclusiveJetJESDefDoublePtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefDoublePtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins, double scaling);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefDoublePtEta->Sumw2();
        hs.hRecoInclusiveJetJESDefNinetyPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefNinetyPtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins, 90% scaling);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefNinetyPtEta->Sumw2();
        hs.hRecoInclusiveJetJESDefPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefPtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefPtEta->Sumw2();
        hs.hRecoInclusiveJetJESUpPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESUpPtEta", 
                                                                    "Reco jet JES up (reco/ref) vs pT vs #eta (standard bins);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESUpPtEta->Sumw2();
        hs.hRecoInclusiveJetJESDownPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDownPtEta", 
                                                                        "Reco jet JES down (reco/ref) vs pT vs #eta (standard bins);p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta",
                                                                        nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                        nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDownPtEta->Sumw2();
    }

    // Reco dijet histograms

    hs.hRecoDijetPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoDijetPtEtaLabUnflipped", 
                                                            "Reco dijet #eta (lab frame, unflipped) vs p_{T}^{ave} (with |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta",
                                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaLabUnflipped->Sumw2();
    hs.hRecoDijetPtEtaLab = std::make_unique<TH2D>("hRecoDijetPtEtaLab", 
                                                   "Reco dijet #eta (lab frame) vs p_{T}^{ave} (with |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta",
                                                   nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                   nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaLab->Sumw2();
    hs.hRecoDijetPtEtaCM = std::make_unique<TH2D>("hRecoDijetPtEtaCM", 
                                                  "Reco dijet #eta_{CM} (CM frame) vs p_{T}^{ave} (with |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta_{CM}",
                                                  nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                  nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaCM->Sumw2();

    hs.hRecoInclusiveJetRawPtEtaLabUnflipped = std::make_unique<TH2D>(
            "hRecoInclusiveJetRawPtEtaLabUnflipped",
            "Selected reco jet #eta (lab frame, unflipped) vs raw p_{T};raw p_{T} (GeV);#eta",
            nJetPtBins, jetPtBins[0], jetPtBins[1],
            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetRawPtEtaLabUnflipped->Sumw2();

    hs.hRecoDijetPtEtaCMRunArr[nDiagnosticRuns] = std::make_unique<TH2D>(
        "hRecoDijetPtEtaCMRunIntegrated",
        "Reco dijet #eta_{CM} vs p_{T}^{ave}, all runs (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta_{CM}^{dijet}",
        nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
        nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaCMRunArr[nDiagnosticRuns]->Sumw2();
    for (int iRun{0}; iRun < nDiagnosticRuns; ++iRun) {
        hs.hRecoDijetPtEtaCMRunArr[iRun] = std::make_unique<TH2D>(
            Form("hRecoDijetPtEtaCMRun_%u", diagnosticRunIds[iRun]),
            Form("Reco dijet #eta_{CM} vs p_{T}^{ave}, run %u (PU %.3g, |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta_{CM}^{dijet}",
                 diagnosticRunIds[iRun], diagnosticRunPileup[iRun]),
            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRecoDijetPtEtaCMRunArr[iRun]->Sumw2();
    }

    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        hs.hRecoDijetPtEtaLabArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaLab_%d", iCut),
                                            Form("Reco dijet #eta (lab frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRecoDijetPtEtaLabArr[iCut]->Sumw2();
        hs.hRecoDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCM_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRecoDijetPtEtaCMArr[iCut]->Sumw2();

        if (isMc) {
            hs.hRecoDijetPtEtaCMControlledSampleArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMControlledSample_%d", iCut),
                                                Form("Reco dijet #eta_{CM}, controlled sample (|#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMControlledSampleArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMStudiedSampleArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMStudiedSample_%d", iCut),
                                                Form("Reco dijet #eta_{CM}, studied sample (|#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMStudiedSampleArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefNoSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER no SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefDoubleSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default, double SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefNinetySF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default, 90%% SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut]->Sumw2();

            hs.hRecoDijetPtEtaCMArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDown[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDef_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDef[iCut]->Sumw2();

            auto makeRecoJerExtraHistogram = [&](const char *variation, const char *title,
                                               std::unique_ptr<TH2D> *fullArr) {
                fullArr[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaCM%s_%d", variation, iCut),
                    Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (%s);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut], title),
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                fullArr[iCut]->Sumw2();
            };
            makeRecoJerExtraHistogram("JerDownExtra", "JER down with eta-dependent scaling",
                hs.hRecoDijetPtEtaCMArrJerDownExtra);
            makeRecoJerExtraHistogram("JerDefExtra", "JER default with eta-dependent scaling",
                hs.hRecoDijetPtEtaCMArrJerDefExtra);
            makeRecoJerExtraHistogram("JerUpExtra", "JER up with eta-dependent scaling",
                hs.hRecoDijetPtEtaCMArrJerUpExtra);
            makeRecoJerExtraHistogram("JerDefExtraUnfold", "JER default with eta-dependent scaling, unfolding",
                hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold);

            hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefUnfold_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default, unfolding);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMJerDefUnfoldControlledSampleArr[iCut] = std::make_unique<TH2D>(
                                            Form("hRecoDijetPtEtaCMJerDefUnfoldControlledSample_%d", iCut),
                                            Form("Reco dijet #eta_{CM}, controlled sample (JER default, unfolding, |#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMJerDefUnfoldControlledSampleArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMJerDefUnfoldStudiedSampleArr[iCut] = std::make_unique<TH2D>(
                                            Form("hRecoDijetPtEtaCMJerDefUnfoldStudiedSample_%d", iCut),
                                            Form("Reco dijet #eta_{CM}, studied sample (JER default, unfolding, |#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMJerDefUnfoldStudiedSampleArr[iCut]->Sumw2();

        } // if (isMc)

        if (!isMc) {
            hs.hRecoDijetPtEtaCMArrJeuUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJeuUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JEU up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJeuUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJeuDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJeuDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JEU down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJeuDown[iCut]->Sumw2();
        } // if (!isMc)
    } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)


    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
} // for (int iCut{0}; iCut < nEtaCuts; ++iCut)

//________________
/// Set up the chains for the different trees
/// @param input: the input file or file list
/// @param hltTree: the HLT tree
/// @param eventTree: the event tree
/// @param skimTree: the skim tree

}
