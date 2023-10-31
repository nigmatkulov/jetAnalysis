/**
 * @file HistoManagerJetESR.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for JES and JER studies
 * @version 0.1
 * @date 2023-10-24
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "HistoManagerJetESR.h"

// ROOT headers
#include "TObject.h"
#include "TIterator.h"
#include "TString.h"
#include "TClass.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

ClassImp(HistoManagerJetESR)

//________________
HistoManagerJetESR::HistoManagerJetESR() :
  fIsMc{kFALSE}, 
  fCentBins{10}, fCentRange{-10., 90.},
  fPtBins{50}, fPtRange{0., 500.}, 
  fEtaBins{50}, fEtaRange{-2.5, 2.5},
  fPhiBins{32}, fPhiRange{-TMath::Pi(), TMath::Pi()},
  fJESBins{500}, fJESRange{0., 5.},
  fPtHatBins{10}, fPtHatRange{15., 215.},
  fFlavorForBBins{14}, fFlavorForBRange{-6.5, 6.5},

  hVz{nullptr}, hVzWeighted{nullptr}, hMult{nullptr},
  hHiBin{nullptr}, hHiBinWeighted{nullptr},
  hPtHat{nullptr}, hPtHatWeighted{nullptr}, hPtHatWeight{nullptr},
  hCentrality{nullptr}, hCentralityWeighted{nullptr},

  hNGenJets{nullptr},
  hGenJetPtEtaPhiCent{nullptr}, hGenJetPtEtaPhiCentWeighted{nullptr},
  hGenJetPtFlavPtHatCent{nullptr}, hGenJetPtFlavPtHatCentWeighted{nullptr},

  hNRecoJets{nullptr},
  hRecoJetRawPtEtaPhiCent{nullptr}, hRecoJetPtEtaPhiCent{nullptr},
  hRecoJetPtEtaPhiCentWeighted{nullptr}, hRecoJetPtFlavPtHatCent{nullptr}, hRecoJetPtFlavPtHatCentWeighted{nullptr},
  hRecoJetPtFlavPtHatCentInclusive{nullptr}, hRecoJetPtFlavPtHatCentInclusiveWeighted{nullptr},
  hRecoJetDeltaRPtCent{nullptr},

  hNRefJets{nullptr},
  hRefJetPtEtaPhiCent{nullptr}, hRefJetPtEtaPhiCentWeighted{nullptr},
  hRefJetPtFlavPtHatCent{nullptr}, hRefJetPtFlavPtHatCentWeighted{nullptr},
  hRecoJetPtRefJetPtPtHatCent{nullptr},
  hRecoJetPtRefJetPtPtHatCentWeighted{nullptr},

  hJESPtEtaPhiCent{nullptr}, hJESPtEtaPhiCentWeighted{nullptr}, 
  hJESRawPtFlavPtHatCent{nullptr}, hJESRawPtFlavPtHatCentWeighted{nullptr},
  hJESPtFlavPtHatCent{nullptr}, hJESPtFlavPtHatCentWeighted{nullptr} { 
    /* Empty */
}

//________________
HistoManagerJetESR::~HistoManagerJetESR() {
    if (hVz)            delete hVz;
    if (hVzWeighted)    delete hVzWeighted;
    if (hMult)          delete hMult;
    if (hHiBin)         delete hHiBin;
    if (hHiBinWeighted) delete hHiBinWeighted;
    if (hPtHat)         delete hPtHat;
    if (hPtHatWeighted) delete hPtHatWeighted;
    if (hPtHatWeight)   delete hPtHatWeight;
    if (hCentrality)    delete hCentrality;
    if (hCentralityWeighted) delete hCentralityWeighted;

    if (hNRecoJets)              delete hNRecoJets;
    if (hRecoJetRawPtEtaPhiCent) delete hRecoJetRawPtEtaPhiCent;
    if (hRecoJetPtEtaPhiCent)    delete hRecoJetPtEtaPhiCent;
    if (hRecoJetPtEtaPhiCentWeighted)    delete hRecoJetPtEtaPhiCentWeighted;
    if (hRecoJetPtFlavPtHatCent)         delete hRecoJetPtFlavPtHatCent;
    if (hRecoJetPtFlavPtHatCentWeighted) delete hRecoJetPtFlavPtHatCentWeighted;
    if (hRecoJetPtFlavPtHatCentInclusive) delete hRecoJetPtFlavPtHatCentInclusive;
    if (hRecoJetPtFlavPtHatCentInclusiveWeighted) delete hRecoJetPtFlavPtHatCentInclusiveWeighted;
    if (hRecoJetDeltaRPtCent) delete hRecoJetDeltaRPtCent;

    if (fIsMc) {

        if (hNGenJets)                   delete hNGenJets;
        if (hGenJetPtEtaPhiCent)         delete hGenJetPtEtaPhiCent;
        if (hGenJetPtEtaPhiCentWeighted) delete hGenJetPtEtaPhiCentWeighted;
        if (hGenJetPtFlavPtHatCent)      delete hGenJetPtFlavPtHatCent;
        if (hGenJetPtFlavPtHatCentWeighted) delete hGenJetPtFlavPtHatCentWeighted;

        if (hNRefJets)                   delete hNRefJets;
        if (hRefJetPtEtaPhiCent)         delete hRefJetPtEtaPhiCent;
        if (hRefJetPtEtaPhiCentWeighted) delete hRefJetPtEtaPhiCentWeighted;
        if (hRefJetPtFlavPtHatCent)      delete hRefJetPtFlavPtHatCent;
        if (hRefJetPtFlavPtHatCentWeighted) delete hRefJetPtFlavPtHatCentWeighted;

        if (hRecoJetPtRefJetPtPtHatCent)         delete hRecoJetPtRefJetPtPtHatCent;
        if (hRecoJetPtRefJetPtPtHatCentWeighted) delete hRecoJetPtRefJetPtPtHatCentWeighted;

        if (hJESPtEtaPhiCent)         delete hJESPtEtaPhiCent;
        if (hJESPtEtaPhiCentWeighted) delete hJESPtEtaPhiCentWeighted;
        if (hJESRawPtFlavPtHatCent)   delete hJESRawPtFlavPtHatCent;
        if (hJESRawPtFlavPtHatCentWeighted) delete hJESRawPtFlavPtHatCentWeighted;
        if (hJESPtFlavPtHatCent)            delete hJESPtFlavPtHatCent;
        if (hJESPtFlavPtHatCentWeighted)    delete hJESPtFlavPtHatCentWeighted;

    } // if (fIsMc)
}

//________________
void HistoManagerJetESR::init(const Bool_t& isMc) {
    
    Int_t    prescale = 1;

    Int_t    bins4D_jet_PtEtaPhiCent[4] = { fPtBins    , fEtaBins    , fPhiBins    , fCentBins     };
    Double_t xmin4D_jet_PtEtaPhiCent[4] = { fPtRange[0], fEtaRange[0], fPhiRange[0], fCentRange[0] };
    Double_t xmax4D_jet_PtEtaPhiCent[4] = { fPtRange[1], fEtaRange[1], fPhiRange[1], fCentRange[1] };

    Int_t    bins4D_jet_PtFlavPthatCent[4] = { fPtBins    , fFlavorForBBins    , fPtHatBins    , fCentBins     };
    Double_t xmin4D_jet_PtFlavPthatCent[4] = { fPtRange[0], fFlavorForBRange[0], fPtHatRange[0], fCentRange[0] };
    Double_t xmax4D_jet_PtFlavPthatCent[4] = { fPtRange[1], fFlavorForBRange[1], fPtHatRange[1], fCentRange[1] };

    Int_t    bins5D_jes_PtEtaPhiCent[5] = { fJESBins    , fPtBins    , fEtaBins    , fPhiBins    , fCentBins  };
    Double_t xmin5D_jes_PtEtaPhiCent[5] = { fJESRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fCentRange[0] };
    Double_t xmax5D_jes_PtEtaPhiCent[5] = { fJESRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fCentRange[1] };

    Int_t    bins5D_jes_PtFlavPthatCent[5] = { fJESBins    , fPtBins    , fFlavorForBBins    , fPtHatBins    , fCentBins     };
    Double_t xmin5D_jes_PtFlavPthatCent[5] = { fJESRange[0], fPtRange[0], fFlavorForBRange[0], fPtHatRange[0], fCentRange[0] };
    Double_t xmax5D_jes_PtFlavPthatCent[5] = { fJESRange[1], fPtRange[1], fFlavorForBRange[1], fPtHatRange[1], fCentRange[1] };

    Int_t    multBins{1800};
    Double_t multRange[2] = {-0.5, 1799.5};
    Int_t    hiBinBins{203};
    Double_t hiBinRange[2] = {-1.5, 201.5};
    Int_t    centralityBins{101};
    Double_t centralityRange[2] = {-0.5, 100.5};
    Int_t    weightBins{110};
    Double_t weightRange[2] = {-0.05, 1.05};
    Int_t    ptHatBins{100};
    Double_t ptHatRange[2] = {0., 1000.};

    Int_t    deltaRBins{80};
    Double_t deltaRRange[2] = {-0.05, 0.35};
    Int_t    bins3D_jet_DeltaRPtCent[3] = { deltaRBins    ,  fPtBins   , fCentBins     };
    Double_t xmin3D_jet_DeltaRPtCent[3] = { deltaRRange[0], fPtRange[0], fCentRange[0] };
    Double_t xmax3D_jet_DeltaRPtCent[3] = { deltaRRange[1], fPtRange[1], fCentRange[1] };

    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", 400, -50., 50.);
    hVz->Sumw2();
    hVzWeighted = new TH1D("hVzWeighted","Vertex z position;vz (cm);Entries", 400, -50., 50.);
    hVzWeighted->Sumw2();
    hMult = new TH1D("hMult","Charged particle multiplicity;Multiplicity;Entries", 
                     multBins, multRange[0], multRange[1]);
    hMult->Sumw2();
    hHiBin = new TH1D("hHiBin","HiBin a.k.a. centrality;HiBin;Entries", 
                      hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBin->Sumw2();
    hHiBinWeighted = new TH1D("hHiBinWeighted","HiBin a.k.a. centrality with #hat{p_{T}} weight;HiBin;Entries", 
                              hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinWeighted->Sumw2();
    hPtHat = new TH1D("hPtHat","#hat{p_{T}};#hat{p_{T}} (GeV/c);Entries", 
                      ptHatBins, ptHatRange[0], ptHatRange[1]);
    hPtHat->Sumw2();
    hPtHatWeighted = new TH1D("hPtHatWeighted","#hat{p_{T}} with #hat{p_{T}} weight;#hat{p_{T}} (GeV/c);Entries", 
                              ptHatBins, ptHatRange[0], ptHatRange[1]);
    hPtHatWeighted->Sumw2();
    hPtHatWeight = new TH1D("hPtHatWeight","#hat{p_{T}} weight;#hat{p_{T}} weight;Entries", 
                            weightBins, weightRange[0], weightRange[1]);
    hPtHatWeight->Sumw2();
    hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries",
                           centralityBins, centralityRange[0], centralityRange[1]);
    hCentrality->Sumw2();
    hCentralityWeighted = new TH1D("hCentralityWeighted","Collision centrality weighted;Centrality (%);Entries",
                                   centralityBins, centralityRange[0], centralityRange[1]);
    hCentralityWeighted->Sumw2();

    //
    // Gen jets
    //

    hNGenJets = new TH1D("hNGenJets","Number of generated jets", 16, -0.5, 15.5);
    hNGenJets->Sumw2();
    hGenJetPtEtaPhiCent = new THnSparseD("hGenJetPtEtaPhiCent","Generated jet;p_{T}^{gen} (GeV/c);#eta;#phi (rad);centrality", 
                                          4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
    hGenJetPtEtaPhiCent->Sumw2();
    hGenJetPtEtaPhiCentWeighted = new THnSparseD("hGenJetPtEtaPhiCentWeighted","Generated jet weighted;p_{T}^{gen} (GeV/c);#eta;#phi (rad);centrality", 
                                                  4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
    hGenJetPtEtaPhiCentWeighted->Sumw2();
    hGenJetPtFlavPtHatCent = new THnSparseD("hGenJetPtFlavPtHatCent","Generated jet;p_{T}^{gen} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                             4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
    hGenJetPtFlavPtHatCent->Sumw2();
    hGenJetPtFlavPtHatCentWeighted = new THnSparseD("hGenJetPtFlavPtHatCentWeighted","Generated jet;p_{T}^{gen} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                                     4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
    hGenJetPtFlavPtHatCentWeighted->Sumw2();

    //
    // Reco jets
    //

    hNRecoJets = new TH1D("hNRecoJets","Number of reconstructed jets", 16, -0.5, 15.5);
    hNRecoJets->Sumw2();
    hRecoJetRawPtEtaPhiCent = new THnSparseD("hRecoJetRawPtEtaPhiCent","Reconstructed jet with raw p_{T};p_{T}^{raw} (GeV/c);#eta;#phi (rad);centrality", 
                                             4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
    hRecoJetRawPtEtaPhiCent->Sumw2();
    hRecoJetPtEtaPhiCent = new THnSparseD("hRecoJetPtEtaPhiCent","Reconstructed jet with p_{T};p_{T}^{corr} (GeV/c);#eta;#phi (rad);centrality", 
                                          4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
    hRecoJetPtEtaPhiCent->Sumw2();
    hRecoJetPtEtaPhiCentWeighted = new THnSparseD("hRecoJetPtEtaPhiCentWeighted","Reconstructed jet with p_{T} with weight;p_{T}^{corr} (GeV/c);#eta;#phi (rad);centrality", 
                                                  4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
    hRecoJetPtEtaPhiCentWeighted->Sumw2();

    hRecoJetPtFlavPtHatCent = new THnSparseD("hRecoJetPtFlavPtHatCent","Reconstructed jet;p_{T}^{corr} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                              4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
    hRecoJetPtFlavPtHatCent->Sumw2();
    hRecoJetPtFlavPtHatCentWeighted = new THnSparseD("hRecoJetPtFlavPtHatCentWeighted","Reconstructed jet weighted;p_{T}^{corr} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                                      4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
    hRecoJetPtFlavPtHatCentWeighted->Sumw2();
    hRecoJetPtFlavPtHatCentInclusive = new THnSparseD("hRecoJetPtFlavPtHatCentInclusive","Reconstructed jet (matched+unmatched);p_{T}^{corr} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                                       4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
    hRecoJetPtFlavPtHatCentInclusive->Sumw2();
    hRecoJetPtFlavPtHatCentInclusiveWeighted = new THnSparseD("hRecoJetPtFlavPtHatCentInclusiveWeighted","Reconstructed jet (matched+unmatched) weighted;p_{T}^{corr} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                                               4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
    hRecoJetPtFlavPtHatCentInclusiveWeighted->Sumw2();
    hRecoJetDeltaRPtCent = new THnSparseD("hRecoJetDeltaRPtCent","Reconstructed jet #Delta R;#DeltaR=#sqrt{(#eta-#eta_{WTA})^{2}+(#phi-#phi_{WTA})^{2}},p_{T}^{corr} (GeV/c);centrality",
                                          3, bins3D_jet_DeltaRPtCent, xmin3D_jet_DeltaRPtCent, xmax3D_jet_DeltaRPtCent);
    hRecoJetDeltaRPtCent->Sumw2();

    if (fIsMc) {

        //
        // Ref jets
        //

        hNRefJets = new TH1D("hNRefJets","Number of reference jets", 16, -0.5, 15.5);
        hNRefJets->Sumw2();
        hRefJetPtEtaPhiCent = new THnSparseD("hRefJetPtEtaPhiCent","Reference jet;p_{T}^{gen} (GeV/c);#eta;#phi (rad);centrality", 
                                              4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
        hRefJetPtEtaPhiCent->Sumw2();
        hRefJetPtEtaPhiCentWeighted = new THnSparseD("hRefJetPtEtaPhiCentWeighted","Reference jet weighted;p_{T}^{gen} (GeV/c);#eta;#phi (rad);centrality", 
                                                    4, bins4D_jet_PtEtaPhiCent, xmin4D_jet_PtEtaPhiCent, xmax4D_jet_PtEtaPhiCent);
        hRefJetPtEtaPhiCentWeighted->Sumw2();
        hRefJetPtFlavPtHatCent = new THnSparseD("hRefJetPtFlavPtHatCent","Reference jet;p_{T}^{gen} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                                 4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
        hRefJetPtFlavPtHatCent->Sumw2();
        hRefJetPtFlavPtHatCentWeighted = new THnSparseD("hRefJetPtFlavPtHatCentWeighted","Reference jet;p_{T}^{gen} (GeV/c);flavorForB;#hat{p_{T}} (GeV/c);centrality",
                                                         4, bins4D_jet_PtFlavPthatCent, xmin4D_jet_PtFlavPthatCent, xmax4D_jet_PtFlavPthatCent);
        hRefJetPtFlavPtHatCentWeighted->Sumw2();


        //
        // Jet energy scale and resolution
        //


        hJESPtEtaPhiCent = new THnSparseD("hJESPtEtaPhiCent", "Jet Energy Scale;p_{T}^{corr} / p_{T}^{gen};p_{T}^{gen} (GeV/c);#eta;#phi (rad);centrality", 
                                          5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
        hJESPtEtaPhiCent->Sumw2();
        hJESPtEtaPhiCentWeighted = new THnSparseD("hJESPtEtaPhiCentWeighted", "Jet Energy Scale weighted;p_{T}^{corr} / p_{T}^{gen};p_{T}^{gen} (GeV/c);#eta;#phi (rad);centrality", 
                                                  5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
        hJESPtEtaPhiCentWeighted->Sumw2();

        hJESRawPtFlavPtHatCent = new THnSparseD("hJESRawPtFlavPtHatCent", "Jet Energy Scale (raw p_{T});p_{T}^{raw} / p_{T}^{gen};p_{T}^{gen} (GeV/c);FlavorForB;#hat{p_{T}} (GeV/c);centrality", 
                                                5, bins5D_jes_PtFlavPthatCent, xmin5D_jes_PtFlavPthatCent, xmax5D_jes_PtFlavPthatCent);
        hJESRawPtFlavPtHatCent->Sumw2();
        hJESRawPtFlavPtHatCentWeighted = new THnSparseD("hJESRawPtFlavPtHatCentWeighted", "Jet Energy Scale (raw p_{T}) weighted;p_{T}^{raw} / p_{T}^{gen};p_{T}^{gen} (GeV/c);FlavorForB;#hat{p_{T}} (GeV/c);centrality", 
                                                        5, bins5D_jes_PtFlavPthatCent, xmin5D_jes_PtFlavPthatCent, xmax5D_jes_PtFlavPthatCent);
        hJESRawPtFlavPtHatCentWeighted->Sumw2();

        hJESPtFlavPtHatCent = new THnSparseD("hJESPtFlavPtHatCent", "Jet Energy Scale;p_{T}^{corr} / p_{T}^{gen};p_{T}^{gen} (GeV/c);FlavorForB;#hat{p_{T}} (GeV/c);centrality", 
                                             5, bins5D_jes_PtFlavPthatCent, xmin5D_jes_PtFlavPthatCent, xmax5D_jes_PtFlavPthatCent);
        hJESPtFlavPtHatCent->Sumw2();
        hJESPtFlavPtHatCentWeighted = new THnSparseD("hJESPtFlavPtHatCentWeighted", "Jet Energy Scale weighted;p_{T}^{corr} / p_{T}^{gen};p_{T}^{gen} (GeV/c);FlavorForB;#hat{p_{T}} (GeV/c);centrality", 
                                                     5, bins5D_jes_PtFlavPthatCent, xmin5D_jes_PtFlavPthatCent, xmax5D_jes_PtFlavPthatCent);
        hJESPtFlavPtHatCentWeighted->Sumw2();

    }
}

//________________
void HistoManagerJetESR::writeOutput() {
    // Event
    hVz->Write();
    hVzWeighted->Write();
    hMult->Write();
    hHiBin->Write();
    hHiBinWeighted->Write();
    hPtHat->Write();
    hPtHatWeighted->Write();
    hPtHatWeight->Write();
    hCentrality->Write();
    hCentralityWeighted->Write();

    // Reco jets
    hNRecoJets->Write();
    hRecoJetRawPtEtaPhiCent->Write();
    hRecoJetPtEtaPhiCent->Write();
    hRecoJetPtEtaPhiCentWeighted->Write();
    hRecoJetPtFlavPtHatCent->Write();
    hRecoJetPtFlavPtHatCentWeighted->Write();
    hRecoJetPtFlavPtHatCentInclusive->Write();
    hRecoJetPtFlavPtHatCentInclusiveWeighted->Write();
    hRecoJetDeltaRPtCent->Write();

    if (fIsMc) {
        // Gen jets
        hNGenJets->Write();
        hGenJetPtEtaPhiCent->Write();
        hGenJetPtEtaPhiCentWeighted->Write();
        hGenJetPtFlavPtHatCent->Write();
        hGenJetPtFlavPtHatCentWeighted->Write();

        // Ref jets
        hNRefJets->Write();
        hRefJetPtEtaPhiCent->Write();
        hRefJetPtEtaPhiCentWeighted->Write();
        hRefJetPtFlavPtHatCent->Write();
        hRefJetPtFlavPtHatCentWeighted->Write();

        if (hRecoJetPtRefJetPtPtHatCent) hRecoJetPtRefJetPtPtHatCent->Write();
        if (hRecoJetPtRefJetPtPtHatCentWeighted) hRecoJetPtRefJetPtPtHatCentWeighted->Write();

        // Jet Energy Scale
        hJESPtEtaPhiCent->Write();
        hJESPtEtaPhiCentWeighted->Write();
        hJESRawPtFlavPtHatCent->Write();
        hJESRawPtFlavPtHatCentWeighted->Write();
        hJESPtFlavPtHatCent->Write();
        hJESPtFlavPtHatCentWeighted->Write();
    }
}