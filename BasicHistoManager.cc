// Jet analysis headers
#include "BasicHistoManager.h"

// ROOT headers
#include "TObject.h"
#include "TIterator.h"
#include "TString.h"
#include "TClass.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

ClassImp(BasicHistoManager)

//________________
BasicHistoManager::BasicHistoManager() :
  fIsMc{kFALSE}, 
  fCentBins{9}, fCentRange{0., 90.},
  fJetPtBins{200}, fJetPtRange{0., 1000.}, 
  fJetEtaBins{10}, fJetEtaRange{-2.5, 2.5},
  fJetPhiBins{32}, fJetPhiRange{-TMath::Pi(), TMath::Pi()},
  fJESBins{500}, fJESRange{0., 5.},
  fJERBins{400}, fJERRange{-2., 2.},
  fPtHatBins{10}, fPtHatRange{15., 215.},
  hVz{nullptr}, hVzWeighted{nullptr}, hMult{nullptr},
  hHiBin{nullptr}, hHiBinWeighted{nullptr},
  hPtHat{nullptr}, hPtHatWeighted{nullptr}, hPtHatWeight{nullptr},
  hCentrality{nullptr}, 
  hNRecoJets{nullptr}, hRecoJetPtRaw{nullptr}, hRecoJetPtCorrVsPtRaw{nullptr},
  hRecoJetPt{nullptr}, hRecoJetPtWeighted{nullptr}, hRecoJetEta{nullptr},
  hRecoJetPhi{nullptr}, hRecoJetPtVsEta{nullptr}, 
  hRecoJetPhiVsPt{nullptr}, hRecoJetEtaVsPhi{nullptr},
  hRecoJet{nullptr}, hRecoJetCorr{nullptr}, hRecoJetCorrWeighted{nullptr},
  hRecoJetCorrEtaPhiPtHat{nullptr}, hRecoJetCorrEtaPhiPtHatWeighted{nullptr},
  hNRefJets{nullptr}, hRefJetPt{nullptr}, hRefJetPtWeighted{nullptr},
  hRefJetEta{nullptr}, hRefJetPhi{nullptr}, hRefJetPtVsEta{nullptr},
  hPtHatVsGenPtWOPtJetGtPtHatCut{nullptr},
  hPtHatVsGenPtWeightedWoPtJetGtPtHatCut{nullptr},
  hPtHatVsGenPtWithPtJetGtPtHatCut{nullptr},
  hPtHatVsGenPtWeightedWithPtJetGtPtHatCut{nullptr},
  hRefJet{nullptr}, hRefJetWeighted{nullptr},
  hRefJetEtaPhiPtHat{nullptr}, hRefJetEtaPhiPtHatWeighted{nullptr},
  hJESRaw{nullptr}, hJESRawWeighted{nullptr}, hJESReco{nullptr}, hJESRecoWeighted{nullptr},
  hJESRecoEtaPhiPtHat{nullptr}, hJESRecoEtaPhiPtHatWeighted{nullptr},
  hJERRaw{nullptr}, hJERRawWeighted{nullptr}, hJERReco{nullptr}, hJERRecoWeighted{nullptr},
  hRecoJetPtCorrVsGenPt{nullptr} { 
    /* Empty */
}

//________________
BasicHistoManager::~BasicHistoManager() {
    if (hVz) delete hVz;
    if (hVzWeighted) delete hVzWeighted;
    if (hMult) delete hMult;
    if (hHiBin) delete hHiBin;
    if (hHiBinWeighted) delete hHiBinWeighted;
    if (hPtHat) delete hPtHat;
    if (hPtHatWeighted) delete hPtHatWeighted;
    if (hPtHatWeight) delete hPtHatWeight;
    if (hCentrality) delete hCentrality;

    if (hNRecoJets) delete hNRecoJets;
    if (hRecoJetPtRaw) delete hRecoJetPtRaw;
    if (hRecoJetPtCorrVsPtRaw) delete hRecoJetPtCorrVsPtRaw;
    if (hRecoJetPt) delete hRecoJetPt;
    if (hRecoJetPtWeighted) delete hRecoJetPtWeighted;
    if (hRecoJetEta) delete hRecoJetEta;
    if (hRecoJetPhi) delete hRecoJetPhi;
    if (hRecoJetPtVsEta) delete hRecoJetPtVsEta;
    if (hRecoJetPhiVsPt) delete hRecoJetPhiVsPt; 
    if (hRecoJetEtaVsPhi) delete hRecoJetEtaVsPhi;

    if (hRecoJet) delete hRecoJet;
    if (hRecoJetCorr) delete hRecoJetCorr;
    if (hRecoJetCorrWeighted) delete hRecoJetCorrWeighted;
    if (hRecoJetCorrEtaPhiPtHat) delete hRecoJetCorrEtaPhiPtHat;
    if (hRecoJetCorrEtaPhiPtHatWeighted) delete hRecoJetCorrEtaPhiPtHatWeighted;

    if (fIsMc) {
        if (hNRefJets) delete hNRefJets;
        if (hRefJetPt) delete hRefJetPt;
        if (hRefJetPtWeighted) delete hRefJetPtWeighted;
        if (hRecoJetEta) delete hRecoJetEta;
        if (hRecoJetPhi) delete hRecoJetPhi;
        if (hRecoJetPtVsEta) delete hRecoJetPtVsEta;

        if (hPtHatVsGenPtWOPtJetGtPtHatCut) delete hPtHatVsGenPtWOPtJetGtPtHatCut;
        if (hPtHatVsGenPtWeightedWoPtJetGtPtHatCut) delete hPtHatVsGenPtWeightedWoPtJetGtPtHatCut;
        if (hPtHatVsGenPtWithPtJetGtPtHatCut) delete hPtHatVsGenPtWithPtJetGtPtHatCut;
        if (hPtHatVsGenPtWeightedWithPtJetGtPtHatCut) delete hPtHatVsGenPtWeightedWithPtJetGtPtHatCut;

        if (hRefJet) delete hRefJet;
        if (hRefJetWeighted) delete hRefJetWeighted;
        if (hRefJetEtaPhiPtHat) delete hRefJetEtaPhiPtHat;
        if (hRefJetEtaPhiPtHatWeighted) delete hRefJetEtaPhiPtHatWeighted;

        if (hJESRaw) delete hJESRaw;
        if (hJESRawWeighted) delete hJESRawWeighted;
        if (hJESReco) delete hJESReco;
        if (hJESRecoWeighted) delete hJESRecoWeighted;
        if (hJESRecoEtaPhiPtHat) delete hJESRecoEtaPhiPtHat;
        if (hJESRecoEtaPhiPtHatWeighted) delete hJESRecoEtaPhiPtHatWeighted;

        if (hJERRaw) delete hJERRaw;
        if (hJERRawWeighted) delete hJERRawWeighted;
        if (hJERReco) delete hJERReco;
        if (hJERRecoWeighted) delete hJERRecoWeighted;
    }

    if (hRecoJetPtCorrVsGenPt) delete hRecoJetPtCorrVsGenPt;
}

//________________
void BasicHistoManager::init(const Bool_t& isMc) {
    
    Int_t    prRescale = 4;
    Int_t    bins4D_jet[4] = { fJetPtBins    , fJetEtaBins    , fJetPhiBins    , fCentBins     };
    Double_t xmin4D_jet[4] = { fJetPtRange[0], fJetEtaRange[0], fJetPhiRange[0], fCentRange[0] };
    Double_t xmax4D_jet[4] = { fJetPtRange[1], fJetEtaRange[1], fJetPhiRange[1], fCentRange[1] };
    
    Int_t    bins4D_jet_pthat[4] = { fJetPtBins    , fJetEtaBins    , fJetPhiBins    , fPtHatBins     };
    Double_t xmin4D_jet_pthat[4] = { fJetPtRange[0], fJetEtaRange[0], fJetPhiRange[0], fPtHatRange[0] };
    Double_t xmax4D_jet_pthat[4] = { fJetPtRange[1], fJetEtaRange[1], fJetPhiRange[1], fPtHatRange[1] };


    Int_t    bins4D_jes[4] = { fJESBins    , fJetPtBins / prRescale, 8, fCentBins  };
    Double_t xmin4D_jes[4] = { fJESRange[0], fJetPtRange[0], 0, fCentRange[0] };
    Double_t xmax4D_jes[4] = { fJESRange[1], fJetPtRange[1], 8, fCentRange[1] };

    Int_t    bins4D_jes_pthat[4] = { fJESBins    , fJetPtBins / prRescale, fJetEtaBins    , fPtHatBins     };
    Double_t xmin4D_jes_pthat[4] = { fJESRange[0], fJetPtRange[0], fJetEtaRange[0], fPtHatRange[0] };
    Double_t xmax4D_jes_pthat[4] = { fJESRange[1], fJetPtRange[1], fJetEtaRange[1], fPtHatRange[1] };

    Int_t    bins4D_jer[4] = { fJERBins    , fJetPtBins / prRescale, 8, fCentBins  };
    Double_t xmin4D_jer[4] = { fJERRange[0], fJetPtRange[0], 0, fCentRange[0] };
    Double_t xmax4D_jer[4] = { fJERRange[1], fJetPtRange[1], 8, fCentRange[1] };

    Int_t multBins = 1800;
    Double_t multRange[2] = {-0.5, 1799.5};
    Int_t hiBinBins = 203;
    Double_t hiBinRange[2] = {-1.5, 201.5};
    Int_t centralityBins = 101;
    Double_t centralityRange[2] = {-0.5, 100.5};
    Int_t weightBins = 110;
    Double_t weightRange[2] = {-0.05, 1.05};
    Int_t ptHatBins = 100;
    Double_t ptHatRange[2] = {0., 1000.};

    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", 320, -40., 40.);
    hVz->Sumw2();
    hVzWeighted = new TH1D("hVzWeighted","Vertex z position;vz (cm);Entries", 320, -40., 40.);
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

    hNRecoJets = new TH1D("hNRecoJets","Number of reconstructed jets", 16, -0.5, 15.5);
    hNRecoJets->Sumw2();
    hRecoJetPtRaw = new TH1D("hRecoJetPtRaw","Reco jet raw p_{T};p_{T}^{raw} (GeV/c)", 
                             fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoJetPtRaw->Sumw2();
    hRecoJetPtCorrVsPtRaw = new TH2F("hRecoJetPtCorrVsPtRaw","Reco jet corrected p_{T} vs. raw p_{T};p_{T}^{raw} (GeV/c);p_{T}^{corr} (GeV/c)", 
                                     fJetPtBins, fJetPtRange[0], fJetPtRange[1],
                                     fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoJetPtCorrVsPtRaw->Sumw2();
    hRecoJetPt = new TH1D("hRecoJetPt","Reco jet corrected p_{T};p_{T}^{corr} (GeV/c)", 
                             fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoJetPt->Sumw2();
    hRecoJetPtWeighted = new TH1D("hRecoJetPtWeighted","Reco jet corrected p_{T} weighted;p_{T}^{corr} (GeV/c)", 
                                  fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoJetPtWeighted->Sumw2();
    hRecoJetEta = new TH1D("hRecoJetEta","Reco jet #eta;#eta;Entries", 
                           fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hRecoJetEta->Sumw2();
    hRecoJetPhi = new TH1D("hRecoJetPhi","Reco jet #phi;#phi (rad);Entries", 
                           fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hRecoJetPhi->Sumw2();

    hRecoJet = new THnSparseD("hRecoJet","Reconstructed jet with raw p_{T};p_{T}^{raw} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
    hRecoJet->Sumw2();
    hRecoJetCorr = new THnSparseD("hRecoJetCorr","Reconstructed jet p_{T}^{corr};p_{T}^{corr} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
    hRecoJetCorr->Sumw2();
    hRecoJetCorrWeighted = new THnSparseD("hRecoJetCorrWeighted","Reconstructed jet p_{T}^{corr} weighted with #hat{p_{T}};p_{T}^{corr} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
    hRecoJetCorrWeighted->Sumw2();
    hRecoJetCorrEtaPhiPtHat = new THnSparseD("hRecoJetCorrEtaPhiPtHat","Reconstructed jet (corrected);p_{T}^{corr} (GeV/c);#eta; #phi (rad);#hat{p_{T}}", 4, bins4D_jet_pthat, xmin4D_jet_pthat, xmax4D_jet_pthat);
    hRecoJetCorrEtaPhiPtHat->Sumw2();
    hRecoJetCorrEtaPhiPtHatWeighted = new THnSparseD("hRecoJetCorrEtaPhiPtHatWeighted","Reconstructed jet (corrected) weighted;p_{T}^{corr}(weighted) (GeV/c);#eta; #phi (rad);#hat{p_{T}}", 4, bins4D_jet_pthat, xmin4D_jet_pthat, xmax4D_jet_pthat);
    hRecoJetCorrEtaPhiPtHatWeighted->Sumw2();

    hRecoJetPtCorrVsGenPt = new TH2D("hRecoJetPtCorrVsGenPt","Reconstructed jet p_{T} vs gen jet p_{T};p_{T}^{gen} (GeV/c);p_{T}^{corr} (GeV/c)", 
                                     fJetPtBins, fJetPtRange[0], fJetPtRange[1],
                                     fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoJetPtCorrVsGenPt->Sumw2();


    if (fIsMc) {

        hNRefJets = new TH1D("hNRefJets","Number of reconstructed jets", 16, -0.5, 15.5);
        hNRefJets->Sumw2();
        hRefJetPt = new TH1D("hRefJetPt","Ref jet corrected p_{T};p_{T} (GeV/c)", 
                              fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
        hRefJetPt->Sumw2();
        hRefJetPtWeighted = new TH1D("hRefJetPtWeighted","Ref jet  p_{T} weighted;p_{T}^{corr} (GeV/c)", 
                                      fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
        hRefJetPtWeighted->Sumw2();
        hRefJetEta = new TH1D("hRefJetEta","Ref jet #eta;#eta;Entries", 
                               fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
        hRefJetEta->Sumw2();
        hRefJetPhi = new TH1D("hRefJetPhi","Ref jet #phi;#phi (rad);Entries", 
                                fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
        hRefJetPhi->Sumw2();


        hRefJet = new THnSparseD("hRefJet","Generated jet p_{T};p_{T} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
        hRefJet->Sumw2();
        hRefJetWeighted = new THnSparseD("hRefJetWeighted","Generated jet p_{T} weighted;p_{T} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
        hRefJetWeighted->Sumw2();
        hRefJetEtaPhiPtHat = new THnSparseD("hRefJetEtaPhiPtHat","Generated jet p_{T};p_{T}^{gen} (GeV/c);#eta;#phi;#hat{p_{T}}", 4, bins4D_jet_pthat, xmin4D_jet_pthat, xmax4D_jet_pthat);
        hRefJetEtaPhiPtHat->Sumw2();
        hRefJetEtaPhiPtHatWeighted = new THnSparseD("hRefJetEtaPhiPtHatWeighted","Generated jet p_{T} weighted;p_{T}^{gen} (GeV/c);#eta;#phi;#hat{p_{T}}", 4, bins4D_jet_pthat, xmin4D_jet_pthat, xmax4D_jet_pthat);
        hRefJetEtaPhiPtHatWeighted->Sumw2();

        hJESRaw = new THnSparseD("hJESRaw", "hJESRaw;#frac{p_{T}^{raw}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        hJESRaw->Sumw2();
        hJESRawWeighted = new THnSparseD("hJESRawWeighted", "hJESRawWeighted;#frac{p_{T}^{raw}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        hJESRawWeighted->Sumw2();
        hJESReco = new THnSparseD("hJESReco", "hJESReco;#frac{p_{T}^{corr}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        hJESReco->Sumw2();
        hJESRecoWeighted = new THnSparseD("hJESRecoWeighted", "hJESRecoWeighted;#frac{p_{T}^{corr}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        hJESRecoWeighted->Sumw2();
        hJESRecoEtaPhiPtHat = new THnSparseD("hJESRecoEtaPhiPtHat", "JES for 0-10\% central collisions;#frac{p_{T}^{corr}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);#eta;#hat{p_{T}}", 4, bins4D_jes_pthat, xmin4D_jes_pthat, xmax4D_jes_pthat);
        hJESRecoEtaPhiPtHat->Sumw2();
        hJESRecoEtaPhiPtHatWeighted = new THnSparseD("hJESRecoEtaPhiPtHatWeighted", "JES for 0-10\% central collisions weighted;#frac{p_{T}^{corr}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);#eta;#hat{p_{T}}", 4, bins4D_jes_pthat, xmin4D_jes_pthat, xmax4D_jes_pthat);
        hJESRecoEtaPhiPtHatWeighted->Sumw2();

        hJERRaw = new THnSparseD("hJERRaw", "hJERRaw;JER^{raw};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        hJERRaw->Sumw2();
        hJERRawWeighted = new THnSparseD("hJERRawWeighted", "hJERRawWeighted;#frac{p_{T}^{raw}-p_{T}^{gen}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
        hJERRawWeighted->Sumw2();
        hJERReco = new THnSparseD("hJERReco", "hJERReco;#frac{p_{T}^{corr}-p_{T}^{gen}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
        hJERReco->Sumw2();
        hJERRecoWeighted = new THnSparseD("hJERRecoWeighted", "hJERRecoWeighted;#frac{p_{T}^{corr}-p_{T}^{gen}}{p_{T}^{gen}};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
        hJERRecoWeighted->Sumw2();
    }
}

//________________
void BasicHistoManager::writeOutput() {
    hVz->Write();
    hVzWeighted->Write();
    hMult->Write();
    hHiBin->Write();
    hHiBinWeighted->Write();
    hPtHat->Write();
    hPtHatWeighted->Write();
    hPtHatWeight->Write();
    hCentrality->Write();
    hNRecoJets->Write();
    hRecoJetPtRaw->Write();
    hRecoJetPtCorrVsPtRaw->Write();
    hRecoJetPt->Write();
    hRecoJetPtWeighted->Write();
    hRecoJetEta->Write();
    hRecoJetPhi->Write();
    hRecoJet->Write();
    hRecoJetCorr->Write();
    hRecoJetCorrWeighted->Write();
    hRecoJetCorrEtaPhiPtHat->Write();
    hRecoJetCorrEtaPhiPtHatWeighted->Write();
    hRecoJetPtCorrVsGenPt->Write();
    if (fIsMc) {
        hNRefJets->Write();
        hRefJetPt->Write();
        hRefJetPtWeighted->Write();
        hRefJetEta->Write();
        hRefJetPhi->Write();
        hRefJet->Write();
        hRefJetWeighted->Write();
        hRefJetEtaPhiPtHat->Write();
        hRefJetEtaPhiPtHatWeighted->Write();
        hJESRaw->Write();
        hJESRawWeighted->Write();
        hJESReco->Write();
        hJESRecoWeighted->Write();
        hJESRecoEtaPhiPtHat->Write();
        hJESRecoEtaPhiPtHatWeighted->Write();
        hJERRaw->Write();
        hJERRawWeighted->Write();
        hJERReco->Write();
        hJERRecoWeighted->Write();
    }
}