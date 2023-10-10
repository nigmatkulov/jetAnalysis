// Jet analysis headers
#include "BasicHistoManager.h"

// ROOT headers
#include "TObject.h"
#include "TIterator.h"
#include "TString.h"
#include "TClass.h"
#include "TKey.h"
#include "TROOT.h"

ClassImp(BasicHistoManager)

//________________
BasicHistoManager::BasicHistoManager() :
  fListOfObjects{nullptr}, fIsMc{kFALSE}, 
  fCentBins{6}, fCentRange{0., 60.},
  fJetPtBins{200}, fJetPtRange{0., 1000.}, 
  fJetEtaBins{50}, fJetEtaRange{-2.5, 2.5},
  fJetPhiBins{64}, fJetPhiRange{-TMath::Pi(), TMath::Pi()},
  fJESBins{500}, fJESRange{0., 5.},
  fJERBins{400}, fJERRange{-2., 2.},
  hVz{nullptr}, hVzWeighted{nullptr}, hMult{nullptr},
  hHiBin{nullptr}, hWeight{nullptr}, hPtHat{nullptr},
  hPtHatWeighted{nullptr}, hCentrality{nullptr}, 
  hNRecoJets{nullptr}, hRecoJetPtRaw{nullptr}, hRecoJetPtCorrVsPtRaw{nullptr},
  hRecoJetPt{nullptr}, hRecoJetPtWeighted{nullptr}, hRecoJetEta{nullptr},
  hRecoJetPhi{nullptr}, hRecoJetPtVsEta{nullptr}, 
  hRecoJetPhiVsPt{nullptr}, hRecoJetEtaVsPhi{nullptr},
  hRecoJet{nullptr}, hRecoJetCorr{nullptr}, hRecoJetCorrWeighted{nullptr},
  hNRefJets{nullptr}, hRefJetPt{nullptr}, hRefJetPtWeighted{nullptr},
  hRefJetEta{nullptr}, hRefJetPhi{nullptr}, hRefJetPtVsEta{nullptr},
  hRefJet{nullptr}, hRefJetWeighted{nullptr},
  hJESRaw{nullptr}, hJESRawWeighted{nullptr}, hJESReco{nullptr}, hJESRecoWeighted{nullptr},
  hJERRaw{nullptr}, hJERRawWeighted{nullptr}, hJERReco{nullptr}, hJERRecoWeighted{nullptr} { 
    fListOfObjects = new TList();
}

//________________
BasicHistoManager::~BasicHistoManager() {
    if (hVz) delete hVz;
    if (hVzWeighted) delete hVzWeighted;
    if (hMult) delete hMult;
    if (hHiBin) delete hHiBin;
    if (hWeight) delete hWeight;
    if (hPtHat) delete hPtHat;
    if (hPtHatWeighted) delete hPtHatWeighted;
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

    if (fIsMc) {
        if (hNRefJets) delete hNRefJets;
        if (hRefJetPt) delete hRefJetPt;
        if (hRefJetPtWeighted) delete hRefJetPtWeighted;
        if (hRecoJetEta) delete hRecoJetEta;
        if (hRecoJetPhi) delete hRecoJetPhi;
        if (hRecoJetPtVsEta) delete hRecoJetPtVsEta;

        if (hRefJet) delete hRefJet;
        if (hRefJetWeighted) delete hRefJetWeighted;

        if (hJESRaw) delete hJESRaw;
        if (hJESRawWeighted) delete hJESRawWeighted;
        if (hJESReco) delete hJESReco;
        if (hJESRecoWeighted) delete hJESRecoWeighted;

        if (hJERRaw) delete hJERRaw;
        if (hJERRawWeighted) delete hJERRawWeighted;
        if (hJERReco) delete hJERReco;
        if (hJERRecoWeighted) delete hJERRecoWeighted;
    }
}

//________________
void BasicHistoManager::init(const Bool_t& isMc) {
    
    Int_t    bins4D_jet[4] = { fJetPtBins    , fJetEtaBins    , fJetPhiBins    , fCentBins     };
    Double_t xmin4D_jet[4] = { fJetPtRange[0], fJetEtaRange[0], fJetPhiRange[0], fCentRange[0] };
    Double_t xmax4D_jet[4] = { fJetPtRange[1], fJetEtaRange[1], fJetPhiRange[1], fCentRange[1] };

    Int_t    bins4D_jes[4] = { fJESBins    , fJetPtBins    , 8, fCentBins     };
    Double_t xmin4D_jes[4] = { fJESRange[0], fJetPtRange[0], 0, fCentRange[0] };
    Double_t xmax4D_jes[4] = { fJESRange[1], fJetPtRange[1], 8, fCentRange[1] };

    Int_t    bins4D_jer[4] = { fJERBins    , fJetPtBins    , 8, fCentBins     };
    Double_t xmin4D_jer[4] = { fJERRange[0], fJetPtRange[0], 0, fCentRange[0] };
    Double_t xmax4D_jer[4] = { fJERRange[1], fJetPtRange[1], 8, fCentRange[1] };

    Int_t multBins = 1800;
    Double_t multRange[2] = {-0.5, 1799.5};
    Int_t hiBinBins = 203;
    Double_t hiBinRange[2] = {-1.5, 201.5};
    Int_t centralityBins = 101;
    Double_t centralityRange[2] = {-0.5, 100.5};
    Int_t weightBins = 220;
    Double_t weightRange[2] = {-1.1, 1.1};
    Int_t ptHatBins = 100;
    Double_t ptHatRange[2] = {0., 1000.};

    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", 320, -40., 40.);
    fListOfObjects->Add(hVz);
    hVzWeighted = new TH1D("hVzWeighted","Vertex z position;vz (cm);Entries", 320, -40., 40.);
    fListOfObjects->Add(hVzWeighted);
    hMult = new TH1D("hMult","Charged particle multiplicity;Multiplicity;Entries", 
                     multBins, multRange[0], multRange[1]);
    fListOfObjects->Add(hMult);
    hHiBin = new TH1D("hHiBin","HiBin a.k.a. centrality;HiBin;Entries", 
                      hiBinBins, hiBinRange[0], hiBinRange[1]);
    fListOfObjects->Add(hHiBin);
    hWeight = new TH1D("hWeight","Event weight;weight;Entries", 
                       weightBins, weightRange[0], weightRange[1]);
    fListOfObjects->Add(hWeight);
    hPtHat = new TH1D("hPtHat","#hat{p_{T}};#hat{p_{T}} (GeV/c);Entries", 
                      ptHatBins, ptHatRange[0], ptHatRange[1]);
    fListOfObjects->Add(hPtHat);
    hPtHatWeighted = new TH1D("hPtHatWeighted","#hat{p_{T}};#hat{p_{T}} (GeV/c);Entries", 
                              ptHatBins, ptHatRange[0], ptHatRange[1]);
    fListOfObjects->Add(hPtHatWeighted);
    hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries",
                           centralityBins, centralityRange[0], centralityRange[1]);
    fListOfObjects->Add(hCentrality);

    hNRecoJets = new TH1D("hNRecoJets","Number of reconstructed jets", 16, -0.5, 15.5);
    fListOfObjects->Add(hNRecoJets);
    hRecoJetPtRaw = new TH1D("hRecoJetPtRaw","Reco jet raw p_{T};p_{T}^{raw} (GeV/c)", 
                             fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    fListOfObjects->Add(hRecoJetPtRaw);
    hRecoJetPtCorrVsPtRaw = new TH2F("hRecoJetPtCorrVsPtRaw","Reco jet corrected p_{T} vs. raw p_{T};p_{T}^{raw} (GeV/c);p_{T}^{corr} (GeV/c)", 
                                     fJetPtBins, fJetPtRange[0], fJetPtRange[1],
                                     fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    fListOfObjects->Add(hRecoJetPtCorrVsPtRaw);
    hRecoJetPt = new TH1D("hRecoJetPt","Reco jet corrected p_{T};p_{T}^{corr} (GeV/c)", 
                             fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    fListOfObjects->Add(hRecoJetPt);
    hRecoJetPtWeighted = new TH1D("hRecoJetPtWeighted","Reco jet corrected p_{T} weighted;p_{T}^{corr} (GeV/c)", 
                                  fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    fListOfObjects->Add(hRecoJetPtWeighted);
    hRecoJetEta = new TH1D("hRecoJetEta","Reco jet #eta;#eta;Entries", 
                           fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    fListOfObjects->Add(hRecoJetEta);
    hRecoJetPhi = new TH1D("hRecoJetPhi","Reco jet #phi;#phi (rad);Entries", 
                           fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    fListOfObjects->Add(hRecoJetPhi);

    hRecoJet = new THnSparseD("hRecoJet","Reconstructed jet with raw p_{T};p_{T}^{raw} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
    fListOfObjects->Add(hRecoJet);
    hRecoJetCorr = new THnSparseD("hRecoJetCorr","Reconstructed jet p_{T}^{corr};p_{T}^{corr} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
    fListOfObjects->Add(hRecoJetCorr);
    hRecoJetCorrWeighted = new THnSparseD("hRecoJetCorrWeighted","Reconstructed jet p_{T}^{corr} weighted;p_{T}^{corr} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
    fListOfObjects->Add(hRecoJetCorrWeighted);


    if (fIsMc) {

        hNRefJets = new TH1D("hNRefJets","Number of reconstructed jets", 16, -0.5, 15.5);
        fListOfObjects->Add(hNRefJets);
        hRefJetPt = new TH1D("hRefJetPt","Ref jet corrected p_{T};p_{T} (GeV/c)", 
                              fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
        fListOfObjects->Add(hRefJetPt);
        hRefJetPtWeighted = new TH1D("hRefJetPtWeighted","Ref jet  p_{T} weighted;p_{T}^{corr} (GeV/c)", 
                                      fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
        fListOfObjects->Add(hRefJetPtWeighted);
        hRefJetEta = new TH1D("hRefJetEta","Ref jet #eta;#eta;Entries", 
                               fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
        fListOfObjects->Add(hRefJetEta);
        hRefJetPhi = new TH1D("hRefJetPhi","Ref jet #phi;#phi (rad);Entries", 
                                fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
        fListOfObjects->Add(hRefJetPhi);


        hRefJet = new THnSparseD("hRefJet","Generated jet p_{T};p_{T} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
        fListOfObjects->Add(hRefJet);
        hRefJetWeighted = new THnSparseD("hRefJetWeighted","Generated jet p_{T} weighted;p_{T} (GeV/c);#eta;#phi;centrality", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
        fListOfObjects->Add(hRefJetWeighted);


        hJESRaw = new THnSparseD("hJESRaw", "hJESRaw;JES^{raw};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJESRaw);
        hJESRawWeighted = new THnSparseD("hJESRawWeighted", "hJESRawWeighted;JES^{raw};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJESRawWeighted);
        hJESReco = new THnSparseD("hJESReco", "hJESReco;JES;p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJESReco);
        hJESRecoWeighted = new THnSparseD("hJESRecoWeighted", "hJESRecoWeighted;JES;p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJESRecoWeighted);

        hJERRaw = new THnSparseD("hJERRaw", "hJERRaw;JER^{raw};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJERRaw);
        hJERRawWeighted = new THnSparseD("hJERRawWeighted", "hJERRawWeighted;JER^{raw};p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJERRawWeighted);
        hJERReco = new THnSparseD("hJERReco", "hJERReco;JER;p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJERReco);
        hJERRecoWeighted = new THnSparseD("hJERRecoWeighted", "hJERRecoWeighted;JER;p_{T}^{gen} (GeV/c);FlavorFromB;centrality", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
        fListOfObjects->Add(hJERRecoWeighted);
    }

    // Iterate over a list and reweight histograms
    TIter next(fListOfObjects);
    TObject *obj{nullptr};
    TKey *key{nullptr};
    while((obj = next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (cl->InheritsFrom("TH1") ) {
        TH1 *h = (TH1*)key->ReadObj();
        h->Sumw2();
      }
      else if (cl->InheritsFrom("TH2") ) {
        TH2 *h = (TH2*)key->ReadObj();
        h->Sumw2();
      }
      else if (cl->InheritsFrom("TH3") ) {
        TH3 *h = (TH3*)key->ReadObj();
        h->Sumw2();
      }
      else if (cl->InheritsFrom("THnSparse") ) {
        THnSparse *h = (THnSparse*)key->ReadObj();
        h->Sumw2();
      }
    }
    delete obj;

    // TH1::Sumw2();
    // TH2::Sumw2();
    // THnSparse::Sumw2();
}

//________________
void BasicHistoManager::writeOutput() {
    // Iterate over a list
    TIter next(fListOfObjects);
    TObject *obj{nullptr};
    while((obj = next())) {
        obj->Write();
    }
    delete obj;
}