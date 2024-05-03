/**
 * @file HistoManagerDiJet.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for dijet studies
 * @version 0.1
 * @date 2023-10-24
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "HistoManagerDiJet.h"

// ROOT headers
#include "TObject.h"
#include "TIterator.h"
#include "TString.h"
#include "TClass.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

ClassImp(HistoManagerDiJet)

//________________
HistoManagerDiJet::HistoManagerDiJet() :
  fIsMc{kFALSE}, 
//   fCentBins{10}, fCentRange{-10., 90.},
  fPtBins{50}, fPtRange{20., 520.}, 
  fEtaBins{50}, fEtaRange{-5.0, 5.0},
  fPhiBins{16}, fPhiRange{-TMath::Pi(), TMath::Pi()},
  fDijetPtBins{120}, fDijetPtRange{20., 620.},
  fDijetEtaBins{50}, fDijetEtaRange{-5., 5.},
  fDijetDphiBins{16}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},
  fPtHatBins{60}, fPtHatRange{15., 615.},
  
  hVz{nullptr}, hVzWeighted{nullptr}, hMult{nullptr},
  hHiBin{nullptr}, hHiBinWeighted{nullptr},
  hPtHat{nullptr}, hPtHatWeighted{nullptr}, hPtHatWeight{nullptr},
//   hCentrality{nullptr}, hCentralityWeighted{nullptr},
  hVzPtHat{nullptr}, hVzPtHatWeighted{nullptr},

  // Gen jets
  hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi{nullptr},
  hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted{nullptr},
  hGenInclusiveJetPt{nullptr},
  hGenInclusiveJetPtEta{nullptr},
  hGenPtLeadPtSublead{nullptr},
  hGenEtaLeadEtaSublead{nullptr},
  hGenDijetEta{nullptr},
  hGenDijetPtEtaDphi{nullptr},

  // Reco jets (single)
  hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
  hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
  hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
  hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
  hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
  hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
  hJESInclusiveJetPtEtaPhiPtHat{nullptr},
  hJESInclusiveJetPtEtaPhiPtHatWeighted{nullptr},
  hRecoMatchedPtEta{nullptr},

  // Dijets (experiment)
  hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi{nullptr},
  hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted{nullptr},
  hRecoInclusiveJetPt{nullptr},
  hRecoPtLeadPtSublead{nullptr},
  hRecoEtaLeadEtaSublead{nullptr},
  hRecoDijetPtEtaDphi{nullptr},
  hRecoInclusiveAllJetPtVsEta{nullptr},
  hRecoInclusiveMatchedJetPtVsEta{nullptr},

  // Dijets (MC)
  hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta{nullptr},
  hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},
  hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},
  // hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat{nullptr},
  // hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted{nullptr},
  hRecoDijetEta{nullptr},
  hRefDijetEta{nullptr},
  hRefDijetEtaVsRecoDijetEta{nullptr},
  hRefDijetPtEtaDphi{nullptr},
  hRefSelDijetEta{nullptr},
  hRefSelDijetPtEtaDphi{nullptr},

  hRefInclusiveJetPt{nullptr},
  hRefInclusiveJetPtEta{nullptr},
  hRefPtLeadPtSublead{nullptr},
  hRefEtaLeadEtaSublead{nullptr}
{ 
    /* Empty */
}

//________________
HistoManagerDiJet::~HistoManagerDiJet() {
    if (hVz)            delete hVz;
    if (hVzWeighted)    delete hVzWeighted;
    if (hMult)          delete hMult;
    if (hHiBin)         delete hHiBin;
    if (hHiBinWeighted) delete hHiBinWeighted;
    if (hPtHat)         delete hPtHat;
    if (hPtHatWeighted) delete hPtHatWeighted;
    if (hPtHatWeight)   delete hPtHatWeight;
    // if (hCentrality)    delete hCentrality;
    // if (hCentralityWeighted) delete hCentralityWeighted;
    if (hVzPtHat) delete hVzPtHat;
    if (hVzPtHatWeighted) delete hVzPtHatWeighted;

    if (fIsMc) {
      // Gen jets
      if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
      if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
      if (hGenInclusiveJetPt) delete hGenInclusiveJetPt;
      if (hGenInclusiveJetPtEta) delete hGenInclusiveJetPtEta;
      if (hGenPtLeadPtSublead) delete hGenPtLeadPtSublead;
      if (hGenEtaLeadEtaSublead) delete hGenEtaLeadEtaSublead;
      if (hGenDijetEta) delete hGenDijetEta;
      if (hGenDijetPtEtaDphi) delete hGenDijetPtEtaDphi;


      if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
      if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
      if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
      if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
      if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
      if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
      if (hJESInclusiveJetPtEtaPhiPtHat) delete hJESInclusiveJetPtEtaPhiPtHat;
      if (hJESInclusiveJetPtEtaPhiPtHatWeighted) delete hJESInclusiveJetPtEtaPhiPtHatWeighted;
      if (hRecoMatchedPtEta) delete hRecoMatchedPtEta;
      if (hRecoInclusiveAllJetPtVsEta) delete hRecoInclusiveAllJetPtVsEta;
      if (hRecoInclusiveMatchedJetPtVsEta) delete hRecoInclusiveMatchedJetPtVsEta;
    } // if (fIsMc)

    //
    // Dijets (experiment)
    //
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    if (hRecoInclusiveJetPt) delete hRecoInclusiveJetPt;
    if (hRecoPtLeadPtSublead) delete hRecoPtLeadPtSublead;
    if (hRecoEtaLeadEtaSublead) delete hRecoEtaLeadEtaSublead;
    if (hRecoDijetEta) delete hRecoDijetEta;
    if (hRecoDijetPtEtaDphi) delete hRecoDijetPtEtaDphi;

    //
    // Dijets exp vs mc
    //
    if ( fIsMc ) {
      if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta;
      if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;
      if (hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;
      // if (hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat) delete hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat;
      // if (hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted) delete hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted;
      if (hRefInclusiveJetPt) delete hRefInclusiveJetPt;
      if (hRefInclusiveJetPtEta) delete hRefInclusiveJetPtEta;
      if (hRefPtLeadPtSublead) delete hRefPtLeadPtSublead;
      if (hRefEtaLeadEtaSublead) delete hRefEtaLeadEtaSublead;
      if (hRefDijetPtEtaDphi) delete hRefDijetPtEtaDphi;
      if (hRefSelDijetEta) delete hRefSelDijetEta;
      if (hRefSelDijetPtEtaDphi) delete hRefSelDijetPtEtaDphi;

      if (hRefDijetEta) delete hRefDijetEta;
      if (hRefDijetEtaVsRecoDijetEta) delete hRefDijetEtaVsRecoDijetEta;

    } // if ( fIsMc )
}

//________________
void HistoManagerDiJet::init(const Bool_t& isMc) {
    
    Int_t    prescale = 2;

    Int_t    vzBins = 320;
    Double_t vzRange[2] {-31., 31.};
    Int_t    multBins{1800};
    Double_t multRange[2] {-0.5, 1799.5};
    Int_t    hiBinBins{203};
    Double_t hiBinRange[2] {-1.5, 201.5};
    // Int_t    centralityBins{101};
    // Double_t centralityRange[2] {-0.5, 100.5};
    Int_t    weightBins{110};
    Double_t weightRange[2] {-0.05, 1.05};
    Int_t    ptHatBins{100};
    Double_t ptHatRange[2] {0., 1000.};
    Int_t    fJESBins{100}; 
    Double_t fJESRange[2] {0., 2.};


    Int_t    bins2D_ev_VzPtHat[2] {     vzBins,     fPtHatBins };
    Double_t xmin2D_ev_VzPtHat[2] { vzRange[0], fPtHatRange[0] };
    Double_t xmax2D_ev_VzPtHat[2] { vzRange[1], fPtHatRange[1] };

    //
    // Gen
    //
    Int_t    bins9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi[9]
    { fDijetPtBins, fDijetEtaBins, fDijetDphiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    Double_t xmin9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi[9]
    { fDijetPtRange[0], fDijetEtaRange[0], fDijetDphiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0]  };
    Double_t xmax9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi[9]
    { fDijetPtRange[1], fDijetEtaRange[1], fDijetDphiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    //
    // Reco single
    //
    Int_t    bins5D_jet_PtPtPtEtaEta[5] { fPtBins, fPtBins, fPtBins, fEtaBins, fEtaBins };
    Double_t xmin5D_jet_PtPtPtEtaEta[5] { fPtRange[0], fPtRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0] };
    Double_t xmax5D_jet_PtPtPtEtaEta[5] { fPtRange[1], fPtRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1] };

    Int_t    bins5D_jet_JESPtEtaPhiPtHat[5] { fJESBins, fPtBins, fEtaBins, fPhiBins, fPtHatBins };
    Double_t xmin5D_jet_JESPtEtaPhiPtHat[5] { fJESRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtHatRange[0] };
    Double_t xmax5D_jet_JESPtEtaPhiPtHat[5] { fJESRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtHatRange[1] };

    //
    // Reco dijets
    //
    Int_t    bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi[9]
    { fDijetPtBins, fDijetEtaBins, fDijetDphiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    Double_t xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi[9]
    { fDijetPtRange[0], fDijetEtaRange[0], fDijetDphiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0] };
    Double_t xmax9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi[9]
    { fDijetPtRange[1], fDijetEtaRange[1], fDijetDphiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    Int_t    bins12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta[12]
    { fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins, fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins };
    Double_t xmin12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta[12]
    { fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0] };
    Double_t xmax12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta[12]
    { fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1] };

    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
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
    // hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries",
    //                        centralityBins, centralityRange[0], centralityRange[1]);
    // hCentrality->Sumw2();
    // hCentralityWeighted = new TH1D("hCentralityWeighted","Collision centrality weighted;Centrality (%);Entries",
    //                                centralityBins, centralityRange[0], centralityRange[1]);
    // hCentralityWeighted->Sumw2();

    hVzPtHat = new THnSparseD( "hVzPtHat","Vertex and #hat{p_{T}};vz (cm);#hat{p_{T}} (GeV/c)",
                               2, bins2D_ev_VzPtHat, xmin2D_ev_VzPtHat, xmax2D_ev_VzPtHat );
    hVzPtHat->Sumw2();
    hVzPtHatWeighted = new THnSparseD( "hVzPtHatWeighted","Vertex and #hat{p_{T}} weighted;vz (cm);#hat{p_{T}} (GeV/c)",
                                       2, bins2D_ev_VzPtHat, xmin2D_ev_VzPtHat, xmax2D_ev_VzPtHat );
    hVzPtHatWeighted->Sumw2();


    //
    // Monte Carlo information
    //
    if (fIsMc) {

      //
      // Gen jets
      //

      hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi = new THnSparseD("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi","Title;Gen p_{T}^{dijet} (GeV/c);Gen #eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV/c);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV/c);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
            9, 
            bins9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
            xmin9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
            xmax9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi);
      hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Sumw2();
      hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted = new THnSparseD("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted","Title;Gen p_{T}^{dijet} (GeV/c);Gen #eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV/c);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV/c);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
            9, 
            bins9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
            xmin9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
            xmax9D_gen_GenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi);
      hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Sumw2();

      hGenInclusiveJetPt = new TH1D("hGenInclusiveJetPt","Inclusive gen jet;Gen p_{T}^{inclusive} (GeV/c)",
                                    fPtBins, fPtRange[0], fPtRange[1] );
      hGenInclusiveJetPt->Sumw2();
      hGenInclusiveJetPtEta = new TH2D("hGenInclusiveJetPtEta","Gen inclusive jet acceptance;Gen #eta;Gen p_{T} (GeV/c)",
                                      fEtaBins, fEtaRange[0], fEtaRange[1],
                                      fPtBins, fPtRange[0], fPtRange[1] );
      hGenInclusiveJetPtEta->Sumw2();
      hGenPtLeadPtSublead = new TH2D("hGenPtLeadPtSublead","Leading gen jet pT vs subleading gen jet pT;Gen p_{T}^{Leading} (GeV/c);Gen p_{T}^{Subleading} (GeV/c)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1] );
      hGenPtLeadPtSublead->Sumw2();
      hGenEtaLeadEtaSublead = new TH2D("hGenEtaLeadEtaSublead","Leading gen jet eta vs subleading gen jet eta;Gen #eta^{Leading};Gen #eta^{Subleading}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1] );
      hGenEtaLeadEtaSublead->Sumw2();
      hGenDijetPtEtaDphi = new TH3D("hGenDijetPtEtaDphi","Gen dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                    fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                    fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
      hGenDijetPtEtaDphi->Sumw2();
      hGenDijetEta = new TH1D("hGenDijetEta", "Gen dijet #eta;#eta^{dijet}",
                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
      hGenDijetEta->Sumw2();

      //
      // Reco single jets
      //
      hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed inclusive jets;Reco p_{T, corr}^{Inclusive} (GeV/c);Reco p_{T, raw}^{Inclusive} (GeV/c);Ref p_{T}^{Inclusive} (GeV/c);Reco #eta^{Inclusive};Ref #eta^{Inclusive}",
            5,
            bins5D_jet_PtPtPtEtaEta,
            xmin5D_jet_PtPtPtEtaEta,
            xmax5D_jet_PtPtPtEtaEta);
      hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
      hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted = new THnSparseD("hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted","Reconstructed inclusive jets weighted;Reco p_{T, corr}^{Inclusive} (GeV/c);Reco p_{T, raw}^{Inclusive} (GeV/c);Ref p_{T}^{Inclusive} (GeV/c);Reco #eta^{Inclusive};Ref #eta^{Inclusive}",
            5,
            bins5D_jet_PtPtPtEtaEta,
            xmin5D_jet_PtPtPtEtaEta,
            xmax5D_jet_PtPtPtEtaEta);
      hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Sumw2();
      hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed leading jets;Reco p_{T, corr}^{Leading} (GeV/c);Reco p_{T, raw}^{Leading} (GeV/c);Ref p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Ref #eta^{Leading}",
            5,
            bins5D_jet_PtPtPtEtaEta,
            xmin5D_jet_PtPtPtEtaEta,
            xmax5D_jet_PtPtPtEtaEta);
      hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
      hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted = new THnSparseD("hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted","Reconstructed leading jets weighted;Reco p_{T, corr}^{Leading} (GeV/c);Reco p_{T, raw}^{Leading} (GeV/c);Ref p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Ref #eta^{Leading}",
            5,
            bins5D_jet_PtPtPtEtaEta,
            xmin5D_jet_PtPtPtEtaEta,
            xmax5D_jet_PtPtPtEtaEta);
      hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Sumw2();
      hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed subleading jets;Reco p_{T, corr}^{Subleading} (GeV/c);Reco p_{T, raw}^{Subleading} (GeV/c);Ref p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref #eta^{Subleading}",
            5,
            bins5D_jet_PtPtPtEtaEta,
            xmin5D_jet_PtPtPtEtaEta,
            xmax5D_jet_PtPtPtEtaEta);
      hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
      hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted = new THnSparseD("hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted","Reconstructed subleading jets weighted;Reco p_{T, corr}^{Subleading} (GeV/c);Reco p_{T, raw}^{Subleading} (GeV/c);Ref p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref #eta^{Subleading}",
            5,
            bins5D_jet_PtPtPtEtaEta,
            xmin5D_jet_PtPtPtEtaEta,
            xmax5D_jet_PtPtPtEtaEta);
      hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Sumw2();

      hJESInclusiveJetPtEtaPhiPtHat = new THnSparseD("hJESInclusiveJetPtEtaPhiPtHat","JES of inclusive jets;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV/c);#eta^{reco};#phi^{reco} (rad);#hat{p_{T}} (GeV/c)",
            5,
            bins5D_jet_JESPtEtaPhiPtHat,
            xmin5D_jet_JESPtEtaPhiPtHat,
            xmax5D_jet_JESPtEtaPhiPtHat);
      hJESInclusiveJetPtEtaPhiPtHat->Sumw2();
      hJESInclusiveJetPtEtaPhiPtHatWeighted = new THnSparseD("hJESInclusiveJetPtEtaPhiPtHatWeighted","JES of inclusive jet weighted;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV/c);#eta^{reco};#phi^{reco} (rad);#hat{p_{T}} (GeV/c)",
            5,
            bins5D_jet_JESPtEtaPhiPtHat,
            xmin5D_jet_JESPtEtaPhiPtHat,
            xmax5D_jet_JESPtEtaPhiPtHat);
      hJESInclusiveJetPtEtaPhiPtHatWeighted->Sumw2();

      hRecoMatchedPtEta = new TH2D("hRecoMatchedPtEta","Reconstructed jets that matched gen;#eta^{reco};p_{T}^{reco} (GeV/c)",
            fEtaBins, fEtaRange[0], fEtaRange[1],
            fPtBins, fPtRange[0], fPtRange[1]);
      hRecoMatchedPtEta->Sumw2();
      hRecoInclusiveMatchedJetPtVsEta = new TH2D("hRecoInclusiveMatchedJetPtVsEta","Inclusive reco jet that has matching pT vs eta;#eta;p_{T} (GeV/c)",
                                                 fEtaBins, fEtaRange[0], fEtaRange[1],
                                                 30, 5., 155.);
      hRecoInclusiveMatchedJetPtVsEta->Sumw2();

      //
      // Reco dijet with MC
      //
      hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta = new THnSparseD("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta",
            "Reco to ref correspondence;Reco p_{T}^{dijet} (GeV/c);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV/c);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV/c);Ref #eta^{Subleading}",
            12,
            bins12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta,
            xmin12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta,
            xmax12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta);
      hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Sumw2();
      hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted = new THnSparseD("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted",
            "Reco to ref correspondence weighted;Reco p_{T}^{dijet} (GeV/c);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV/c);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV/c);Ref #eta^{Subleading}",
            12,
            bins12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta,
            xmin12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta,
            xmax12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta);
      hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Sumw2();

      hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted = new THnSparseD("hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted",
            "Reco to ref correspondence (via ref selection) weighted;Reco p_{T}^{dijet} (GeV/c);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV/c);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV/c);Ref #eta^{Subleading}",
            12,
            bins12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta,
            xmin12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta,
            xmax12D_dijet_PtEtPtEtaPtEtaPtEtaPtEtaPtEta);
      hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Sumw2();
      hRefSelDijetEta = new TH1D("hRefSelDijetEta","Ref selected dijets;#eta^{dijet}",
                                 fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
      hRefSelDijetEta->Sumw2();

      hRefInclusiveJetPt = new TH1D("hRefInclusiveJetPt","Ref jet p_{T};Ref p_{T} (GeV/c);Entries",
                                    fPtBins, fPtRange[0], fPtRange[1]);
      hRefInclusiveJetPt->Sumw2();
      hRefInclusiveJetPtEta = new TH2D("hRefInclusiveJetPtEta","Ref jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV/c)",
                                    fEtaBins, fEtaRange[0], fEtaRange[1],
                                    fPtBins, fPtRange[0], fPtRange[1]);
      hRefPtLeadPtSublead = new TH2D("hRefPtLeadPtSublead","Ref leading vs subleading p_{T};Ref p_{T}^{Leading} (GeV/c);Ref p_{T}^{Subleading} (GeV/c)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
      hRefPtLeadPtSublead->Sumw2();
      hRefEtaLeadEtaSublead = new TH2D("hRefEtaLeadEtaSublead","Ref leading vs subleading #eta;Ref #eta^{Leading};Ref #eta^{Subleading}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
      hRefEtaLeadEtaSublead->Sumw2();
      hRefDijetEta = new TH1D("hRefDijetEta","Ref dijet #eta;Ref #eta^{dijet};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
      hRefDijetEta->Sumw2();
      hRefDijetEtaVsRecoDijetEta = new TH2D("hRefDijetEtaVsRecoDijetEta","Ref dijet #eta vs reco dijet #eta;Reco #eta^{dijet};Ref #eta^{dijet}",
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
      hRefDijetEtaVsRecoDijetEta->Sumw2();
      hRefDijetPtEtaDphi = new TH3D("hRefDijetPtEtaDphi","Ref dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                    fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                    fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
      hRefDijetPtEtaDphi->Sumw2();
      hRefSelDijetPtEtaDphi = new TH3D("hRefSelDijetPtEtaDphi","RefSel dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                       fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                       fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                       fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
      hRefSelDijetPtEtaDphi->Sumw2();
    } // if (fIsMc)

    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi = new THnSparseD("hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi",
            "Reconstructed dijet and jet info;p_{T}^{dijet} (GeV/c);#eta^{dijet};#Delta#phi^{dijet} (rad);p_{T}^{Leading} (GeV/c);#eta^{Leading};#phi^{Leading} (rad);p_{T}^{Subleading} (GeV/c);#eta^{Subleading};#phi^{Subleading} (rad)",
            9,
            bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Sumw2();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted = new THnSparseD("hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted",
            "Reconstructed dijet and jet info weighted;p_{T}^{dijet} (GeV/c);#eta^{dijet};#Delta#phi^{dijet} (rad);p_{T}^{Leading} (GeV/c);#eta^{Leading};#phi^{Leading} (rad);p_{T}^{Subleading} (GeV/c);#eta^{Subleading};#phi^{Subleading} (rad)",
            9,
            bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhiPtEtaPhi);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Sumw2();
    
    hRecoInclusiveJetPt = new TH1D("hRecoInclusiveJetPt","Reco jet p_{T};Reco p_{T} (GeV/c);Entries",
                                   fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveJetPt->Sumw2();
    hRecoInclusiveAllJetPtVsEta = new TH2D("hRecoInclusiveAllJetPtVsEta", "Inclusive reco jet pT vs eta;#eta;p^{T} (GeV/c)",
                                           fEtaBins, fEtaRange[0], fEtaRange[1],
                                           30, 5., 155.);
    hRecoInclusiveAllJetPtVsEta->Sumw2();
    hRecoPtLeadPtSublead = new TH2D("hRecoPtLeadPtSublead","Reco leading vs subleading p_{T};Reco p_{T}^{Leading} (GeV/c);Reco p_{T}^{Subleading} (GeV/c)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSublead->Sumw2();
    hRecoEtaLeadEtaSublead = new TH2D("hRecoEtaLeadEtaSublead","Reco leading vs subleading #eta;Reco #eta^{Leading};Reco #eta^{Subleading}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSublead->Sumw2();
    hRecoDijetEta = new TH1D("hRecoDijetEta","Reco dijet #eta;Reco #eta^{dijet};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEta->Sumw2();
    hRecoDijetPtEtaDphi = new TH3D("hRecoDijetPtEtaDphi","Reco dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphi->Sumw2();

    const Int_t dijetEtaBins{30};
    Double_t dijetEtaVals[dijetEtaBins+1] { -5.0, -4.0, -3.0, -2.4, -2.2, 
                                            -2.0, -1.8, -1.6, -1.4, -1.2, 
                                            -1.0, -0.8, -0.6, -0.4, -0.2,  
                                             0.0,  0.2,  0.4,  0.6,  0.8,  
                                             1.0,  1.2,  1.4,  1.6,  1.8,  
                                             2.0,  2.2,  2.4,  3.0,  4.0,  
                                             5.0 };
    const Int_t dijetPtBins{17};
    Double_t dijetPtVals[dijetPtBins+1] {  40.,  50.,   60.,  70.,  80.,
                                           90., 100.,  110., 120., 130.,
                                          140., 150.,  160., 180., 200., 
                                          240., 300., 1000.};
    // Modify bins of gen dijets
    hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hGenDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hGenDijetPtEtaDphi->GetXaxis()->Set(dijetPtBins, dijetPtVals);
    hGenDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

    // Modify bins of reco dijets
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphi->GetXaxis()->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

    // Modify bins of reco <-> ref dijets
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->GetAxis(6)->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(6)->Set(dijetPtBins, dijetPtVals);
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);

    hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(0)->Set(dijetPtBins, dijetPtVals);
    hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(6)->Set(dijetPtBins, dijetPtVals);
    hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);

    hRefDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRefDijetEtaVsRecoDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRefDijetEtaVsRecoDijetEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRefDijetPtEtaDphi->GetXaxis()->Set(dijetPtBins, dijetPtVals);
    hRefDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRefSelDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRefSelDijetPtEtaDphi->GetXaxis()->Set(dijetPtBins, dijetPtVals);
    hRefSelDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
}

//________________
void HistoManagerDiJet::writeOutput() {
    // Event
    hVz->Write();
    hVzWeighted->Write();
    hMult->Write();
    hHiBin->Write();
    hHiBinWeighted->Write();
    hPtHat->Write();
    hPtHatWeighted->Write();
    hPtHatWeight->Write();
    // hCentrality->Write();
    // hCentralityWeighted->Write();
    hVzPtHat->Write();
    hVzPtHatWeighted->Write();

    // Reco jets

    if (fIsMc) {
        // Gen jets
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Write();
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Write();
        hGenInclusiveJetPt->Write();
        hGenInclusiveJetPtEta->Write();
        hGenPtLeadPtSublead->Write();
        hGenEtaLeadEtaSublead->Write();
        hGenDijetPtEtaDphi->Write();
        hGenDijetEta->Write();

        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hJESInclusiveJetPtEtaPhiPtHat->Write();
        hJESInclusiveJetPtEtaPhiPtHatWeighted->Write();
        hRecoMatchedPtEta->Write();
        hRecoInclusiveMatchedJetPtVsEta->Write();

        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Write();
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Write();
        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Write();
        hRefDijetEta->Write();
        hRefDijetEtaVsRecoDijetEta->Write();
        hRefDijetPtEtaDphi->Write();
        hRefSelDijetPtEtaDphi->Write();
        hRefSelDijetEta->Write();
        hRefInclusiveJetPt->Write();
        hRefInclusiveJetPtEta->Write();
        hRefPtLeadPtSublead->Write();
        hRefEtaLeadEtaSublead->Write();
    } // if ( fIsMc )

    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Write();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Write();
    hRecoInclusiveJetPt->Write();
    hRecoInclusiveAllJetPtVsEta->Write();
    hRecoPtLeadPtSublead->Write();
    hRecoEtaLeadEtaSublead->Write();
    hRecoDijetEta->Write();
    hRecoDijetPtEtaDphi->Write();
}