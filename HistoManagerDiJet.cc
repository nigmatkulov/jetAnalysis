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

  // Dijets (MC)
  hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta{nullptr},
  hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},
  // hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat{nullptr},
  // hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted{nullptr},
  hRefInclusiveJetPt{nullptr},
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
      if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiPtHat;
      if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiPtHatWeighted;
      if (hGenInclusiveJetPt) delete hGenInclusiveJetPt;
      if (hGenInclusiveJetPtEta) delete hGenInclusiveJetPtEta;
      if (hGenPtLeadPtSublead) delete hGenPtLeadPtSublead;
      if (hGenEtaLeadEtaSublead) delete hGenEtaLeadEtaSublead;


      if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
      if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
      if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
      if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
      if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
      if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
      if (hJESInclusiveJetPtEtaPhiPtHat) delete hJESInclusiveJetPtEtaPhiPtHat;
      if (hJESInclusiveJetPtEtaPhiPtHatWeighted) delete hJESInclusiveJetPtEtaPhiPtHatWeighted;
      if (hRecoMatchedPtEta) delete hRecoMatchedPtEta;
    } // if (fIsMc)

    //
    // Dijets (experiment)
    //
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    if (hRecoInclusiveJetPt) delete hRecoInclusiveJetPt;
    if (hRecoPtLeadPtSublead) delete hRecoPtLeadPtSublead;
    if (hRecoEtaLeadEtaSublead) delete hRecoEtaLeadEtaSublead;

    //
    // Dijets exp vs mc
    //
    if ( fIsMc ) {
      if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta;
      if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;
      // if (hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat) delete hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat;
      // if (hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted) delete hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted;
      if (hRefInclusiveJetPt) delete hRefInclusiveJetPt;
      if (hRefPtLeadPtSublead) delete hRefPtLeadPtSublead;
      if (hRefEtaLeadEtaSublead) delete hRefEtaLeadEtaSublead;
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
    Double_t fJESRange[2] {0., 2.},


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

      //
      // Reco single jets
      //
      hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenPtHat = new THnSparseD("hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenPtHat","Reconstructed inclusive jets;Reco p_{T, corr}^{Inclusive} (GeV/c);Reco p_{T, raw}^{Inclusive} (GeV/c);Ref p_{T}^{Inclusive} (GeV/c);Reco #eta^{Inclusive};Ref #eta^{Inclusive};#hat{p_{T}} (GeV/c)",
      )

    } // if (fIsMc)
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
    hCentrality->Write();
    hCentralityWeighted->Write();
    hVzPtHatCent->Write();
    hVzPtHatCentWeighted->Write();

    // Reco jets

    if (fIsMc) {
        // Gen jets
        hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent->Write();
        hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted->Write();
        hGenDiJetPtEtaDeltaPhiCent->Write();
        hGenDiJetPtEtaDeltaPhiCentWeighted->Write();
    }
}