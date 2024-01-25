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
  fCentBins{10}, fCentRange{-10., 90.},
  fPtBins{50}, fPtRange{20., 520.}, 
  fEtaBins{50}, fEtaRange{-2.5, 2.5},
  fPhiBins{32}, fPhiRange{-TMath::Pi(), TMath::Pi()},
  fPtHatBins{120}, fPtHatRange{15., 615.},
  hVzPtHatCent{nullptr}, hVzPtHatCentWeighted{nullptr},

  hVz{nullptr}, hVzWeighted{nullptr}, hMult{nullptr},
  hHiBin{nullptr}, hHiBinWeighted{nullptr},
  hPtHat{nullptr}, hPtHatWeighted{nullptr}, hPtHatWeight{nullptr},
  hCentrality{nullptr}, hCentralityWeighted{nullptr},

  hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent{nullptr},
  hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted{nullptr},
  hGenDiJetPtEtaDeltaPhiCent{nullptr},
  hGenDiJetPtEtaDeltaPhiCentWeighted{nullptr}
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
    if (hCentrality)    delete hCentrality;
    if (hCentralityWeighted) delete hCentralityWeighted;
    if (hVzPtHatCent) delete hVzPtHatCent;
    if (hVzPtHatCentWeighted) delete hVzPtHatCentWeighted;

    if (fIsMc) {
        // Gen jets
        if (hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent) delete hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent;
        if (hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted) delete hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted;
        if (hGenDiJetPtEtaDeltaPhiCent) delete hGenDiJetPtEtaDeltaPhiCent;
        if (hGenDiJetPtEtaDeltaPhiCentWeighted) delete hGenDiJetPtEtaDeltaPhiCentWeighted;
    } // if (fIsMc)
}

//________________
void HistoManagerDiJet::init(const Bool_t& isMc) {
    
    Int_t    prescale = 2;

    Int_t    vzBins = 360;
    Double_t vzRange[2] {-45., 45.};
    Int_t    multBins{1800};
    Double_t multRange[2] {-0.5, 1799.5};
    Int_t    hiBinBins{203};
    Double_t hiBinRange[2] {-1.5, 201.5};
    Int_t    centralityBins{101};
    Double_t centralityRange[2] {-0.5, 100.5};
    Int_t    weightBins{110};
    Double_t weightRange[2] {-0.05, 1.05};
    Int_t    ptHatBins{100};
    Double_t ptHatRange[2] {0., 1000.};
    Int_t    deltaRBins{80};
    Double_t deltaRRange[2] {-0.05, 0.35};

    Int_t    diJetPtBins{80};
    Double_t diJetPtRange[2] { 50., 450. };
    Int_t    diJetEtaBins{21};
    Double_t diJetEtaRange[2] { -2.922, 3.0 };
    // Int_t    diJetEtaBins{70};
    // Double_t diJetEtaRange[2] { -3.5, 3.5 };
    Int_t    diJetPhiBins{62};
    Double_t diJetPhiRange[2] { -TMath::Pi(), TMath::Pi() };

    Int_t    bins3D_ev_VzCent[3] {     vzBins,     fPtHatBins, fCentBins     };
    Double_t xmin3D_ev_VzCent[3] { vzRange[0], fPtHatRange[0], fCentRange[0] };
    Double_t xmax3D_ev_VzCent[3] { vzRange[1], fPtHatRange[1], fCentRange[1] };

    Int_t    bins7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent[7] { fPtBins,     fEtaBins,     fPhiBins, fPtBins, fEtaBins, fPhiBins, fCentBins };
    Double_t xmin7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent[7] { fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fCentRange[0] };
    Double_t xmax7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent[7] { fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fCentRange[1] };

    Int_t    bins4D_diJet_PtEtaPhiCent[4] { diJetPtBins,  diJetEtaBins, diJetPhiBins, fCentBins };
    Double_t xmin4D_diJet_PtEtaPhiCent[4] { diJetPtRange[0], diJetEtaRange[0], diJetPhiRange[0], fCentRange[0] };
    Double_t xmax4D_diJet_PtEtaPhiCent[4] { diJetPtRange[1], diJetEtaRange[1], diJetPhiRange[1], fCentRange[1] };

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
    hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries",
                           centralityBins, centralityRange[0], centralityRange[1]);
    hCentrality->Sumw2();
    hCentralityWeighted = new TH1D("hCentralityWeighted","Collision centrality weighted;Centrality (%);Entries",
                                   centralityBins, centralityRange[0], centralityRange[1]);
    hCentralityWeighted->Sumw2();

    hVzPtHatCent = new THnSparseD( "hVzPtHatCent","Vertex and centrality;vz (cm);#hat{p_{T}} (GeV/c);centrality",
                                   3, bins3D_ev_VzCent, xmin3D_ev_VzCent, xmax3D_ev_VzCent );
    hVzPtHatCent->Sumw2();
    hVzPtHatCentWeighted = new THnSparseD( "hVzPtHatCentWeighted","Vertex and centrality weighted;vz (cm);#hat{p_{T}} (GeV/c);centrality",
                                           3, bins3D_ev_VzCent, xmin3D_ev_VzCent, xmax3D_ev_VzCent );
    hVzPtHatCentWeighted->Sumw2();


    //
    // Monte Carlo information
    //
    if (fIsMc) {

        //
        // Gen jets
        //

        hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent = new THnSparseD("hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent","Generated leading and subleading jet parameters;p_{T}^{Lead} (GeV/c);#eta^{Lead};#phi^{Lead} (rad);p_{T}^{SubLead} (GeV/c);#eta^{SubLead};#phi^{SubLead} (rad);Centrality (%)",
                                                                7, bins7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent, xmin7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent, xmax7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent);
        hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent->Sumw2();
        hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted = new THnSparseD("hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted","Generated leading and subleading jet parameters weighted with #hat{p_{T}} * centrality;p_{T}^{Lead} (GeV/c);#eta^{Lead};#phi^{Lead} (rad);p_{T}^{SubLead} (GeV/c);#eta^{SubLead};#phi^{SubLead} (rad);Centrality (%)",
                                                                        7, bins7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent, xmin7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent, xmax7D_jet_LeadPtEtaPhiSubLeadPtEtaPhiCent);
        hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted->Sumw2();

        hGenDiJetPtEtaDeltaPhiCent = new THnSparseD("hGenDiJetPtEtaDeltaPhiCent","Generated dijets;p_{T}^{dijet} (GeV/c);#eta^{dijet};#Delta#phi^{dijet} (rad); Centrality (%)",
                                                    4, bins4D_diJet_PtEtaPhiCent, xmin4D_diJet_PtEtaPhiCent, xmax4D_diJet_PtEtaPhiCent);
        hGenDiJetPtEtaDeltaPhiCent->Sumw2();
        hGenDiJetPtEtaDeltaPhiCentWeighted = new THnSparseD("hGenDiJetPtEtaDeltaPhiCentWeighted","Generated dijets weighted with #hat{p_{T}} * centrality;p_{T}^{dijet} (GeV/c);#eta^{dijet};#Delta#phi^{dijet} (rad); Centrality (%)",
                                                            4, bins4D_diJet_PtEtaPhiCent, xmin4D_diJet_PtEtaPhiCent, xmax4D_diJet_PtEtaPhiCent);
        hGenDiJetPtEtaDeltaPhiCentWeighted->Sumw2();

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