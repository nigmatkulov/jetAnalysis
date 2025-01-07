/**
 * @file HistoManagerJetESR.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for JES and JER studies
 * @version 1.01
 * @date 2025-01-06
 * 
 * @copyright Copyright (c) 2025
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
    fVerbose{false},
    fPtHatBins{100}, fPtHatRange{15., 1015.},
    fCentBins{10}, fCentRange{-10., 90.},
    fJetPtBins{10}, fJetPtRange{15., 1015.}, 
    fJetEtaBins{52}, fJetEtaRange{-5.2, 5.2},
    fJetPhiBins{32}, fJetPhiRange{-TMath::Pi(), TMath::Pi()},
    fJESBins{100}, fJESRange{0., 2.},
    fJetFlavorForBBins{14}, fJetFlavorForBRange{-6.5, 6.5},
    fDijetPtBins{194}, fDijetPtRange{30., 100.},
    fDijetEtaBins{48}, fDijetEtaRange{-4.8, 4.8},
    fDijetDphiBins{32}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},

    hVz{nullptr}, hVzCentWeighted{nullptr}, hVzPtHatWeighted{nullptr}, hVzWeighted{nullptr},
    hHiBin{nullptr}, hHiBinPtHatWeighted{nullptr}, hHiBinWeighted{nullptr},
    hPtHat{nullptr}, hPtHatPtHatWeighted{nullptr}, hPtHatCentWeighted{nullptr}, hPtHatWeighted{nullptr},
    hPtHatWeight{nullptr}, hPtHatWeightWeighted{nullptr},
    hCentrality{nullptr}, hCentralityPtHatWeighted{nullptr}, hCentralityWeighted{nullptr},
    hVzPtHatCent{nullptr}, hVzPtHatCentPtHatWeighted{nullptr}, hVzPtHatCentWeighted{nullptr},

    // Gen histograms
    hGenInclusiveJetPt{nullptr}, hGenInclusiveJetPtWeighted{nullptr},
    hGenInclusiveJetEtaPt{nullptr}, hGenInclusiveJetEtaPtWeighted{nullptr},
    hGenInclusiveJetPtEtaPhiFlavPtHatCent{nullptr}, 
    hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hGenInclusiveLeadJetPt{nullptr}, hGenInclusiveLeadJetPtWeighted{nullptr},
    hGenInclusiveLeadJetEtaPt{nullptr}, hGenInclusiveLeadJetEtaPtWeighted{nullptr},
    hGenInclusiveLeadJetPtEtaPhiFlavPtHatCent{nullptr},
    hGenInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hGenInclusiveSubLeadJetPt{nullptr}, hGenInclusiveSubLeadJetPtWeighted{nullptr},
    hGenInclusiveSubLeadJetEtaPt{nullptr}, hGenInclusiveSubLeadJetEtaPtWeighted{nullptr},
    hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCent{nullptr},
    hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hGenInclusiveDijetDphi{nullptr}, hGenInclusiveDijetDphiWeighted{nullptr},
    hGenInclusiveDijetEtaPt{nullptr}, hGenInclusiveDijetEtaPtWeighted{nullptr},
    hGenInclusiveDijetEtaPtDphi{nullptr}, hGenInclusiveDijetEtaPtDphiWeighted{nullptr},
    hGenInclusiveDijetDetaCM{nullptr}, hGenInclusiveDijetDetaCMWeighted{nullptr},
    hGenInclusiveDijetDetaCMPt{nullptr}, hGenInclusiveDijetDetaCMPtWeighted{nullptr},
    hGenInclusiveDijetEtaDetaCMPt{nullptr}, hGenInclusiveDijetEtaDetaCMPtWeighted{nullptr},
    hGenInclusiveDijetXPb{nullptr}, hGenInclusiveDijetXPbWeighted{nullptr},
    hGenInclusiveDijetXp{nullptr}, hGenInclusiveDijetXpWeighted{nullptr},
    hGenInclusiveDijetXPbOverXp{nullptr}, hGenInclusiveDijetXPbOverXpWeighted{nullptr},
    hGenInclusiveDijetXPbOverXpEta{nullptr}, hGenInclusiveDijetXPbOverXpEtaWeighted{nullptr},

    hGenSelectedLeadJetPt{nullptr}, hGenSelectedLeadJetPtWeighted{nullptr},
    hGenSelectedLeadJetEtaPt{nullptr}, hGenSelectedLeadJetEtaPtWeighted{nullptr},
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCent{nullptr},
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hGenSelectedSubLeadJetPt{nullptr}, hGenSelectedSubLeadJetPtWeighted{nullptr},
    hGenSelectedSubLeadJetEtaPt{nullptr}, hGenSelectedSubLeadJetEtaPtWeighted{nullptr},
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent{nullptr},
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hGenSelectedDijetPt{nullptr}, hGenSelectedDijetPtWeighted{nullptr},
    hGenSelectedDijetEta{nullptr}, hGenSelectedDijetEtaWeighted{nullptr},
    hGenSelectedDijetEtaPt{nullptr}, hGenSelectedDijetEtaPtWeighted{nullptr},
    hGenSelectedDijetEtaPtDphi{nullptr}, hGenSelectedDijetEtaPtDphiWeighted{nullptr},
    hGenSelectedDijetDetaCM{nullptr}, hGenSelectedDijetDetaCMWeighted{nullptr},
    hGenSelectedDijetDetaCMPt{nullptr}, hGenSelectedDijetDetaCMPtWeighted{nullptr},
    hGenSelectedDijetEtaDetaCMPt{nullptr}, hGenSelectedDijetEtaDetaCMPtWeighted{nullptr},
    hGenSelectedDijetXPb{nullptr}, hGenSelectedDijetXPbWeighted{nullptr},
    hGenSelectedDijetXp{nullptr}, hGenSelectedDijetXpWeighted{nullptr},
    hGenSelectedDijetXPbOverXp{nullptr}, hGenSelectedDijetXPbOverXpWeighted{nullptr},
    hGenSelectedDijetXPbOverXpEta{nullptr}, hGenSelectedDijetXPbOverXpEtaWeighted{nullptr},

    // Reco histograms
    hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent{nullptr},

    hRecoInclusiveJetPt{nullptr}, hRecoInclusiveJetPtWeighted{nullptr},
    hRecoInclusiveJetEtaPt{nullptr}, hRecoInclusiveJetEtaPtWeighted{nullptr},
    hRecoInclusiveJetPtRawPtCorrEta{nullptr}, hRecoInclusiveJetPtRawPtCorrEtaWeighted{nullptr},

    hRecoInclusiveUnmatchedJetPt{nullptr}, hRecoInclusiveUnmatchedJetPtWeighted{nullptr},
    hRecoInclusiveUnmatchedJetEtaPt{nullptr}, hRecoInclusiveUnmatchedJetEtaPtWeighted{nullptr},
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent{nullptr}, hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRecoInclusiveMatchedJetPt{nullptr}, hRecoInclusiveMatchedJetPtWeighted{nullptr},
    hRecoInclusiveMatchedJetEtaPt{nullptr}, hRecoInclusiveMatchedJetEtaPtWeighted{nullptr},
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent{nullptr}, hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted{nullptr},
    hRecoInclusiveJetJESPtEtaPhiCent{nullptr}, hRecoInclusiveJetJESPtEtaPhiCentWeighted{nullptr},

    hReco2RefInclusiveJetPt{nullptr}, hReco2RefInclusiveJetPtWeighted{nullptr},
    hReco2RefInclusiveJetEta{nullptr}, hReco2RefInclusiveJetEtaWeighted{nullptr},
    hReco2RefInclusiveJetPhi{nullptr}, hReco2RefInclusiveJetPhiWeighted{nullptr},
    hReco2RefInclusiveJetPtEtaPhiPtHatCentrality{nullptr}, hReco2RefInclusiveJetPtEtaPhiPtHatCentralityWeighted{nullptr},

    hRecoInclusiveLeadJetPt{nullptr}, hRecoInclusiveLeadJetPtWeighted{nullptr},
    hRecoInclusiveLeadJetEtaPt{nullptr}, hRecoInclusiveLeadJetEtaPtWeighted{nullptr},
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},
    hRecoInclusiveLeadJetJESPtEtaPhiCent{nullptr}, hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted{nullptr},

    hReco2RefInclusiveLeadJetPt{nullptr}, hReco2RefInclusiveLeadJetPtWeighted{nullptr},
    hReco2RefInclusiveLeadJetEta{nullptr}, hReco2RefInclusiveLeadJetEtaWeighted{nullptr},
    hReco2RefInclusiveLeadJetPhi{nullptr}, hReco2RefInclusiveLeadJetPhiWeighted{nullptr},
    hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality{nullptr}, hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentralityWeighted{nullptr},

    hRecoInclusiveSubLeadJetPt{nullptr}, hRecoInclusiveSubLeadJetPtWeighted{nullptr},
    hRecoInclusiveSubLeadJetEtaPt{nullptr}, hRecoInclusiveSubLeadJetEtaPtWeighted{nullptr},
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},
    hRecoInclusiveSubLeadJetJESPtEtaPhiCent{nullptr}, hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted{nullptr},

    hReco2RefInclusiveSubLeadJetPt{nullptr}, hReco2RefInclusiveSubLeadJetPtWeighted{nullptr},
    hReco2RefInclusiveSubLeadJetEta{nullptr}, hReco2RefInclusiveSubLeadJetEtaWeighted{nullptr},
    hReco2RefInclusiveSubLeadJetPhi{nullptr}, hReco2RefInclusiveSubLeadJetPhiWeighted{nullptr},
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality{nullptr}, hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted{nullptr},

    hRecoInclusiveDijetDphi{nullptr}, hRecoInclusiveDijetDphiWeighted{nullptr},
    hRecoInclusiveDijetEtaPt{nullptr}, hRecoInclusiveDijetEtaPtWeighted{nullptr},
    hRecoInclusiveDijetEtaPtDphi{nullptr}, hRecoInclusiveDijetEtaPtDphiWeighted{nullptr},
    hRecoInclusiveDijetDetaCM{nullptr}, hRecoInclusiveDijetDetaCMWeighted{nullptr},
    hRecoInclusiveDijetDetaCMPt{nullptr}, hRecoInclusiveDijetDetaCMPtWeighted{nullptr},
    hRecoInclusiveDijetEtaDetaCMPt{nullptr}, hRecoInclusiveDijetEtaDetaCMPtWeighted{nullptr},
    hRecoInclusiveDijetJESPtEtaDphiCent{nullptr}, hRecoInclusiveDijetJESPtEtaDphiCentWeighted{nullptr},
    
    hRecoSelectedLeadJetPt{nullptr}, hRecoSelectedLeadJetPtWeighted{nullptr},
    hRecoSelectedLeadJetEtaPt{nullptr}, hRecoSelectedLeadJetEtaPtWeighted{nullptr},
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRecoSelectedSubLeadJetPt{nullptr}, hRecoSelectedSubLeadJetPtWeighted{nullptr},
    hRecoSelectedSubLeadJetEtaPt{nullptr}, hRecoSelectedSubLeadJetEtaPtWeighted{nullptr},
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRecoSelectedDijetPt{nullptr}, hRecoSelectedDijetPtWeighted{nullptr},
    hRecoSelectedDijetEta{nullptr}, hRecoSelectedDijetEtaWeighted{nullptr},
    hRecoSelectedDijetEtaPt{nullptr}, hRecoSelectedDijetEtaPtWeighted{nullptr},
    hRecoSelectedDijetEtaPtDphi{nullptr}, hRecoSelectedDijetEtaPtDphiWeighted{nullptr},
    hRecoSelectedDijetDetaCM{nullptr}, hRecoSelectedDijetDetaCMWeighted{nullptr},
    hRecoSelectedDijetDetaCMPt{nullptr}, hRecoSelectedDijetDetaCMPtWeighted{nullptr},
    hRecoSelectedDijetEtaDetaCMPt{nullptr}, hRecoSelectedDijetEtaDetaCMPtWeighted{nullptr},

    hReco2RefSelectedDijetPtEtaFull{nullptr}, hReco2RefSelectedDijetPtEtaFullWeighted{nullptr},

    // Ref histograms
    hRefInclusiveJetPt{nullptr}, hRefInclusiveJetPtWeighted{nullptr},
    hRefInclusiveJetEtaPt{nullptr}, hRefInclusiveJetEtaPtWeighted{nullptr},
    hRefInclusiveJetPtEtaPhiFlavPtHatCent{nullptr}, hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRefInclusiveLeadJetPt{nullptr}, hRefInclusiveLeadJetPtWeighted{nullptr},
    hRefInclusiveLeadJetEtaPt{nullptr}, hRefInclusiveLeadJetEtaPtWeighted{nullptr},
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRefInclusiveSubLeadJetPt{nullptr}, hRefInclusiveSubLeadJetPtWeighted{nullptr},
    hRefInclusiveSubLeadJetEtaPt{nullptr}, hRefInclusiveSubLeadJetEtaPtWeighted{nullptr},
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRefInclusiveDijetDphi{nullptr}, hRefInclusiveDijetDphiWeighted{nullptr},
    hRefInclusiveDijetPt{nullptr}, hRefInclusiveDijetPtWeighted{nullptr},
    hRefInclusiveDijetEta{nullptr}, hRefInclusiveDijetEtaWeighted{nullptr},
    hRefInclusiveDijetEtaPt{nullptr}, hRefInclusiveDijetEtaPtWeighted{nullptr},
    hRefInclusiveDijetEtaPtDphi{nullptr}, hRefInclusiveDijetEtaPtDphiWeighted{nullptr},

    hRefSelectedLeadJetPt{nullptr}, hRefSelectedLeadJetPtWeighted{nullptr},
    hRefSelectedLeadJetEtaPt{nullptr}, hRefSelectedLeadJetEtaPtWeighted{nullptr},
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRefSelectedSubLeadJetPt{nullptr}, hRefSelectedSubLeadJetPtWeighted{nullptr},
    hRefSelectedSubLeadJetEtaPt{nullptr}, hRefSelectedSubLeadJetEtaPtWeighted{nullptr},
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent{nullptr}, hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted{nullptr},

    hRefSelectedDijetPt{nullptr}, hRefSelectedDijetPtWeighted{nullptr},
    hRefSelectedDijetEta{nullptr}, hRefSelectedDijetEtaWeighted{nullptr},
    hRefSelectedDijetEtaPt{nullptr}, hRefSelectedDijetEtaPtWeighted{nullptr},
    hRefSelectedDijetEtaPtDphi{nullptr}, hRefSelectedDijetEtaPtDphiWeighted{nullptr} { 
    
    // Initialize jetId histograms
    for (int iEtaJetId{0}; iEtaJetId<4; iEtaJetId++) {
        hRecoInclusiveJetNHF[iEtaJetId] = nullptr;
        hRecoInclusiveJetNEmF[iEtaJetId] = nullptr;
        hRecoInclusiveJetNumOfConst[iEtaJetId] = nullptr;
        hRecoInclusiveJetMUF[iEtaJetId] = nullptr;
        hRecoInclusiveJetCHF[iEtaJetId] = nullptr;
        hRecoInclusiveJetChargedMult[iEtaJetId] = nullptr;
        hRecoInclusiveJetCEmF[iEtaJetId] = nullptr;
        hRecoInclusiveJetNumOfNeutPart[iEtaJetId] = nullptr;
    } // for (int iEtaJetId{0}; iEtaJetId<4; iEtaJetId++)
}

//________________
HistoManagerJetESR::~HistoManagerJetESR() {
    // Destructor
    delete hVz;
    delete hVzCentWeighted;
    delete hVzPtHatWeighted;
    delete hVzWeighted;
    delete hHiBin;
    delete hHiBinPtHatWeighted;
    delete hHiBinWeighted;
    delete hPtHat;
    delete hPtHatPtHatWeighted;
    delete hPtHatCentWeighted;
    delete hPtHatWeighted;
    delete hPtHatWeight;
    delete hPtHatWeightWeighted;
    delete hCentrality;
    delete hCentralityPtHatWeighted;
    delete hCentralityWeighted;
    delete hVzPtHatCent;
    delete hVzPtHatCentPtHatWeighted;
    delete hVzPtHatCentWeighted;

    // Gen histograms
    if (hGenInclusiveJetPt) delete hGenInclusiveJetPt;
    if (hGenInclusiveJetPtWeighted) delete hGenInclusiveJetPtWeighted;
    if (hGenInclusiveJetEtaPt) delete hGenInclusiveJetEtaPt;
    if (hGenInclusiveJetEtaPtWeighted) delete hGenInclusiveJetEtaPtWeighted;
    if (hGenInclusiveJetPtEtaPhiFlavPtHatCent) delete hGenInclusiveJetPtEtaPhiFlavPtHatCent;
    if (hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted) delete hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted;

    if (hGenInclusiveLeadJetPt) delete hGenInclusiveLeadJetPt;
    if (hGenInclusiveLeadJetPtWeighted) delete hGenInclusiveLeadJetPtWeighted;
    if (hGenInclusiveLeadJetEtaPt) delete hGenInclusiveLeadJetEtaPt;
    if (hGenInclusiveLeadJetEtaPtWeighted) delete hGenInclusiveLeadJetEtaPtWeighted;
    if (hGenInclusiveLeadJetPtEtaPhiFlavPtHatCent) delete hGenInclusiveLeadJetPtEtaPhiFlavPtHatCent;
    if (hGenInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hGenInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hGenInclusiveSubLeadJetPt) delete hGenInclusiveSubLeadJetPt;
    if (hGenInclusiveSubLeadJetPtWeighted) delete hGenInclusiveSubLeadJetPtWeighted;
    if (hGenInclusiveSubLeadJetEtaPt) delete hGenInclusiveSubLeadJetEtaPt;
    if (hGenInclusiveSubLeadJetEtaPtWeighted) delete hGenInclusiveSubLeadJetEtaPtWeighted;
    if (hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCent) delete hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCent;
    if (hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hGenInclusiveDijetDphi) delete hGenInclusiveDijetDphi;
    if (hGenInclusiveDijetDphiWeighted) delete hGenInclusiveDijetDphiWeighted;
    if (hGenInclusiveDijetEtaPt) delete hGenInclusiveDijetEtaPt;
    if (hGenInclusiveDijetEtaPtWeighted) delete hGenInclusiveDijetEtaPtWeighted;
    if (hGenInclusiveDijetEtaPtDphi) delete hGenInclusiveDijetEtaPtDphi;
    if (hGenInclusiveDijetEtaPtDphiWeighted) delete hGenInclusiveDijetEtaPtDphiWeighted;
    if (hGenInclusiveDijetDetaCM) delete hGenInclusiveDijetDetaCM;
    if (hGenInclusiveDijetDetaCMWeighted) delete hGenInclusiveDijetDetaCMWeighted;
    if (hGenInclusiveDijetDetaCMPt) delete hGenInclusiveDijetDetaCMPt;
    if (hGenInclusiveDijetDetaCMPtWeighted) delete hGenInclusiveDijetDetaCMPtWeighted;
    if (hGenInclusiveDijetEtaDetaCMPt) delete hGenInclusiveDijetEtaDetaCMPt;
    if (hGenInclusiveDijetEtaDetaCMPtWeighted) delete hGenInclusiveDijetEtaDetaCMPtWeighted;
    if (hGenInclusiveDijetXPb) delete hGenInclusiveDijetXPb;
    if (hGenInclusiveDijetXPbWeighted) delete hGenInclusiveDijetXPbWeighted;
    if (hGenInclusiveDijetXp) delete hGenInclusiveDijetXp;
    if (hGenInclusiveDijetXpWeighted) delete hGenInclusiveDijetXpWeighted;
    if (hGenInclusiveDijetXPbOverXp) delete hGenInclusiveDijetXPbOverXp;
    if (hGenInclusiveDijetXPbOverXpWeighted) delete hGenInclusiveDijetXPbOverXpWeighted;
    if (hGenInclusiveDijetXPbOverXpEta) delete hGenInclusiveDijetXPbOverXpEta;
    if (hGenInclusiveDijetXPbOverXpEtaWeighted) delete hGenInclusiveDijetXPbOverXpEtaWeighted;

    if (hGenSelectedLeadJetPt) delete hGenSelectedLeadJetPt;
    if (hGenSelectedLeadJetPtWeighted) delete hGenSelectedLeadJetPtWeighted;
    if (hGenSelectedLeadJetEtaPt) delete hGenSelectedLeadJetEtaPt;
    if (hGenSelectedLeadJetEtaPtWeighted) delete hGenSelectedLeadJetEtaPtWeighted;
    if (hGenSelectedLeadJetPtEtaPhiFlavPtHatCent) delete hGenSelectedLeadJetPtEtaPhiFlavPtHatCent;
    if (hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hGenSelectedSubLeadJetPt) delete hGenSelectedSubLeadJetPt;
    if (hGenSelectedSubLeadJetPtWeighted) delete hGenSelectedSubLeadJetPtWeighted;
    if (hGenSelectedSubLeadJetEtaPt) delete hGenSelectedSubLeadJetEtaPt;
    if (hGenSelectedSubLeadJetEtaPtWeighted) delete hGenSelectedSubLeadJetEtaPtWeighted;
    if (hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent) delete hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent;
    if (hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hGenSelectedDijetPt) delete hGenSelectedDijetPt;
    if (hGenSelectedDijetPtWeighted) delete hGenSelectedDijetPtWeighted;
    if (hGenSelectedDijetEta) delete hGenSelectedDijetEta;
    if (hGenSelectedDijetEtaWeighted) delete hGenSelectedDijetEtaWeighted;
    if (hGenSelectedDijetEtaPt) delete hGenSelectedDijetEtaPt;
    if (hGenSelectedDijetEtaPtWeighted) delete hGenSelectedDijetEtaPtWeighted;
    if (hGenSelectedDijetEtaPtDphi) delete hGenSelectedDijetEtaPtDphi;
    if (hGenSelectedDijetEtaPtDphiWeighted) delete hGenSelectedDijetEtaPtDphiWeighted;
    if (hGenSelectedDijetDetaCM) delete hGenSelectedDijetDetaCM;
    if (hGenSelectedDijetDetaCMWeighted) delete hGenSelectedDijetDetaCMWeighted;
    if (hGenSelectedDijetDetaCMPt) delete hGenSelectedDijetDetaCMPt;
    if (hGenSelectedDijetDetaCMPtWeighted) delete hGenSelectedDijetDetaCMPtWeighted;
    if (hGenSelectedDijetEtaDetaCMPt) delete hGenSelectedDijetEtaDetaCMPt;
    if (hGenSelectedDijetEtaDetaCMPtWeighted) delete hGenSelectedDijetEtaDetaCMPtWeighted;
    if (hGenSelectedDijetXPb) delete hGenSelectedDijetXPb;
    if (hGenSelectedDijetXPbWeighted) delete hGenSelectedDijetXPbWeighted;
    if (hGenSelectedDijetXp) delete hGenSelectedDijetXp;
    if (hGenSelectedDijetXpWeighted) delete hGenSelectedDijetXpWeighted;
    if (hGenSelectedDijetXPbOverXp) delete hGenSelectedDijetXPbOverXp;
    if (hGenSelectedDijetXPbOverXpWeighted) delete hGenSelectedDijetXPbOverXpWeighted;
    if (hGenSelectedDijetXPbOverXpEta) delete hGenSelectedDijetXPbOverXpEta;
    if (hGenSelectedDijetXPbOverXpEtaWeighted) delete hGenSelectedDijetXPbOverXpEtaWeighted;

    // Reco histograms
    if (hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent) delete hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent;

    if (hRecoInclusiveJetPt) delete hRecoInclusiveJetPt;
    if (hRecoInclusiveJetPtWeighted) delete hRecoInclusiveJetPtWeighted;
    if (hRecoInclusiveJetEtaPt) delete hRecoInclusiveJetEtaPt;
    if (hRecoInclusiveJetEtaPtWeighted) delete hRecoInclusiveJetEtaPtWeighted;
    if (hRecoInclusiveJetPtRawPtCorrEta) delete hRecoInclusiveJetPtRawPtCorrEta;
    if (hRecoInclusiveJetPtRawPtCorrEtaWeighted) delete hRecoInclusiveJetPtRawPtCorrEtaWeighted;

    if (hRecoInclusiveUnmatchedJetPt) delete hRecoInclusiveUnmatchedJetPt;
    if (hRecoInclusiveUnmatchedJetPtWeighted) delete hRecoInclusiveUnmatchedJetPtWeighted;
    if (hRecoInclusiveUnmatchedJetEtaPt) delete hRecoInclusiveUnmatchedJetEtaPt;
    if (hRecoInclusiveUnmatchedJetEtaPtWeighted) delete hRecoInclusiveUnmatchedJetEtaPtWeighted;
    if (hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent) delete hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent;
    if (hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted) delete hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRecoInclusiveMatchedJetPt) delete hRecoInclusiveMatchedJetPt;
    if (hRecoInclusiveMatchedJetPtWeighted) delete hRecoInclusiveMatchedJetPtWeighted;
    if (hRecoInclusiveMatchedJetEtaPt) delete hRecoInclusiveMatchedJetEtaPt;
    if (hRecoInclusiveMatchedJetEtaPtWeighted) delete hRecoInclusiveMatchedJetEtaPtWeighted;
    if (hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent) delete hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent;
    if (hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted) delete hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted;
    if (hRecoInclusiveJetJESPtEtaPhiCent) delete hRecoInclusiveJetJESPtEtaPhiCent;
    if (hRecoInclusiveJetJESPtEtaPhiCentWeighted) delete hRecoInclusiveJetJESPtEtaPhiCentWeighted;

    if (hReco2RefInclusiveJetPt) delete hReco2RefInclusiveJetPt;
    if (hReco2RefInclusiveJetPtWeighted) delete hReco2RefInclusiveJetPtWeighted;
    if (hReco2RefInclusiveJetEta) delete hReco2RefInclusiveJetEta;
    if (hReco2RefInclusiveJetEtaWeighted) delete hReco2RefInclusiveJetEtaWeighted;
    if (hReco2RefInclusiveJetPhi) delete hReco2RefInclusiveJetPhi;
    if (hReco2RefInclusiveJetPhiWeighted) delete hReco2RefInclusiveJetPhiWeighted;
    if (hReco2RefInclusiveJetPtEtaPhiPtHatCentrality) delete hReco2RefInclusiveJetPtEtaPhiPtHatCentrality;
    if (hReco2RefInclusiveJetPtEtaPhiPtHatCentralityWeighted) delete hReco2RefInclusiveJetPtEtaPhiPtHatCentralityWeighted;

    if (hRecoInclusiveLeadJetPt) delete hRecoInclusiveLeadJetPt;
    if (hRecoInclusiveLeadJetPtWeighted) delete hRecoInclusiveLeadJetPtWeighted;
    if (hRecoInclusiveLeadJetEtaPt) delete hRecoInclusiveLeadJetEtaPt;
    if (hRecoInclusiveLeadJetEtaPtWeighted) delete hRecoInclusiveLeadJetEtaPtWeighted;
    if (hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent) delete hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent;
    if (hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted;
    if (hRecoInclusiveLeadJetJESPtEtaPhiCent) delete hRecoInclusiveLeadJetJESPtEtaPhiCent;
    if (hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted) delete hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted;

    if (hReco2RefInclusiveLeadJetPt) delete hReco2RefInclusiveLeadJetPt;
    if (hReco2RefInclusiveLeadJetPtWeighted) delete hReco2RefInclusiveLeadJetPtWeighted;
    if (hReco2RefInclusiveLeadJetEta) delete hReco2RefInclusiveLeadJetEta;
    if (hReco2RefInclusiveLeadJetEtaWeighted) delete hReco2RefInclusiveLeadJetEtaWeighted;
    if (hReco2RefInclusiveLeadJetPhi) delete hReco2RefInclusiveLeadJetPhi;
    if (hReco2RefInclusiveLeadJetPhiWeighted) delete hReco2RefInclusiveLeadJetPhiWeighted;
    if (hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality) delete hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality;
    if (hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentralityWeighted) delete hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentralityWeighted;

    if (hRecoInclusiveSubLeadJetPt) delete hRecoInclusiveSubLeadJetPt;
    if (hRecoInclusiveSubLeadJetPtWeighted) delete hRecoInclusiveSubLeadJetPtWeighted;
    if (hRecoInclusiveSubLeadJetEtaPt) delete hRecoInclusiveSubLeadJetEtaPt;
    if (hRecoInclusiveSubLeadJetEtaPtWeighted) delete hRecoInclusiveSubLeadJetEtaPtWeighted;
    if (hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent) delete hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent;
    if (hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted;
    if (hRecoInclusiveSubLeadJetJESPtEtaPhiCent) delete hRecoInclusiveSubLeadJetJESPtEtaPhiCent;
    if (hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted) delete hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted;

    if (hReco2RefInclusiveSubLeadJetPt) delete hReco2RefInclusiveSubLeadJetPt;
    if (hReco2RefInclusiveSubLeadJetPtWeighted) delete hReco2RefInclusiveSubLeadJetPtWeighted;
    if (hReco2RefInclusiveSubLeadJetEta) delete hReco2RefInclusiveSubLeadJetEta;
    if (hReco2RefInclusiveSubLeadJetEtaWeighted) delete hReco2RefInclusiveSubLeadJetEtaWeighted;
    if (hReco2RefInclusiveSubLeadJetPhi) delete hReco2RefInclusiveSubLeadJetPhi;
    if (hReco2RefInclusiveSubLeadJetPhiWeighted) delete hReco2RefInclusiveSubLeadJetPhiWeighted;
    if (hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality) delete hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality;
    if (hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted) delete hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted;

    if (hRecoInclusiveDijetDphi) delete hRecoInclusiveDijetDphi;
    if (hRecoInclusiveDijetDphiWeighted) delete hRecoInclusiveDijetDphiWeighted;
    if (hRecoInclusiveDijetEtaPt) delete hRecoInclusiveDijetEtaPt;
    if (hRecoInclusiveDijetEtaPtWeighted) delete hRecoInclusiveDijetEtaPtWeighted;
    if (hRecoInclusiveDijetEtaPtDphi) delete hRecoInclusiveDijetEtaPtDphi;
    if (hRecoInclusiveDijetEtaPtDphiWeighted) delete hRecoInclusiveDijetEtaPtDphiWeighted;
    if (hRecoInclusiveDijetDetaCM) delete hRecoInclusiveDijetDetaCM;
    if (hRecoInclusiveDijetDetaCMWeighted) delete hRecoInclusiveDijetDetaCMWeighted;
    if (hRecoInclusiveDijetDetaCMPt) delete hRecoInclusiveDijetDetaCMPt;
    if (hRecoInclusiveDijetDetaCMPtWeighted) delete hRecoInclusiveDijetDetaCMPtWeighted;
    if (hRecoInclusiveDijetEtaDetaCMPt) delete hRecoInclusiveDijetEtaDetaCMPt;
    if (hRecoInclusiveDijetEtaDetaCMPtWeighted) delete hRecoInclusiveDijetEtaDetaCMPtWeighted;
    if (hRecoInclusiveDijetJESPtEtaDphiCent) delete hRecoInclusiveDijetJESPtEtaDphiCent;
    if (hRecoInclusiveDijetJESPtEtaDphiCentWeighted) delete hRecoInclusiveDijetJESPtEtaDphiCentWeighted;

    if (hRecoSelectedLeadJetPt) delete hRecoSelectedLeadJetPt;
    if (hRecoSelectedLeadJetPtWeighted) delete hRecoSelectedLeadJetPtWeighted;
    if (hRecoSelectedLeadJetEtaPt) delete hRecoSelectedLeadJetEtaPt;
    if (hRecoSelectedLeadJetEtaPtWeighted) delete hRecoSelectedLeadJetEtaPtWeighted;
    if (hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent) delete hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent;
    if (hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRecoSelectedSubLeadJetPt) delete hRecoSelectedSubLeadJetPt;
    if (hRecoSelectedSubLeadJetPtWeighted) delete hRecoSelectedSubLeadJetPtWeighted;
    if (hRecoSelectedSubLeadJetEtaPt) delete hRecoSelectedSubLeadJetEtaPt;
    if (hRecoSelectedSubLeadJetEtaPtWeighted) delete hRecoSelectedSubLeadJetEtaPtWeighted;
    if (hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent) delete hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent;
    if (hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRecoSelectedDijetPt) delete hRecoSelectedDijetPt;
    if (hRecoSelectedDijetPtWeighted) delete hRecoSelectedDijetPtWeighted;
    if (hRecoSelectedDijetEta) delete hRecoSelectedDijetEta;
    if (hRecoSelectedDijetEtaWeighted) delete hRecoSelectedDijetEtaWeighted;
    if (hRecoSelectedDijetEtaPt) delete hRecoSelectedDijetEtaPt;
    if (hRecoSelectedDijetEtaPtWeighted) delete hRecoSelectedDijetEtaPtWeighted;
    if (hRecoSelectedDijetEtaPtDphi) delete hRecoSelectedDijetEtaPtDphi;
    if (hRecoSelectedDijetEtaPtDphiWeighted) delete hRecoSelectedDijetEtaPtDphiWeighted;
    if (hRecoSelectedDijetDetaCM) delete hRecoSelectedDijetDetaCM;
    if (hRecoSelectedDijetDetaCMWeighted) delete hRecoSelectedDijetDetaCMWeighted;
    if (hRecoSelectedDijetDetaCMPt) delete hRecoSelectedDijetDetaCMPt;
    if (hRecoSelectedDijetDetaCMPtWeighted) delete hRecoSelectedDijetDetaCMPtWeighted;
    if (hRecoSelectedDijetEtaDetaCMPt) delete hRecoSelectedDijetEtaDetaCMPt;
    if (hRecoSelectedDijetEtaDetaCMPtWeighted) delete hRecoSelectedDijetEtaDetaCMPtWeighted;

    if (hReco2RefSelectedDijetPtEtaFull) delete hReco2RefSelectedDijetPtEtaFull;
    if (hReco2RefSelectedDijetPtEtaFullWeighted) delete hReco2RefSelectedDijetPtEtaFullWeighted;

    // Ref histograms
    if (hRefInclusiveJetPt) delete hRefInclusiveJetPt;
    if (hRefInclusiveJetPtWeighted) delete hRefInclusiveJetPtWeighted;
    if (hRefInclusiveJetEtaPt) delete hRefInclusiveJetEtaPt;
    if (hRefInclusiveJetEtaPtWeighted) delete hRefInclusiveJetEtaPtWeighted;
    if (hRefInclusiveJetPtEtaPhiFlavPtHatCent) delete hRefInclusiveJetPtEtaPhiFlavPtHatCent;
    if (hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted) delete hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRefInclusiveLeadJetPt) delete hRefInclusiveLeadJetPt;
    if (hRefInclusiveLeadJetPtWeighted) delete hRefInclusiveLeadJetPtWeighted;
    if (hRefInclusiveLeadJetEtaPt) delete hRefInclusiveLeadJetEtaPt;
    if (hRefInclusiveLeadJetEtaPtWeighted) delete hRefInclusiveLeadJetEtaPtWeighted;
    if (hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent) delete hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent;
    if (hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRefInclusiveSubLeadJetPt) delete hRefInclusiveSubLeadJetPt;
    if (hRefInclusiveSubLeadJetPtWeighted) delete hRefInclusiveSubLeadJetPtWeighted;
    if (hRefInclusiveSubLeadJetEtaPt) delete hRefInclusiveSubLeadJetEtaPt;
    if (hRefInclusiveSubLeadJetEtaPtWeighted) delete hRefInclusiveSubLeadJetEtaPtWeighted;
    if (hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent) delete hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent;
    if (hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRefInclusiveDijetDphi) delete hRefInclusiveDijetDphi;
    if (hRefInclusiveDijetDphiWeighted) delete hRefInclusiveDijetDphiWeighted;
    if (hRefInclusiveDijetPt) delete hRefInclusiveDijetPt;
    if (hRefInclusiveDijetPtWeighted) delete hRefInclusiveDijetPtWeighted;
    if (hRefInclusiveDijetEta) delete hRefInclusiveDijetEta;
    if (hRefInclusiveDijetEtaWeighted) delete hRefInclusiveDijetEtaWeighted;
    if (hRefInclusiveDijetEtaPt) delete hRefInclusiveDijetEtaPt;
    if (hRefInclusiveDijetEtaPtWeighted) delete hRefInclusiveDijetEtaPtWeighted;
    if (hRefInclusiveDijetEtaPtDphi) delete hRefInclusiveDijetEtaPtDphi;
    if (hRefInclusiveDijetEtaPtDphiWeighted) delete hRefInclusiveDijetEtaPtDphiWeighted;

    if (hRefSelectedLeadJetPt) delete hRefSelectedLeadJetPt;
    if (hRefSelectedLeadJetPtWeighted) delete hRefSelectedLeadJetPtWeighted;
    if (hRefSelectedLeadJetEtaPt) delete hRefSelectedLeadJetEtaPt;
    if (hRefSelectedLeadJetEtaPtWeighted) delete hRefSelectedLeadJetEtaPtWeighted;
    if (hRefSelectedLeadJetPtEtaPhiFlavPtHatCent) delete hRefSelectedLeadJetPtEtaPhiFlavPtHatCent;
    if (hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRefSelectedSubLeadJetPt) delete hRefSelectedSubLeadJetPt;
    if (hRefSelectedSubLeadJetPtWeighted) delete hRefSelectedSubLeadJetPtWeighted;
    if (hRefSelectedSubLeadJetEtaPt) delete hRefSelectedSubLeadJetEtaPt;
    if (hRefSelectedSubLeadJetEtaPtWeighted) delete hRefSelectedSubLeadJetEtaPtWeighted;
    if (hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent) delete hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent;
    if (hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted) delete hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    if (hRefSelectedDijetPt) delete hRefSelectedDijetPt;
    if (hRefSelectedDijetPtWeighted) delete hRefSelectedDijetPtWeighted;
    if (hRefSelectedDijetEta) delete hRefSelectedDijetEta;
    if (hRefSelectedDijetEtaWeighted) delete hRefSelectedDijetEtaWeighted;
    if (hRefSelectedDijetEtaPt) delete hRefSelectedDijetEtaPt;
    if (hRefSelectedDijetEtaPtWeighted) delete hRefSelectedDijetEtaPtWeighted;
    if (hRefSelectedDijetEtaPtDphi) delete hRefSelectedDijetEtaPtDphi;
    if (hRefSelectedDijetEtaPtDphiWeighted) delete hRefSelectedDijetEtaPtDphiWeighted;

    // Delete jetId histograms
    for (int iEtaJetId{0}; iEtaJetId<4; iEtaJetId++) {
        if (hRecoInclusiveJetNHF[iEtaJetId]) delete hRecoInclusiveJetNHF[iEtaJetId];
        if (hRecoInclusiveJetNEmF[iEtaJetId]) delete hRecoInclusiveJetNEmF[iEtaJetId];
        if (hRecoInclusiveJetNumOfConst[iEtaJetId]) delete hRecoInclusiveJetNumOfConst[iEtaJetId];
        if (hRecoInclusiveJetMUF[iEtaJetId]) delete hRecoInclusiveJetMUF[iEtaJetId];
        if (hRecoInclusiveJetCHF[iEtaJetId]) delete hRecoInclusiveJetCHF[iEtaJetId];
        if (hRecoInclusiveJetChargedMult[iEtaJetId]) delete hRecoInclusiveJetChargedMult[iEtaJetId];
        if (hRecoInclusiveJetCEmF[iEtaJetId]) delete hRecoInclusiveJetCEmF[iEtaJetId];
        if (hRecoInclusiveJetNumOfNeutPart[iEtaJetId]) delete hRecoInclusiveJetNumOfNeutPart[iEtaJetId];
    }
}

//________________
void HistoManagerJetESR::init(const Bool_t& isMc) {
    
    int    vzBins = 320;
    double vzRange[2] = {-31., 31.};
    int    hiBinBins{203};
    double hiBinRange[2] = {-1.5, 201.5};
    int    centralityBins{101};
    double centralityRange[2] = {-0.5, 100.5};
    int    weightBins{110};
    double weightRange[2] = {-0.05, 1.05};
    int    ptHatBins{100};
    double ptHatRange[2] = {0., 1000.};
    int fFracBins{100}; 
    double fFracRange[2] = {0., 1.};
    int fMultBins{32}; 
    double fMultRange[2] = {-0.5, 31.5};
    int xBins = 10000;
    double xRange[2] = {0.0001, 1.};
    int xPbOverXpBins = 10000;
    double xPbOverXpRange[2] = {0.0001, 1000.};

    int bins6D_jet_PtEtaPhiFlavPtHatCent[6] = {fJetPtBins, fJetEtaBins, fJetPhiBins, fJetFlavorForBBins, fPtHatBins, fCentBins};
    double xmin6D_jet_PtEtaPhiFlavPtHatCent[6] = {fJetPtRange[0], fJetEtaRange[0], fJetPhiRange[0], fJetFlavorForBRange[0], fPtHatRange[0], fCentRange[0]};
    double xmax6D_jet_PtEtaPhiFlavPtHatCent[6] = {fJetPtRange[1], fJetEtaRange[1], fJetPhiRange[1], fJetFlavorForBRange[1], fPtHatRange[1], fCentRange[1]};

    int bins14D_jet_PtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent[14] = {fJetPtBins, fJetEtaBins, fJetPhiBins, 1, fFracBins, fFracBins, fMultBins, fFracBins, fFracBins, fMultBins, fFracBins, fMultBins, fPtHatBins, fCentBins};
    double xmin14D_jet_PtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent[14] = {fJetPtRange[0], fJetEtaRange[0], fJetPhiRange[0], 0, fFracRange[0], fFracRange[0], fMultRange[0], fFracRange[0], fFracRange[0], fMultRange[0], fFracRange[0], fMultRange[0], fPtHatRange[0], fCentRange[0]};
    double xmax14D_jet_PtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent[14] = {fJetPtRange[1], fJetEtaRange[1], fJetPhiRange[1], 1, fFracRange[1], fFracRange[1], fMultRange[1], fFracRange[1], fFracRange[1], fMultRange[1], fFracRange[1], fMultRange[1], fPtHatRange[1], fCentRange[1]};

    int bins5D_jes_PtEtaPhiCent[5] = {fJetPtBins, fJetEtaBins, fJetPhiBins, fJESBins, fCentBins};
    double xmin5D_jes_PtEtaPhiCent[5] = {fJetPtRange[0], fJetEtaRange[0], fJetPhiRange[0], fJESRange[0], fCentRange[0]};
    double xmax5D_jes_PtEtaPhiCent[5] = {fJetPtRange[1], fJetEtaRange[1], fJetPhiRange[1], fJESRange[1], fCentRange[1]};

    int bins5D_jes_PtEtaDphiCent[5] = {fDijetPtBins, fDijetEtaBins, fDijetDphiBins, fJESBins, fCentBins};
    double xmin5D_jes_PtEtaDphiCent[5] = {fDijetPtRange[0], fDijetEtaRange[0], fDijetDphiRange[0], fJESRange[0], fCentRange[0]};
    double xmax5D_jes_PtEtaDphiCent[5] = {fDijetPtRange[1], fDijetEtaRange[1], fDijetDphiRange[1], fJESRange[1], fCentRange[1]};
    int bins8D_jet_PtEtaPhiPtHatCentrality[8] = {fJetPtBins, fJetPtBins, fJetEtaBins, fJetEtaBins, fJetPhiBins, fJetPhiBins, fPtHatBins, fCentBins};
    double xmin8D_jet_PtEtaPhiPtHatCentrality[8] = {fJetPtRange[0], fJetPtRange[0], fJetEtaRange[0], fJetEtaRange[0], fJetPhiRange[0], fJetPhiRange[0], fPtHatRange[0], fCentRange[0]};
    double xmax8D_jet_PtEtaPhiPtHatCentrality[8] = {fJetPtRange[1], fJetPtRange[1], fJetEtaRange[1], fJetEtaRange[1], fJetPhiRange[1], fJetPhiRange[1], fPtHatRange[1], fCentRange[1]};

    int bins12D_jet_PtEtaFull[12] = {fDijetPtBins, fDijetEtaBins, fJetPtBins, fJetEtaBins, fJetPtBins, fJetEtaBins, fDijetPtBins, fDijetEtaBins, fJetPtBins, fJetEtaBins, fJetPtBins, fJetEtaBins};
    double xmin12D_jet_PtEtaFull[12] = {fDijetPtRange[0], fDijetEtaRange[0], fJetPtRange[0], fJetEtaRange[0], fJetPtRange[0], fJetEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0], fJetPtRange[0], fJetEtaRange[0], fJetPtRange[0], fJetEtaRange[0]};
    double xmax12D_jet_PtEtaFull[12] = {fDijetPtRange[1], fDijetEtaRange[1], fJetPtRange[1], fJetEtaRange[1], fJetPtRange[1], fJetEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1], fJetPtRange[1], fJetEtaRange[1], fJetPtRange[1], fJetEtaRange[1]};

    int ptHatWeightBins = 100;
    double ptHatWeightRange[2] = {0., 10.};


    //
    // Event histograms
    //

    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVz->Sumw2();
    hVzCentWeighted = new TH1D("hVzCentWeighted", "Vertex z position (centrality weighted);vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzCentWeighted->Sumw2();
    hVzPtHatWeighted = new TH1D("hVzPtHatWeighted", "Vertex z position (#hat{p_{T}} weighted);vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzPtHatWeighted->Sumw2();
    hVzWeighted = new TH1D("hVzWeighted","Vertex z position (weighted);vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzWeighted->Sumw2();

    hHiBin = new TH1D("hHiBin","HiBin a.k.a. centrality;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBin->Sumw2();
    hHiBinPtHatWeighted = new TH1D("hHiBinPtHatWeighted", "HiBin a.k.a. centrality (#hat{p_{T}} weighted);HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinPtHatWeighted->Sumw2();
    hHiBinWeighted = new TH1D("hHiBinWeighted","HiBin a.k.a. centrality (weighted);HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinWeighted->Sumw2();

    hPtHat = new TH1D("hPtHat","#hat{p_{T}};#hat{p_{T}} (GeV);Entries", fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hPtHat->Sumw2();
    hPtHatPtHatWeighted = new TH1D("hPtHatPtHatWeighted", "#hat{p_{T}} (#hat{p_{T}} weighted);#hat{p_{T}} (GeV);Entries", fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hPtHatPtHatWeighted->Sumw2();
    hPtHatWeighted = new TH1D("hPtHatWeighted","#hat{p_{T}} (weighted);#hat{p_{T}} (GeV);Entries", fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hPtHatWeighted->Sumw2();

    hPtHatWeight = new TH1D("hPtHatWeight","#hat{p_{T}} weight;#hat{p_{T}}^{weight};Entries", ptHatWeightBins, ptHatWeightRange[0], ptHatWeightRange[1]);
    hPtHatWeight->Sumw2();
    hPtHatWeightWeighted = new TH1D("hPtHatWeightWeighted","#hat{p_{T}} weight (#hat{p_{T}} weighted);#hat{p_{T}}^{weight};Entries", ptHatWeightBins, ptHatWeightRange[0], ptHatWeightRange[1]);
    hPtHatWeightWeighted->Sumw2();

    hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries", fCentBins, fCentRange[0], fCentRange[1]);
    hCentrality->Sumw2();
    hCentralityPtHatWeighted = new TH1D("hCentralityPtHatWeighted", "Collision centrality (#hat{p_{T}} weighted);Centrality (%);Entries", fCentBins, fCentRange[0], fCentRange[1]);
    hCentralityPtHatWeighted->Sumw2();
    hCentralityWeighted = new TH1D("hCentralityWeighted","Collision centrality (weighted);Centrality (%);Entries", fCentBins, fCentRange[0], fCentRange[1]);
    hCentralityWeighted->Sumw2();

    hVzPtHatCent = new TH3D("hVzPtHatCent","Vertex and centrality;vz (cm);#hat{p_{T}} (GeV);centrality (%)",
                            vzBins, vzRange[0], vzRange[1], fPtHatBins, fPtHatRange[0], fPtHatRange[1], fCentBins, fCentRange[0], fCentRange[1]);
    hVzPtHatCent->Sumw2();
    hVzPtHatCentPtHatWeighted = new TH3D("hVzPtHatCentPtHatWeighted","Vertex, #hat{p_{T}} and centrality (#hat{p_{T}} weighted);vz (cm);#hat{p_{T}} (GeV);centrality (%)",
                                         vzBins, vzRange[0], vzRange[1], fPtHatBins, fPtHatRange[0], fPtHatRange[1], fCentBins, fCentRange[0], fCentRange[1]);
    hVzPtHatCentPtHatWeighted->Sumw2();
    hVzPtHatCentWeighted = new TH3D("hVzPtHatCentWeighted","Vertex, #hat{p_{T}} and centrality (weighted);vz (cm);#hat{p_{T}} (GeV);centrality (%)",
                                    vzBins, vzRange[0], vzRange[1], fPtHatBins, fPtHatRange[0], fPtHatRange[1], fCentBins, fCentRange[0], fCentRange[1]);
                            

    //
    // Gen jets
    //

    // Inclusive jets
    hGenInclusiveJetPt = new TH1D("hGenInclusiveJetPt", "Inclusive Gen Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveJetPt->Sumw2();
    hGenInclusiveJetPtWeighted = new TH1D("hGenInclusiveJetPtWeighted", "Inclusive Gen Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveJetPtWeighted->Sumw2();
    hGenInclusiveJetEtaPt = new TH2D("hGenInclusiveJetEtaPt", "Inclusive Gen Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveJetEtaPt->Sumw2();
    hGenInclusiveJetEtaPtWeighted = new TH2D("hGenInclusiveJetEtaPtWeighted", "Inclusive Gen Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveJetEtaPtWeighted->Sumw2();
    // 0 - pT, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    hGenInclusiveJetPtEtaPhiFlavPtHatCent = new THnSparseD("hGenInclusiveJetPtEtaPhiFlavPtHatCent", "Inclusive Gen Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hGenInclusiveJetPtEtaPhiFlavPtHatCent->Sumw2();
    hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Gen Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Inclusive leading jet
    hGenInclusiveLeadJetPt = new TH1D("hGenInclusiveLeadJetPt", "Inclusive Leading Gen Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveLeadJetPt->Sumw2();
    hGenInclusiveLeadJetPtWeighted = new TH1D("hGenInclusiveLeadJetPtWeighted", "Inclusive Leading Gen Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveLeadJetPtWeighted->Sumw2();
    hGenInclusiveLeadJetEtaPt = new TH2D("hGenInclusiveLeadJetEtaPt", "Inclusive Leading Gen Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveLeadJetEtaPt->Sumw2();
    hGenInclusiveLeadJetEtaPtWeighted = new TH2D("hGenInclusiveLeadJetEtaPtWeighted", "Inclusive Leading Gen Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveLeadJetEtaPtWeighted->Sumw2();

    // Inclusive subleading jet
    hGenInclusiveSubLeadJetPt = new TH1D("hGenInclusiveSubLeadJetPt", "Inclusive Sub-Leading Gen Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveSubLeadJetPt->Sumw2();
    hGenInclusiveSubLeadJetPtWeighted = new TH1D("hGenInclusiveSubLeadJetPtWeighted", "Inclusive Sub-Leading Gen Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveSubLeadJetPtWeighted->Sumw2();
    hGenInclusiveSubLeadJetEtaPt = new TH2D("hGenInclusiveSubLeadJetEtaPt", "Inclusive Sub-Leading Gen Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveSubLeadJetEtaPt->Sumw2();
    hGenInclusiveSubLeadJetEtaPtWeighted = new TH2D("hGenInclusiveSubLeadJetEtaPtWeighted", "Inclusive Sub-Leading Gen Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenInclusiveSubLeadJetEtaPtWeighted->Sumw2();

    // Inclusive dijet
    hGenInclusiveDijetDphi = new TH1D("hGenInclusiveDijetDphi", "Inclusive Gen Dijet #Delta#phi^{dijet};#Delta#phi^{dijet} (rad);Entries", fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hGenInclusiveDijetDphi->Sumw2();
    hGenInclusiveDijetDphiWeighted = new TH1D("hGenInclusiveDijetDphiWeighted", "Inclusive Gen Dijet #Delta#phi^{dijet} (weighted);#Delta#phi^{dijet} (rad);Entries", fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hGenInclusiveDijetDphiWeighted->Sumw2();
    hGenInclusiveDijetEtaPt = new TH2D("hGenInclusiveDijetEtaPt", "Inclusive Gen Dijet #eta^{dijet} vs p_{T}^{ave};#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenInclusiveDijetEtaPt->Sumw2();
    hGenInclusiveDijetEtaPtWeighted = new TH2D("hGenInclusiveDijetEtaPtWeighted", "Inclusive Gen Dijet #eta^{dijet} vs p_{T}^{ave} (weighted);#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenInclusiveDijetEtaPtWeighted->Sumw2();
    hGenInclusiveDijetEtaPtDphi = new TH3D("hGenInclusiveDijetEtaPtDphi", "Inclusive Gen Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet};#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hGenInclusiveDijetEtaPtDphi->Sumw2();
    hGenInclusiveDijetEtaPtDphiWeighted = new TH3D("hGenInclusiveDijetEtaPtDphiWeighted", "Inclusive Gen Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet} (weighted);#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hGenInclusiveDijetEtaPtDphiWeighted->Sumw2();
    hGenInclusiveDijetDetaCM = new TH1D("hGenInclusiveDijetDetaCM", "Inclusive Gen Dijet #Delta#eta_{CM};#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenInclusiveDijetDetaCM->Sumw2();
    hGenInclusiveDijetDetaCMWeighted = new TH1D("hGenInclusiveDijetDetaCMWeighted", "Inclusive Gen Dijet #Delta#eta_{CM} (weighted);#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenInclusiveDijetDetaCMWeighted->Sumw2();
    hGenInclusiveDijetDetaCMPt = new TH2D("hGenInclusiveDijetDetaCMPt", "Inclusive Gen Dijet #Delta#eta_{CM} vs p_{T}^{ave};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenInclusiveDijetDetaCMPt->Sumw2();
    hGenInclusiveDijetDetaCMPtWeighted = new TH2D("hGenInclusiveDijetDetaCMPtWeighted", "Inclusive Gen Dijet #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenInclusiveDijetDetaCMPtWeighted->Sumw2();
    hGenInclusiveDijetEtaDetaCMPt = new TH3D("hGenInclusiveDijetEtaDetaCMPt", "Inclusive Gen Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave};#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenInclusiveDijetEtaDetaCMPt->Sumw2();
    hGenInclusiveDijetEtaDetaCMPtWeighted = new TH3D("hGenInclusiveDijetEtaDetaCMPtWeighted", "Inclusive Gen Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenInclusiveDijetEtaDetaCMPtWeighted->Sumw2();
    hGenInclusiveDijetXPb = new TH1D("hGenInclusiveDijetXPb", "Inclusive Gen Dijet x_{Pb};x_{Pb};Entries", xBins, xRange[0], xRange[1]);
    hGenInclusiveDijetXPb->Sumw2();
    hGenInclusiveDijetXPbWeighted = new TH1D("hGenInclusiveDijetXPbWeighted", "Inclusive Gen Dijet x_{Pb} (weighted);x_{Pb};Entries", xBins, xRange[0], xRange[1]);
    hGenInclusiveDijetXPbWeighted->Sumw2();
    hGenInclusiveDijetXp = new TH1D("hGenInclusiveDijetXp", "Inclusive Gen Dijet x_{p};x_{p};Entries", xBins, xRange[0], xRange[1]);
    hGenInclusiveDijetXp->Sumw2();
    hGenInclusiveDijetXpWeighted = new TH1D("hGenInclusiveDijetXpWeighted", "Inclusive Gen Dijet x_{p} (weighted);x_{p};Entries", xBins, xRange[0], xRange[1]);
    hGenInclusiveDijetXpWeighted->Sumw2();
    hGenInclusiveDijetXPbOverXp = new TH1D("hGenInclusiveDijetXPbOverXp", "Inclusive Gen Dijet x_{Pb}/x_{p};x_{Pb}/x_{p};Entries", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1]);
    hGenInclusiveDijetXPbOverXp->Sumw2();
    hGenInclusiveDijetXPbOverXpWeighted = new TH1D("hGenInclusiveDijetXPbOverXpWeighted", "Inclusive Gen Dijet x_{Pb}/x_{p} (weighted);x_{Pb}/x_{p};Entries", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1]);
    hGenInclusiveDijetXPbOverXpWeighted->Sumw2();
    hGenInclusiveDijetXPbOverXpEta = new TH2D("hGenInclusiveDijetXPbOverXpEta", "Inclusive Gen Dijet x_{Pb}/x_{p} vs #eta^{dijet};x_{Pb}/x_{p};#eta^{dijet}", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenInclusiveDijetXPbOverXpEta->Sumw2();
    hGenInclusiveDijetXPbOverXpEtaWeighted = new TH2D("hGenInclusiveDijetXPbOverXpEtaWeighted", "Inclusive Gen Dijet x_{Pb}/x_{p} vs #eta^{dijet} (weighted);x_{Pb}/x_{p};#eta^{dijet}", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenInclusiveDijetXPbOverXpEtaWeighted->Sumw2();

    // Selected leading jet
    hGenSelectedLeadJetPt = new TH1D("hGenSelectedLeadJetPt", "Selected Leading Gen Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedLeadJetPt->Sumw2();
    hGenSelectedLeadJetPtWeighted = new TH1D("hGenSelectedLeadJetPtWeighted", "Selected Leading Gen Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedLeadJetPtWeighted->Sumw2();
    hGenSelectedLeadJetEtaPt = new TH2D("hGenSelectedLeadJetEtaPt", "Selected Leading Gen Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedLeadJetEtaPt->Sumw2();
    hGenSelectedLeadJetEtaPtWeighted = new TH2D("hGenSelectedLeadJetEtaPtWeighted", "Selected Leading Gen Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedLeadJetEtaPtWeighted->Sumw2();
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hGenSelectedLeadJetPtEtaPhiFlavPtHatCent", "Selected Leading Gen Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted", "Selected Leading Gen Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Selected subleading jet
    hGenSelectedSubLeadJetPt = new TH1D("hGenSelectedSubLeadJetPt", "Selected Sub-Leading Gen Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedSubLeadJetPt->Sumw2();
    hGenSelectedSubLeadJetPtWeighted = new TH1D("hGenSelectedSubLeadJetPtWeighted", "Selected Sub-Leading Gen Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedSubLeadJetPtWeighted->Sumw2();
    hGenSelectedSubLeadJetEtaPt = new TH2D("hGenSelectedSubLeadJetEtaPt", "Selected Sub-Leading Gen Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedSubLeadJetEtaPt->Sumw2();
    hGenSelectedSubLeadJetEtaPtWeighted = new TH2D("hGenSelectedSubLeadJetEtaPtWeighted", "Selected Sub-Leading Gen Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hGenSelectedSubLeadJetEtaPtWeighted->Sumw2();
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent", "Selected Sub-Leading Gen Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted", "Selected Sub-Leading Gen Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Selected dijet
    hGenSelectedDijetPt = new TH1D("hGenSelectedDijetPt", "Selected Gen Dijet p_{T}^{ave};p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetPt->Sumw2();
    hGenSelectedDijetPtWeighted = new TH1D("hGenSelectedDijetPtWeighted", "Selected Gen Dijet p_{T}^{ave} (weighted);p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetPtWeighted->Sumw2();
    hGenSelectedDijetEta = new TH1D("hGenSelectedDijetEta", "Selected Gen Dijet #eta^{dijet};#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenSelectedDijetEta->Sumw2();
    hGenSelectedDijetEtaWeighted = new TH1D("hGenSelectedDijetEtaWeighted", "Selected Gen Dijet #eta^{dijet} (weighted);#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenSelectedDijetEtaWeighted->Sumw2();
    hGenSelectedDijetEtaPt = new TH2D("hGenSelectedDijetEtaPt", "Selected Gen Dijet #eta^{dijet} vs p_{T}^{ave};#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetEtaPt->Sumw2();
    hGenSelectedDijetEtaPtWeighted = new TH2D("hGenSelectedDijetEtaPtWeighted", "Selected Gen Dijet #eta^{dijet} vs p_{T}^{ave} (weighted);#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetEtaPtWeighted->Sumw2();
    hGenSelectedDijetEtaPtDphi = new TH3D("hGenSelectedDijetEtaPtDphi", "Selected Gen Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet};#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hGenSelectedDijetEtaPtDphi->Sumw2();
    hGenSelectedDijetEtaPtDphiWeighted = new TH3D("hGenSelectedDijetEtaPtDphiWeighted", "Selected Gen Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet} (weighted);#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hGenSelectedDijetEtaPtDphiWeighted->Sumw2();
    hGenSelectedDijetDetaCM = new TH1D("hGenSelectedDijetDetaCM", "Selected Gen Dijet #Delta#eta_{CM};#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenSelectedDijetDetaCM->Sumw2();
    hGenSelectedDijetDetaCMWeighted = new TH1D("hGenSelectedDijetDetaCMWeighted", "Selected Gen Dijet #Delta#eta_{CM} (weighted);#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hGenSelectedDijetDetaCMWeighted->Sumw2();
    hGenSelectedDijetDetaCMPt = new TH2D("hGenSelectedDijetDetaCMPt", "Selected Gen Dijet #Delta#eta_{CM} vs p_{T}^{ave};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetDetaCMPt->Sumw2();
    hGenSelectedDijetDetaCMPtWeighted = new TH2D("hGenSelectedDijetDetaCMPtWeighted", "Selected Gen Dijet #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetDetaCMPtWeighted->Sumw2();
    hGenSelectedDijetEtaDetaCMPt = new TH3D("hGenSelectedDijetEtaDetaCMPt", "Selected Gen Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave};#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetEtaDetaCMPt->Sumw2();
    hGenSelectedDijetEtaDetaCMPtWeighted = new TH3D("hGenSelectedDijetEtaDetaCMPtWeighted", "Selected Gen Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hGenSelectedDijetEtaDetaCMPtWeighted->Sumw2();
    hGenSelectedDijetXPb = new TH1D("hGenSelectedDijetXPb", "Selected Gen Dijet x_{Pb};x_{Pb};Entries", xBins, xRange[0], xRange[1]);
    hGenSelectedDijetXPb->Sumw2();
    hGenSelectedDijetXPbWeighted = new TH1D("hGenSelectedDijetXPbWeighted", "Selected Gen Dijet x_{Pb} (weighted);x_{Pb};Entries", xBins, xRange[0], xRange[1]);
    hGenSelectedDijetXPbWeighted->Sumw2();
    hGenSelectedDijetXp = new TH1D("hGenSelectedDijetXp", "Selected Gen Dijet x_{p};x_{p};Entries", xBins, xRange[0], xRange[1]);
    hGenSelectedDijetXp->Sumw2();
    hGenSelectedDijetXpWeighted = new TH1D("hGenSelectedDijetXpWeighted", "Selected Gen Dijet x_{p} (weighted);x_{p};Entries", xBins, xRange[0], xRange[1]);
    hGenSelectedDijetXpWeighted->Sumw2();
    hGenSelectedDijetXPbOverXp = new TH1D("hGenSelectedDijetXPbOverXp", "Selected Gen Dijet x_{Pb}/x_{p};x_{Pb}/x_{p};Entries", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1]);
    hGenSelectedDijetXPbOverXp->Sumw2();
    hGenSelectedDijetXPbOverXpWeighted = new TH1D("hGenSelectedDijetXPbOverXpWeighted", "Selected Gen Dijet x_{Pb}/x_{p} (weighted);x_{Pb}/x_{p};Entries", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1]);
    hGenSelectedDijetXPbOverXpWeighted->Sumw2();

    //
    // Reco jets
    //

    // Initialize jetId histograms
    for (int iEtaJetId{0}; iEtaJetId<4; iEtaJetId++) {
        float low{0}, hi{0};
        if      ( iEtaJetId == 0 ) { low = {-2.4f}; hi = {2.4f}; }
        else if ( iEtaJetId == 1 ) { low = {2.4f}; hi = {2.7f}; }
        else if ( iEtaJetId == 2 ) { low = {2.7f}; hi = {3.0f}; }
        else if ( iEtaJetId == 3 ) { low = {3.f}; hi = {100.0f}; }
        hRecoInclusiveJetNHF[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetNHF_%d", iEtaJetId), Form("Neutral hadron fraction for %3.1f< #eta<=%3.1f;Neutral hadron fraction;1/N dN/dNHF", low, hi), 
                                                   fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetNHF[iEtaJetId]->Sumw2();
        hRecoInclusiveJetNEmF[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetNEmF_%d", iEtaJetId), Form("Neutral EM fraction for %3.1f< #eta<=%3.1f;Neutral EM fraction;1/N dN/dNEF", low, hi), 
                                                    fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetNEmF[iEtaJetId]->Sumw2();
        hRecoInclusiveJetNumOfConst[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetNumOfConst_%d", iEtaJetId), Form("Number of constituents for %3.1f< #eta<=%3.1f;Number of constituents;1/N dN/dNconst", low, hi), 
                                                          fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetNumOfConst[iEtaJetId]->Sumw2();
        hRecoInclusiveJetMUF[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetMUF_%d", iEtaJetId), Form("Muon fraction for %3.1f< #eta<=%3.1f;Muon fraction;1/N dN/dMUF", low, hi), 
                                                   fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetMUF[iEtaJetId]->Sumw2();
        hRecoInclusiveJetCHF[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetCHF_%d", iEtaJetId), Form("Charged hadron fraction for %3.1f< #eta<=%3.1f;Charged hadron fraction;1/N dN/dCHF", low, hi), 
                                                   fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetCHF[iEtaJetId]->Sumw2();
        hRecoInclusiveJetChargedMult[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetChargedMult_%d", iEtaJetId), Form("Charged particle multiplicity for %3.1f< #eta<=%3.1f;Charged multiplicity;1/N dN/dChMult", low, hi), 
                                                           fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetChargedMult[iEtaJetId]->Sumw2();
        hRecoInclusiveJetCEmF[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetCEmF_%d", iEtaJetId), Form("Charged EM fraction for %3.1f< #eta<=%3.1f;Charged EM fraction;1/N dN/dChEmF", low, hi), 
                                                    fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetCEmF[iEtaJetId]->Sumw2();
        hRecoInclusiveJetNumOfNeutPart[iEtaJetId] = new TH1D(Form("hRecoInclusiveJetNumOfNeutPart_%d", iEtaJetId), Form("Number of neutral particles for %3.1f< #eta<=%3.1f;Number of neutral particles;1/N dN/dNumOfNeutrals", low, hi), 
                                                             fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetNumOfNeutPart[iEtaJetId]->Sumw2();
    }
    // 0 - pt, 1 - eta, 2 - phi, 3 - dummyIter, 4 - NHF, 5 - NEF, 6 - num of constituents, 
    // 7 - MUF, 8 - CHF, 9 - charged mult, 10 - CEF, 11 - neutral mult, 12 - ptHat, 13 - centrality
    hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent = new THnSparseD("hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent", "Inclusive Reco Jet;p_{T} (GeV);#eta;#phi;dummyIter;NHF;NEF;NumOfConst;MUF;CHF;ChargedMult;CEF;NeutralMult;#hat{p_{T}} (GeV);centrality (%)", 14, bins14D_jet_PtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent, xmin14D_jet_PtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent, xmax14D_jet_PtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent);
    hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent->Sumw2();

    hRecoInclusiveJetPt = new TH1D("hRecoInclusiveJetPt", "Inclusive Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveJetPt->Sumw2();
    hRecoInclusiveJetPtWeighted = new TH1D("hRecoInclusiveJetPtWeighted", "Inclusive Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveJetPtWeighted->Sumw2();
    hRecoInclusiveJetEtaPt = new TH2D("hRecoInclusiveJetEtaPt", "Inclusive Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveJetEtaPt->Sumw2();
    hRecoInclusiveJetEtaPtWeighted = new TH2D("hRecoInclusiveJetEtaPtWeighted", "Inclusive Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveJetEtaPtWeighted->Sumw2();
    hRecoInclusiveJetPtRawPtCorrEta = new TH3D("hRecoInclusiveJetPtRawPtCorrEta", "Inclusive Reco Jet p_{T} vs raw p_{T} vs corrected p_{T};p_{T} (GeV);Raw p_{T} (GeV);Corrected p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveJetPtRawPtCorrEta->Sumw2();
    hRecoInclusiveJetPtRawPtCorrEtaWeighted = new TH3D("hRecoInclusiveJetPtRawPtCorrEtaWeighted", "Inclusive Reco Jet p_{T} vs raw p_{T} vs corrected p_{T} (weighted);p_{T} (GeV);Raw p_{T} (GeV);Corrected p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveJetPtRawPtCorrEtaWeighted->Sumw2();

    // Inclusive unmatched jets
    hRecoInclusiveUnmatchedJetPt = new TH1D("hRecoInclusiveUnmatchedJetPt", "Inclusive Unmatched Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveUnmatchedJetPt->Sumw2();
    hRecoInclusiveUnmatchedJetPtWeighted = new TH1D("hRecoInclusiveUnmatchedJetPtWeighted", "Inclusive Unmatched Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveUnmatchedJetPtWeighted->Sumw2();
    hRecoInclusiveUnmatchedJetEtaPt = new TH2D("hRecoInclusiveUnmatchedJetEtaPt", "Inclusive Unmatched Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveUnmatchedJetEtaPt->Sumw2();
    hRecoInclusiveUnmatchedJetEtaPtWeighted = new TH2D("hRecoInclusiveUnmatchedJetEtaPtWeighted", "Inclusive Unmatched Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveUnmatchedJetEtaPtWeighted->Sumw2();
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent", "Inclusive Unmatched Reco Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Unmatched Reco Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Inclusive matched jets
    hRecoInclusiveMatchedJetPt = new TH1D("hRecoInclusiveMatchedJetPt", "Inclusive Matched Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveMatchedJetPt->Sumw2();
    hRecoInclusiveMatchedJetPtWeighted = new TH1D("hRecoInclusiveMatchedJetPtWeighted", "Inclusive Matched Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveMatchedJetPtWeighted->Sumw2();
    hRecoInclusiveMatchedJetEtaPt = new TH2D("hRecoInclusiveMatchedJetEtaPt", "Inclusive Matched Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveMatchedJetEtaPt->Sumw2();
    hRecoInclusiveMatchedJetEtaPtWeighted = new TH2D("hRecoInclusiveMatchedJetEtaPtWeighted", "Inclusive Matched Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveMatchedJetEtaPtWeighted->Sumw2();
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent", "Inclusive Matched Reco Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Matched Reco Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();
    hRecoInclusiveJetJESPtEtaPhiCent = new THnSparseD("hRecoInclusiveJetJESPtEtaPhiCent", "Inclusive Reco Jet JES;p_{T} (GeV);#eta;#phi;JES;centrality", 5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
    hRecoInclusiveJetJESPtEtaPhiCent->Sumw2();
    hRecoInclusiveJetJESPtEtaPhiCentWeighted = new THnSparseD("hRecoInclusiveJetJESPtEtaPhiCentWeighted", "Inclusive Reco Jet JES (weighted);p_{T} (GeV);#eta;#phi;JES;centrality", 5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
    hRecoInclusiveJetJESPtEtaPhiCentWeighted->Sumw2();

    // Correlation between reco and ref jets for inclusive jets
    hReco2RefInclusiveJetPt = new TH2D("hReco2RefInclusiveJetPt", "Reco vs Ref Inclusive Jet p_{T};Reco Jet p_{T} (GeV);Ref Jet p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hReco2RefInclusiveJetPt->Sumw2();
    hReco2RefInclusiveJetPtWeighted = new TH2D("hReco2RefInclusiveJetPtWeighted", "Reco vs Ref Inclusive Jet p_{T} (weighted);Reco Jet p_{T} (GeV);Ref Jet p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hReco2RefInclusiveJetPtWeighted->Sumw2();
    hReco2RefInclusiveJetEta = new TH2D("hReco2RefInclusiveJetEta", "Reco vs Ref Inclusive Jet #eta;Reco Jet #eta;Ref Jet #eta", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hReco2RefInclusiveJetEta->Sumw2();
    hReco2RefInclusiveJetEtaWeighted = new TH2D("hReco2RefInclusiveJetEtaWeighted", "Reco vs Ref Inclusive Jet #eta (weighted);Reco Jet #eta;Ref Jet #eta", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hReco2RefInclusiveJetEtaWeighted->Sumw2();
    hReco2RefInclusiveJetPhi = new TH2D("hReco2RefInclusiveJetPhi", "Reco vs Ref Inclusive Jet #phi;Reco Jet #phi;Ref Jet #phi", fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1], fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hReco2RefInclusiveJetPhi->Sumw2();
    hReco2RefInclusiveJetPhiWeighted = new TH2D("hReco2RefInclusiveJetPhiWeighted", "Reco vs Ref Inclusive Jet #phi (weighted);Reco Jet #phi;Ref Jet #phi", fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1], fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hReco2RefInclusiveJetPhiWeighted->Sumw2();
    // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
    hReco2RefInclusiveJetPtEtaPhiPtHatCentrality = new THnSparseD("hReco2RefInclusiveJetPtEtaPhiPtHatCentrality", "Reco vs Ref Inclusive Jet;Reco p_{T} (GeV);Ref p_{T} (GeV);Reco #eta;Ref #eta;Reco #phi;Ref #phi;#hat{p_{T}} (GeV);centrality (%)", 8, bins8D_jet_PtEtaPhiPtHatCentrality, xmin8D_jet_PtEtaPhiPtHatCentrality, xmax8D_jet_PtEtaPhiPtHatCentrality);

    // Inclusive leading jets
    hRecoInclusiveLeadJetPt = new TH1D("hRecoInclusiveLeadJetPt", "Inclusive Leading Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveLeadJetPt->Sumw2();
    hRecoInclusiveLeadJetPtWeighted = new TH1D("hRecoInclusiveLeadJetPtWeighted", "Inclusive Leading Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveLeadJetPtWeighted->Sumw2();
    hRecoInclusiveLeadJetEtaPt = new TH2D("hRecoInclusiveLeadJetEtaPt", "Inclusive Leading Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveLeadJetEtaPt->Sumw2();
    hRecoInclusiveLeadJetEtaPtWeighted = new TH2D("hRecoInclusiveLeadJetEtaPtWeighted", "Inclusive Leading Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveLeadJetEtaPtWeighted->Sumw2();
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent", "Inclusive Leading Reco Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Leading Reco Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();
    hRecoInclusiveLeadJetJESPtEtaPhiCent = new THnSparseD("hRecoInclusiveLeadJetJESPtEtaPhiCent", "Inclusive Leading Reco Jet JES;p_{T} (GeV);#eta;#phi;JES;centrality", 5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
    hRecoInclusiveLeadJetJESPtEtaPhiCent->Sumw2();
    hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted = new THnSparseD("hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted", "Inclusive Leading Reco Jet JES (weighted);p_{T} (GeV);#eta;#phi;JES;centrality", 5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
    hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted->Sumw2();

    // Correlation between reco and ref jets for inclusive leading jets
    hReco2RefInclusiveLeadJetPt = new TH2D("hReco2RefInclusiveLeadJetPt", "Reco vs Ref Inclusive Leading Jet p_{T};Reco Jet p_{T} (GeV);Ref Jet p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hReco2RefInclusiveLeadJetPt->Sumw2();
    hReco2RefInclusiveLeadJetPtWeighted = new TH2D("hReco2RefInclusiveLeadJetPtWeighted", "Reco vs Ref Inclusive Leading Jet p_{T} (weighted);Reco Jet p_{T} (GeV);Ref Jet p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hReco2RefInclusiveLeadJetPtWeighted->Sumw2();
    hReco2RefInclusiveLeadJetEta = new TH2D("hReco2RefInclusiveLeadJetEta", "Reco vs Ref Inclusive Leading Jet #eta;Reco Jet #eta;Ref Jet #eta", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hReco2RefInclusiveLeadJetEta->Sumw2();
    hReco2RefInclusiveLeadJetEtaWeighted = new TH2D("hReco2RefInclusiveLeadJetEtaWeighted", "Reco vs Ref Inclusive Leading Jet #eta (weighted);Reco Jet #eta;Ref Jet #eta", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hReco2RefInclusiveLeadJetEtaWeighted->Sumw2();
    hReco2RefInclusiveLeadJetPhi = new TH2D("hReco2RefInclusiveLeadJetPhi", "Reco vs Ref Inclusive Leading Jet #phi;Reco Jet #phi;Ref Jet #phi", fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1], fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hReco2RefInclusiveLeadJetPhi->Sumw2();
    hReco2RefInclusiveLeadJetPhiWeighted = new TH2D("hReco2RefInclusiveLeadJetPhiWeighted", "Reco vs Ref Inclusive Leading Jet #phi (weighted);Reco Jet #phi;Ref Jet #phi", fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1], fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hReco2RefInclusiveLeadJetPhiWeighted->Sumw2();
    // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
    hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality = new THnSparseD("hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality", "Reco vs Ref Inclusive Leading Jet;Reco p_{T} (GeV);Ref p_{T} (GeV);Reco #eta;Ref #eta;Reco #phi;Ref #phi;#hat{p_{T}} (GeV);centrality (%)", 8, bins8D_jet_PtEtaPhiPtHatCentrality, xmin8D_jet_PtEtaPhiPtHatCentrality, xmax8D_jet_PtEtaPhiPtHatCentrality);
    hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality->Sumw2();

    // Inclusive subleading jets
    hRecoInclusiveSubLeadJetPt = new TH1D("hRecoInclusiveSubLeadJetPt", "Inclusive Sub-Leading Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveSubLeadJetPt->Sumw2();
    hRecoInclusiveSubLeadJetPtWeighted = new TH1D("hRecoInclusiveSubLeadJetPtWeighted", "Inclusive Sub-Leading Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveSubLeadJetPtWeighted->Sumw2();
    hRecoInclusiveSubLeadJetEtaPt = new TH2D("hRecoInclusiveSubLeadJetEtaPt", "Inclusive Sub-Leading Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveSubLeadJetEtaPt->Sumw2();
    hRecoInclusiveSubLeadJetEtaPtWeighted = new TH2D("hRecoInclusiveSubLeadJetEtaPtWeighted", "Inclusive Sub-Leading Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoInclusiveSubLeadJetEtaPtWeighted->Sumw2();
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent", "Inclusive Sub-Leading Reco Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Sub-Leading Reco Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();
    hRecoInclusiveSubLeadJetJESPtEtaPhiCent = new THnSparseD("hRecoInclusiveSubLeadJetJESPtEtaPhiCent", "Inclusive Sub-Leading Reco Jet JES;p_{T} (GeV);#eta;#phi;JES;centrality", 5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
    hRecoInclusiveSubLeadJetJESPtEtaPhiCent->Sumw2();
    hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted = new THnSparseD("hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted", "Inclusive Sub-Leading Reco Jet JES (weighted);p_{T} (GeV);#eta;#phi;JES;centrality", 5, bins5D_jes_PtEtaPhiCent, xmin5D_jes_PtEtaPhiCent, xmax5D_jes_PtEtaPhiCent);
    hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted->Sumw2();

    // Correlation between reco and ref jets for inclusive subleading jets
    hReco2RefInclusiveSubLeadJetPt = new TH2D("hReco2RefInclusiveSubLeadJetPt", "Reco vs Ref Inclusive Sub-Leading Jet p_{T};Reco Jet p_{T} (GeV);Ref Jet p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hReco2RefInclusiveSubLeadJetPt->Sumw2();
    hReco2RefInclusiveSubLeadJetPtWeighted = new TH2D("hReco2RefInclusiveSubLeadJetPtWeighted", "Reco vs Ref Inclusive Sub-Leading Jet p_{T} (weighted);Reco Jet p_{T} (GeV);Ref Jet p_{T} (GeV)", fJetPtBins, fJetPtRange[0], fJetPtRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hReco2RefInclusiveSubLeadJetPtWeighted->Sumw2();
    hReco2RefInclusiveSubLeadJetEta = new TH2D("hReco2RefInclusiveSubLeadJetEta", "Reco vs Ref Inclusive Sub-Leading Jet #eta;Reco Jet #eta;Ref Jet #eta", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hReco2RefInclusiveSubLeadJetEta->Sumw2();
    hReco2RefInclusiveSubLeadJetEtaWeighted = new TH2D("hReco2RefInclusiveSubLeadJetEtaWeighted", "Reco vs Ref Inclusive Sub-Leading Jet #eta (weighted);Reco Jet #eta;Ref Jet #eta", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1]);
    hReco2RefInclusiveSubLeadJetEtaWeighted->Sumw2();
    hReco2RefInclusiveSubLeadJetPhi = new TH2D("hReco2RefInclusiveSubLeadJetPhi", "Reco vs Ref Inclusive Sub-Leading Jet #phi;Reco Jet #phi;Ref Jet #phi", fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1], fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hReco2RefInclusiveSubLeadJetPhi->Sumw2();
    hReco2RefInclusiveSubLeadJetPhiWeighted = new TH2D("hReco2RefInclusiveSubLeadJetPhiWeighted", "Reco vs Ref Inclusive Sub-Leading Jet #phi (weighted);Reco Jet #phi;Ref Jet #phi", fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1], fJetPhiBins, fJetPhiRange[0], fJetPhiRange[1]);
    hReco2RefInclusiveSubLeadJetPhiWeighted->Sumw2();
    // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality = new THnSparseD("hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality", "Reco vs Ref Inclusive Sub-Leading Jet;Reco p_{T} (GeV);Ref p_{T} (GeV);Reco #eta;Ref #eta;Reco #phi;Ref #phi;#hat{p_{T}} (GeV);centrality (%)", 8, bins8D_jet_PtEtaPhiPtHatCentrality, xmin8D_jet_PtEtaPhiPtHatCentrality, xmax8D_jet_PtEtaPhiPtHatCentrality);
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality->Sumw2();
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted = new THnSparseD("hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted", "Reco vs Ref Inclusive Sub-Leading Jet (weighted);Reco p_{T} (GeV);Ref p_{T} (GeV);Reco #eta;Ref #eta;Reco #phi;Ref #phi;#hat{p_{T}} (GeV);centrality (%)", 8, bins8D_jet_PtEtaPhiPtHatCentrality, xmin8D_jet_PtEtaPhiPtHatCentrality, xmax8D_jet_PtEtaPhiPtHatCentrality);
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted->Sumw2();

    // Inclusive dijets
    hRecoInclusiveDijetDphi = new TH1D("hRecoInclusiveDijetDphi", "Inclusive Reco Dijet #Delta#phi^{dijet};#Delta#phi^{dijet} (rad);Entries", fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRecoInclusiveDijetDphi->Sumw2();
    hRecoInclusiveDijetDphiWeighted = new TH1D("hRecoInclusiveDijetDphiWeighted", "Inclusive Reco Dijet #Delta#phi^{dijet} (weighted);#Delta#phi^{dijet} (rad);Entries", fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRecoInclusiveDijetDphiWeighted->Sumw2();
    hRecoInclusiveDijetEtaPt = new TH2D("hRecoInclusiveDijetEtaPt", "Inclusive Reco Dijet #eta^{dijet} vs p_{T}^{ave};#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoInclusiveDijetEtaPt->Sumw2();
    hRecoInclusiveDijetEtaPtWeighted = new TH2D("hRecoInclusiveDijetEtaPtWeighted", "Inclusive Reco Dijet #eta^{dijet} vs p_{T}^{ave} (weighted);#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoInclusiveDijetEtaPtWeighted->Sumw2();
    hRecoInclusiveDijetEtaPtDphi = new TH3D("hRecoInclusiveDijetEtaPtDphi", "Inclusive Reco Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet};#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRecoInclusiveDijetEtaPtDphi->Sumw2();
    hRecoInclusiveDijetEtaPtDphiWeighted = new TH3D("hRecoInclusiveDijetEtaPtDphiWeighted", "Inclusive Reco Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet} (weighted);#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRecoInclusiveDijetEtaPtDphiWeighted->Sumw2();
    hRecoInclusiveDijetDetaCM = new TH1D("hRecoInclusiveDijetDetaCM", "Inclusive Reco Dijet #Delta#eta_{CM};#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoInclusiveDijetDetaCM->Sumw2();
    hRecoInclusiveDijetDetaCMWeighted = new TH1D("hRecoInclusiveDijetDetaCMWeighted", "Inclusive Reco Dijet #Delta#eta_{CM} (weighted);#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoInclusiveDijetDetaCMWeighted->Sumw2();
    hRecoInclusiveDijetDetaCMPt = new TH2D("hRecoInclusiveDijetDetaCMPt", "Inclusive Reco Dijet #Delta#eta_{CM} vs p_{T}^{ave};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoInclusiveDijetDetaCMPt->Sumw2();
    hRecoInclusiveDijetDetaCMPtWeighted = new TH2D("hRecoInclusiveDijetDetaCMPtWeighted", "Inclusive Reco Dijet #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoInclusiveDijetDetaCMPtWeighted->Sumw2();
    hRecoInclusiveDijetEtaDetaCMPt = new TH3D("hRecoInclusiveDijetEtaDetaCMPt", "Inclusive Reco Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave};#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoInclusiveDijetEtaDetaCMPt->Sumw2();
    hRecoInclusiveDijetEtaDetaCMPtWeighted = new TH3D("hRecoInclusiveDijetEtaDetaCMPtWeighted", "Inclusive Reco Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoInclusiveDijetEtaDetaCMPtWeighted->Sumw2();
    hRecoInclusiveDijetJESPtEtaDphiCent = new THnSparseD("hRecoInclusiveDijetJESPtEtaDphiCent", "Inclusive Reco Dijet JES;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi^{dijet};JES;centrality", 5, bins5D_jes_PtEtaDphiCent, xmin5D_jes_PtEtaDphiCent, xmax5D_jes_PtEtaDphiCent);
    hRecoInclusiveDijetJESPtEtaDphiCent->Sumw2();
    hRecoInclusiveDijetJESPtEtaDphiCentWeighted = new THnSparseD("hRecoInclusiveDijetJESPtEtaDphiCentWeighted", "Inclusive Reco Dijet JES (weighted);p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi^{dijet};JES;centrality", 5, bins5D_jes_PtEtaDphiCent, xmin5D_jes_PtEtaDphiCent, xmax5D_jes_PtEtaDphiCent);
    hRecoInclusiveDijetJESPtEtaDphiCentWeighted->Sumw2();

    // Selected leading jets
    hRecoSelectedLeadJetPt = new TH1D("hRecoSelectedLeadJetPt", "Selected Leading Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedLeadJetPt->Sumw2();
    hRecoSelectedLeadJetPtWeighted = new TH1D("hRecoSelectedLeadJetPtWeighted", "Selected Leading Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedLeadJetPtWeighted->Sumw2();
    hRecoSelectedLeadJetEtaPt = new TH2D("hRecoSelectedLeadJetEtaPt", "Selected Leading Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedLeadJetEtaPt->Sumw2();
    hRecoSelectedLeadJetEtaPtWeighted = new TH2D("hRecoSelectedLeadJetEtaPtWeighted", "Selected Leading Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedLeadJetEtaPtWeighted->Sumw2();
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent", "Selected Leading Reco Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted", "Selected Leading Reco Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Selected subleading jets
    hRecoSelectedSubLeadJetPt = new TH1D("hRecoSelectedSubLeadJetPt", "Selected Sub-Leading Reco Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedSubLeadJetPt->Sumw2();
    hRecoSelectedSubLeadJetPtWeighted = new TH1D("hRecoSelectedSubLeadJetPtWeighted", "Selected Sub-Leading Reco Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedSubLeadJetPtWeighted->Sumw2();
    hRecoSelectedSubLeadJetEtaPt = new TH2D("hRecoSelectedSubLeadJetEtaPt", "Selected Sub-Leading Reco Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedSubLeadJetEtaPt->Sumw2();
    hRecoSelectedSubLeadJetEtaPtWeighted = new TH2D("hRecoSelectedSubLeadJetEtaPtWeighted", "Selected Sub-Leading Reco Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRecoSelectedSubLeadJetEtaPtWeighted->Sumw2();
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent", "Selected Sub-Leading Reco Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted", "Selected Sub-Leading Reco Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Selected dijets
    hRecoSelectedDijetPt = new TH1D("hRecoSelectedDijetPt", "Selected Reco Dijet p_{T}^{ave};p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetPt->Sumw2();
    hRecoSelectedDijetPtWeighted = new TH1D("hRecoSelectedDijetPtWeighted", "Selected Reco Dijet p_{T}^{ave} (weighted);p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetPtWeighted->Sumw2();
    hRecoSelectedDijetEta = new TH1D("hRecoSelectedDijetEta", "Selected Reco Dijet #eta^{dijet};#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoSelectedDijetEta->Sumw2();
    hRecoSelectedDijetEtaWeighted = new TH1D("hRecoSelectedDijetEtaWeighted", "Selected Reco Dijet #eta^{dijet} (weighted);#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoSelectedDijetEtaWeighted->Sumw2();
    hRecoSelectedDijetEtaPt = new TH2D("hRecoSelectedDijetEtaPt", "Selected Reco Dijet #eta^{dijet} vs p_{T}^{ave};#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetEtaPt->Sumw2();
    hRecoSelectedDijetEtaPtWeighted = new TH2D("hRecoSelectedDijetEtaPtWeighted", "Selected Reco Dijet #eta^{dijet} vs p_{T}^{ave} (weighted);#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetEtaPtWeighted->Sumw2();
    hRecoSelectedDijetEtaPtDphi = new TH3D("hRecoSelectedDijetEtaPtDphi", "Selected Reco Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet};#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRecoSelectedDijetEtaPtDphi->Sumw2();
    hRecoSelectedDijetEtaPtDphiWeighted = new TH3D("hRecoSelectedDijetEtaPtDphiWeighted", "Selected Reco Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet} (weighted);#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRecoSelectedDijetEtaPtDphiWeighted->Sumw2();
    hRecoSelectedDijetDetaCM = new TH1D("hRecoSelectedDijetDetaCM", "Selected Reco Dijet #Delta#eta_{CM};#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoSelectedDijetDetaCM->Sumw2();
    hRecoSelectedDijetDetaCMWeighted = new TH1D("hRecoSelectedDijetDetaCMWeighted", "Selected Reco Dijet #Delta#eta_{CM} (weighted);#Delta#eta_{CM};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoSelectedDijetDetaCMWeighted->Sumw2();
    hRecoSelectedDijetDetaCMPt = new TH2D("hRecoSelectedDijetDetaCMPt", "Selected Reco Dijet #Delta#eta_{CM} vs p_{T}^{ave};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetDetaCMPt->Sumw2();
    hRecoSelectedDijetDetaCMPtWeighted = new TH2D("hRecoSelectedDijetDetaCMPtWeighted", "Selected Reco Dijet #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetDetaCMPtWeighted->Sumw2();
    hRecoSelectedDijetEtaDetaCMPt = new TH3D("hRecoSelectedDijetEtaDetaCMPt", "Selected Reco Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave};#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetEtaDetaCMPt->Sumw2();
    hRecoSelectedDijetEtaDetaCMPtWeighted = new TH3D("hRecoSelectedDijetEtaDetaCMPtWeighted", "Selected Reco Dijet #eta^{dijet} vs #Delta#eta_{CM} vs p_{T}^{ave} (weighted);#eta^{dijet};#Delta#eta_{CM};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRecoSelectedDijetEtaDetaCMPtWeighted->Sumw2();

    // Correlation between reco and ref for selected dijets
    // 0 - reco dijetpt, 1 - reco dijeteta, 
    // 2 - reco leading jet pt, 3 - reco leading jet eta,
    // 4 - reco subleading jet pt, 5 - reco subleading jet eta, 
    // 6 - ref dijet pt, 7 - ref dijet eta,
    // 8 - ref leading jet pt, 9 - ref leading jet eta,
    // 10 - ref subleading jet pt, 11 - ref subleading jet eta
    hReco2RefSelectedDijetPtEtaFull = new THnSparseD("hReco2RefSelectedDijetPtEtaFull", "Reco vs Ref Selected Dijet;Reco dijet p_{T}^{ave} (GeV);Reco dijet #eta^{dijet};Reco leading jet p_{T} (GeV);Reco leading jet #eta;Reco subleading jet p_{T} (GeV);Reco subleading jet #eta;Ref dijet p_{T}^{ave} (GeV);Ref dijet #eta^{dijet};Ref leading jet p_{T} (GeV);Ref leading jet #eta;Ref subleading jet p_{T} (GeV);Ref subleading jet #eta", 12, bins12D_jet_PtEtaFull, xmin12D_jet_PtEtaFull, xmax12D_jet_PtEtaFull);
    hReco2RefSelectedDijetPtEtaFull->Sumw2();
    hReco2RefSelectedDijetPtEtaFullWeighted = new THnSparseD("hReco2RefSelectedDijetPtEtaFullWeighted", "Reco vs Ref Selected Dijet (weighted);Reco dijet p_{T}^{ave} (GeV);Reco dijet #eta^{dijet};Reco leading jet p_{T} (GeV);Reco leading jet #eta;Reco subleading jet p_{T} (GeV);Reco subleading jet #eta;Ref dijet p_{T}^{ave} (GeV);Ref dijet #eta^{dijet};Ref leading jet p_{T} (GeV);Ref leading jet #eta;Ref subleading jet p_{T} (GeV);Ref subleading jet #eta", 12, bins12D_jet_PtEtaFull, xmin12D_jet_PtEtaFull, xmax12D_jet_PtEtaFull);
    hReco2RefSelectedDijetPtEtaFullWeighted->Sumw2();


    //
    // Ref jets
    //

    // Inclusive jets
    hRefInclusiveJetPt = new TH1D("hRefInclusiveJetPt", "Inclusive Ref Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveJetPt->Sumw2();
    hRefInclusiveJetPtWeighted = new TH1D("hRefInclusiveJetPtWeighted", "Inclusive Ref Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveJetPtWeighted->Sumw2();
    hRefInclusiveJetEtaPt = new TH2D("hRefInclusiveJetEtaPt", "Inclusive Ref Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveJetEtaPt->Sumw2();
    hRefInclusiveJetEtaPtWeighted = new TH2D("hRefInclusiveJetEtaPtWeighted", "Inclusive Ref Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveJetEtaPtWeighted->Sumw2();
    hRefInclusiveJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRefInclusiveJetPtEtaPhiFlavPtHatCent", "Inclusive Ref Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefInclusiveJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Ref Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Inclusive leading jet
    hRefInclusiveLeadJetPt = new TH1D("hRefInclusiveLeadJetPt", "Inclusive Leading Ref Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveLeadJetPt->Sumw2();
    hRefInclusiveLeadJetPtWeighted = new TH1D("hRefInclusiveLeadJetPtWeighted", "Inclusive Leading Ref Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveLeadJetPtWeighted->Sumw2();
    hRefInclusiveLeadJetEtaPt = new TH2D("hRefInclusiveLeadJetEtaPt", "Inclusive Leading Ref Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveLeadJetEtaPt->Sumw2();
    hRefInclusiveLeadJetEtaPtWeighted = new TH2D("hRefInclusiveLeadJetEtaPtWeighted", "Inclusive Leading Ref Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveLeadJetEtaPtWeighted->Sumw2();
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent", "Inclusive Leading Ref Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Leading Ref Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Inclusive subleading jet
    hRefInclusiveSubLeadJetPt = new TH1D("hRefInclusiveSubLeadJetPt", "Inclusive Sub-Leading Ref Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveSubLeadJetPt->Sumw2();
    hRefInclusiveSubLeadJetPtWeighted = new TH1D("hRefInclusiveSubLeadJetPtWeighted", "Inclusive Sub-Leading Ref Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveSubLeadJetPtWeighted->Sumw2();
    hRefInclusiveSubLeadJetEtaPt = new TH2D("hRefInclusiveSubLeadJetEtaPt", "Inclusive Sub-Leading Ref Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveSubLeadJetEtaPt->Sumw2();
    hRefInclusiveSubLeadJetEtaPtWeighted = new TH2D("hRefInclusiveSubLeadJetEtaPtWeighted", "Inclusive Sub-Leading Ref Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefInclusiveSubLeadJetEtaPtWeighted->Sumw2();
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent", "Inclusive Sub-Leading Ref Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted", "Inclusive Sub-Leading Ref Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Inclusive dijets
    hRefInclusiveDijetDphi = new TH1D("hRefInclusiveDijetDphi", "Inclusive Ref Dijet #Delta#phi^{dijet};#Delta#phi^{dijet} (rad);Entries", fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRefInclusiveDijetDphi->Sumw2();
    hRefInclusiveDijetDphiWeighted = new TH1D("hRefInclusiveDijetDphiWeighted", "Inclusive Ref Dijet #Delta#phi^{dijet} (weighted);#Delta#phi^{dijet} (rad);Entries", fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRefInclusiveDijetDphiWeighted->Sumw2();
    hRefInclusiveDijetPt = new TH1D("hRefInclusiveDijetPt", "Inclusive Ref Dijet p_{T}^{ave};p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefInclusiveDijetPt->Sumw2();
    hRefInclusiveDijetPtWeighted = new TH1D("hRefInclusiveDijetPtWeighted", "Inclusive Ref Dijet p_{T}^{ave} (weighted);p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefInclusiveDijetPtWeighted->Sumw2();
    hRefInclusiveDijetEta = new TH1D("hRefInclusiveDijetEta", "Inclusive Ref Dijet #eta^{dijet};#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRefInclusiveDijetEta->Sumw2();
    hRefInclusiveDijetEtaWeighted = new TH1D("hRefInclusiveDijetEtaWeighted", "Inclusive Ref Dijet #eta^{dijet} (weighted);#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRefInclusiveDijetEtaWeighted->Sumw2();
    hRefInclusiveDijetEtaPt = new TH2D("hRefInclusiveDijetEtaPt", "Inclusive Ref Dijet #eta^{dijet} vs p_{T}^{ave};#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefInclusiveDijetEtaPt->Sumw2();
    hRefInclusiveDijetEtaPtWeighted = new TH2D("hRefInclusiveDijetEtaPtWeighted", "Inclusive Ref Dijet #eta^{dijet} vs p_{T}^{ave} (weighted);#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefInclusiveDijetEtaPtWeighted->Sumw2();
    hRefInclusiveDijetEtaPtDphi = new TH3D("hRefInclusiveDijetEtaPtDphi", "Inclusive Ref Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet};#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRefInclusiveDijetEtaPtDphi->Sumw2();
    hRefInclusiveDijetEtaPtDphiWeighted = new TH3D("hRefInclusiveDijetEtaPtDphiWeighted", "Inclusive Ref Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet} (weighted);#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRefInclusiveDijetEtaPtDphiWeighted->Sumw2();

    // Selected leading jet
    hRefSelectedLeadJetPt = new TH1D("hRefSelectedLeadJetPt", "Selected Leading Ref Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedLeadJetPt->Sumw2();
    hRefSelectedLeadJetPtWeighted = new TH1D("hRefSelectedLeadJetPtWeighted", "Selected Leading Ref Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedLeadJetPtWeighted->Sumw2();
    hRefSelectedLeadJetEtaPt = new TH2D("hRefSelectedLeadJetEtaPt", "Selected Leading Ref Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedLeadJetEtaPt->Sumw2();
    hRefSelectedLeadJetEtaPtWeighted = new TH2D("hRefSelectedLeadJetEtaPtWeighted", "Selected Leading Ref Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedLeadJetEtaPtWeighted->Sumw2();
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRefSelectedLeadJetPtEtaPhiFlavPtHatCent", "Selected Leading Ref Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted", "Selected Leading Ref Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Selected subleading jet
    hRefSelectedSubLeadJetPt = new TH1D("hRefSelectedSubLeadJetPt", "Selected Sub-Leading Ref Jet p_{T};p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedSubLeadJetPt->Sumw2();
    hRefSelectedSubLeadJetPtWeighted = new TH1D("hRefSelectedSubLeadJetPtWeighted", "Selected Sub-Leading Ref Jet p_{T} (weighted);p_{T} (GeV);Entries", fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedSubLeadJetPtWeighted->Sumw2();
    hRefSelectedSubLeadJetEtaPt = new TH2D("hRefSelectedSubLeadJetEtaPt", "Selected Sub-Leading Ref Jet #eta vs p_{T};#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedSubLeadJetEtaPt->Sumw2();
    hRefSelectedSubLeadJetEtaPtWeighted = new TH2D("hRefSelectedSubLeadJetEtaPtWeighted", "Selected Sub-Leading Ref Jet #eta vs p_{T} (weighted);#eta;p_{T} (GeV)", fJetEtaBins, fJetEtaRange[0], fJetEtaRange[1], fJetPtBins, fJetPtRange[0], fJetPtRange[1]);
    hRefSelectedSubLeadJetEtaPtWeighted->Sumw2();
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent = new THnSparseD("hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent", "Selected Sub-Leading Ref Jet;p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Sumw2();
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted = new THnSparseD("hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted", "Selected Sub-Leading Ref Jet (weighted);p_{T} (GeV);#eta;#phi;flavorForB;#hat{p_{T}} (GeV);centrality", 6, bins6D_jet_PtEtaPhiFlavPtHatCent, xmin6D_jet_PtEtaPhiFlavPtHatCent, xmax6D_jet_PtEtaPhiFlavPtHatCent);
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Sumw2();

    // Selected dijets
    hRefSelectedDijetPt = new TH1D("hRefSelectedDijetPt", "Selected Ref Dijet p_{T}^{ave};p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefSelectedDijetPt->Sumw2();
    hRefSelectedDijetPtWeighted = new TH1D("hRefSelectedDijetPtWeighted", "Selected Ref Dijet p_{T}^{ave} (weighted);p_{T}^{ave} (GeV);Entries", fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefSelectedDijetPtWeighted->Sumw2();
    hRefSelectedDijetEta = new TH1D("hRefSelectedDijetEta", "Selected Ref Dijet #eta^{dijet};#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRefSelectedDijetEta->Sumw2();
    hRefSelectedDijetEtaWeighted = new TH1D("hRefSelectedDijetEtaWeighted", "Selected Ref Dijet #eta^{dijet} (weighted);#eta^{dijet};Entries", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRefSelectedDijetEtaWeighted->Sumw2();
    hRefSelectedDijetEtaPt = new TH2D("hRefSelectedDijetEtaPt", "Selected Ref Dijet #eta^{dijet} vs p_{T}^{ave};#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefSelectedDijetEtaPt->Sumw2();
    hRefSelectedDijetEtaPtWeighted = new TH2D("hRefSelectedDijetEtaPtWeighted", "Selected Ref Dijet #eta^{dijet} vs p_{T}^{ave} (weighted);#eta^{dijet};p_{T}^{ave} (GeV)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
    hRefSelectedDijetEtaPtWeighted->Sumw2();
    hRefSelectedDijetEtaPtDphi = new TH3D("hRefSelectedDijetEtaPtDphi", "Selected Ref Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet};#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRefSelectedDijetEtaPtDphi->Sumw2();
    hRefSelectedDijetEtaPtDphiWeighted = new TH3D("hRefSelectedDijetEtaPtDphiWeighted", "Selected Ref Dijet #eta^{dijet} vs p_{T}^{ave} vs #Delta#phi^{dijet} (weighted);#eta^{dijet};p_{T}^{ave} (GeV);#Delta#phi^{dijet} (rad)", fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1], fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1], fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1]);
    hRefSelectedDijetEtaPtDphiWeighted->Sumw2();

}

//________________
void HistoManagerJetESR::writeOutput() {

    //
    // Event information
    //
    hVz->Write();
    hVzCentWeighted->Write();
    hVzPtHatWeighted->Write();
    hVzWeighted->Write();
    hHiBin->Write();
    hHiBinPtHatWeighted->Write();
    hHiBinWeighted->Write();
    hPtHat->Write();
    hPtHatPtHatWeighted->Write();
    hPtHatWeighted->Write();
    hPtHatWeight->Write();
    hPtHatWeightWeighted->Write();
    hCentrality->Write();
    hCentralityPtHatWeighted->Write();
    hCentralityWeighted->Write();
    hVzPtHatCent->Write();
    hVzPtHatCentPtHatWeighted->Write();
    hVzPtHatCentWeighted->Write();

    //
    // Gen jets
    //
    hGenInclusiveJetPt->Write();
    hGenInclusiveJetPtWeighted->Write();
    hGenInclusiveJetEtaPt->Write();
    hGenInclusiveJetEtaPtWeighted->Write();
    hGenInclusiveJetPtEtaPhiFlavPtHatCent->Write();
    hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hGenInclusiveLeadJetPt->Write();
    hGenInclusiveLeadJetPtWeighted->Write();
    hGenInclusiveLeadJetEtaPt->Write();
    hGenInclusiveLeadJetEtaPtWeighted->Write();
    hGenInclusiveLeadJetPtEtaPhiFlavPtHatCent->Write();
    hGenInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hGenInclusiveSubLeadJetPt->Write();
    hGenInclusiveSubLeadJetPtWeighted->Write();
    hGenInclusiveSubLeadJetEtaPt->Write();
    hGenInclusiveSubLeadJetEtaPtWeighted->Write();
    hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Write();
    hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hGenInclusiveDijetDphi->Write();
    hGenInclusiveDijetDphiWeighted->Write();
    hGenInclusiveDijetEtaPt->Write();
    hGenInclusiveDijetEtaPtWeighted->Write();
    hGenInclusiveDijetEtaPtDphi->Write();
    hGenInclusiveDijetEtaPtDphiWeighted->Write();
    hGenInclusiveDijetDetaCM->Write();
    hGenInclusiveDijetDetaCMWeighted->Write();
    hGenInclusiveDijetDetaCMPt->Write();
    hGenInclusiveDijetDetaCMPtWeighted->Write();
    hGenInclusiveDijetEtaDetaCMPt->Write();
    hGenInclusiveDijetEtaDetaCMPtWeighted->Write();
    hGenInclusiveDijetXPb->Write();
    hGenInclusiveDijetXPbWeighted->Write();
    hGenInclusiveDijetXp->Write();
    hGenInclusiveDijetXpWeighted->Write();
    hGenInclusiveDijetXPbOverXp->Write();
    hGenInclusiveDijetXPbOverXpWeighted->Write();
    hGenInclusiveDijetXPbOverXpEta->Write();
    hGenInclusiveDijetXPbOverXpEtaWeighted->Write();

    hGenSelectedLeadJetPt->Write();
    hGenSelectedLeadJetPtWeighted->Write();
    hGenSelectedLeadJetEtaPt->Write();
    hGenSelectedLeadJetEtaPtWeighted->Write();
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCent->Write();
    hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hGenSelectedSubLeadJetPt->Write();
    hGenSelectedSubLeadJetPtWeighted->Write();
    hGenSelectedSubLeadJetEtaPt->Write();
    hGenSelectedSubLeadJetEtaPtWeighted->Write();
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Write();
    hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hGenSelectedDijetPt->Write();
    hGenSelectedDijetPtWeighted->Write();
    hGenSelectedDijetEta->Write();
    hGenSelectedDijetEtaWeighted->Write();
    hGenSelectedDijetEtaPt->Write();
    hGenSelectedDijetEtaPtWeighted->Write();
    hGenSelectedDijetEtaPtDphi->Write();
    hGenSelectedDijetEtaPtDphiWeighted->Write();
    hGenSelectedDijetDetaCM->Write();
    hGenSelectedDijetDetaCMWeighted->Write();
    hGenSelectedDijetDetaCMPt->Write();
    hGenSelectedDijetDetaCMPtWeighted->Write();
    hGenSelectedDijetEtaDetaCMPt->Write();
    hGenSelectedDijetEtaDetaCMPtWeighted->Write();
    hGenSelectedDijetXPb->Write();
    hGenSelectedDijetXPbWeighted->Write();
    hGenSelectedDijetXp->Write();
    hGenSelectedDijetXpWeighted->Write();
    hGenSelectedDijetXPbOverXp->Write();
    hGenSelectedDijetXPbOverXpWeighted->Write();


    //
    // Reco jets
    //

    for (int iEtaJetId{0}; iEtaJetId<4; iEtaJetId++) {
        if (hRecoInclusiveJetNHF[iEtaJetId]) hRecoInclusiveJetNHF[iEtaJetId]->Write();
        if (hRecoInclusiveJetNEmF[iEtaJetId]) hRecoInclusiveJetNEmF[iEtaJetId]->Write();
        if (hRecoInclusiveJetNumOfConst[iEtaJetId]) hRecoInclusiveJetNumOfConst[iEtaJetId]->Write();
        if (hRecoInclusiveJetMUF[iEtaJetId]) hRecoInclusiveJetMUF[iEtaJetId]->Write();
        if (hRecoInclusiveJetCHF[iEtaJetId]) hRecoInclusiveJetCHF[iEtaJetId]->Write();
        if (hRecoInclusiveJetChargedMult[iEtaJetId]) hRecoInclusiveJetChargedMult[iEtaJetId]->Write();
        if (hRecoInclusiveJetCEmF[iEtaJetId]) hRecoInclusiveJetCEmF[iEtaJetId]->Write();
        if (hRecoInclusiveJetNumOfNeutPart[iEtaJetId]) hRecoInclusiveJetNumOfNeutPart[iEtaJetId]->Write();
    }

    hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent->Write();

    hRecoInclusiveJetPt->Write();
    hRecoInclusiveJetPtWeighted->Write();
    hRecoInclusiveJetEtaPt->Write();
    hRecoInclusiveJetEtaPtWeighted->Write();
    hRecoInclusiveJetPtRawPtCorrEta->Write();
    hRecoInclusiveJetPtRawPtCorrEtaWeighted->Write();

    hRecoInclusiveUnmatchedJetPt->Write();
    hRecoInclusiveUnmatchedJetPtWeighted->Write();
    hRecoInclusiveUnmatchedJetEtaPt->Write();
    hRecoInclusiveUnmatchedJetEtaPtWeighted->Write();
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent->Write();
    hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRecoInclusiveMatchedJetPt->Write();
    hRecoInclusiveMatchedJetPtWeighted->Write();
    hRecoInclusiveMatchedJetEtaPt->Write();
    hRecoInclusiveMatchedJetEtaPtWeighted->Write();
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent->Write();
    hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted->Write();
    hRecoInclusiveJetJESPtEtaPhiCent->Write();
    hRecoInclusiveJetJESPtEtaPhiCentWeighted->Write();

    hReco2RefInclusiveJetPt->Write();
    hReco2RefInclusiveJetPtWeighted->Write();
    hReco2RefInclusiveJetEta->Write();
    hReco2RefInclusiveJetEtaWeighted->Write();
    hReco2RefInclusiveJetPhi->Write();
    hReco2RefInclusiveJetPhiWeighted->Write();
    hReco2RefInclusiveJetPtEtaPhiPtHatCentrality->Write();

    hRecoInclusiveLeadJetPt->Write();
    hRecoInclusiveLeadJetPtWeighted->Write();
    hRecoInclusiveLeadJetEtaPt->Write();
    hRecoInclusiveLeadJetEtaPtWeighted->Write();
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();
    hRecoInclusiveLeadJetJESPtEtaPhiCent->Write();
    hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted->Write();

    hReco2RefInclusiveLeadJetPt->Write();
    hReco2RefInclusiveLeadJetPtWeighted->Write();
    hReco2RefInclusiveLeadJetEta->Write();
    hReco2RefInclusiveLeadJetEtaWeighted->Write();
    hReco2RefInclusiveLeadJetPhi->Write();
    hReco2RefInclusiveLeadJetPhiWeighted->Write();
    hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality->Write();

    hRecoInclusiveSubLeadJetPt->Write();
    hRecoInclusiveSubLeadJetPtWeighted->Write();
    hRecoInclusiveSubLeadJetEtaPt->Write();
    hRecoInclusiveSubLeadJetEtaPtWeighted->Write();
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();
    hRecoInclusiveSubLeadJetJESPtEtaPhiCent->Write();
    hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted->Write();

    hReco2RefInclusiveSubLeadJetPt->Write();
    hReco2RefInclusiveSubLeadJetPtWeighted->Write();
    hReco2RefInclusiveSubLeadJetEta->Write();
    hReco2RefInclusiveSubLeadJetEtaWeighted->Write();
    hReco2RefInclusiveSubLeadJetPhi->Write();
    hReco2RefInclusiveSubLeadJetPhiWeighted->Write();
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality->Write();
    hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted->Write();

    hRecoInclusiveDijetDphi->Write();
    hRecoInclusiveDijetDphiWeighted->Write();
    hRecoInclusiveDijetEtaPt->Write();
    hRecoInclusiveDijetEtaPtWeighted->Write();
    hRecoInclusiveDijetEtaPtDphi->Write();
    hRecoInclusiveDijetEtaPtDphiWeighted->Write();
    hRecoInclusiveDijetDetaCM->Write();
    hRecoInclusiveDijetDetaCMWeighted->Write();
    hRecoInclusiveDijetDetaCMPt->Write();
    hRecoInclusiveDijetDetaCMPtWeighted->Write();
    hRecoInclusiveDijetEtaDetaCMPt->Write();
    hRecoInclusiveDijetEtaDetaCMPtWeighted->Write();
    hRecoInclusiveDijetJESPtEtaDphiCent->Write();
    hRecoInclusiveDijetJESPtEtaDphiCentWeighted->Write();

    hRecoSelectedLeadJetPt->Write();
    hRecoSelectedLeadJetPtWeighted->Write();
    hRecoSelectedLeadJetEtaPt->Write();
    hRecoSelectedLeadJetEtaPtWeighted->Write();
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRecoSelectedSubLeadJetPt->Write();
    hRecoSelectedSubLeadJetPtWeighted->Write();
    hRecoSelectedSubLeadJetEtaPt->Write();
    hRecoSelectedSubLeadJetEtaPtWeighted->Write();
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRecoSelectedDijetPt->Write();
    hRecoSelectedDijetPtWeighted->Write();
    hRecoSelectedDijetEta->Write();
    hRecoSelectedDijetEtaWeighted->Write();
    hRecoSelectedDijetEtaPt->Write();
    hRecoSelectedDijetEtaPtWeighted->Write();
    hRecoSelectedDijetEtaPtDphi->Write();
    hRecoSelectedDijetEtaPtDphiWeighted->Write();
    hRecoSelectedDijetDetaCM->Write();
    hRecoSelectedDijetDetaCMWeighted->Write();
    hRecoSelectedDijetDetaCMPt->Write();
    hRecoSelectedDijetDetaCMPtWeighted->Write();
    hRecoSelectedDijetEtaDetaCMPt->Write();
    hRecoSelectedDijetEtaDetaCMPtWeighted->Write();

    hReco2RefSelectedDijetPtEtaFull->Write();
    hReco2RefSelectedDijetPtEtaFullWeighted->Write();

    //
    // Ref jets
    //

    hRefInclusiveJetPt->Write();
    hRefInclusiveJetPtWeighted->Write();
    hRefInclusiveJetEtaPt->Write();
    hRefInclusiveJetEtaPtWeighted->Write();
    hRefInclusiveJetPtEtaPhiFlavPtHatCent->Write();
    hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRefInclusiveLeadJetPt->Write();
    hRefInclusiveLeadJetPtWeighted->Write();
    hRefInclusiveLeadJetEtaPt->Write();
    hRefInclusiveLeadJetEtaPtWeighted->Write();
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRefInclusiveSubLeadJetPt->Write();
    hRefInclusiveSubLeadJetPtWeighted->Write();
    hRefInclusiveSubLeadJetEtaPt->Write();
    hRefInclusiveSubLeadJetEtaPtWeighted->Write();
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();


    hRefInclusiveDijetDphi->Write();
    hRefInclusiveDijetDphiWeighted->Write();
    hRefInclusiveDijetPt->Write();
    hRefInclusiveDijetPtWeighted->Write();
    hRefInclusiveDijetEta->Write();
    hRefInclusiveDijetEtaWeighted->Write();
    hRefInclusiveDijetEtaPt->Write();
    hRefInclusiveDijetEtaPtWeighted->Write();
    hRefInclusiveDijetEtaPtDphi->Write();
    hRefInclusiveDijetEtaPtDphiWeighted->Write();

    hRefSelectedLeadJetPt->Write();
    hRefSelectedLeadJetPtWeighted->Write();
    hRefSelectedLeadJetEtaPt->Write();
    hRefSelectedLeadJetEtaPtWeighted->Write();
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRefSelectedSubLeadJetPt->Write();
    hRefSelectedSubLeadJetPtWeighted->Write();
    hRefSelectedSubLeadJetEtaPt->Write();
    hRefSelectedSubLeadJetEtaPtWeighted->Write();
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Write();
    hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Write();

    hRefSelectedDijetPt->Write();
    hRefSelectedDijetPtWeighted->Write();
    hRefSelectedDijetEta->Write();
    hRefSelectedDijetEtaWeighted->Write();
    hRefSelectedDijetEtaPt->Write();
    hRefSelectedDijetEtaPtWeighted->Write();
    hRefSelectedDijetEtaPtDphi->Write();
    hRefSelectedDijetEtaPtDphiWeighted->Write();
}