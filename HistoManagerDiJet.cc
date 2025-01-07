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
  fPtBins{150}, fPtRange{20., 1520.}, 
  fEtaBins{52}, fEtaRange{-5.2, 5.2},
  fPhiBins{32}, fPhiRange{-TMath::Pi(), TMath::Pi()},
  fDijetPtBins{194}, fDijetPtRange{30., 1000.},
  fDijetEtaBins{48}, fDijetEtaRange{-4.8, 4.8},
  fDijetDphiBins{16}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},
  fPtHatBins{100}, fPtHatRange{15., 1015.},
  fFracBins{100}, fFracRange{0., 1.},
  fMultBins{32}, fMultRange{-0.5, 31.5},
  
  hVz{nullptr}, hVzWeighted{nullptr}, hMult{nullptr},
  hHiBin{nullptr}, hHiBinWeighted{nullptr},
  hPtHat{nullptr}, hPtHatWeighted{nullptr}, hPtHatWeight{nullptr},
//   hCentrality{nullptr}, hCentralityWeighted{nullptr},
  hVzPtHat{nullptr}, hVzPtHatWeighted{nullptr},

  hNHF{nullptr}, hNEmF{nullptr}, hNumOfConst{nullptr}, hMUF{nullptr},
  hCHF{nullptr}, hChargedMult{nullptr}, hCEmF{nullptr}, hNumOfNeutPart{nullptr},

  // Gen jets
  hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi{nullptr},
  hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted{nullptr},
  hGenInclusiveJetPt{nullptr},
  hGenInclusiveJetPtEta{nullptr},
  hGenPtLeadPtSublead{nullptr},
  hGenEtaLeadEtaSublead{nullptr},
  hGenEtaCMLeadEtaCMSublead{nullptr},
  hGenPtLeadPtSubleadMcReweight{nullptr},
  hGenEtaLeadEtaSubleadMcReweight{nullptr},
  hGenDijetEta{nullptr},
  hGenDijetEta1D{nullptr},
  hGenDijetEta1DOldPt{nullptr},
  hGenDijetEta1DOldPtBinning{nullptr},
  hGenDijetPtEtaDphi{nullptr},
  hGenDijetPtEtaDphiWeighted{nullptr},
  hGenDijetEtaCM{nullptr},
  hGenDijetPtEtaDphiCM{nullptr},
  hGenDijetPtEtaDphiCMWeighted{nullptr},
  hGenDijetPtEtaForward{nullptr},
  hGenDijetPtEtaBackward{nullptr},
  hGenDijetPtEtaCMForward{nullptr},
  hGenDijetPtEtaCMBackward{nullptr},
  hGenDijetPtEtaForwardWeighted{nullptr},
  hGenDijetPtEtaBackwardWeighted{nullptr},
  hGenDijetPtEtaCMForwardWeighted{nullptr},
  hGenDijetPtEtaCMBackwardWeighted{nullptr},

  hGenGoodInclusiveJetEtaLabFrame{nullptr},
  hGenGoodInclusiveJetEtaCMFrame{nullptr},

  // Reco jets (single)
  hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
  hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
  hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
  hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
  hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
  hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
  hJESInclusiveJetPtEtaPhi{nullptr},
  hJESInclusiveJetPtEtaPhiWeighted{nullptr},
  hRecoMatchedPtEta{nullptr},

  // Dijets (experiment)
  hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi{nullptr},
  hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted{nullptr},
  hRecoInclusiveJetPt{nullptr},
  hRecoPtLeadPtSublead{nullptr},
  hRecoEtaLeadEtaSublead{nullptr},
  hRecoEtaCMLeadEtaCMSublead{nullptr},
  hRecoPtLeadPtSubleadMcReweight{nullptr},
  hRecoEtaLeadEtaSubleadMcReweight{nullptr},
  hRecoDijetPtEta{nullptr},
  hRecoDijetPtEtaDphi{nullptr},
  hRecoDijetPtEtaDphiWeighted{nullptr},
  hRecoDijetPtEtaDphiCM{nullptr},
  hRecoDijetPtEtaDphiCMWeighted{nullptr},
  hRecoDijetPtEtaForward{nullptr},
  hRecoDijetPtEtaBackward{nullptr},
  hRecoDijetPtEtaCMForward{nullptr},
  hRecoDijetPtEtaCMBackward{nullptr},
  hRecoDijetPtEtaForwardWeighted{nullptr},
  hRecoDijetPtEtaBackwardWeighted{nullptr},
  hRecoDijetPtEtaCMForwardWeighted{nullptr},
  hRecoDijetPtEtaCMBackwardWeighted{nullptr},

  hRecoDijetPtEtaDphiJetId{nullptr},
  hRecoInclusiveAllJetPtVsEta{nullptr},
  hRecoInclusiveMatchedJetPtVsEta{nullptr},
  hRecoInclusiveUnmatchedJetPtVsEta{nullptr},

  // Jet selection algo checks
  hRecoInclusiveJetPtVsEtaKineCut{nullptr},
  hRecoInclusiveJetJESPtEtaPhiKineCut{nullptr},
  hRecoInclusiveJetDEtaPtEtaKineCut{nullptr},
  hRecoInclusiveMatchedJetPtVsEtaKineCut{nullptr},
  hRecoInclusiveUnmatchedJetPtVsEtaKineCut{nullptr},
  hRecoInclusiveJetRefPtVsEtaKineCut{nullptr},

  hRecoInclusiveJetPtVsEtaTrkMaxCut{nullptr},
  hRecoInclusiveJetJESPtEtaPhiTrkMaxCut{nullptr},
  hRecoInclusiveJetDEtaPtEtaTrkMaxCut{nullptr},
  hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut{nullptr},
  hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut{nullptr},
  hRecoInclusiveJetRefPtVsEtaTrkMaxCut{nullptr},

  hRecoInclusiveJetPtVsEtaJetIdCut{nullptr},
  hRecoInclusiveJetJESPtEtaPhiJetIdCut{nullptr},
  hRecoInclusiveJetDEtaPtEtaJetIdCut{nullptr},
  hRecoInclusiveMatchedJetPtVsEtaJetIdCut{nullptr},
  hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut{nullptr},
  hRecoInclusiveJetRefPtVsEtaJetIdCut{nullptr},

  hRecoLeadJetAllPtVsEta{nullptr},
  hRecoLeadJetMatchedPtVsEta{nullptr},
  hRecoLeadJetUnmatchedPtVsEta{nullptr},
  hRecoSubLeadJetAllPtVsEta{nullptr},
  hRecoSubLeadJetMatchedPtVsEta{nullptr},
  hRecoSubLeadJetUnmatchedPtVsEta{nullptr},

  hRecoLeadJetAllPtVsEtaJetIdCut{nullptr},
  hRecoLeadJetMatchedPtVsEtaJetIdCut{nullptr},
  hRecoLeadJetUnmatchedPtVsEtaJetIdCut{nullptr},
  hRecoSubLeadJetAllPtVsEtaJetIdCut{nullptr},
  hRecoSubLeadJetMatchedPtVsEtaJetIdCut{nullptr},
  hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut{nullptr},

  hRecoTrkMaxToJetIdDijetMatching{nullptr},

  // Dijets (MC)
  hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta{nullptr},
  hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},
  hRecoDijetPtEtaRefDijetPtEta{nullptr},
  hRecoDijetPtEtaRefDijetPtEtaWeighted{nullptr},
  hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},
  // hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHat{nullptr},
  // hJESDijetPtDijetEtaDijetDeltaPhiGenDijetPtEtaDeltaPhiPtHatWeighted{nullptr},

  hRecoGoodInclusiveJetEtaLabFrame{nullptr},
  hRecoGoodInclusiveJetEtaCMFrame{nullptr},

  hRecoDijetEta{nullptr},
  hRecoDijetEtaCM{nullptr},
  hRecoDijetEta1D{nullptr},
  hRecoDijetEta1DOldPt{nullptr},
  hRecoDijetEta1DOldPtBinning{nullptr},
  hRefDijetEta{nullptr},
  hRefDijetEta1D{nullptr},
  hRefDijetEta1DOldPt{nullptr},
  hRefDijetEta1DOldPtBinning{nullptr},
  hRefDijetEtaVsRecoDijetEta{nullptr},
  hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt{nullptr},
  hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted{nullptr},
  hRefDijetPtEtaDphi{nullptr},
  hRefDijetPtEtaDphiWeighted{nullptr},
  hRefDijetEtaCM{nullptr},
  hRefDijetPtEtaDphiCM{nullptr},
  hRefDijetPtEtaDphiCMWeighted{nullptr},
  hRefDijetPtEtaForward{nullptr},
  hRefDijetPtEtaBackward{nullptr},
  hRefDijetPtEtaCMForward{nullptr},
  hRefDijetPtEtaCMBackward{nullptr},
  hRefDijetPtEtaForwardWeighted{nullptr},
  hRefDijetPtEtaBackwardWeighted{nullptr},
  hRefDijetPtEtaCMForwardWeighted{nullptr},
  hRefDijetPtEtaCMBackwardWeighted{nullptr},

  hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM{nullptr},
  hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted{nullptr},

  hRefSelDijetEta{nullptr},
  hRefSelDijetPtEtaDphi{nullptr},
  hRefSelDijetPtEtaDphiWeighted{nullptr},
  hRefSelDijetEtaCM{nullptr},
  hRefSelDijetPtEtaDphiCM{nullptr},
  hRefSelDijetPtEtaDphiCMWeighted{nullptr},

  hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId{nullptr},
  hRefDijetPtEtaDphiJetId{nullptr},

  hRefInclusiveJetPt{nullptr},
  hRefInclusiveJetPtEta{nullptr},
  hRefPtLeadPtSublead{nullptr},
  hRefEtaLeadEtaSublead{nullptr},
  hRefPtLeadPtSubleadMcReweight{nullptr},
  hRefEtaLeadEtaSubleadMcReweight{nullptr}
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
    for (int i=0; i<4; i++) {
        if (hNHF[i]) delete hNHF[i];
        if (hNEmF[i]) delete hNEmF[i];
        if (hNumOfConst[i]) delete hNumOfConst[i];
        if (hMUF[i]) delete hMUF[i];
        if (hCHF[i]) delete hCHF[i];
        if (hChargedMult[i]) delete hChargedMult[i];
        if (hCEmF[i]) delete hCEmF[i];
        if (hNumOfNeutPart[i]) delete hNumOfNeutPart[i];
    }

    if (hRecoLeadJetAllPtVsEta) delete hRecoLeadJetAllPtVsEta;
    if (hRecoSubLeadJetAllPtVsEta) delete hRecoSubLeadJetAllPtVsEta;
    if (hRecoLeadJetAllPtVsEtaJetIdCut) delete hRecoLeadJetAllPtVsEtaJetIdCut;
    if (hRecoSubLeadJetAllPtVsEtaJetIdCut) delete hRecoSubLeadJetAllPtVsEtaJetIdCut;

    if (hRecoInclusiveAllJetPtVsEta) delete hRecoInclusiveAllJetPtVsEta;
    if (hRecoInclusiveJetPtVsEtaKineCut) delete hRecoInclusiveJetPtVsEtaKineCut;
    if (hRecoInclusiveJetPtVsEtaTrkMaxCut) delete hRecoInclusiveJetPtVsEtaTrkMaxCut;
    if (hRecoInclusiveJetPtVsEtaJetIdCut) delete hRecoInclusiveJetPtVsEtaJetIdCut;

    if (fIsMc) {
        // Gen jets
        if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
        if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
        if (hGenInclusiveJetPt) delete hGenInclusiveJetPt;
        if (hGenInclusiveJetPtEta) delete hGenInclusiveJetPtEta;
        if (hGenPtLeadPtSublead) delete hGenPtLeadPtSublead;
        if (hGenEtaLeadEtaSublead) delete hGenEtaLeadEtaSublead;
        if (hGenEtaCMLeadEtaCMSublead) delete hGenEtaCMLeadEtaCMSublead;
        if (hGenPtLeadPtSubleadMcReweight) delete hGenPtLeadPtSubleadMcReweight;
        if (hGenEtaLeadEtaSubleadMcReweight) delete hGenEtaLeadEtaSubleadMcReweight;
        if (hGenDijetEta) delete hGenDijetEta;
        for (int i=0; i<16; i++) {
            if (hGenDijetEta1D[i]) delete hGenDijetEta1D[i];
        }
        for (int i=0; i<5; i++) {
            if (hGenDijetEta1DOldPt[i]) delete hGenDijetEta1DOldPt[i];
            if (hGenDijetEta1DOldPtBinning[i]) delete hGenDijetEta1DOldPtBinning[i];
        }
        if (hGenDijetPtEtaDphi) delete hGenDijetPtEtaDphi;
        if (hGenDijetPtEtaDphiWeighted) delete hGenDijetPtEtaDphiWeighted;
        if (hGenDijetEtaCM) delete hGenDijetEtaCM;
        if (hGenDijetPtEtaDphiCM) delete hGenDijetPtEtaDphiCM;
        if (hGenDijetPtEtaDphiCMWeighted) delete hGenDijetPtEtaDphiCMWeighted;
        if (hGenDijetPtEtaForward) delete hGenDijetPtEtaForward;
        if (hGenDijetPtEtaBackward) delete hGenDijetPtEtaBackward;
        if (hGenDijetPtEtaCMForward) delete hGenDijetPtEtaCMForward;
        if (hGenDijetPtEtaCMBackward) delete hGenDijetPtEtaCMBackward;
        if (hGenDijetPtEtaForwardWeighted) delete hGenDijetPtEtaForwardWeighted;
        if (hGenDijetPtEtaBackwardWeighted) delete hGenDijetPtEtaBackwardWeighted;
        if (hGenDijetPtEtaCMForwardWeighted) delete hGenDijetPtEtaCMForwardWeighted;
        if (hGenDijetPtEtaCMBackwardWeighted) delete hGenDijetPtEtaCMBackwardWeighted;

        if (hGenGoodInclusiveJetEtaLabFrame) delete hGenGoodInclusiveJetEtaLabFrame;
        if (hGenGoodInclusiveJetEtaCMFrame) delete hGenGoodInclusiveJetEtaCMFrame;


        if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
        if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
        if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
        if (hJESInclusiveJetPtEtaPhi) delete hJESInclusiveJetPtEtaPhi;
        if (hJESInclusiveJetPtEtaPhiWeighted) delete hJESInclusiveJetPtEtaPhiWeighted;
        if (hRecoMatchedPtEta) delete hRecoMatchedPtEta;
        if (hRecoInclusiveMatchedJetPtVsEta) delete hRecoInclusiveMatchedJetPtVsEta;
        if (hRecoInclusiveUnmatchedJetPtVsEta) delete hRecoInclusiveUnmatchedJetPtVsEta;


        if (hRecoLeadJetMatchedPtVsEta) delete hRecoLeadJetMatchedPtVsEta;
        if (hRecoLeadJetUnmatchedPtVsEta) delete hRecoLeadJetUnmatchedPtVsEta;
        if (hRecoSubLeadJetMatchedPtVsEta) delete hRecoSubLeadJetMatchedPtVsEta;
        if (hRecoSubLeadJetUnmatchedPtVsEta) delete hRecoSubLeadJetUnmatchedPtVsEta;

        if (hRecoLeadJetMatchedPtVsEtaJetIdCut) delete hRecoLeadJetMatchedPtVsEtaJetIdCut;
        if (hRecoLeadJetUnmatchedPtVsEtaJetIdCut) delete hRecoLeadJetUnmatchedPtVsEtaJetIdCut;
        if (hRecoSubLeadJetMatchedPtVsEtaJetIdCut) delete hRecoSubLeadJetMatchedPtVsEtaJetIdCut;
        if (hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut) delete hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut;

        
        if (hRecoInclusiveJetJESPtEtaPhiKineCut) delete hRecoInclusiveJetJESPtEtaPhiKineCut;
        if (hRecoInclusiveJetDEtaPtEtaKineCut) delete hRecoInclusiveJetDEtaPtEtaKineCut;
        if (hRecoInclusiveMatchedJetPtVsEtaKineCut) delete hRecoInclusiveMatchedJetPtVsEtaKineCut;
        if (hRecoInclusiveUnmatchedJetPtVsEtaKineCut) delete hRecoInclusiveUnmatchedJetPtVsEtaKineCut;
        if (hRecoInclusiveJetRefPtVsEtaKineCut) delete hRecoInclusiveJetRefPtVsEtaKineCut;

        if (hRecoInclusiveJetJESPtEtaPhiTrkMaxCut) delete hRecoInclusiveJetJESPtEtaPhiTrkMaxCut;
        if (hRecoInclusiveJetDEtaPtEtaTrkMaxCut) delete hRecoInclusiveJetDEtaPtEtaTrkMaxCut;
        if (hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut) delete hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut;
        if (hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut) delete hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut;
        if (hRecoInclusiveJetRefPtVsEtaTrkMaxCut) delete hRecoInclusiveJetRefPtVsEtaTrkMaxCut;

        if (hRecoInclusiveJetJESPtEtaPhiJetIdCut) delete hRecoInclusiveJetJESPtEtaPhiJetIdCut;
        if (hRecoInclusiveJetDEtaPtEtaJetIdCut) delete hRecoInclusiveJetDEtaPtEtaJetIdCut;
        if (hRecoInclusiveMatchedJetPtVsEtaJetIdCut) delete hRecoInclusiveMatchedJetPtVsEtaJetIdCut;
        if (hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut) delete hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut;
        if (hRecoInclusiveJetRefPtVsEtaJetIdCut) delete hRecoInclusiveJetRefPtVsEtaJetIdCut;
    } // if (fIsMc)

    //
    // Dijets (experiment)
    //
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    if (hRecoInclusiveJetPt) delete hRecoInclusiveJetPt;
    if (hRecoPtLeadPtSublead) delete hRecoPtLeadPtSublead;
    if (hRecoEtaLeadEtaSublead) delete hRecoEtaLeadEtaSublead;
    if (hRecoEtaCMLeadEtaCMSublead) delete hRecoEtaCMLeadEtaCMSublead;
    if (hRecoPtLeadPtSubleadMcReweight) delete hRecoPtLeadPtSubleadMcReweight;
    if (hRecoEtaLeadEtaSubleadMcReweight) delete hRecoEtaLeadEtaSubleadMcReweight;
    if (hRecoGoodInclusiveJetEtaLabFrame) delete hRecoGoodInclusiveJetEtaLabFrame;
    if (hRecoGoodInclusiveJetEtaCMFrame) delete hRecoGoodInclusiveJetEtaCMFrame;
    if (hRecoDijetEta) delete hRecoDijetEta;
    for (int i=0; i<16; i++) {
        if (hRecoDijetEta1D[i]) delete hRecoDijetEta1D[i];
    }
    for (int i=0; i<5; i++) {
            if (hRecoDijetEta1DOldPt[i]) delete hRecoDijetEta1DOldPt[i];
            if (hRecoDijetEta1DOldPtBinning[i]) delete hRecoDijetEta1DOldPtBinning[i];
    }
    if (hRecoDijetPtEta) delete hRecoDijetPtEta;
    if (hRecoDijetPtEtaDphi) delete hRecoDijetPtEtaDphi;
    if (hRecoDijetPtEtaDphiWeighted) delete hRecoDijetPtEtaDphiWeighted;
    if (hRecoDijetEtaCM) delete hRecoDijetEtaCM;
    if (hRecoDijetPtEtaDphiCM) delete hRecoDijetPtEtaDphiCM;
    if (hRecoDijetPtEtaDphiCMWeighted) delete hRecoDijetPtEtaDphiCMWeighted;
    if (hRecoDijetPtEtaForward) delete hRecoDijetPtEtaForward;
    if (hRecoDijetPtEtaBackward) delete hRecoDijetPtEtaBackward;
    if (hRecoDijetPtEtaCMForward) delete hRecoDijetPtEtaCMForward;
    if (hRecoDijetPtEtaCMBackward) delete hRecoDijetPtEtaCMBackward;
    if (hRecoDijetPtEtaForwardWeighted) delete hRecoDijetPtEtaForwardWeighted;
    if (hRecoDijetPtEtaBackwardWeighted) delete hRecoDijetPtEtaBackwardWeighted;
    if (hRecoDijetPtEtaCMForwardWeighted) delete hRecoDijetPtEtaCMForwardWeighted;
    if (hRecoDijetPtEtaCMBackwardWeighted) delete hRecoDijetPtEtaCMBackwardWeighted;

    if (hRecoDijetPtEtaDphiJetId) delete hRecoDijetPtEtaDphiJetId;

    if (hRecoTrkMaxToJetIdDijetMatching) delete hRecoTrkMaxToJetIdDijetMatching;

    //
    // Dijets exp vs mc
    //
    if ( fIsMc ) {
        if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta;
        if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;
        if (hRecoDijetPtEtaRefDijetPtEta) delete hRecoDijetPtEtaRefDijetPtEta;
        if (hRecoDijetPtEtaRefDijetPtEtaWeighted) delete hRecoDijetPtEtaRefDijetPtEtaWeighted;
        if (hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;
    
        if (hRefInclusiveJetPt) delete hRefInclusiveJetPt;
        if (hRefInclusiveJetPtEta) delete hRefInclusiveJetPtEta;
        if (hRefPtLeadPtSublead) delete hRefPtLeadPtSublead;
        if (hRefEtaLeadEtaSublead) delete hRefEtaLeadEtaSublead;
        if (hRefEtaCMLeadEtaCMSublead) delete hRefEtaCMLeadEtaCMSublead;
        if (hRefPtLeadPtSubleadMcReweight) delete hRefPtLeadPtSubleadMcReweight;
        if (hRefEtaLeadEtaSubleadMcReweight) delete hRefEtaLeadEtaSubleadMcReweight;
        if (hRefDijetEta) delete hRefDijetEta;
        if (hRefDijetPtEtaDphi) delete hRefDijetPtEtaDphi;
        if (hRefDijetPtEtaDphiWeighted) delete hRefDijetPtEtaDphiWeighted;

        if (hRefDijetEtaVsRecoDijetEta) delete hRefDijetEtaVsRecoDijetEta;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted;

        if (hRefDijetEtaCM) delete hRefDijetEtaCM;
        for (int i=0; i<16; i++) {
            if (hRefDijetEta1D[i]) delete hRefDijetEta1D[i];
        }
        for (int i=0; i<5; i++) {
            if (hRefDijetEta1DOldPt[i]) delete hRefDijetEta1DOldPt[i];
            if (hRefDijetEta1DOldPtBinning[i]) delete hRefDijetEta1DOldPtBinning[i];
        }
        if (hRefDijetPtEtaDphiCM) delete hRefDijetPtEtaDphiCM;
        if (hRefDijetPtEtaDphiCMWeighted) delete hRefDijetPtEtaDphiCMWeighted;

        if (hRefDijetPtEtaForward) delete hRefDijetPtEtaForward;
        if (hRefDijetPtEtaBackward) delete hRefDijetPtEtaBackward;
        if (hRefDijetPtEtaCMForward) delete hRefDijetPtEtaCMForward;
        if (hRefDijetPtEtaCMBackward) delete hRefDijetPtEtaCMBackward;
        if (hRefDijetPtEtaForwardWeighted) delete hRefDijetPtEtaForwardWeighted;
        if (hRefDijetPtEtaBackwardWeighted) delete hRefDijetPtEtaBackwardWeighted;
        if (hRefDijetPtEtaCMForwardWeighted) delete hRefDijetPtEtaCMForwardWeighted;
        if (hRefDijetPtEtaCMBackwardWeighted) delete hRefDijetPtEtaCMBackwardWeighted;

        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted;

        if (hRefSelDijetEta) delete hRefSelDijetEta;
        if (hRefSelDijetPtEtaDphi) delete hRefSelDijetPtEtaDphi;
        if (hRefSelDijetPtEtaDphiWeighted) delete hRefSelDijetPtEtaDphiWeighted;
        if (hRefSelDijetEtaCM) delete hRefSelDijetEtaCM;
        if (hRefSelDijetPtEtaDphiCM) delete hRefSelDijetPtEtaDphiCM;
        if (hRefSelDijetPtEtaDphiCMWeighted) delete hRefSelDijetPtEtaDphiCMWeighted;

        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId;
        if (hRefDijetPtEtaDphiJetId) delete hRefDijetPtEtaDphiJetId;
    } // if ( fIsMc )
}

//________________
void HistoManagerDiJet::init(const Bool_t& isMc) {

    const int dijetEtaBins{30};
    double dijetEtaVals[dijetEtaBins+1] { -5.0, -4.0, -3.0, -2.4, -2.2, 
                                            -2.0, -1.8, -1.6, -1.4, -1.2, 
                                            -1.0, -0.8, -0.6, -0.4, -0.2,  
                                             0.0,  0.2,  0.4,  0.6,  0.8,  
                                             1.0,  1.2,  1.4,  1.6,  1.8,  
                                             2.0,  2.2,  2.4,  3.0,  4.0,  
                                             5.0 };

    // Old binning convention
    const int dijetEtaOldBins{18};
    double dijetEtaOldVals[dijetEtaBins+1] = { -2.915, -2.63333333333, -2.07, -1.78833333333, -1.50666666667,
                                                -1.225, -0.94333333333, -0.66166666666, -0.38, -0.09833333333,
                                                0.18333333333, 0.465, 0.74666666666, 1.02833333333, 1.31,
                                                1.59166666667, 1.87333333333, 2.43666666667, 3.};

    const int dijetEtaFBBins{13};
    double dijetEtaFBVals[dijetEtaFBBins+1] { 0.0,  0.2,  0.4,  0.6,  0.8,  
                                                1.0,  1.2,  1.4,  1.6,  1.8,  
                                                2.0,  2.2,  2.4,  3.0 };


    const int dijetPtBins{16};
    double dijetPtVals[dijetPtBins+1] {  50.,  60.,   70.,  80.,  90.,
                                         100., 110.,  120., 130., 140.,
                                         150., 160.,  180., 200., 250., 
                                         300., 500.};

    // Old binning convention
    const int dijetPtOldBins{6};
    double dijetPtOldVals[dijetPtBins+1] {25., 55., 75., 95., 115., 150., 400.}; // 6 bins
    
    int    prescale = 2;

    int    vzBins = 320;
    double vzRange[2] {-31., 31.};
    int    multBins{1800};
    double multRange[2] {-0.5, 1799.5};
    int    hiBinBins{203};
    double hiBinRange[2] {-1.5, 201.5};
    // int    centralityBins{101};
    // double centralityRange[2] {-0.5, 100.5};
    int    weightBins{110};
    double weightRange[2] {-0.05, 1.05};
    int    ptHatBins{100};
    double ptHatRange[2] {0., 1000.};
    int    fJESBins{100}; 
    double fJESRange[2] {0., 2.};


    int    bins2D_ev_VzPtHat[2] {     vzBins,     fPtHatBins };
    double xmin2D_ev_VzPtHat[2] { vzRange[0], fPtHatRange[0] };
    double xmax2D_ev_VzPtHat[2] { vzRange[1], fPtHatRange[1] };

    int  ptShortBins{40}; 
    double ptShortRange[2] {15., 215.};

    int    dEtaBins{50}; 
    double dEtaRange[2] {-0.05, 0.05};

    //
    // Gen
    //
    int    bins9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi[9]
    { fDijetPtBins, fDijetEtaBins, fDijetDphiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    double xmin9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi[9]
    { fDijetPtRange[0], fDijetEtaRange[0], fDijetDphiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0]  };
    double xmax9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi[9]
    { fDijetPtRange[1], fDijetEtaRange[1], fDijetDphiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    //
    // Reco single
    //
    int    bins5D_jet_PtPtPtEtaEta[5] { fPtBins, fPtBins, fPtBins, fEtaBins, fEtaBins };
    double xmin5D_jet_PtPtPtEtaEta[5] { fPtRange[0], fPtRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0] };
    double xmax5D_jet_PtPtPtEtaEta[5] { fPtRange[1], fPtRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1] };

    int    bins4D_jet_JESPtEtaPhi[4] { fJESBins,     fPtBins,     fEtaBins,     fPhiBins };
    double xmin4D_jet_JESPtEtaPhi[4] { fJESRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0] };
    double xmax4D_jet_JESPtEtaPhi[4] { fJESRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    //
    // Reco dijets
    //
    int    bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi[9]
    { fDijetPtBins, fDijetEtaBins, fDijetDphiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    double xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi[9]
    { fDijetPtRange[0], fDijetEtaRange[0], fDijetDphiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0] };
    double xmax9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi[9]
    { fDijetPtRange[1], fDijetEtaRange[1], fDijetDphiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    int    bins12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta[12]
    { fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins, fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins };
    double xmin12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta[12]
    { fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0] };
    double xmax12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta[12]
    { fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1] };

    int    bins4D_dijet_PtEtaPtEta[4] { fDijetPtBins, fDijetEtaBins, fDijetPtBins, fDijetEtaBins };
    double xmin4D_dijet_PtEtaPtEta[4] { fDijetPtRange[0], fDijetEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0] };
    double xmax4D_dijet_PtEtaPtEta[4] { fDijetPtRange[1], fDijetEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1] };

    //
    // Event histograms
    //
    
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
    hVzPtHat = new THnSparseD( "hVzPtHat","Vertex and #hat{p_{T}};vz (cm);#hat{p_{T}} (GeV/c)",
                               2, bins2D_ev_VzPtHat, xmin2D_ev_VzPtHat, xmax2D_ev_VzPtHat );
    hVzPtHat->Sumw2();
    hVzPtHatWeighted = new THnSparseD( "hVzPtHatWeighted","Vertex and #hat{p_{T}} weighted;vz (cm);#hat{p_{T}} (GeV/c)",
                                       2, bins2D_ev_VzPtHat, xmin2D_ev_VzPtHat, xmax2D_ev_VzPtHat );
    hVzPtHatWeighted->Sumw2();

    // For 4 eta ranges: <=2.4, <=2.7, <=3, >3
    for (int i{0}; i<4; i++) {
        float low{0}, hi{0};
        if      ( i == 0 ) { low = {-2.4f}; hi = {2.4f}; }
        else if ( i == 1 ) { low = {2.4f}; hi = {2.7f}; }
        else if ( i == 2 ) { low = {2.7f}; hi = {3.0f}; }
        else if ( i == 3 ) { low = {3.f}; hi = {100.0f}; }
        hNHF[i] = new TH1D(Form("hNHF_%d",i), Form("Neutral hadron fraction for %3.1f< #eta<=%3.1f;Neutral hadron fraction;1/N dN/dNHF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hNHF[i]->Sumw2();
        hNEmF[i] = new TH1D(Form("hNEmF_%d",i), Form("Neutral EM fraction for %3.1f< #eta<=%3.1f;Neutral EM fraction;1/N dN/dNEF", low, hi), 
                            fFracBins, fFracRange[0], fFracRange[1]);
        hNEmF[i]->Sumw2();
        hNumOfConst[i] = new TH1D(Form("hNumOfConst_%d",i), Form("Number of constituents for %3.1f< #eta<=%3.1f;Number of constituents;1/N dN/dNconst", low, hi), 
                                  fMultBins, fMultRange[0], fMultRange[1]);
        hNumOfConst[i]->Sumw2();
        hMUF[i] = new TH1D(Form("hMUF_%d",i), Form("Muon fraction for %3.1f< #eta<=%3.1f;Muon fraction;1/N dN/dMUF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hMUF[i]->Sumw2();
        hCHF[i] = new TH1D(Form("hCHF_%d",i), Form("Charged hadron fraction for %3.1f< #eta<=%3.1f;Charged hadron fraction;1/N dN/dCHF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hCHF[i]->Sumw2();
        hChargedMult[i] = new TH1D(Form("hChargedMult_%d",i), Form("Charged particle multiplicity for %3.1f< #eta<=%3.1f;Charged multiplicity;1/N dN/dChMult", low, hi), 
                                   fMultBins, fMultRange[0], fMultRange[1]);
        hChargedMult[i]->Sumw2();
        hCEmF[i] = new TH1D(Form("hCEmF_%d",i), Form("Charged EM fraction for %3.1f< #eta<=%3.1f;Charged EM fraction;1/N dN/dChEmF", low, hi), 
                            fFracBins, fFracRange[0], fFracRange[1]);
        hCEmF[i]->Sumw2();
        hNumOfNeutPart[i] = new TH1D(Form("hNumOfNeutPart_%d",i), Form("Number of neutral particles for %3.1f< #eta<=%3.1f;Number of neutral particles;1/N dN/dNumOfNeutrals", low, hi), 
                                     fMultBins, fMultRange[0], fMultRange[1]);
        hNumOfNeutPart[i]->Sumw2();
    } // for (int i{0}; i<4; i++)


    //
    // Monte Carlo information
    //
    if (fIsMc) {

        //
        // Gen inclusive jets
        //

        hGenInclusiveJetPt = new TH1D("hGenInclusiveJetPt","Inclusive gen jet;Gen p_{T}^{inclusive} (GeV/c)",
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPt->Sumw2();
        hGenInclusiveJetPtEta = new TH2D("hGenInclusiveJetPtEta","Gen inclusive jet acceptance;Gen #eta;Gen p_{T} (GeV/c)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        ptShortBins, ptShortRange[0], ptShortRange[1] );
        hGenInclusiveJetPtEta->Sumw2();

        //
        // Gen dijets
        //
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi = new THnSparseD("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi","Title;Gen p_{T}^{dijet} (GeV/c);Gen #eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV/c);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV/c);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
                9, 
                bins9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmin9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmax9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi);
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Sumw2();
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted = new THnSparseD("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted","Title;Gen p_{T}^{dijet} (GeV/c);Gen #eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV/c);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV/c);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
                9, 
                bins9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmin9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmax9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi);
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Sumw2();


        hGenPtLeadPtSublead = new TH2D("hGenPtLeadPtSublead","Leading gen jet pT vs subleading gen jet pT;Gen p_{T}^{Leading} (GeV/c);Gen p_{T}^{Subleading} (GeV/c)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenPtLeadPtSublead->Sumw2();
        hGenEtaLeadEtaSublead = new TH2D("hGenEtaLeadEtaSublead","Leading gen jet eta vs subleading gen jet eta;Gen #eta^{Leading};Gen #eta^{Subleading}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaLeadEtaSublead->Sumw2();
        hGenEtaCMLeadEtaCMSublead = new TH2D("hGenEtaCMLeadEtaCMSublead","Leading gen jet eta in CM vs subleading gen jet eta in CM;Gen #eta^{Leading}_{CM};Gen #eta^{Subleading}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaCMLeadEtaCMSublead->Sumw2();
        hGenPtLeadPtSubleadMcReweight = new TH2D("hGenPtLeadPtSubleadMcReweight","Leading gen jet pT vs subleading gen jet pT (MC reweighted to data);Gen p_{T}^{Leading} (GeV/c);Gen p_{T}^{Subleading} (GeV/c)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenPtLeadPtSubleadMcReweight->Sumw2();
        hGenEtaLeadEtaSubleadMcReweight = new TH2D("hGenEtaLeadEtaSubleadMcReweight","Leading gen jet eta vs subleading gen jet eta (MC reweighted to data);Gen #eta^{Leading};Gen #eta^{Subleading}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaLeadEtaSubleadMcReweight->Sumw2();
        hGenDijetEta = new TH1D("hGenDijetEta", "Gen dijet #eta;#eta^{dijet}",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEta->Sumw2();
        for (int i=0; i<16; i++) {
            hGenDijetEta1D[i] = new TH1D(Form("hGenDijetEta1D_%d",i), Form("Gen dijet #eta bin %d;#eta^{dijet}",i),
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            hGenDijetEta1D[i]->Sumw2();
        }
        for (int i=0; i<5; i++) {
            hGenDijetEta1DOldPt[i] = new TH1D(Form("hGenDijetEta1DOldPt_%d",i), Form("Gen dijet #eta old p_{T} bin %d;#eta^{dijet}",i),
                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            hGenDijetEta1DOldPt[i]->Sumw2();
            hGenDijetEta1DOldPtBinning[i] = new TH1D(Form("hGenDijetEta1DOldPtBinning_%d",i), Form("Gen dijet #eta old p_{T} bin %d old #eta binning;#eta^{dijet}",i),
                                                     dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinning[i]->Sumw2();
        }

        hGenDijetPtEtaDphi = new TH3D("hGenDijetPtEtaDphi","Gen dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                      fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                      fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                      fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphi->Sumw2();
        hGenDijetPtEtaDphiWeighted = new TH3D("hGenDijetPtEtaDphiWeighted","Gen dijet info weighted;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphiWeighted->Sumw2();

        hGenDijetEtaCM = new TH1D("hGenDijetEtaCM", "Gen dijet #eta in CM;#eta^{dijet}_{CM}",
                                  fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEtaCM->Sumw2();
        hGenDijetPtEtaDphiCM = new TH3D("hGenDijetPtEtaDphiCM","Gen dijet info in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphiCM->Sumw2();
        hGenDijetPtEtaDphiCMWeighted = new TH3D("hGenDijetPtEtaDphiCMWeighted","Gen dijet info weighted in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphiCMWeighted->Sumw2();


        hGenDijetPtEtaForward = new TH2D("hGenDijetPtEtaForward", "Gen dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaForward->Sumw2();
        hGenDijetPtEtaBackward = new TH2D("hGenDijetPtEtaBackward", "Gen dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaBackward->Sumw2();
        hGenDijetPtEtaCMForward = new TH2D("hGenDijetPtEtaCMForward", "Gen dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaCMForward->Sumw2();
        hGenDijetPtEtaCMBackward = new TH2D("hGenDijetPtEtaCMBackward", "Gen dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaCMBackward->Sumw2();

        hGenDijetPtEtaForwardWeighted = new TH2D("hGenDijetPtEtaForwardWeighted", "Gen dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaForwardWeighted->Sumw2();
        hGenDijetPtEtaBackwardWeighted = new TH2D("hGenDijetPtEtaBackwardWeighted", "Gen dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaBackwardWeighted->Sumw2();
        hGenDijetPtEtaCMForwardWeighted = new TH2D("hGenDijetPtEtaCMForwardWeighted", "Gen dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaCMForwardWeighted->Sumw2();

        hGenDijetPtEtaCMBackwardWeighted = new TH2D("hGenDijetPtEtaCMBackwardWeighted", "Gen dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaCMBackwardWeighted->Sumw2();



        hGenGoodInclusiveJetEtaLabFrame = new TH1D("hGenGoodInclusiveJetEtaLabFrame","Gen good inclusive jet #eta in lab frame;#eta^{Inclusive}",
                                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenGoodInclusiveJetEtaLabFrame->Sumw2();
        hGenGoodInclusiveJetEtaCMFrame = new TH1D("hGenGoodInclusiveJetEtaCMFrame","Gen good inclusive jet #eta in CM frame;#eta^{Inclusive}_{CM}",
                                                  fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenGoodInclusiveJetEtaCMFrame->Sumw2();


        //
        // Reco single jets
        //
        hRecoInclusiveUnmatchedJetPtVsEta = new TH2D("hRecoInclusiveUnmatchedJetPtVsEta", "Inclusive reco jet unmatched gen pT vs eta;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveUnmatchedJetPtVsEta->Sumw2();

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

        hJESInclusiveJetPtEtaPhi = new THnSparseD("hJESInclusiveJetPtEtaPhi","JES of inclusive jets;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV/c);#eta^{reco};#phi^{reco} (rad)",
                4,
                bins4D_jet_JESPtEtaPhi,
                xmin4D_jet_JESPtEtaPhi,
                xmax4D_jet_JESPtEtaPhi);
        hJESInclusiveJetPtEtaPhi->Sumw2();
        hJESInclusiveJetPtEtaPhiWeighted = new THnSparseD("hJESInclusiveJetPtEtaPhiWeighted","JES of inclusive jet weighted;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV/c);#eta^{reco};#phi^{reco} (rad)",
                4,
                bins4D_jet_JESPtEtaPhi,
                xmin4D_jet_JESPtEtaPhi,
                xmax4D_jet_JESPtEtaPhi);
        hJESInclusiveJetPtEtaPhiWeighted->Sumw2();

        hRecoMatchedPtEta = new TH2D("hRecoMatchedPtEta","Reconstructed jets that matched gen;#eta^{reco};p_{T}^{reco} (GeV/c)",
                fEtaBins, fEtaRange[0], fEtaRange[1],
                ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoMatchedPtEta->Sumw2();
        hRecoInclusiveMatchedJetPtVsEta = new TH2D("hRecoInclusiveMatchedJetPtVsEta","Inclusive reco jet that has matching pT vs eta;#eta;p_{T} (GeV/c)",
                                                    fEtaBins, fEtaRange[0], fEtaRange[1],
                                                    ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveMatchedJetPtVsEta->Sumw2();



        hRecoLeadJetMatchedPtVsEta = new TH2D("hRecoLeadJetMatchedPtVsEta","Leading jet matched p_{T} vs #eta;#eta;p_{T} (GeV/c)", 
                                              fEtaBins, fEtaRange[0], fEtaRange[1], 
                                              ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoLeadJetMatchedPtVsEta->Sumw2();
        hRecoLeadJetUnmatchedPtVsEta = new TH2D("hRecoLeadJetUnmatchedPtVsEta","Leading jet unmatched p_{T} vs #eta;#eta;p_{T} (GeV/c)", 
                                                fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoLeadJetUnmatchedPtVsEta->Sumw2();
        hRecoSubLeadJetMatchedPtVsEta = new TH2D("hRecoSubLeadJetMatchedPtVsEta","Subleading jet matched p_{T} vs #eta;#eta;p_{T} (GeV/c)", 
                                                 fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                 ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoSubLeadJetMatchedPtVsEta->Sumw2();
        hRecoSubLeadJetUnmatchedPtVsEta = new TH2D("hRecoSubLeadJetUnmatchedPtVsEta","Subleading jet unmatched p_{T} vs #eta;#eta;p_{T} (GeV/c)", 
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                   ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoSubLeadJetUnmatchedPtVsEta->Sumw2();


        hRecoLeadJetMatchedPtVsEtaJetIdCut = new TH2D("hRecoLeadJetMatchedPtVsEtaJetIdCut","Leading jet matched p_{T} vs #eta with jetId selection;#eta;p_{T} (GeV/c)", 
                                                      fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                      ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoLeadJetMatchedPtVsEtaJetIdCut->Sumw2();
        hRecoLeadJetUnmatchedPtVsEtaJetIdCut = new TH2D("hRecoLeadJetUnmatchedPtVsEtaJetIdCut","Leading jet unmatched p_{T} vs #eta with jetId selection;#eta;p_{T} (GeV/c)", 
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                        ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoLeadJetUnmatchedPtVsEtaJetIdCut->Sumw2();
        hRecoSubLeadJetMatchedPtVsEtaJetIdCut = new TH2D("hRecoSubLeadJetMatchedPtVsEtaJetIdCut","Subleading jet matched p_{T} vs #eta with jetId selection;#eta;p_{T} (GeV/c)", 
                                                         fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                         ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoSubLeadJetMatchedPtVsEtaJetIdCut->Sumw2();
        hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut = new TH2D("hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut","Subleading jet unmatched p_{T} vs #eta with jetId selection;#eta;p_{T} (GeV/c)", 
                                                           fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                           ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut->Sumw2();

        //
        // Reco jets selection algo
        //
        hRecoInclusiveJetJESPtEtaPhiKineCut = new THnSparseD("hRecoInclusiveJetJESPtEtaPhiKineCut","Inclusive jet JES for kine cut;p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta;#phi (rad)",
                        4,
                        bins4D_jet_JESPtEtaPhi,
                        xmin4D_jet_JESPtEtaPhi,
                        xmax4D_jet_JESPtEtaPhi);
        hRecoInclusiveJetJESPtEtaPhiKineCut->Sumw2();
        hRecoInclusiveJetDEtaPtEtaKineCut = new TH3D("hRecoInclusiveJetDEtaPtEtaKineCut", "#Delta#eta for kine cut;#eta^{reco}-#eta^{ref};Ref p_{T} (GeV/c);Ref #eta", 
                                                     dEtaBins, dEtaRange[0], dEtaRange[1],
                                                     fPtBins, fPtRange[0], fPtRange[1],
                                                     fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetDEtaPtEtaKineCut->Sumw2();
        hRecoInclusiveMatchedJetPtVsEtaKineCut = new TH2D("hRecoInclusiveMatchedJetPtVsEtaKineCut", "Inclusive reco matched jet pT vs eta kine cut;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveMatchedJetPtVsEtaKineCut->Sumw2();
        hRecoInclusiveUnmatchedJetPtVsEtaKineCut = new TH2D("hRecoInclusiveUnmatchedJetPtVsEtaKineCut", "Inclusive reco unmatched jet pT vs eta kine cut;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveUnmatchedJetPtVsEtaKineCut->Sumw2();
        hRecoInclusiveJetRefPtVsEtaKineCut = new TH2D("hRecoInclusiveJetRefPtVsEtaKineCut", "Inclusive reco jet ref pT vs ref eta kine cut;#eta^{ref};p_{T}^{ref} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveJetRefPtVsEtaKineCut->Sumw2();


        hRecoInclusiveJetJESPtEtaPhiTrkMaxCut = new THnSparseD("hRecoInclusiveJetJESPtEtaPhiTrkMaxCut","Inclusive jet JES for trkMax cut;p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta;#phi (rad)",
                        4,
                        bins4D_jet_JESPtEtaPhi,
                        xmin4D_jet_JESPtEtaPhi,
                        xmax4D_jet_JESPtEtaPhi);
        hRecoInclusiveJetJESPtEtaPhiTrkMaxCut->Sumw2();
        hRecoInclusiveJetDEtaPtEtaTrkMaxCut = new TH3D("hRecoInclusiveJetDEtaPtEtaTrkMaxCut", "#Delta#eta for trkMax cut;#eta^{reco}-#eta^{ref};Ref p_{T} (GeV/c);Ref #eta", 
                                                       dEtaBins, dEtaRange[0], dEtaRange[1],
                                                       fPtBins, fPtRange[0], fPtRange[1],
                                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetDEtaPtEtaTrkMaxCut->Sumw2();
        hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut = new TH2D("hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut", "Inclusive reco matched jet pT vs eta trkMax cut;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut->Sumw2();
        hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut = new TH2D("hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut", "Inclusive reco unmatched jet pT vs eta trkMax cut;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut->Sumw2();
        hRecoInclusiveJetRefPtVsEtaTrkMaxCut = new TH2D("hRecoInclusiveJetRefPtVsEtaTrkMaxCut", "Inclusive reco matched jet ref pT vs ref eta trkMax cut;#eta^{ref};p_{T}^{ref} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveJetRefPtVsEtaTrkMaxCut->Sumw2();


        hRecoInclusiveJetJESPtEtaPhiJetIdCut = new THnSparseD("hRecoInclusiveJetJESPtEtaPhiJetIdCut","Inclusive jet JES for jetId cut;p_{T}^{reco}/p_{T}^{ref};p_{T} (GeV);#eta;#phi (rad)",
                        4,
                        bins4D_jet_JESPtEtaPhi,
                        xmin4D_jet_JESPtEtaPhi,
                        xmax4D_jet_JESPtEtaPhi);
        hRecoInclusiveJetJESPtEtaPhiJetIdCut->Sumw2();
        hRecoInclusiveJetDEtaPtEtaJetIdCut = new TH3D("hRecoInclusiveJetDEtaPtEtaJetIdCut", "#Delta#eta for jetId cut;#eta^{reco}-#eta^{ref};Ref p_{T} (GeV/c);Ref #eta", 
                                                       dEtaBins, dEtaRange[0], dEtaRange[1],
                                                       fPtBins, fPtRange[0], fPtRange[1],
                                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetDEtaPtEtaJetIdCut->Sumw2();
        hRecoInclusiveMatchedJetPtVsEtaJetIdCut = new TH2D("hRecoInclusiveMatchedJetPtVsEtaJetIdCut", "Inclusive reco matched jet pT vs eta jetId cut;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveMatchedJetPtVsEtaJetIdCut->Sumw2();
        hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut = new TH2D("hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut", "Inclusive reco unmatched jet pT vs eta jetId cut;#eta;p_{T} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut->Sumw2();
        hRecoInclusiveJetRefPtVsEtaJetIdCut = new TH2D("hRecoInclusiveJetRefPtVsEtaJetIdCut", "Inclusive reco matched jet ref pT vs ref eta jetId cut;#eta^{ref};p_{T}^{ref} (GeV/c)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     ptShortBins, ptShortRange[0], ptShortRange[1]);
        hRecoInclusiveJetRefPtVsEtaJetIdCut->Sumw2();


        //
        // Reco dijet with MC
        //
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta = new THnSparseD("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta",
                "Reco to ref correspondence;Reco p_{T}^{dijet} (GeV/c);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV/c);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV/c);Ref #eta^{Subleading}",
                12,
                bins12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmin12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmax12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta);
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Sumw2();
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted = new THnSparseD("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted",
                "Reco to ref correspondence weighted;Reco p_{T}^{dijet} (GeV/c);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV/c);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV/c);Ref #eta^{Subleading}",
                12,
                bins12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmin12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmax12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta);
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Sumw2();

        hRecoDijetPtEtaRefDijetPtEta = new THnSparseD("hRecoDijetPtEtaRefDijetPtEta","Reco dijet vs Ref dijet;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}; Ref p_{T}^{ave} (GeV); Ref #eta^{dijet}",
                                                      4,
                                                      bins4D_dijet_PtEtaPtEta,
                                                      xmin4D_dijet_PtEtaPtEta,
                                                      xmax4D_dijet_PtEtaPtEta);
        hRecoDijetPtEtaRefDijetPtEta->Sumw2();
        hRecoDijetPtEtaRefDijetPtEtaWeighted = new THnSparseD("hRecoDijetPtEtaRefDijetPtEtaWeighted","Reco dijet vs Ref dijet;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}; Ref p_{T}^{ave} (GeV); Ref #eta^{dijet}",
                                                      4,
                                                      bins4D_dijet_PtEtaPtEta,
                                                      xmin4D_dijet_PtEtaPtEta,
                                                      xmax4D_dijet_PtEtaPtEta);
        hRecoDijetPtEtaRefDijetPtEtaWeighted->Sumw2();

        //
        // Ref selected dijets
        //
        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted = new THnSparseD("hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted",
                "Reco to ref correspondence (via ref selection) weighted;Reco p_{T}^{dijet} (GeV/c);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV/c);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV/c);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV/c);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV/c);Ref #eta^{Subleading}",
                12,
                bins12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmin12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmax12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta);
        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Sumw2();
        hRefSelDijetEta = new TH1D("hRefSelDijetEta","Ref selected dijets;#eta^{dijet}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEta->Sumw2();
        hRefSelDijetPtEtaDphi = new TH3D("hRefSelDijetPtEtaDphi","RefSel dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphi->Sumw2();
        hRefSelDijetPtEtaDphiWeighted = new TH3D("hRefSelDijetPtEtaDphiWeighted","RefSel dijet info weighted;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphiWeighted->Sumw2();
        hRefSelDijetEtaCM = new TH1D("hRefSelDijetEtaCM","Ref selected dijets in CM;#eta^{dijet}_{CM}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEtaCM->Sumw2();
        hRefSelDijetPtEtaDphiCM = new TH3D("hRefSelDijetPtEtaDphiCM","RefSel dijet info in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphiCM->Sumw2();
        hRefSelDijetPtEtaDphiCMWeighted = new TH3D("hRefSelDijetPtEtaDphiCMWeighted","RefSel dijet info weighted in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphiCMWeighted->Sumw2();

        //
        // Ref inclusive jets
        //
        hRefInclusiveJetPt = new TH1D("hRefInclusiveJetPt","Ref jet p_{T};Ref p_{T} (GeV/c);Entries",
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPt->Sumw2();
        hRefInclusiveJetPtEta = new TH2D("hRefInclusiveJetPtEta","Ref jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV/c)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        //
        // Ref dijets
        //
        hRefPtLeadPtSublead = new TH2D("hRefPtLeadPtSublead","Ref leading vs subleading p_{T};Ref p_{T}^{Leading} (GeV/c);Ref p_{T}^{Subleading} (GeV/c)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSublead->Sumw2();
        hRefEtaLeadEtaSublead = new TH2D("hRefEtaLeadEtaSublead","Ref leading vs subleading #eta;Ref #eta^{Leading};Ref #eta^{Subleading}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaLeadEtaSublead->Sumw2();
        hRefEtaCMLeadEtaCMSublead = new TH2D("hRefEtaCMLeadEtaCMSublead","Ref leading vs subleading #eta in CM;Ref #eta^{Leading}_{CM};Ref #eta^{Subleading}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaCMLeadEtaCMSublead->Sumw2();
        hRefPtLeadPtSubleadMcReweight = new TH2D("hRefPtLeadPtSubleadMcReweight","Ref leading vs subleading p_{T} (MC reweighted to data);Ref p_{T}^{Leading} (GeV/c);Ref p_{T}^{Subleading} (GeV/c)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSubleadMcReweight->Sumw2();
        hRefEtaLeadEtaSubleadMcReweight = new TH2D("hRefEtaLeadEtaSubleadMcReweight","Ref leading vs subleading #eta (MC reweighted to data);Ref #eta^{Leading};Ref #eta^{Subleading}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaLeadEtaSubleadMcReweight->Sumw2();
        hRefDijetEta = new TH1D("hRefDijetEta","Ref dijet #eta;Ref #eta^{dijet};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEta->Sumw2();

        for (int i{0}; i<16; i++) {
            hRefDijetEta1D[i] = new TH1D(Form("hRefDijetEta1D_%d",i),Form("Ref dijet #eta for p_{T} bin %d;Ref #eta^{dijet};Entries",i),
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            hRefDijetEta1D[i]->Sumw2();
        }
        for (int i{0}; i<5; i++) {
            hRefDijetEta1DOldPt[i] = new TH1D(Form("hRefDijetEta1DOldPt_%d",i),Form("Ref dijet #eta for old p_{T} bin %d;Ref #eta^{dijet};Entries",i),
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            hRefDijetEta1DOldPt[i]->Sumw2();
            hRefDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefDijetEta1DOldPtBinning_%d",i),Form("Ref dijet #eta for old p_{T} bin %d old #eta binning;Ref #eta^{dijet};Entries",i),
                                                     dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinning[i]->Sumw2();
                                                     
        }

        hRefDijetEtaVsRecoDijetEta = new TH2D("hRefDijetEtaVsRecoDijetEta","Ref dijet #eta vs reco dijet #eta;Reco #eta^{dijet};Ref #eta^{dijet}",
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaVsRecoDijetEta->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T};Reco #eta^{dijet};Ref #eta^{dijet}; Reco p_{T}^{dijet} (GeV/c)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} weighted;Reco #eta^{dijet};Ref #eta^{dijet}; Reco p_{T}^{dijet} (GeV/c)",
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Sumw2();
        hRefDijetPtEtaDphi = new TH3D("hRefDijetPtEtaDphi","Ref dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphi->Sumw2();
        hRefDijetPtEtaDphiWeighted = new TH3D("hRefDijetPtEtaDphiWeighted","Ref dijet info weighted;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphiWeighted->Sumw2();

        hRefDijetEtaCM = new TH1D("hRefDijetEtaCM","Ref dijet #eta in CM;Ref #eta^{dijet}_{CM};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV/c)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} weighted in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV/c)",
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Sumw2();
        hRefDijetPtEtaDphiCM = new TH3D("hRefDijetPtEtaDphiCM","Ref dijet info in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphiCM->Sumw2();
        hRefDijetPtEtaDphiCMWeighted = new TH3D("hRefDijetPtEtaDphiCMWeighted","Ref dijet info weighted in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphiCMWeighted->Sumw2();

        hRefDijetPtEtaForward = new TH2D("hRefDijetPtEtaForward", "Ref dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaForward->Sumw2();
        hRefDijetPtEtaBackward = new TH2D("hRefDijetPtEtaBackward", "Ref dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaBackward->Sumw2();
        hRefDijetPtEtaCMForward = new TH2D("hRefDijetPtEtaCMForward", "Ref dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaCMForward->Sumw2();
        hRefDijetPtEtaCMBackward = new TH2D("hRefDijetPtEtaCMBackward", "Ref dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaCMBackward->Sumw2();
        
        hRefDijetPtEtaForwardWeighted = new TH2D("hRefDijetPtEtaForwardWeighted", "Ref dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaForwardWeighted->Sumw2();
        hRefDijetPtEtaBackwardWeighted = new TH2D("hRefDijetPtEtaBackwardWeighted", "Ref dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaBackwardWeighted->Sumw2();
        hRefDijetPtEtaCMForwardWeighted = new TH2D("hRefDijetPtEtaCMForwardWeighted", "Ref dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaCMForwardWeighted->Sumw2();

        hRefDijetPtEtaCMBackwardWeighted = new TH2D("hRefDijetPtEtaCMBackwardWeighted", "Ref dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetPtEtaCMBackwardWeighted->Sumw2();       


        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} jetId cut;Reco #eta^{dijet};Ref #eta^{dijet}; Reco p_{T}^{dijet} (GeV/c)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId->Sumw2();
        hRefDijetPtEtaDphiJetId = new TH3D("hRefDijetPtEtaDphiJetId","Ref dijet info jetId;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphiJetId->Sumw2();


        //
        // Modify bins of reco <-> ref dijets
        //
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        
        hRecoDijetPtEtaRefDijetPtEta->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaRefDijetPtEta->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaRefDijetPtEtaWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaRefDijetPtEtaWeighted->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);

        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);

        hRefDijetEtaVsRecoDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        for (int i{0}; i<16; i++) {
            hRefDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        }
        for (int i{0}; i<5; i++) {
            hRefDijetEta1DOldPt[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPtBinning[i]->GetXaxis()->Set(dijetEtaOldBins, dijetEtaOldVals);
        }
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetPtEtaDphiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

        hRefDijetEtaCM->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetPtEtaDphiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetPtEtaDphiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

        hRefSelDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaDphiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetEtaCM->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaDphiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaDphiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

        hRefDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);

        hRefDijetPtEtaDphiJetId->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

        //
        // Modify bins of gen dijets
        //
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        for (int i{0}; i<16; i++) {
            hGenDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        }
        for (int i{0}; i<5; i++) {
            hGenDijetEta1DOldPt[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPtBinning[i]->GetXaxis()->Set(dijetEtaOldBins, dijetEtaOldVals);
        }
        hGenDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaDphiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

        hGenDijetEtaCM->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaDphiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaDphiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

        hGenDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);

    } // if (fIsMc)

    
    //
    // Reco inclusive jets
    //

    hRecoLeadJetAllPtVsEta = new TH2D("hRecoLeadJetAllPtVsEta","Leading jet all p_{T} vs #eta;#eta;p_{T} (GeV/c)", 
                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                        ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoLeadJetAllPtVsEta->Sumw2();
    hRecoSubLeadJetAllPtVsEta = new TH2D("hRecoSubLeadJetAllPtVsEta","Subleading jet all p_{T} vs #eta;#eta;p_{T} (GeV/c)", 
                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                        ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoSubLeadJetAllPtVsEta->Sumw2();
    hRecoLeadJetAllPtVsEtaJetIdCut = new TH2D("hRecoLeadJetAllPtVsEtaJetIdCut","Leading jet all p_{T} vs #eta with jetId selection;#eta;p_{T} (GeV/c)", 
                                            fEtaBins, fEtaRange[0], fEtaRange[1], 
                                            ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoLeadJetAllPtVsEtaJetIdCut->Sumw2();
    hRecoSubLeadJetAllPtVsEtaJetIdCut = new TH2D("hRecoSubLeadJetAllPtVsEtaJetIdCut","Subleading jet all p_{T} vs #eta with jetId selection;#eta;p_{T} (GeV/c)", 
                                                fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoSubLeadJetAllPtVsEtaJetIdCut->Sumw2();


    hRecoInclusiveJetPt = new TH1D("hRecoInclusiveJetPt","Reco jet p_{T};Reco p_{T} (GeV/c);Entries",
                                   fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveJetPt->Sumw2();
    hRecoInclusiveAllJetPtVsEta = new TH2D("hRecoInclusiveAllJetPtVsEta", "Inclusive reco jet pT vs eta;#eta;p_{T} (GeV/c)",
                                           fEtaBins, fEtaRange[0], fEtaRange[1],
                                           fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtVsEta->Sumw2();

    hRecoGoodInclusiveJetEtaLabFrame = new TH1D("hRecoGoodInclusiveJetEtaLabFrame","Reco good jet #eta in lab frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaLabFrame->Sumw2();
    hRecoGoodInclusiveJetEtaCMFrame = new TH1D("hRecoGoodInclusiveJetEtaCMFrame","Reco good jet #eta in CM frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaCMFrame->Sumw2();

    //
    // Reco jet selection algo
    //
    hRecoInclusiveJetPtVsEtaKineCut = new TH2D("hRecoInclusiveJetPtVsEtaKineCut", "Reco inclusive jet passed kine cut;#eta;p_{T} (GeV/c);Entries", 
                                               fEtaBins, fEtaRange[0], fEtaRange[1], 
                                               ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoInclusiveJetPtVsEtaKineCut->Sumw2();
    hRecoInclusiveJetPtVsEtaTrkMaxCut = new TH2D("hRecoInclusiveJetPtVsEtaTrkMaxCut", "Reco inclusive jet passed trkMax cut;#eta;p_{T} (GeV/c);Entries", 
                                                 fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                 ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoInclusiveJetPtVsEtaTrkMaxCut->Sumw2();
    hRecoInclusiveJetPtVsEtaJetIdCut = new TH2D("hRecoInclusiveJetPtVsEtaJetIdCut", "Reco inclusive jet passed jetId cut;#eta;p_{T} (GeV/c);Entries", 
                                                fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                ptShortBins, ptShortRange[0], ptShortRange[1]);
    hRecoInclusiveJetPtVsEtaJetIdCut->Sumw2();


    //
    // Reco dijets
    //
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi = new THnSparseD("hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi",
            "Reconstructed dijet and jet info;p_{T}^{dijet} (GeV/c);#eta^{dijet};#Delta#phi^{dijet} (rad);p_{T}^{Leading} (GeV/c);#eta^{Leading};#phi^{Leading} (rad);p_{T}^{Subleading} (GeV/c);#eta^{Subleading};#phi^{Subleading} (rad)",
            9,
            bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmax9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Sumw2();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted = new THnSparseD("hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted",
            "Reconstructed dijet and jet info weighted;p_{T}^{dijet} (GeV/c);#eta^{dijet};#Delta#phi^{dijet} (rad);p_{T}^{Leading} (GeV/c);#eta^{Leading};#phi^{Leading} (rad);p_{T}^{Subleading} (GeV/c);#eta^{Subleading};#phi^{Subleading} (rad)",
            9,
            bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmax9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Sumw2();
    

    hRecoPtLeadPtSublead = new TH2D("hRecoPtLeadPtSublead","Reco leading vs subleading p_{T};Reco p_{T}^{Leading} (GeV/c);Reco p_{T}^{Subleading} (GeV/c)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSublead->Sumw2();
    hRecoEtaLeadEtaSublead = new TH2D("hRecoEtaLeadEtaSublead","Reco leading vs subleading #eta;Reco #eta^{Leading};Reco #eta^{Subleading}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSublead->Sumw2();
    hRecoEtaCMLeadEtaCMSublead = new TH2D("hRecoEtaCMLeadEtaCMSublead","Reco leading vs subleading #eta in CM;Reco #eta^{Leading}_{CM};Reco #eta^{Subleading}_{CM}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaCMLeadEtaCMSublead->Sumw2();
    hRecoPtLeadPtSubleadMcReweight = new TH2D("hRecoPtLeadPtSubleadMcReweight","Reco leading vs subleading p_{T} (MC reweighted to data);Reco p_{T}^{Leading} (GeV/c);Reco p_{T}^{Subleading} (GeV/c)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSubleadMcReweight->Sumw2();
    hRecoEtaLeadEtaSubleadMcReweight = new TH2D("hRecoEtaLeadEtaSubleadMcReweight","Reco leading vs subleading #eta (MC reweighted to data);Reco #eta^{Leading};Reco #eta^{Subleading}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSubleadMcReweight->Sumw2();
    hRecoDijetEta = new TH1D("hRecoDijetEta","Reco dijet #eta;Reco #eta^{dijet};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEta->Sumw2();
    for (int i{0}; i<16; ++i) {
        hRecoDijetEta1D[i] = new TH1D(Form("hRecoDijetEta1D_%d",i),Form("Reco dijet #eta bin %d;Reco #eta^{dijet};Entries",i),
                                      fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRecoDijetEta1D[i]->Sumw2();
    }
    for (int i{0}; i<5; ++i) {
        hRecoDijetEta1DOldPt[i] = new TH1D(Form("hRecoDijetEta1DOldPt_%d",i),Form("Reco dijet #eta old p_{T} bin %d;Reco #eta^{dijet};Entries",i),
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRecoDijetEta1DOldPt[i]->Sumw2();
        hRecoDijetEta1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinning_%d",i),Form("Reco dijet #eta old p_{T} bin %d old #eta binning;Reco #eta^{dijet};Entries",i),
                                                  dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinning[i]->Sumw2();
    }

    hRecoDijetPtEta = new TH2D("hRecoDijetPtEta", "Reco dijet #eta vs p_{T};p_{T}^{ave} (GeV/c);#eta^{dijet}", 
                               fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                               fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEta->Sumw2();
    hRecoDijetPtEtaDphi = new TH3D("hRecoDijetPtEtaDphi","Reco dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphi->Sumw2();
    hRecoDijetPtEtaDphiWeighted = new TH3D("hRecoDijetPtEtaDphiWeighted","Reco dijet info;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphiWeighted->Sumw2();
    hRecoDijetEtaCM = new TH1D("hRecoDijetEtaCM","Reco dijet #eta in CM;Reco #eta^{dijet}_{CM};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEtaCM->Sumw2();
    hRecoDijetPtEtaDphiCM = new TH3D("hRecoDijetPtEtaDphiCM","Reco dijet info in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphiCM->Sumw2();
    hRecoDijetPtEtaDphiCMWeighted = new TH3D("hRecoDijetPtEtaDphiCMWeighted","Reco dijet info in CM;p_{T}^{ave} (GeV/c);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphiCMWeighted->Sumw2();


    hRecoDijetPtEtaForward = new TH2D("hRecoDijetPtEtaForward", "Reco dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaForward->Sumw2();
    hRecoDijetPtEtaBackward = new TH2D("hRecoDijetPtEtaBackward", "Reco dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaBackward->Sumw2();
    hRecoDijetPtEtaCMForward = new TH2D("hRecoDijetPtEtaCMForward", "Reco dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaCMForward->Sumw2();
    hRecoDijetPtEtaCMBackward = new TH2D("hRecoDijetPtEtaCMBackward", "Reco dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaCMBackward->Sumw2();
    
    hRecoDijetPtEtaForwardWeighted = new TH2D("hRecoDijetPtEtaForwardWeighted", "Reco dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaForwardWeighted->Sumw2();
    hRecoDijetPtEtaBackwardWeighted = new TH2D("hRecoDijetPtEtaBackwardWeighted", "Reco dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaBackwardWeighted->Sumw2();
    hRecoDijetPtEtaCMForwardWeighted = new TH2D("hRecoDijetPtEtaCMForwardWeighted", "Reco dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaCMForwardWeighted->Sumw2();

    hRecoDijetPtEtaCMBackwardWeighted = new TH2D("hRecoDijetPtEtaCMBackwardWeighted", "Reco dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEtaCMBackwardWeighted->Sumw2();


    hRecoDijetPtEtaDphiJetId  = new TH3D("hRecoDijetPtEtaDphiJetId","Reco dijet info jetId cut;p_{T}^{ave} (GeV/c);#eta^{dijet};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphiJetId->Sumw2();


    hRecoTrkMaxToJetIdDijetMatching = new TH1D("hRecoTrkMaxToJetIdDijetMatching", "Matching of dijets between trkMax and jetId selection;;Entries", 
                                               9, -0.5, 8.5);
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(7, "None found");
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(1, "Both match");
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(2, "Lead match");
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(3, "Sublead match");
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(4, "Both diff");
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(5, "TrkMax && !JetId");
    hRecoTrkMaxToJetIdDijetMatching->GetXaxis()->SetBinLabel(6, "!TrkMax && JetId");

    //
    // Modify bins of reco dijets
    //
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    for (int i{0}; i<16; i++) {
        hRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    }
    for (int i{0}; i<5; i++) {
        hRecoDijetEta1DOldPt[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPtBinning[i]->GetXaxis()->Set(dijetEtaOldBins, dijetEtaOldVals);
    }
    hRecoDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetEtaCM->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiJetId->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

    hRecoDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
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

   for (int i=0; i<4; i++) {
      hNHF[i]->Write();
      hNEmF[i]->Write();
      hNumOfConst[i]->Write();
      hMUF[i]->Write();
      hCHF[i]->Write();
      hChargedMult[i]->Write();
      hCEmF[i]->Write();
      hNumOfNeutPart[i]->Write();
    }

    // Reco jets

    if (fIsMc) {
        // Gen jets
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Write();
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Write();
        hGenInclusiveJetPt->Write();
        hGenInclusiveJetPtEta->Write();
        hGenPtLeadPtSublead->Write();
        hGenEtaLeadEtaSublead->Write();
        hGenEtaCMLeadEtaCMSublead->Write();
        hGenPtLeadPtSubleadMcReweight->Write();
        hGenEtaLeadEtaSubleadMcReweight->Write();
        hGenDijetEta->Write();
        for (int i{0}; i<16; i++) {
            hGenDijetEta1D[i]->Write();
        }
        for (int i{0}; i<5; i++) {
            hGenDijetEta1DOldPt[i]->Write();
            hGenDijetEta1DOldPtBinning[i]->Write();
        }
        hGenDijetPtEtaDphi->Write();
        hGenDijetPtEtaDphiWeighted->Write();
        hGenDijetEtaCM->Write();
        hGenDijetPtEtaDphiCM->Write();
        hGenDijetPtEtaDphiCMWeighted->Write();
        hGenDijetPtEtaForward->Write();
        hGenDijetPtEtaBackward->Write();
        hGenDijetPtEtaCMForward->Write();
        hGenDijetPtEtaCMBackward->Write();
        hGenDijetPtEtaForwardWeighted->Write();
        hGenDijetPtEtaBackwardWeighted->Write();
        hGenDijetPtEtaCMForwardWeighted->Write();
        hGenDijetPtEtaCMBackwardWeighted->Write();

        hGenGoodInclusiveJetEtaLabFrame->Write();
        hGenGoodInclusiveJetEtaCMFrame->Write();
        

        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hJESInclusiveJetPtEtaPhi->Write();
        hJESInclusiveJetPtEtaPhiWeighted->Write();
        hRecoMatchedPtEta->Write();
        hRecoInclusiveMatchedJetPtVsEta->Write();
        hRecoInclusiveUnmatchedJetPtVsEta->Write();

        
        hRecoLeadJetMatchedPtVsEta->Write();
        hRecoLeadJetUnmatchedPtVsEta->Write();
        hRecoSubLeadJetMatchedPtVsEta->Write();
        hRecoSubLeadJetUnmatchedPtVsEta->Write();

        
        hRecoLeadJetMatchedPtVsEtaJetIdCut->Write();
        hRecoLeadJetUnmatchedPtVsEtaJetIdCut->Write();
        hRecoSubLeadJetMatchedPtVsEtaJetIdCut->Write();
        hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut->Write();

        hRecoInclusiveJetJESPtEtaPhiKineCut->Write();
        hRecoInclusiveJetDEtaPtEtaKineCut->Write();
        hRecoInclusiveMatchedJetPtVsEtaKineCut->Write();
        hRecoInclusiveUnmatchedJetPtVsEtaKineCut->Write();
        hRecoInclusiveJetRefPtVsEtaKineCut->Write();

        hRecoInclusiveJetJESPtEtaPhiTrkMaxCut->Write();
        hRecoInclusiveJetDEtaPtEtaTrkMaxCut->Write();
        hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut->Write();
        hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut->Write();
        hRecoInclusiveJetRefPtVsEtaTrkMaxCut->Write();

        hRecoInclusiveJetJESPtEtaPhiJetIdCut->Write();
        hRecoInclusiveJetDEtaPtEtaJetIdCut->Write();
        hRecoInclusiveMatchedJetPtVsEtaJetIdCut->Write();
        hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut->Write();
        hRecoInclusiveJetRefPtVsEtaJetIdCut->Write();

        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Write();
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Write();
        hRecoDijetPtEtaRefDijetPtEta->Write();
        hRecoDijetPtEtaRefDijetPtEtaWeighted->Write();
        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Write();

        hRefDijetEta->Write();
        for (int i{0}; i<16; i++) {
            hRefDijetEta1D[i]->Write();
        }
        for (int i{0}; i<5; i++) {
            hRefDijetEta1DOldPt[i]->Write();
            hRefDijetEta1DOldPtBinning[i]->Write();
        }
        hRefDijetPtEtaDphi->Write();
        hRefDijetPtEtaDphiWeighted->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Write();
        hRefDijetEtaVsRecoDijetEta->Write();
        hRefDijetEtaCM->Write();
        hRefDijetPtEtaDphiCM->Write();
        hRefDijetPtEtaDphiCMWeighted->Write();
        hRefDijetPtEtaForward->Write();
        hRefDijetPtEtaBackward->Write();
        hRefDijetPtEtaCMForward->Write();
        hRefDijetPtEtaCMBackward->Write();
        hRefDijetPtEtaForwardWeighted->Write();
        hRefDijetPtEtaBackwardWeighted->Write();
        hRefDijetPtEtaCMForwardWeighted->Write();
        hRefDijetPtEtaCMBackwardWeighted->Write();

        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Write();

        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId->Write();
        hRefDijetPtEtaDphiJetId->Write();

        hRefSelDijetEta->Write();
        hRefSelDijetPtEtaDphi->Write();
        hRefSelDijetPtEtaDphiWeighted->Write();
        hRefSelDijetEtaCM->Write();
        hRefSelDijetPtEtaDphiCM->Write();
        hRefSelDijetPtEtaDphiCMWeighted->Write();

        hRefInclusiveJetPt->Write();
        hRefInclusiveJetPtEta->Write();
        hRefPtLeadPtSublead->Write();
        hRefEtaLeadEtaSublead->Write();
        hRefEtaCMLeadEtaCMSublead->Write();
        hRefPtLeadPtSubleadMcReweight->Write();
        hRefEtaLeadEtaSubleadMcReweight->Write();
    } // if ( fIsMc )

    hRecoLeadJetAllPtVsEta->Write();
    hRecoSubLeadJetAllPtVsEta->Write();
    hRecoLeadJetAllPtVsEtaJetIdCut->Write();
    hRecoSubLeadJetAllPtVsEtaJetIdCut->Write();

    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Write();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Write();
    hRecoInclusiveJetPt->Write();
    hRecoPtLeadPtSublead->Write();
    hRecoEtaLeadEtaSublead->Write();
    hRecoEtaCMLeadEtaCMSublead->Write();
    hRecoPtLeadPtSubleadMcReweight->Write();
    hRecoEtaLeadEtaSubleadMcReweight->Write();
    hRecoDijetPtEta->Write();
    hRecoDijetEta->Write();
    for (int i{0}; i<16; i++) {
        hRecoDijetEta1D[i]->Write();
    }
    for (int i{0}; i<5; i++) {
        hRecoDijetEta1DOldPt[i]->Write();
        hRecoDijetEta1DOldPtBinning[i]->Write();
    }
    hRecoDijetPtEtaDphi->Write();
    hRecoDijetPtEtaDphiWeighted->Write();
    hRecoDijetEtaCM->Write();
    hRecoDijetPtEtaDphiCM->Write();
    hRecoDijetPtEtaDphiCMWeighted->Write();
    hRecoDijetPtEtaForward->Write();
    hRecoDijetPtEtaBackward->Write();
    hRecoDijetPtEtaCMForward->Write();
    hRecoDijetPtEtaCMBackward->Write();
    hRecoDijetPtEtaForwardWeighted->Write();
    hRecoDijetPtEtaBackwardWeighted->Write();
    hRecoDijetPtEtaCMForwardWeighted->Write();
    hRecoDijetPtEtaCMBackwardWeighted->Write();


    hRecoGoodInclusiveJetEtaLabFrame->Write();
    hRecoGoodInclusiveJetEtaCMFrame->Write();

    hRecoDijetPtEtaDphiJetId->Write();
    hRecoTrkMaxToJetIdDijetMatching->Write();

    hRecoInclusiveAllJetPtVsEta->Write();
    hRecoInclusiveJetPtVsEtaKineCut->Write();
    hRecoInclusiveJetPtVsEtaTrkMaxCut->Write();
    hRecoInclusiveJetPtVsEtaJetIdCut->Write();
}