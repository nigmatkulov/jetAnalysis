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

//________________
HistoManagerDiJet::HistoManagerDiJet() :
    BaseHistoManager(),
    //
    // Event histograms
    // 
    hVz{nullptr},
    hVzWeighted{nullptr},
    hPtHat{nullptr},
    hPtHatWeighted{nullptr},
    hHiBin{nullptr},
    hHiBinWeighted{nullptr},
    hVzGenDijetLab{nullptr},
    hVzGenDijetLabWeighted{nullptr},
    hHiBinGenDijetLab{nullptr},
    hHiBinGenDijetLabWeighted{nullptr},
    hVzGenDijetCM{nullptr},
    hVzGenDijetCMWeighted{nullptr},
    hHiBinGenDijetCM{nullptr},
    hHiBinGenDijetCMWeighted{nullptr},
    hVzRecoDijetLab{nullptr},
    hVzRecoDijetLabWeighted{nullptr},
    hHiBinRecoDijetLab{nullptr},
    hHiBinRecoDijetLabWeighted{nullptr},
    hVzRecoDijetCM{nullptr},
    hVzRecoDijetCMWeighted{nullptr},
    hHiBinRecoDijetCM{nullptr},
    hHiBinRecoDijetCMWeighted{nullptr},
    hVzRefSelDijetLab{nullptr},
    hVzRefSelDijetLabWeighted{nullptr},
    hHiBinRefSelDijetLab{nullptr},
    hHiBinRefSelDijetLabWeighted{nullptr},
    hVzRefSelDijetCM{nullptr},
    hVzRefSelDijetCMWeighted{nullptr},
    hHiBinRefSelDijetCM{nullptr},
    hHiBinRefSelDijetCMWeighted{nullptr},

    //
    // Gen jets
    //

    hGenJetCollectionSize{nullptr},
    hGenVsRecoJetCollectionSize{nullptr},

    hGenLeadJetPtOverPtHatVsLeadJetPt{nullptr},
    hGenLeadJetPtOverPtHatVsLeadJetPtWeighted{nullptr},
    hGenDijetPtOverPtHatVsDijetPt{nullptr},
    hGenDijetPtOverPtHatVsDijetPtWeighted{nullptr},
    hGenDijetPtAveOverPtHatVsDijetPtAve{nullptr},
    hGenDijetPtAveOverPtHatVsDijetPtAveWeighted{nullptr},

    // hGenDijetInfo{nullptr},
    hGenInclusiveJetPt{nullptr},
    hGenInclusiveJetEta{nullptr},
    hGenInclusiveJetEtaUnweighted{nullptr},
    hGenInclusiveJetPtEta{nullptr},
    hGenInclusiveJetPtEtaCM{nullptr},
    hGenInclusiveJetPtEtaPtHat{nullptr},
    hGenLeadJetPtEta{nullptr},
    hGenLeadJetPtEtaCM{nullptr},
    hGenLeadJetPtEtaPtHat{nullptr},
    hGenSubLeadJetPtEta{nullptr},
    hGenSubLeadJetPtEtaCM{nullptr},
    hGenSubLeadJetPtEtaPtHat{nullptr},
    hGenPtLeadPtSublead{nullptr},
    hGenEtaLeadEtaSublead{nullptr},
    hGenEtaCMLeadEtaCMSublead{nullptr},
    hGenPtLeadPtSubleadMcReweight{nullptr},
    hGenEtaLeadEtaSubleadMcReweight{nullptr},

    hGenDijetEta{nullptr},
    hGenDijetPtEta{nullptr},
    hGenDijetPtEtaWeighted{nullptr},
    hGenDijetPtEtaCMInLab{nullptr},
    hGenDijetEtaCM{nullptr},
    hGenDijetPtEtaCM{nullptr},
    hGenDijetPtEtaCMWeighted{nullptr},
    hGenDijetPtEtaLabInCM{nullptr},
    hGenDijetPtEtaForward{nullptr},
    hGenDijetPtEtaBackward{nullptr},
    hGenDijetPtEtaForwardCMInLab{nullptr},
    hGenDijetPtEtaBackwardCMInLab{nullptr},
    hGenDijetPtEtaCMForward{nullptr},
    hGenDijetPtEtaCMBackward{nullptr},
    hGenDijetPtEtaForwardLabInCM{nullptr},
    hGenDijetPtEtaBackwardLabInCM{nullptr},
    hGenDijetPtEtaForwardWeighted{nullptr},
    hGenDijetPtEtaBackwardWeighted{nullptr},
    hGenDijetPtEtaCMForwardWeighted{nullptr},
    hGenDijetPtEtaCMBackwardWeighted{nullptr},
    hGenDijetPtEtaForwardArr{nullptr},
    hGenDijetPtEtaBackwardArr{nullptr},

    hGenDijetPtAveLeadPtSubLeadPt{nullptr},
    hGenDijetPtAveLeadPtSubLeadPtCM{nullptr},
    hGenDijetPtAveLeadEtaSubLeadEta{nullptr},
    hGenDijetPtAveLeadEtaSubLeadEtaCM{nullptr},
    hGenDijetEtaLeadEtaSubLeadEta{nullptr},
    hGenDijetEtaLeadEtaSubLeadEtaCM{nullptr},

    hGenGoodInclusiveJetEtaLabFrame{nullptr},
    hGenGoodInclusiveJetEtaCMFrame{nullptr},
    hGenInclusiveDijetDetaCM{nullptr},
    hGenInclusiveDijetDetaCMWeighted{nullptr},
    hGenInclusiveDijetDetaCMPt{nullptr},
    hGenInclusiveDijetDetaCMPtWeighted{nullptr},
    hGenInclusiveDijetEtaDetaCMPt{nullptr},
    hGenInclusiveDijetEtaDetaCMPtWeighted{nullptr},
    hGenInclusiveDijetXPb{nullptr},
    hGenInclusiveDijetXPbWeighted{nullptr},
    hGenInclusiveDijetXp{nullptr},
    hGenInclusiveDijetXpWeighted{nullptr},
    hGenInclusiveDijetXPbOverXp{nullptr},
    hGenInclusiveDijetXPbOverXpWeighted{nullptr},
    hGenInclusiveDijetXPbOverXpEta{nullptr},
    hGenInclusiveDijetXPbOverXpEtaWeighted{nullptr},

    hGenSelectedDijetDetaCM{nullptr},
    hGenSelectedDijetDetaCMWeighted{nullptr},
    hGenSelectedDijetDetaCMPt{nullptr},
    hGenSelectedDijetDetaCMPtWeighted{nullptr},
    hGenSelectedDijetEtaDetaCMPt{nullptr},
    hGenSelectedDijetEtaDetaCMPtWeighted{nullptr},
    hGenSelectedDijetXPb{nullptr},
    hGenSelectedDijetXPbWeighted{nullptr},
    hGenSelectedDijetXp{nullptr},
    hGenSelectedDijetXpWeighted{nullptr},
    hGenSelectedDijetXPbOverXp{nullptr},
    hGenSelectedDijetXPbOverXpWeighted{nullptr},
    hGenSelectedDijetXPbOverXpEta{nullptr},
    hGenSelectedDijetXPbOverXpEtaWeighted{nullptr},

    //
    // Reco jets
    //

    hRecoJetCollectionSize{nullptr},
    hRecoInclusiveJetNHF{nullptr}, 
    hRecoInclusiveJetNEmF{nullptr}, 
    hRecoInclusiveJetNumOfConst{nullptr}, 
    hRecoInclusiveJetMUF{nullptr},
    hRecoInclusiveJetCHF{nullptr}, 
    hRecoInclusiveJetChargedMult{nullptr}, 
    hRecoInclusiveJetCEmF{nullptr}, 
    hRecoInclusiveJetNumOfNeutPart{nullptr},

    hRecoInclusiveAllJetPtRawEta{nullptr},
    // hRecoDijetInfo{nullptr},

    hRecoInclusiveAllJetPt{nullptr},
    hRecoInclusiveAllJetEta{nullptr},
    hRecoInclusiveAllJetEtaUnweighted{nullptr},
    hRecoInclusiveAllJetPtEta{nullptr},
    hRecoInclusiveAllJetPtEtaCM{nullptr},
    hRecoInclusiveAllJetPtRawEtaStdBins{nullptr},
    hRecoInclusiveAllJetPtEtaStdBins{nullptr},
    hRecoInclusiveAllJetPtEtaPtHat{nullptr},
    hRecoInclusiveAllJetPtRawEtaPtHat{nullptr},
    hRecoInclusiveMatchedJetPtEtaPtHat{nullptr},
    hRecoInclusiveUnmatchedJetPtEtaPtHat{nullptr},

    hRecoInclusiveJetEtaRun{nullptr},
    hRecoLeadJetEtaRun{nullptr},
    hRecoSubLeadJetEtaRun{nullptr},

    hRecoPtLeadPtSublead{nullptr},
    hRecoEtaLeadEtaSublead{nullptr},
    hRecoEtaCMLeadEtaCMSublead{nullptr},
    hRecoPtLeadPtSubleadMcReweight{nullptr},
    hRecoEtaLeadEtaSubleadMcReweight{nullptr},
    hRecoDijetEtaCMRun{nullptr},

    hRecoDijetPtEtaForward{nullptr},
    hRecoDijetPtEtaBackward{nullptr},
    hRecoDijetPtEtaForwardCMInLab{nullptr},
    hRecoDijetPtEtaBackwardCMInLab{nullptr},
    hRecoDijetPtEtaCMForward{nullptr},
    hRecoDijetPtEtaCMBackward{nullptr},
    hRecoDijetPtEtaForwardLabInCM{nullptr},
    hRecoDijetPtEtaBackwardLabInCM{nullptr},
    hRecoDijetPtEtaForwardWeighted{nullptr},
    hRecoDijetPtEtaBackwardWeighted{nullptr},
    hRecoDijetPtEtaCMForwardWeighted{nullptr},
    hRecoDijetPtEtaCMBackwardWeighted{nullptr},
    hRecoDijetPtRawEtaForwardArr{nullptr},
    hRecoDijetPtRawEtaBackwardArr{nullptr},
    hRecoDijetPtEtaForwardArr{nullptr},
    hRecoDijetPtEtaBackwardArr{nullptr},
    hRecoDijetPtEta{nullptr},
    hRecoDijetPtEtaWeighted{nullptr},
    hRecoDijetPtEtaCMInLab{nullptr},
    hRecoDijetPtEtaCM{nullptr},
    hRecoDijetPtEtaCMWeighted{nullptr},
    hRecoDijetPtEtaLabInCM{nullptr},
    hRecoDijetPtEtaMatched{nullptr},
    hRecoDijetPtEtaCMMatched{nullptr},

    hRecoDijetPtAveLeadPtSubLeadPt{nullptr},
    hRecoDijetPtAveLeadPtSubLeadPtCM{nullptr},
    hRecoDijetPtAveLeadEtaSubLeadEta{nullptr},
    hRecoDijetPtAveLeadEtaSubLeadEtaCM{nullptr},
    hRecoDijetEtaLeadEtaSubLeadEta{nullptr},
    hRecoDijetEtaLeadEtaSubLeadEtaCM{nullptr},

    hRecoDijetLeadPtEta{nullptr},
    hRecoDijetLeadPtEtaStdBins{nullptr},
    hRecoDijetSubLeadPtEta{nullptr},
    hRecoDijetSubLeadPtEtaStdBins{nullptr},
    hRecoDijetXj{nullptr},
    hRecoDijetXjCM{nullptr},

    hRecoLeadAllJetPtEta{nullptr},
    hRecoLeadAllJetPtEtaCM{nullptr},
    hRecoLeadAllJetPtRawEtaStdBins{nullptr},
    hRecoLeadAllJetPtEtaStdBins{nullptr},
    hRecoLeadAllJetPtEtaPtHat{nullptr},
    hRecoSubLeadAllJetPtEta{nullptr},
    hRecoSubLeadAllJetPtEtaCM{nullptr},
    hRecoSubLeadAllJetPtRawEtaStdBins{nullptr},
    hRecoSubLeadAllJetPtEtaStdBins{nullptr},
    hRecoSubLeadAllJetPtEtaPtHat{nullptr},
    hRecoGoodInclusiveJetEtaLabFrame{nullptr},
    hRecoGoodInclusiveJetEtaCMFrame{nullptr},
    hRecoDijetEta{nullptr},
    hRecoDijetEtaCM{nullptr},

    //
    // Ref jet histograms
    //

    hRecoInclusiveJetReco2Ref{nullptr},
    hRecoLeadJetReco2Ref{nullptr},
    hRecoSubLeadJetReco2Ref{nullptr},
    hJESInclusiveJetPtEtaPhi{nullptr},

    hRecoLeadJetPtOverPtHatVsLeadJetPt{nullptr},
    hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted{nullptr},
    hRecoDijetPtOverPtHatVsDijetPt{nullptr},
    hRecoDijetPtOverPtHatVsDijetPtWeighted{nullptr},
    hRecoDijetPtAveOverPtHatVsDijetPtAve{nullptr},
    hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted{nullptr},

    // hRecoInclusiveJetJECFactorVsPtEta{nullptr},
    // hRecoInclusiveJetJEC2FactorVsPtEta{nullptr},
    // hRecoInclusiveJetPtRawOverPtRefVsPtEta{nullptr},
    // hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning{nullptr},
    // hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning{nullptr},

    hInclusiveJetJESVsPtGen{nullptr},
    hInclusiveJetJESGenPtGenEtaPtHatWeighted{nullptr},
    hInclusiveJetJESGenPtGenEtaCMPtHatWeighted{nullptr},
    hInclusiveJetJESRecoPtRecoEtaPtHatWeighted{nullptr},

    hLeadJetJESGenPtEtaPtHatWeighted{nullptr},
    hLeadJetJESGenPtEtaCMPtHatWeighted{nullptr},
    hSubLeadJetJESGenPtEtaPtHatWeighted{nullptr},
    hSubLeadJetJESGenPtEtaCMPtHatWeighted{nullptr},

    hRecoInclusiveMatchedJetPt{nullptr},
    hRecoInclusiveMatchedJetPtEta{nullptr},
    hRecoInclusiveUnmatchedJetPtEta{nullptr},
    hRecoLeadMatchedJetPtEta{nullptr},
    hRecoLeadMatchedJetPtEtaPtHat{nullptr},
    hRecoLeadUnmatchedJetPtEta{nullptr},
    hRecoLeadUnmatchedJetPtEtaPtHat{nullptr},
    hRecoSubLeadMatchedJetPtEta{nullptr},
    hRecoSubLeadMatchedJetPtEtaPtHat{nullptr},
    hRecoSubLeadUnmatchedJetPtEta{nullptr},
    hRecoSubLeadUnmatchedJetPtEtaPtHat{nullptr},

    hRefInclusiveJetPt{nullptr},
    hRefInclusiveJetEta{nullptr},
    hRefInclusiveJetEtaUnweighted{nullptr},
    hRefInclusiveJetPtEta{nullptr},
    hRefInclusiveJetPtEtaCM{nullptr},
    hRefInclusiveJetPtEtaPtHat{nullptr},
    hRefLeadJetPtEta{nullptr},
    hRefLeadJetPtEtaCM{nullptr},
    hRefLeadJetPtEtaPtHat{nullptr},
    hRefLeadUnswappedJetPtEta{nullptr},
    hRefLeadUnswappedJetPtEtaPtHat{nullptr},
    hRefSubLeadJetPtEta{nullptr},
    hRefSubLeadJetPtEtaCM{nullptr},
    hRefSubLeadJetPtEtaPtHat{nullptr},
    hRefSubLeadUnswappedJetPtEta{nullptr},
    hRefSubLeadUnswappedJetPtEtaPtHat{nullptr},

    // hReco2RefFull{nullptr},
    hRecoDijetPtEtaRefDijetPtEta{nullptr},
    hRecoDijetPtEtaRefDijetPtEtaWeighted{nullptr},
    
    // hRefSel2RecoFull{nullptr},

    hRefDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted{nullptr},
    hRefDijetPtEta{nullptr},
    hRefDijetPtEtaWeighted{nullptr},
    hRefDijetPtEtaCMInLab{nullptr},

    hRefDijetPtEtaForward{nullptr},
    hRefDijetPtEtaBackward{nullptr},
    hRefDijetPtEtaForwardCMInLab{nullptr},
    hRefDijetPtEtaBackwardCMInLab{nullptr},
    hRefDijetPtEtaCMForward{nullptr},
    hRefDijetPtEtaCMBackward{nullptr},
    hRefDijetPtEtaForwardLabInCM{nullptr},
    hRefDijetPtEtaBackwardLabInCM{nullptr},
    hRefDijetPtEtaForwardWeighted{nullptr},
    hRefDijetPtEtaBackwardWeighted{nullptr},
    hRefDijetPtEtaCMForwardWeighted{nullptr},
    hRefDijetPtEtaCMBackwardWeighted{nullptr},

    hRefDijetEtaCM{nullptr},
    hRefDijetPtEtaCM{nullptr},
    hRefDijetPtEtaCMWeighted{nullptr},
    hRefDijetPtEtaLabInCM{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted{nullptr},

    hRefPtLeadPtSublead{nullptr},
    hRefEtaLeadEtaSublead{nullptr},
    hRefEtaCMLeadEtaCMSublead{nullptr},
    hRefPtLeadPtSubleadMcReweight{nullptr},
    hRefEtaLeadEtaSubleadMcReweight{nullptr},

    //
    // Ref-selected jet histograms
    //

    hRefSelInclusiveJetPt{nullptr},
    hRefSelInclusiveJetPtEta{nullptr},
    hRefSelInclusiveJetPtEtaPtHat{nullptr},
    hRefSelLeadJetPtEta{nullptr},
    hRefSelLeadJetPtEtaPtHat{nullptr},
    hRefSelSubLeadJetPtEta{nullptr},
    hRefSelSubLeadJetPtEtaPtHat{nullptr},

    hRefSelDijetEta{nullptr},
    hRefSelDijetPtEta{nullptr},
    hRefSelDijetPtEtaWeighted{nullptr},
    hRefSelDijetPtEtaCMInLab{nullptr},
    hRefSelDijetEtaCM{nullptr},
    hRefSelDijetPtEtaCM{nullptr},
    hRefSelDijetPtEtaCMWeighted{nullptr},
    hRefSelDijetPtEtaLabInCM{nullptr},

    //
    // Variables
    //
    fIsMc{false}, fUseVariableBinning{true},
    fPtBins{50}, fPtRange{5., 505.}, 
    fEtaBins{36}, fEtaRange{-3.6, 3.6},
    fPhiBins{16}, fPhiRange{-TMath::Pi(), TMath::Pi()},
    fDijetPtBins{45}, fDijetPtRange{50., 500.},
    fDijetEtaBins{32}, fDijetEtaRange{-3.2, 3.2},
    fDijetEtaFBBins{32}, fDijetEtaFBRange{0., 3.2},
    // fDijetDphiBins{8}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},
    fPtHatBins{100}, fPtHatRange{15., 1015.},
    fFracBins{100}, fFracRange{0., 1.},
    fMultBins{32}, fMultRange{-0.5, 31.5} { 
    
    // Empty
}

//________________
HistoManagerDiJet::~HistoManagerDiJet() {

    // Event histograms
    if (hVz) delete hVz;
    if (hVzWeighted) delete hVzWeighted;
    if (hPtHat) delete hPtHat;
    if (hPtHatWeighted) delete hPtHatWeighted;
    if (hHiBin) delete hHiBin;
    if (hHiBinWeighted) delete hHiBinWeighted;

    if (hVzRecoDijetLab) delete hVzRecoDijetLab;
    if (hVzRecoDijetLabWeighted) delete hVzRecoDijetLabWeighted;
    if (hHiBinRecoDijetLab) delete hHiBinRecoDijetLab;
    if (hHiBinRecoDijetLabWeighted) delete hHiBinRecoDijetLabWeighted;
    if (hVzRecoDijetCM) delete hVzRecoDijetCM;
    if (hVzRecoDijetCMWeighted) delete hVzRecoDijetCMWeighted;
    if (hHiBinRecoDijetCM) delete hHiBinRecoDijetCM;
    if (hHiBinRecoDijetCMWeighted) delete hHiBinRecoDijetCMWeighted;

    if ( fIsMc ) {
        if (hVzGenDijetLab) delete hVzGenDijetLab;
        if (hVzGenDijetLabWeighted) delete hVzGenDijetLabWeighted;
        if (hHiBinGenDijetLab) delete hHiBinGenDijetLab;
        if (hHiBinGenDijetLabWeighted) delete hHiBinGenDijetLabWeighted;
        if (hVzGenDijetCM) delete hVzGenDijetCM;
        if (hVzGenDijetCMWeighted) delete hVzGenDijetCMWeighted;
        if (hHiBinGenDijetCM) delete hHiBinGenDijetCM;
        if (hHiBinGenDijetCMWeighted) delete hHiBinGenDijetCMWeighted;

        if (hVzRefSelDijetLab) delete hVzRefSelDijetLab;
        if (hVzRefSelDijetLabWeighted) delete hVzRefSelDijetLabWeighted;
        if (hHiBinRefSelDijetLab) delete hHiBinRefSelDijetLab;
        if (hHiBinRefSelDijetLabWeighted) delete hHiBinRefSelDijetLabWeighted;
        if (hVzRefSelDijetCM) delete hVzRefSelDijetCM;
        if (hVzRefSelDijetCMWeighted) delete hVzRefSelDijetCMWeighted;
        if (hHiBinRefSelDijetCM) delete hHiBinRefSelDijetCM;
        if (hHiBinRefSelDijetCMWeighted) delete hHiBinRefSelDijetCMWeighted;
    }

    if ( fIsMc ) {
        // Gen histograms
        if (hGenJetCollectionSize) delete hGenJetCollectionSize;
        if (hGenVsRecoJetCollectionSize) delete hGenVsRecoJetCollectionSize;
        if (hGenLeadJetPtOverPtHatVsLeadJetPt) delete hGenLeadJetPtOverPtHatVsLeadJetPt;
        if (hGenLeadJetPtOverPtHatVsLeadJetPtWeighted) delete hGenLeadJetPtOverPtHatVsLeadJetPtWeighted;
        if (hGenDijetPtOverPtHatVsDijetPt) delete hGenDijetPtOverPtHatVsDijetPt;
        if (hGenDijetPtOverPtHatVsDijetPtWeighted) delete hGenDijetPtOverPtHatVsDijetPtWeighted;
        if (hGenDijetPtAveOverPtHatVsDijetPtAve) delete hGenDijetPtAveOverPtHatVsDijetPtAve;
        if (hGenDijetPtAveOverPtHatVsDijetPtAveWeighted) delete hGenDijetPtAveOverPtHatVsDijetPtAveWeighted;
        // if (hGenDijetInfo) delete hGenDijetInfo;
        if (hGenInclusiveJetPt) delete hGenInclusiveJetPt;
        if (hGenInclusiveJetEta) delete hGenInclusiveJetEta;
        if (hGenInclusiveJetEtaUnweighted) delete hGenInclusiveJetEtaUnweighted;
        if (hGenInclusiveJetPtEta) delete hGenInclusiveJetPtEta;
        if (hGenInclusiveJetPtEtaPtHat) delete hGenInclusiveJetPtEtaPtHat;
        if (hGenLeadJetPtEta) delete hGenLeadJetPtEta;
        if (hGenLeadJetPtEtaPtHat) delete hGenLeadJetPtEtaPtHat;
        if (hGenSubLeadJetPtEta) delete hGenSubLeadJetPtEta;
        if (hGenSubLeadJetPtEtaPtHat) delete hGenSubLeadJetPtEtaPtHat;
        if (hGenPtLeadPtSublead) delete hGenPtLeadPtSublead;
        if (hGenEtaLeadEtaSublead) delete hGenEtaLeadEtaSublead;
        if (hGenEtaCMLeadEtaCMSublead) delete hGenEtaCMLeadEtaCMSublead;
        if (hGenPtLeadPtSubleadMcReweight) delete hGenPtLeadPtSubleadMcReweight;
        if (hGenEtaLeadEtaSubleadMcReweight) delete hGenEtaLeadEtaSubleadMcReweight;
        if (hGenDijetEta) delete hGenDijetEta;
        if (hGenDijetPtEta) delete hGenDijetPtEta;
        if (hGenDijetPtEtaWeighted) delete hGenDijetPtEtaWeighted;
        if (hGenDijetPtEtaCMInLab) delete hGenDijetPtEtaCMInLab;
        if (hGenDijetEtaCM) delete hGenDijetEtaCM;
        if (hGenDijetPtEtaCM) delete hGenDijetPtEtaCM;
        if (hGenDijetPtEtaCMWeighted) delete hGenDijetPtEtaCMWeighted;
        if (hGenDijetPtEtaLabInCM) delete hGenDijetPtEtaLabInCM;
        if (hGenDijetPtEtaForward) delete hGenDijetPtEtaForward;
        if (hGenDijetPtEtaBackward) delete hGenDijetPtEtaBackward;
        if (hGenDijetPtEtaForwardCMInLab) delete hGenDijetPtEtaForwardCMInLab;
        if (hGenDijetPtEtaBackwardCMInLab) delete hGenDijetPtEtaBackwardCMInLab;
        if (hGenDijetPtEtaCMForward) delete hGenDijetPtEtaCMForward;
        if (hGenDijetPtEtaCMBackward) delete hGenDijetPtEtaCMBackward;
        if (hGenDijetPtEtaForwardLabInCM) delete hGenDijetPtEtaForwardLabInCM;
        if (hGenDijetPtEtaBackwardLabInCM) delete hGenDijetPtEtaBackwardLabInCM;
        if (hGenDijetPtEtaForwardWeighted) delete hGenDijetPtEtaForwardWeighted;
        if (hGenDijetPtEtaBackwardWeighted) delete hGenDijetPtEtaBackwardWeighted;
        if (hGenDijetPtEtaCMForwardWeighted) delete hGenDijetPtEtaCMForwardWeighted;
        if (hGenDijetPtEtaCMBackwardWeighted) delete hGenDijetPtEtaCMBackwardWeighted;
        for (int i{0}; i<6; i++) {
            if (hGenDijetPtEtaForwardArr[i]) delete hGenDijetPtEtaForwardArr[i];
            if (hGenDijetPtEtaBackwardArr[i]) delete hGenDijetPtEtaBackwardArr[i];
        }
        if (hGenDijetPtAveLeadPtSubLeadPt) delete hGenDijetPtAveLeadPtSubLeadPt;
        if (hGenDijetPtAveLeadPtSubLeadPtCM) delete hGenDijetPtAveLeadPtSubLeadPtCM;
        if (hGenDijetPtAveLeadEtaSubLeadEta) delete hGenDijetPtAveLeadEtaSubLeadEta;
        if (hGenDijetPtAveLeadEtaSubLeadEtaCM) delete hGenDijetPtAveLeadEtaSubLeadEtaCM;
        if (hGenDijetEtaLeadEtaSubLeadEta) delete hGenDijetEtaLeadEtaSubLeadEta;
        if (hGenDijetEtaLeadEtaSubLeadEtaCM) delete hGenDijetEtaLeadEtaSubLeadEtaCM;
        if (hGenGoodInclusiveJetEtaLabFrame) delete hGenGoodInclusiveJetEtaLabFrame;
        if (hGenGoodInclusiveJetEtaCMFrame) delete hGenGoodInclusiveJetEtaCMFrame;
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
    } // if ( fIsMc )

    // Reco jet histograms
    if (hRecoJetCollectionSize) delete hRecoJetCollectionSize;
    for (int i{0}; i<4; i++) {
        if (hRecoInclusiveJetNHF[i]) delete hRecoInclusiveJetNHF[i];
        if (hRecoInclusiveJetNEmF[i]) delete hRecoInclusiveJetNEmF[i];
        if (hRecoInclusiveJetNumOfConst[i]) delete hRecoInclusiveJetNumOfConst[i];
        if (hRecoInclusiveJetMUF[i]) delete hRecoInclusiveJetMUF[i];
        if (hRecoInclusiveJetCHF[i]) delete hRecoInclusiveJetCHF[i];
        if (hRecoInclusiveJetChargedMult[i]) delete hRecoInclusiveJetChargedMult[i];
        if (hRecoInclusiveJetCEmF[i]) delete hRecoInclusiveJetCEmF[i];
        if (hRecoInclusiveJetNumOfNeutPart[i]) delete hRecoInclusiveJetNumOfNeutPart[i];
    }

    if (hRecoInclusiveAllJetPtRawEta) delete hRecoInclusiveAllJetPtRawEta;
    // if (hRecoDijetInfo) delete hRecoDijetInfo;
    if (hRecoInclusiveAllJetPt) delete hRecoInclusiveAllJetPt;
    if (hRecoInclusiveAllJetEta) delete hRecoInclusiveAllJetEta;
    if (hRecoInclusiveAllJetEtaUnweighted) delete hRecoInclusiveAllJetEtaUnweighted;
    if (hRecoInclusiveAllJetPtEta) delete hRecoInclusiveAllJetPtEta;
    if (hRecoInclusiveAllJetPtEtaCM) delete hRecoInclusiveAllJetPtEtaCM;
    if (hRecoInclusiveAllJetPtRawEtaStdBins) delete hRecoInclusiveAllJetPtRawEtaStdBins;
    if (hRecoInclusiveAllJetPtEtaStdBins) delete hRecoInclusiveAllJetPtEtaStdBins;
    if (hRecoInclusiveAllJetPtEtaPtHat) delete hRecoInclusiveAllJetPtEtaPtHat;
    if (hRecoInclusiveAllJetPtRawEtaPtHat) delete hRecoInclusiveAllJetPtRawEtaPtHat;
    if (hRecoInclusiveMatchedJetPtEtaPtHat) delete hRecoInclusiveMatchedJetPtEtaPtHat;
    if (hRecoInclusiveUnmatchedJetPtEtaPtHat) delete hRecoInclusiveUnmatchedJetPtEtaPtHat;
    for (int i{0}; i<6; ++i) {
        if (hRecoInclusiveJetEtaRun[i]) delete hRecoInclusiveJetEtaRun[i];
        if (hRecoLeadJetEtaRun[i]) delete hRecoLeadJetEtaRun[i];
        if (hRecoSubLeadJetEtaRun[i]) delete hRecoSubLeadJetEtaRun[i];
    }
    if (hRecoPtLeadPtSublead) delete hRecoPtLeadPtSublead;
    if (hRecoEtaLeadEtaSublead) delete hRecoEtaLeadEtaSublead;
    if (hRecoEtaCMLeadEtaCMSublead) delete hRecoEtaCMLeadEtaCMSublead;
    if (hRecoPtLeadPtSubleadMcReweight) delete hRecoPtLeadPtSubleadMcReweight;
    if (hRecoEtaLeadEtaSubleadMcReweight) delete hRecoEtaLeadEtaSubleadMcReweight;
    for (int iRun{0}; iRun<6; ++iRun) {
        if (hRecoDijetEtaCMRun[iRun]) delete hRecoDijetEtaCMRun[iRun];
    }
    if (hRecoDijetPtEtaForward) delete hRecoDijetPtEtaForward;
    if (hRecoDijetPtEtaBackward) delete hRecoDijetPtEtaBackward;
    if (hRecoDijetPtEtaForwardCMInLab) delete hRecoDijetPtEtaForwardCMInLab;
    if (hRecoDijetPtEtaBackwardCMInLab) delete hRecoDijetPtEtaBackwardCMInLab;
    if (hRecoDijetPtEtaCMForward) delete hRecoDijetPtEtaCMForward;
    if (hRecoDijetPtEtaCMBackward) delete hRecoDijetPtEtaCMBackward;
    if (hRecoDijetPtEtaForwardLabInCM) delete hRecoDijetPtEtaForwardLabInCM;
    if (hRecoDijetPtEtaBackwardLabInCM) delete hRecoDijetPtEtaBackwardLabInCM;
    if (hRecoDijetPtEtaForwardWeighted) delete hRecoDijetPtEtaForwardWeighted;
    if (hRecoDijetPtEtaBackwardWeighted) delete hRecoDijetPtEtaBackwardWeighted;
    if (hRecoDijetPtEtaCMForwardWeighted) delete hRecoDijetPtEtaCMForwardWeighted;
    if (hRecoDijetPtEtaCMBackwardWeighted) delete hRecoDijetPtEtaCMBackwardWeighted;
    for (int i{0}; i<6; i++) {
        if (hRecoDijetPtRawEtaForwardArr[i]) delete hRecoDijetPtRawEtaForwardArr[i];
        if (hRecoDijetPtRawEtaBackwardArr[i]) delete hRecoDijetPtRawEtaBackwardArr[i];
        if (hRecoDijetPtEtaForwardArr[i]) delete hRecoDijetPtEtaForwardArr[i];
        if (hRecoDijetPtEtaBackwardArr[i]) delete hRecoDijetPtEtaBackwardArr[i];
    }
    if (hRecoDijetPtEta) delete hRecoDijetPtEta;
    if (hRecoDijetPtEtaWeighted) delete hRecoDijetPtEtaWeighted;
    if (hRecoDijetPtEtaCMInLab) delete hRecoDijetPtEtaCMInLab;
    if (hRecoDijetPtEtaCM) delete hRecoDijetPtEtaCM;
    if (hRecoDijetPtEtaCMWeighted) delete hRecoDijetPtEtaCMWeighted;
    if (hRecoDijetPtEtaLabInCM) delete hRecoDijetPtEtaLabInCM;
    if (hRecoDijetPtEtaMatched) delete hRecoDijetPtEtaMatched;
    if (hRecoDijetPtEtaCMMatched) delete hRecoDijetPtEtaCMMatched;
    if (hRecoDijetPtAveLeadPtSubLeadPt) delete hRecoDijetPtAveLeadPtSubLeadPt;
    if (hRecoDijetPtAveLeadPtSubLeadPtCM) delete hRecoDijetPtAveLeadPtSubLeadPtCM;
    if (hRecoDijetPtAveLeadEtaSubLeadEta) delete hRecoDijetPtAveLeadEtaSubLeadEta;
    if (hRecoDijetPtAveLeadEtaSubLeadEtaCM) delete hRecoDijetPtAveLeadEtaSubLeadEtaCM;
    if (hRecoDijetEtaLeadEtaSubLeadEta) delete hRecoDijetEtaLeadEtaSubLeadEta;
    if (hRecoDijetEtaLeadEtaSubLeadEtaCM) delete hRecoDijetEtaLeadEtaSubLeadEtaCM;

    if (hRecoDijetLeadPtEta) delete hRecoDijetLeadPtEta;
    if (hRecoDijetLeadPtEtaStdBins) delete hRecoDijetLeadPtEtaStdBins;
    if (hRecoDijetSubLeadPtEta) delete hRecoDijetSubLeadPtEta;
    if (hRecoDijetSubLeadPtEtaStdBins) delete hRecoDijetSubLeadPtEtaStdBins;
    for (int i{0}; i<3; i++) {
        if (hRecoDijetXj[i]) delete hRecoDijetXj[i];
        if (hRecoDijetXjCM[i]) delete hRecoDijetXjCM[i];
    }

    if (hRecoLeadAllJetPtEta) delete hRecoLeadAllJetPtEta;
    if (hRecoLeadAllJetPtEtaCM) delete hRecoLeadAllJetPtEtaCM;
    if (hRecoLeadAllJetPtRawEtaStdBins) delete hRecoLeadAllJetPtRawEtaStdBins;
    if (hRecoLeadAllJetPtEtaStdBins) delete hRecoLeadAllJetPtEtaStdBins;
    if (hRecoLeadAllJetPtEtaPtHat) delete hRecoLeadAllJetPtEtaPtHat;
    if (hRecoSubLeadAllJetPtEta) delete hRecoSubLeadAllJetPtEta;
    if (hRecoSubLeadAllJetPtEtaCM) delete hRecoSubLeadAllJetPtEtaCM;
    if (hRecoSubLeadAllJetPtRawEtaStdBins) delete hRecoSubLeadAllJetPtRawEtaStdBins;
    if (hRecoSubLeadAllJetPtEtaStdBins) delete hRecoSubLeadAllJetPtEtaStdBins;
    if (hRecoSubLeadAllJetPtEtaPtHat) delete hRecoSubLeadAllJetPtEtaPtHat;
    if (hRecoGoodInclusiveJetEtaLabFrame) delete hRecoGoodInclusiveJetEtaLabFrame;
    if (hRecoGoodInclusiveJetEtaCMFrame) delete hRecoGoodInclusiveJetEtaCMFrame;
    if (hRecoDijetEta) delete hRecoDijetEta;
    if (hRecoDijetEtaCM) delete hRecoDijetEtaCM;

    if ( fIsMc ) {

        //
        // Ref jet histograms
        //
        if (hRecoInclusiveJetReco2Ref) delete hRecoInclusiveJetReco2Ref;
        if (hRecoLeadJetReco2Ref) delete hRecoLeadJetReco2Ref;
        if (hRecoSubLeadJetReco2Ref) delete hRecoSubLeadJetReco2Ref;
        if (hJESInclusiveJetPtEtaPhi) delete hJESInclusiveJetPtEtaPhi;

        if (hRecoLeadJetPtOverPtHatVsLeadJetPt) delete hRecoLeadJetPtOverPtHatVsLeadJetPt;
        if (hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted) delete hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted;
        if (hRecoDijetPtOverPtHatVsDijetPt) delete hRecoDijetPtOverPtHatVsDijetPt;
        if (hRecoDijetPtOverPtHatVsDijetPtWeighted) delete hRecoDijetPtOverPtHatVsDijetPtWeighted;
        if (hRecoDijetPtAveOverPtHatVsDijetPtAve) delete hRecoDijetPtAveOverPtHatVsDijetPtAve;
        if (hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted) delete hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted;

        // if (hRecoInclusiveJetJECFactorVsPtEta) delete hRecoInclusiveJetJECFactorVsPtEta;
        // if (hRecoInclusiveJetJEC2FactorVsPtEta) delete hRecoInclusiveJetJEC2FactorVsPtEta;
        // if (hRecoInclusiveJetPtRawOverPtRefVsPtEta) delete hRecoInclusiveJetPtRawOverPtRefVsPtEta;
        // if (hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning) delete hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning;
        // if (hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning) delete hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning;;

        if ( hInclusiveJetJESVsPtGen ) delete hInclusiveJetJESVsPtGen;
        if ( hInclusiveJetJESGenPtGenEtaPtHatWeighted ) delete hInclusiveJetJESGenPtGenEtaPtHatWeighted;
        if ( hInclusiveJetJESGenPtGenEtaCMPtHatWeighted ) delete hInclusiveJetJESGenPtGenEtaCMPtHatWeighted;
        if ( hInclusiveJetJESRecoPtRecoEtaPtHatWeighted ) delete hInclusiveJetJESRecoPtRecoEtaPtHatWeighted;
        if ( hLeadJetJESGenPtEtaPtHatWeighted ) delete hLeadJetJESGenPtEtaPtHatWeighted;
        if ( hLeadJetJESGenPtEtaCMPtHatWeighted ) delete hLeadJetJESGenPtEtaCMPtHatWeighted;
        if ( hSubLeadJetJESGenPtEtaPtHatWeighted ) delete hSubLeadJetJESGenPtEtaPtHatWeighted;
        if ( hSubLeadJetJESGenPtEtaCMPtHatWeighted ) delete hSubLeadJetJESGenPtEtaCMPtHatWeighted;

        if (hRecoInclusiveMatchedJetPt) delete hRecoInclusiveMatchedJetPt;
        if (hRecoInclusiveMatchedJetPtEta) delete hRecoInclusiveMatchedJetPtEta;
        if (hRecoInclusiveUnmatchedJetPtEta) delete hRecoInclusiveUnmatchedJetPtEta;
        if (hRecoLeadMatchedJetPtEta) delete hRecoLeadMatchedJetPtEta;
        if (hRecoLeadMatchedJetPtEtaPtHat) delete hRecoLeadMatchedJetPtEtaPtHat;
        if (hRecoLeadUnmatchedJetPtEta) delete hRecoLeadUnmatchedJetPtEta;
        if (hRecoLeadUnmatchedJetPtEtaPtHat) delete hRecoLeadUnmatchedJetPtEtaPtHat;
        if (hRecoSubLeadMatchedJetPtEta) delete hRecoSubLeadMatchedJetPtEta;
        if (hRecoSubLeadMatchedJetPtEtaPtHat) delete hRecoSubLeadMatchedJetPtEtaPtHat;
        if (hRecoSubLeadUnmatchedJetPtEta) delete hRecoSubLeadUnmatchedJetPtEta;
        if (hRecoSubLeadUnmatchedJetPtEtaPtHat) delete hRecoSubLeadUnmatchedJetPtEtaPtHat;

        if (hRefInclusiveJetPt) delete hRefInclusiveJetPt;
        if (hRefInclusiveJetEta) delete hRefInclusiveJetEta;
        if (hRefInclusiveJetEtaUnweighted) delete hRefInclusiveJetEtaUnweighted;
        if (hRefInclusiveJetPtEta) delete hRefInclusiveJetPtEta;
        if (hRefInclusiveJetPtEtaPtHat) delete hRefInclusiveJetPtEtaPtHat;
        if (hRefLeadJetPtEta) delete hRefLeadJetPtEta;
        if (hRefLeadJetPtEtaPtHat) delete hRefLeadJetPtEtaPtHat;
        if (hRefLeadUnswappedJetPtEta) delete hRefLeadUnswappedJetPtEta;
        if (hRefLeadUnswappedJetPtEtaPtHat) delete hRefLeadUnswappedJetPtEtaPtHat;
        if (hRefSubLeadJetPtEta) delete hRefSubLeadJetPtEta;
        if (hRefSubLeadJetPtEtaPtHat) delete hRefSubLeadJetPtEtaPtHat;
        if (hRefSubLeadUnswappedJetPtEta) delete hRefSubLeadUnswappedJetPtEta;
        if (hRefSubLeadUnswappedJetPtEtaPtHat) delete hRefSubLeadUnswappedJetPtEtaPtHat;

        // if (hReco2RefFull) delete hReco2RefFull;
        if (hRecoDijetPtEtaRefDijetPtEta) delete hRecoDijetPtEtaRefDijetPtEta;
        if (hRecoDijetPtEtaRefDijetPtEtaWeighted) delete hRecoDijetPtEtaRefDijetPtEtaWeighted;
        // if (hRefSel2RecoFull) delete hRefSel2RecoFull;

        if (hRefDijetEta) delete hRefDijetEta;
        if (hRefDijetEtaVsRecoDijetEta) delete hRefDijetEtaVsRecoDijetEta;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted;
        if (hRefDijetPtEta) delete hRefDijetPtEta;
        if (hRefDijetPtEtaWeighted) delete hRefDijetPtEtaWeighted;
        if (hRefDijetPtEtaCMInLab) delete hRefDijetPtEtaCMInLab;

        if (hRefDijetPtEtaForward) delete hRefDijetPtEtaForward;
        if (hRefDijetPtEtaBackward) delete hRefDijetPtEtaBackward;
        if (hRefDijetPtEtaForwardCMInLab) delete hRefDijetPtEtaForwardCMInLab;
        if (hRefDijetPtEtaBackwardCMInLab) delete hRefDijetPtEtaBackwardCMInLab;
        if (hRefDijetPtEtaCMForward) delete hRefDijetPtEtaCMForward;
        if (hRefDijetPtEtaCMBackward) delete hRefDijetPtEtaCMBackward;
        if (hRefDijetPtEtaForwardLabInCM) delete hRefDijetPtEtaForwardLabInCM;
        if (hRefDijetPtEtaBackwardLabInCM) delete hRefDijetPtEtaBackwardLabInCM;
        if (hRefDijetPtEtaForwardWeighted) delete hRefDijetPtEtaForwardWeighted;
        if (hRefDijetPtEtaBackwardWeighted) delete hRefDijetPtEtaBackwardWeighted;
        if (hRefDijetPtEtaCMForwardWeighted) delete hRefDijetPtEtaCMForwardWeighted;
        if (hRefDijetPtEtaCMBackwardWeighted) delete hRefDijetPtEtaCMBackwardWeighted;

        if (hRefDijetEtaCM) delete hRefDijetEtaCM;
        if (hRefDijetPtEtaCM) delete hRefDijetPtEtaCM;
        if (hRefDijetPtEtaCMWeighted) delete hRefDijetPtEtaCMWeighted;
        if (hRefDijetPtEtaLabInCM) delete hRefDijetPtEtaLabInCM;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted;

        if (hRefPtLeadPtSublead) delete hRefPtLeadPtSublead;
        if (hRefEtaLeadEtaSublead) delete hRefEtaLeadEtaSublead;
        if (hRefEtaCMLeadEtaCMSublead) delete hRefEtaCMLeadEtaCMSublead;
        if (hRefPtLeadPtSubleadMcReweight) delete hRefPtLeadPtSubleadMcReweight;
        if (hRefEtaLeadEtaSubleadMcReweight) delete hRefEtaLeadEtaSubleadMcReweight;

        // Ref-selected jet histograms

        if (hRefSelInclusiveJetPt) delete hRefSelInclusiveJetPt;
        if (hRefSelInclusiveJetEta) delete hRefSelInclusiveJetEta;
        if (hRefSelInclusiveJetEtaUnweighted) delete hRefSelInclusiveJetEtaUnweighted;
        if (hRefSelInclusiveJetPtEta) delete hRefSelInclusiveJetPtEta;
        if (hRefSelInclusiveJetPtEtaPtHat) delete hRefSelInclusiveJetPtEtaPtHat;
        if (hRefSelLeadJetPtEtaPtHat) delete hRefSelLeadJetPtEtaPtHat;
        if (hRefSelSubLeadJetPtEtaPtHat) delete hRefSelSubLeadJetPtEtaPtHat;

        if (hRefSelDijetEta) delete hRefSelDijetEta;
        if (hRefSelDijetPtEta) delete hRefSelDijetPtEta;
        if (hRefSelDijetPtEtaWeighted) delete hRefSelDijetPtEtaWeighted;
        if (hRefDijetPtEtaCMInLab) delete hRefDijetPtEtaCMInLab;
        if (hRefSelDijetEtaCM) delete hRefSelDijetEtaCM;
        if (hRefSelDijetPtEtaCM) delete hRefSelDijetPtEtaCM;
        if (hRefSelDijetPtEtaCMWeighted) delete hRefSelDijetPtEtaCMWeighted;
        if (hRefSelDijetPtEtaLabInCM) delete hRefSelDijetPtEtaLabInCM;

    } // if (fIsMc)
}

//________________
void HistoManagerDiJet::init() {

    double dijetEtaVals[] { -3.0, -2.4, -2.2, 
                            -2.0, -1.8, -1.6, -1.4, -1.2, 
                            -1.0, -0.8, -0.6, -0.4, -0.2,  
                            0.0,  0.2,  0.4,  0.6,  0.8,  
                            1.0,  1.2,  1.4,  1.6,  1.8,  
                            2.0,  2.2,  2.4,  3.0 };
    int dijetEtaBins = sizeof(dijetEtaVals)/sizeof(double)-1;

    // Old binning convention
    // double dijetEtaOldVals[] = { -2.915, -2.63333333333, -2.07, -1.78833333333, -1.50666666667,
    //                             -1.225, -0.94333333333, -0.66166666666, -0.38, -0.09833333333,
    //                             0.18333333333, 0.465, 0.74666666666, 1.02833333333, 1.31,
    //                             1.59166666667, 1.87333333333, 2.43666666667, 3.};
    // int dijetEtaOldBins = sizeof(dijetEtaOldVals)/sizeof(double) - 1;

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

    double dijetEtaFBVals[] { 0.0,  0.2,  0.4,  0.6,  0.8,  
                              1.0,  1.2,  1.4,  1.6,  1.8,  
                              2.0,  2.2,  2.4,  3.0 };
    int dijetEtaFBBins = sizeof(dijetEtaFBVals)/sizeof(double) - 1;
    
    // int    prescale = 1;

    int xBins = 10000;
    double xRange[2] = {0.0001, 1.};
    int xPbOverXpBins = 10000;
    double xPbOverXpRange[2] = {0.0001, 1000.};

    int    vzBins = 320;
    double vzRange[2] {-31., 31.};
    // int    multBins{1800};
    // double multRange[2] {-0.5, 1799.5};
    int    hiBinBins{203};
    double hiBinRange[2] {-1.5, 201.5};
    // int    centralityBins{101};
    // double centralityRange[2] {-0.5, 100.5};
    // int    weightBins{110};
    // double weightRange[2] {-0.05, 1.05};
    int    ptHatBins{100};
    double ptHatRange[2] {0., 1000.};
    int    fJESBins{100}; 
    double fJESRange[2] {0., 2.};


    // int    bins2D_ev_VzPtHat[2] {     vzBins,     fPtHatBins };
    // double xmin2D_ev_VzPtHat[2] { vzRange[0], fPtHatRange[0] };
    // double xmax2D_ev_VzPtHat[2] { vzRange[1], fPtHatRange[1] };

    // int    dEtaBins{50}; 
    // double dEtaRange[2] {-0.05, 0.05};

    int dijetDphiBins{16}; 
    double dijetDphiRange[2] {-TMath::Pi(), TMath::Pi()};

    //
    // Gen
    //
    int    bins10D_gen_dijetInfo[10]
    { fDijetPtBins, fDijetEtaBins, fDijetEtaBins, dijetDphiBins, fPtBins, fEtaBins, fEtaBins, fPtBins, fEtaBins, fEtaBins };
    double xmin10D_gen_dijetInfo[10]
    { fDijetPtRange[0], fDijetEtaRange[0], fDijetEtaRange[0], dijetDphiRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0]  };
    double xmax10D_gen_dijetInfo[10]
    { fDijetPtRange[1], fDijetEtaRange[1], fDijetEtaRange[1], dijetDphiRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1] };

    //
    // Reco single
    //
    int    bins5D_jet_PtPtPtEtaEta[5] { fPtBins, fPtBins, fPtBins, fEtaBins, fEtaBins };
    double xmin5D_jet_PtPtPtEtaEta[5] { fPtRange[0], fPtRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0] };
    double xmax5D_jet_PtPtPtEtaEta[5] { fPtRange[1], fPtRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1] };

    int    bins4D_jet_JESPtEtaPhi[4] { fJESBins,     fPtBins,     fEtaBins,     fPhiBins };
    double xmin4D_jet_JESPtEtaPhi[4] { fJESRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0] };
    double xmax4D_jet_JESPtEtaPhi[4] { fJESRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    int    bins4D_jet_JESPtEtaPtHat[4] { fJESBins, fPtBins, fEtaBins, fPtHatBins };
    double xmin4D_jet_JESPtEtaPtHat[4] { fJESRange[0], fPtRange[0], fEtaRange[0], fPtHatRange[0] };
    double xmax4D_jet_JESPtEtaPtHat[4] { fJESRange[1], fPtRange[1], fEtaRange[1], fPtHatRange[1] };

    //
    // Reco dijets
    //
    int    bins10D_dijet_info[10]
    { fDijetPtBins, fDijetEtaBins, fDijetEtaBins, dijetDphiBins, fPtBins, fEtaBins, fEtaBins, fPtBins, fEtaBins, fEtaBins };
    double xmin10D_dijet_info[10]
    { fDijetPtRange[0], fDijetEtaRange[0], fDijetEtaRange[0], dijetDphiRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0], fEtaRange[0] };
    double xmax10D_dijet_info[10]
    { fDijetPtRange[1], fDijetEtaRange[1], fDijetEtaRange[1], dijetDphiRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1], fEtaRange[1] };

    int    bins12D_ref2reco_info[12]
    { fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins, fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins };
    double xmin12D_ref2reco_info[12]
    { fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0] };
    double xmax12D_ref2reco_info[12]
    { fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1] };

    int    bins4D_ref2reco_info[4] { fDijetPtBins, fDijetEtaBins, fDijetPtBins, fDijetEtaBins };
    double xmin4D_ref2reco_info[4] { fDijetPtRange[0], fDijetEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0] };
    double xmax4D_ref2reco_info[4] { fDijetPtRange[1], fDijetEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1] };


    //
    // Event histograms
    //
    
    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVz->Sumw2();
    hVzWeighted = new TH1D("hVzWeighted","Vertex z position;vz (cm);Entries", 400, -50., 50.);
    hVzWeighted->Sumw2();
    hPtHat = new TH1D("hPtHat","#hat{p_{T}};#hat{p_{T}} (GeV);Entries", 
                      ptHatBins, ptHatRange[0], ptHatRange[1]);
    hPtHat->Sumw2();
    hPtHatWeighted = new TH1D("hPtHatWeighted","#hat{p_{T}} with #hat{p_{T}} weight;#hat{p_{T}} (GeV);Entries", 
                              ptHatBins, ptHatRange[0], ptHatRange[1]);
    hPtHatWeighted->Sumw2();
    hHiBin = new TH1D("hHiBin","HiBin a.k.a. centrality;HiBin;Entries", 
                      hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBin->Sumw2();
    hHiBinWeighted = new TH1D("hHiBinWeighted","HiBin a.k.a. centrality with #hat{p_{T}} weight;HiBin;Entries", 
                              hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinWeighted->Sumw2();

    hVzRecoDijetLab = new TH1D("hVzRecoDijetLab","Reco Dijet Lab Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzRecoDijetLab->Sumw2();
    hVzRecoDijetLabWeighted = new TH1D("hVzRecoDijetLabWeighted","Reco Dijet Lab Vertex z position weighted;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzRecoDijetLabWeighted->Sumw2();
    hHiBinRecoDijetLab = new TH1D("hHiBinRecoDijetLab","Reco Dijet Lab HiBin;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinRecoDijetLab->Sumw2();
    hHiBinRecoDijetLabWeighted = new TH1D("hHiBinRecoDijetLabWeighted","Reco Dijet Lab HiBin weighted;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinRecoDijetLabWeighted->Sumw2();

    hVzRecoDijetCM = new TH1D("hVzRecoDijetCM","Reco Dijet CM Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzRecoDijetCM->Sumw2();
    hVzRecoDijetCMWeighted = new TH1D("hVzRecoDijetCMWeighted","Reco Dijet CM Vertex z position weighted;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzRecoDijetCMWeighted->Sumw2();
    hHiBinRecoDijetCM = new TH1D("hHiBinRecoDijetCM","Reco Dijet CM HiBin;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinRecoDijetCM->Sumw2();
    hHiBinRecoDijetCMWeighted = new TH1D("hHiBinRecoDijetCMWeighted","Reco Dijet CM HiBin weighted;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinRecoDijetCMWeighted->Sumw2();

    if ( fIsMc ) {
        hVzGenDijetLab = new TH1D("hVzGenDijetLab","Gen Dijet Lab Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzGenDijetLab->Sumw2();
        hVzGenDijetLabWeighted = new TH1D("hVzGenDijetLabWeighted","Gen Dijet Lab Vertex z position weighted;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzGenDijetLabWeighted->Sumw2();
        hHiBinGenDijetLab = new TH1D("hHiBinGenDijetLab","Gen Dijet Lab HiBin;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinGenDijetLab->Sumw2();
        hHiBinGenDijetLabWeighted = new TH1D("hHiBinGenDijetLabWeighted","Gen Dijet Lab HiBin weighted;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinGenDijetLabWeighted->Sumw2();

        hVzGenDijetCM = new TH1D("hVzGenDijetCM","Gen Dijet CM Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzGenDijetCM->Sumw2();
        hVzGenDijetCMWeighted = new TH1D("hVzGenDijetCMWeighted","Gen Dijet CM Vertex z position weighted;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzGenDijetCMWeighted->Sumw2();
        hHiBinGenDijetCM = new TH1D("hHiBinGenDijetCM","Gen Dijet CM HiBin;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinGenDijetCM->Sumw2();
        hHiBinGenDijetCMWeighted = new TH1D("hHiBinGenDijetCMWeighted","Gen Dijet CM HiBin weighted;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinGenDijetCMWeighted->Sumw2();

        hVzRefSelDijetLab = new TH1D("hVzRefSelDijetLab","RefSel Dijet Lab Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzRefSelDijetLab->Sumw2();
        hVzRefSelDijetLabWeighted = new TH1D("hVzRefSelDijetLabWeighted","RefSel Dijet Lab Vertex z position weighted;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzRefSelDijetLabWeighted->Sumw2();
        hHiBinRefSelDijetLab = new TH1D("hHiBinRefSelDijetLab","RefSel Dijet Lab HiBin;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinRefSelDijetLab->Sumw2();
        hHiBinRefSelDijetLabWeighted = new TH1D("hHiBinRefSelDijetLabWeighted","RefSel Dijet Lab HiBin weighted;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinRefSelDijetLabWeighted->Sumw2();

        hVzRefSelDijetCM = new TH1D("hVzRefSelDijetCM","RefSel Dijet CM Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzRefSelDijetCM->Sumw2();
        hVzRefSelDijetCMWeighted = new TH1D("hVzRefSelDijetCMWeighted","RefSel Dijet CM Vertex z position weighted;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
        hVzRefSelDijetCMWeighted->Sumw2();
        hHiBinRefSelDijetCM = new TH1D("hHiBinRefSelDijetCM","RefSel Dijet CM HiBin;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinRefSelDijetCM->Sumw2();
        hHiBinRefSelDijetCMWeighted = new TH1D("hHiBinRefSelDijetCMWeighted","RefSel Dijet CM HiBin weighted;HiBin;Entries", hiBinBins, hiBinRange[0], hiBinRange[1]);
        hHiBinRefSelDijetCMWeighted->Sumw2();
    } // if ( fIsMc )

    //
    // Reco histograms
    // 

    // Inclusive jet histograms

    // For 4 eta ranges: <=2.4, <=2.7, <=3, >3
    for (int i{0}; i<4; i++) {
        float low{0}, hi{0};
        if      ( i == 0 ) { low = {-2.4f}; hi = {2.4f}; }
        else if ( i == 1 ) { low = {2.4f}; hi = {2.7f}; }
        else if ( i == 2 ) { low = {2.7f}; hi = {3.0f}; }
        else if ( i == 3 ) { low = {3.f}; hi = {100.0f}; }
        hRecoInclusiveJetNHF[i] = new TH1D(Form("hRecoInclusiveJetNHF_%d",i), Form("Neutral hadron fraction for %3.1f< #eta<=%3.1f;Neutral hadron fraction;1/N dN/dNHF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetNHF[i]->Sumw2();
        hRecoInclusiveJetNEmF[i] = new TH1D(Form("hRecoInclusiveJetNEmF_%d",i), Form("Neutral EM fraction for %3.1f< #eta<=%3.1f;Neutral EM fraction;1/N dN/dNEF", low, hi), 
                            fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetNEmF[i]->Sumw2();
        hRecoInclusiveJetNumOfConst[i] = new TH1D(Form("hRecoInclusiveJetNumOfConst_%d",i), Form("Number of constituents for %3.1f< #eta<=%3.1f;Number of constituents;1/N dN/dNconst", low, hi), 
                                  fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetNumOfConst[i]->Sumw2();
        hRecoInclusiveJetMUF[i] = new TH1D(Form("hRecoInclusiveJetMUF_%d",i), Form("Muon fraction for %3.1f< #eta<=%3.1f;Muon fraction;1/N dN/dMUF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetMUF[i]->Sumw2();
        hRecoInclusiveJetCHF[i] = new TH1D(Form("hRecoInclusiveJetCHF_%d",i), Form("Charged hadron fraction for %3.1f< #eta<=%3.1f;Charged hadron fraction;1/N dN/dCHF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetCHF[i]->Sumw2();
        hRecoInclusiveJetChargedMult[i] = new TH1D(Form("hRecoInclusiveJetChargedMult_%d",i), Form("Charged particle multiplicity for %3.1f< #eta<=%3.1f;Charged multiplicity;1/N dN/dChMult", low, hi), 
                                   fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetChargedMult[i]->Sumw2();
        hRecoInclusiveJetCEmF[i] = new TH1D(Form("hRecoInclusiveJetCEmF_%d",i), Form("Charged EM fraction for %3.1f< #eta<=%3.1f;Charged EM fraction;1/N dN/dChEmF", low, hi), 
                            fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetCEmF[i]->Sumw2();
        hRecoInclusiveJetNumOfNeutPart[i] = new TH1D(Form("hRecoInclusiveJetNumOfNeutPart_%d",i), Form("Number of neutral particles for %3.1f< #eta<=%3.1f;Number of neutral particles;1/N dN/dNumOfNeutrals", low, hi), 
                                     fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetNumOfNeutPart[i]->Sumw2();
    } // for (int i{0}; i<4; i++)

    hRecoJetCollectionSize = new TH1D("hRecoJetCollectionSize","Reco jet collection size;Number of jets;Entries", 100, -0.5, 99.5);
    hRecoJetCollectionSize->Sumw2();

    hRecoLeadAllJetPtEta = new TH2D("hRecoLeadAllJetPtEta","Lead jet all p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoLeadAllJetPtEta->Sumw2();
    hRecoLeadAllJetPtEtaCM = new TH2D("hRecoLeadAllJetPtEtaCM","Lead jet all p_{T} vs #eta CM frame;#eta;p_{T} (GeV)", 
                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoLeadAllJetPtEtaCM->Sumw2();
    hRecoLeadAllJetPtRawEtaStdBins = new TH2D("hRecoLeadAllJetPtRawEtaStdBins","Lead jet p_{T} (raw) vs #eta (std bins);#eta;p_{T}^{raw} (GeV)",
                                              208, -5.2, 5.2,
                                              fPtBins, fPtRange[0], fPtRange[1]);
    hRecoLeadAllJetPtRawEtaStdBins->Sumw2();
    hRecoLeadAllJetPtRawEtaStdBins->GetYaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoLeadAllJetPtRawEtaStdBins->SetBinsLength(-1);
    hRecoLeadAllJetPtEtaStdBins = new TH2D("hRecoLeadAllJetPtEtaStdBins","Lead jet p_{T} vs #eta (std bins);#eta;p_{T} (GeV)",
                                              208, -5.2, 5.2,
                                              fPtBins, fPtRange[0], fPtRange[1]);
    hRecoLeadAllJetPtEtaStdBins->Sumw2();
    hRecoLeadAllJetPtEtaStdBins->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoLeadAllJetPtEtaStdBins->SetBinsLength(-1);

    hRecoLeadAllJetPtEtaPtHat = new TH3D("hRecoLeadAllJetPtEtaPtHat","Lead jet (matched+unmatched) p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                         fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoLeadAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoLeadAllJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoLeadAllJetPtEtaPtHat->Sumw2();

    hRecoInclusiveAllJetPtRawEta = new TH2D("hRecoInclusiveAllJetPtRawEta","Reco jet raw p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtRawEta->Sumw2();

    hRecoSubLeadAllJetPtEta = new TH2D("hRecoSubLeadAllJetPtEta","SubLead jet all p_{T} vs #eta;#eta;p_{T} (GeV)",
                                         fEtaBins, fEtaRange[0], fEtaRange[1], 
                                            fPtBins, fPtRange[0], fPtRange[1]);
    hRecoSubLeadAllJetPtEta->Sumw2();
    hRecoSubLeadAllJetPtEtaCM = new TH2D("hRecoSubLeadAllJetPtEtaCM","SubLead jet all p_{T} vs #eta CM frame;#eta;p_{T} (GeV)",
                                         fEtaBins, fEtaRange[0], fEtaRange[1], 
                                         fPtBins, fPtRange[0], fPtRange[1]);
    hRecoSubLeadAllJetPtEtaCM->Sumw2();
    hRecoSubLeadAllJetPtRawEtaStdBins = new TH2D("hRecoSubLeadAllJetPtRawEtaStdBins","SubLead jet p_{T} (raw) vs #eta (std bins);#eta;p_{T}^{raw} (GeV)",
                                              208, -5.2, 5.2,
                                              fPtBins, fPtRange[0], fPtRange[1]);
    hRecoSubLeadAllJetPtRawEtaStdBins->Sumw2();
    hRecoSubLeadAllJetPtRawEtaStdBins->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoSubLeadAllJetPtRawEtaStdBins->SetBinsLength(-1);
    hRecoSubLeadAllJetPtEtaStdBins = new TH2D("hRecoSubLeadAllJetPtEtaStdBins","SubLead jet p_{T} vs #eta (std bins);#eta;p_{T} (GeV)",
                                              208, -5.2, 5.2,
                                              fPtBins, fPtRange[0], fPtRange[1]);
    hRecoSubLeadAllJetPtEtaStdBins->Sumw2();
    hRecoSubLeadAllJetPtEtaStdBins->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoSubLeadAllJetPtEtaStdBins->SetBinsLength(-1);

    hRecoSubLeadAllJetPtEtaPtHat = new TH3D("hRecoSubLeadAllJetPtEtaPtHat","SubLead jet (matched+unmatched) p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoSubLeadAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoSubLeadAllJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoSubLeadAllJetPtEtaPtHat->Sumw2();


    hRecoInclusiveAllJetPt = new TH1D("hRecoInclusiveAllJetPt","Reco jet p_{T} (matched+unmatched);Reco p_{T} (GeV);Entries",
                                   fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPt->Sumw2();
    hRecoInclusiveAllJetEta = new TH1D("hRecoInclusiveAllJetEta","Reco jet #eta (matched+unmatched);#eta;Entries",
                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoInclusiveAllJetEta->Sumw2();
    hRecoInclusiveAllJetEtaUnweighted = new TH1D("hRecoInclusiveAllJetEtaUnweighted","Reco jet #eta (matched+unmatched) unweighted;#eta;Entries",
                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoInclusiveAllJetEtaUnweighted->Sumw2();
    hRecoInclusiveAllJetPtEta = new TH2D("hRecoInclusiveAllJetPtEta","Reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                      fEtaBins, fEtaRange[0], fEtaRange[1],
                                      fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtEta->Sumw2();
    hRecoInclusiveAllJetPtEtaCM = new TH2D("hRecoInclusiveAllJetPtEtaCM","Reco jet p_{T} vs #eta CM frame;#eta;p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtEtaCM->Sumw2();

    hRecoInclusiveAllJetPtRawEtaStdBins = new TH2D("hRecoInclusiveAllJetPtRawEtaStdBins","Reco jet p_{T} (raw) vs #eta (std bins);#eta;p_{T}^{raw} (GeV)",
                                              208, -5.2, 5.2,
                                              fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtRawEtaStdBins->Sumw2();
    hRecoInclusiveAllJetPtRawEtaStdBins->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoInclusiveAllJetPtRawEtaStdBins->SetBinsLength(-1);
    hRecoInclusiveAllJetPtEtaStdBins = new TH2D("hRecoInclusiveAllJetPtEtaStdBins","Reco jet p_{T} vs #eta (std bins);#eta;p_{T} (GeV)",
                                              208, -5.2, 5.2,
                                              fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtEtaStdBins->Sumw2();
    hRecoInclusiveAllJetPtEtaStdBins->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoInclusiveAllJetPtEtaStdBins->SetBinsLength(-1);

    hRecoInclusiveAllJetPtEtaPtHat = new TH3D("hRecoInclusiveAllJetPtEtaPtHat","Reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                              fEtaBins, fEtaRange[0], fEtaRange[1],
                                           fPtBins, fPtRange[0], fPtRange[1],
                                           fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveAllJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveAllJetPtEtaPtHat->Sumw2();
    hRecoInclusiveAllJetPtRawEtaPtHat = new TH3D("hRecoInclusiveAllJetPtRawEtaPtHat","Reco jet p_{T} (raw) vs #eta vs #hat{p_{T}};#eta;p_{T}^{raw} (GeV);#hat{p_{T}} (GeV)",
                                         fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveAllJetPtRawEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveAllJetPtRawEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveAllJetPtRawEtaPtHat->Sumw2();

    hRecoInclusiveMatchedJetPtEtaPtHat = new TH3D("hRecoInclusiveMatchedJetPtEtaPtHat","Matched reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                         fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveMatchedJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveMatchedJetPtEtaPtHat->Sumw2();
    hRecoInclusiveUnmatchedJetPtEtaPtHat = new TH3D("hRecoInclusiveUnmatchedJetPtEtaPtHat","Unmatched reco jet p_{T} vs #eta vs #hat{p}_{T}};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                         fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveUnmatchedJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveUnmatchedJetPtEtaPtHat->Sumw2();

    // Eta distributions of good jets for specific runIds
    for (int i{0}; i<6; ++i) {
        int runId{0};
        if (i == 1) { runId = 285480; }
        else if (i == 2) { runId = 285505; }
        else if (i == 3) { runId = 285517; }
        else if (i == 4) { runId = 285832; }
        else if (i == 5) { runId = 285993; }
        else { runId = 0; }

        hRecoInclusiveJetEtaRun[i] = new TH1D(Form("hRecoInclusiveJetEtaRun_%d",i), Form("Reco good jet #eta (40<p_{T}<90 GeV) for runId %d;#eta^{Inclusive};Entries", runId),
                                              fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetEtaRun[i]->Sumw2();
        hRecoLeadJetEtaRun[i] = new TH1D(Form("hRecoLeadJetEtaRun_%d",i), Form("Reco lead jet #eta (40<p_{T}<90 GeV) for runId %d;#eta^{Lead};Entries", runId),
                                              fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoLeadJetEtaRun[i]->Sumw2();
        hRecoSubLeadJetEtaRun[i] = new TH1D(Form("hRecoSubLeadJetEtaRun_%d",i), Form("Reco sublead jet #eta (40<p_{T}<90 GeV) for runId %d;#eta^{SubLead};Entries", runId),
                                              fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoSubLeadJetEtaRun[i]->Sumw2();
    }

    hRecoGoodInclusiveJetEtaLabFrame = new TH1D("hRecoGoodInclusiveJetEtaLabFrame","Reco good jet #eta in lab frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaLabFrame->Sumw2();
    hRecoGoodInclusiveJetEtaCMFrame = new TH1D("hRecoGoodInclusiveJetEtaCMFrame","Reco good jet #eta in CM frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaCMFrame->Sumw2();

    // Dijet histograms

    // Reco dijet info [10 dimensions]
    // 0 - dijet pt ave, 1 - dijet eta lab, 2 - dijet eta cm, 3 - dijet delta phi
    // 4 - lead pt, 5 - lead eta, 6 - lead eta cm,
    // 7 - sublead pt, 8 - sublead eta, 9 - sublead eta cm 
    // hRecoDijetInfo = new THnSparseD("hRecoDijetInfo",
    //         "Reconstructed dijet and jet info;p_{T}^{ave} (GeV);#eta^{dijet};;#eta^{dijet}_{CM};#Delta#phi^{dijet} (rad);p_{T}^{Lead} (GeV);#eta^{Lead};#eta^{Lead}_{CM};p_{T}^{SubLead} (GeV);#eta^{SubLead};#eta^{SubLead}_{CM}",
    //         10,
    //         bins10D_dijet_info,
    //         xmin10D_dijet_info,
    //         xmax10D_dijet_info);
    // if ( fUseVariableBinning ) {
    //     hRecoDijetInfo->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    //     hRecoDijetInfo->GetAxis(2)->Set(dijetEtaBins, dijetEtaVals);
    //     hRecoDijetInfo->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
    //     hRecoDijetInfo->GetAxis(6)->Set(dijetEtaBins, dijetEtaVals);
    //     hRecoDijetInfo->GetAxis(8)->Set(dijetEtaBins, dijetEtaVals);
    //     hRecoDijetInfo->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
    // }
    // hRecoDijetInfo->Sumw2();
    

    hRecoPtLeadPtSublead = new TH2D("hRecoPtLeadPtSublead","Reco Lead vs SubLead p_{T};Reco p_{T}^{Lead} (GeV);Reco p_{T}^{SubLead} (GeV)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSublead->Sumw2();
    hRecoEtaLeadEtaSublead = new TH2D("hRecoEtaLeadEtaSublead","Reco Lead vs SubLead #eta;Reco #eta^{Lead};Reco #eta^{SubLead}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSublead->Sumw2();
    hRecoEtaCMLeadEtaCMSublead = new TH2D("hRecoEtaCMLeadEtaCMSublead","Reco Lead vs SubLead #eta in CM;Reco #eta^{Lead}_{CM};Reco #eta^{SubLead}_{CM}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaCMLeadEtaCMSublead->Sumw2();
    hRecoPtLeadPtSubleadMcReweight = new TH2D("hRecoPtLeadPtSubleadMcReweight","Reco Lead vs SubLead p_{T} (MC reweighted to data);Reco p_{T}^{Lead} (GeV);Reco p_{T}^{SubLead} (GeV)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSubleadMcReweight->Sumw2();
    hRecoEtaLeadEtaSubleadMcReweight = new TH2D("hRecoEtaLeadEtaSubleadMcReweight","Reco Lead vs SubLead #eta (MC reweighted to data);Reco #eta^{Lead};Reco #eta^{SubLead}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSubleadMcReweight->Sumw2();
    for (int iRun{0}; iRun<6; iRun++) {
        int runId{0};
        if (iRun == 1) { runId = 285480; }
        else if (iRun == 2) { runId = 285505; }
        else if (iRun == 3) { runId = 285517; }
        else if (iRun == 4) { runId = 285832; }
        else if (iRun == 5) { runId = 285993; }
        else { runId = 0; }
        hRecoDijetEtaCMRun[iRun] = new TH1D(Form("hRecoDijetEtaCMRun_%d", iRun), Form("Reco dijet #eta in CM 50<p_{T}^{ave}<90 GeV for runId %d;Reco #eta^{dijet}_{CM};Entries", runId),
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRecoDijetEtaCMRun[iRun]->Sumw2();
    }
    hRecoDijetEta = new TH1D("hRecoDijetEta","Reco dijet #eta;Reco #eta^{dijet};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEta->Sumw2();

    hRecoDijetPtEta = new TH2D("hRecoDijetPtEta","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet}",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEta->SetBinsLength(-1);
    }
    hRecoDijetPtEta->Sumw2();
    hRecoDijetPtEtaWeighted = new TH2D("hRecoDijetPtEtaWeighted","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet}",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaWeighted->Sumw2();
    hRecoDijetPtEtaCMInLab = new TH2D("hRecoDijetPtEtaCMInLab","Reco dijet info in CM with lab frame selection;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaCMInLab->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMInLab->Sumw2();
    hRecoDijetEtaCM = new TH1D("hRecoDijetEtaCM","Reco dijet #eta in CM;Reco #eta^{dijet}_{CM};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEtaCM->Sumw2();
    hRecoDijetPtEtaCM = new TH2D("hRecoDijetPtEtaCM","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCM->Sumw2();
    hRecoDijetPtEtaCMWeighted = new TH2D("hRecoDijetPtEtaCMWeighted","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaCMWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMWeighted->Sumw2();
    hRecoDijetPtEtaLabInCM = new TH2D("hRecoDijetPtEtaLabInCM","Reco dijet info in lab with CM selection;p_{T}^{ave} (GeV);#eta^{dijet}",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaLabInCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaLabInCM->Sumw2();
    hRecoDijetPtEtaMatched = new TH2D("hRecoDijetPtEtaMatched","Matched reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet}",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaMatched->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaMatched->SetBinsLength(-1);
    }
    hRecoDijetPtEtaMatched->Sumw2();
    hRecoDijetPtEtaCMMatched = new TH2D("hRecoDijetPtEtaCMMatched","Matched reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMMatched->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaCMMatched->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMMatched->Sumw2();

    hRecoDijetPtAveLeadPtSubLeadPt = new TH3D("hRecoDijetPtAveLeadPtSubLeadPt","Reco dijet info;p_{T}^{ave} (GeV);p_{T}^{Lead} (GeV);p_{T}^{SubLead} (GeV)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoDijetPtAveLeadPtSubLeadPt->Sumw2();
    hRecoDijetPtAveLeadPtSubLeadPtCM = new TH3D("hRecoDijetPtAveLeadPtSubLeadPtCM","Reco dijet info in CM;p_{T}^{ave} (GeV);p_{T}^{Lead} (GeV);p_{T}^{SubLead} (GeV)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoDijetPtAveLeadPtSubLeadPtCM->Sumw2();
    hRecoDijetPtAveLeadEtaSubLeadEta = new TH3D("hRecoDijetPtAveLeadEtaSubLeadEta","Reco dijet info;p_{T}^{ave} (GeV);#eta^{Lead};#eta^{SubLead}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoDijetPtAveLeadEtaSubLeadEta->Sumw2();
    hRecoDijetPtAveLeadEtaSubLeadEtaCM = new TH3D("hRecoDijetPtAveLeadEtaSubLeadEtaCM","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoDijetPtAveLeadEtaSubLeadEtaCM->Sumw2();
    hRecoDijetEtaLeadEtaSubLeadEta = new TH3D("hRecoDijetEtaLeadEtaSubLeadEta","Reco dijet eta info;;#eta^{dijet};#eta^{Lead};#eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoDijetEtaLeadEtaSubLeadEta->Sumw2();
    hRecoDijetEtaLeadEtaSubLeadEtaCM = new TH3D("hRecoDijetEtaLeadEtaSubLeadEtaCM","Reco dijet eta info in CM;#eta^{dijet}_{CM};#eta^{Lead}_{CM};#eta^{SubLead}_{CM};",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoDijetEtaLeadEtaSubLeadEtaCM->Sumw2();

    hRecoDijetLeadPtEta = new TH2D("hRecoDijetLeadPtEta","Reco lead jet after dijet selection;p_{T}^{Lead} (GeV);#eta^{Lead}",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoDijetLeadPtEta->Sumw2();
    hRecoDijetLeadPtEtaStdBins = new TH2D("hRecoDijetLeadPtEtaStdBins","Reco lead jet after dijet selection (std bins);p_{T}^{Lead} (GeV);#eta^{Lead}",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        208, -5.2, 5.2);
    hRecoDijetLeadPtEtaStdBins->GetYaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoDijetLeadPtEtaStdBins->SetBinsLength(-1);
    hRecoDijetLeadPtEtaStdBins->Sumw2();
    hRecoDijetSubLeadPtEta = new TH2D("hRecoDijetSubLeadPtEta","Reco sublead jet after dijet selection;p_{T}^{SubLead} (GeV);#eta^{SubLead}",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoDijetSubLeadPtEta->Sumw2();
    hRecoDijetSubLeadPtEtaStdBins = new TH2D("hRecoDijetSubLeadPtEtaStdBins","Reco sublead jet after dijet selection (std bins);p_{T}^{SubLead} (GeV);#eta^{SubLead}",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        208, -5.2, 5.2);
    hRecoDijetSubLeadPtEtaStdBins->GetYaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoDijetSubLeadPtEtaStdBins->SetBinsLength(-1);
    hRecoDijetSubLeadPtEtaStdBins->Sumw2();

    for (int i{0}; i<3; i++) {
        hRecoDijetXj[i] = new TH1D( Form("hRecoDijetXj_%d", i), 
                                    Form("Reco dijet x_{J} the #eta region %d (backward, midrange, forward);x_{J};Entries", i),
                                    22, 0., 1.1);
        hRecoDijetXj[i]->Sumw2();
        hRecoDijetXjCM[i] = new TH1D( Form("hRecoDijetXjCM_%d", i), 
                                    Form("Reco dijet x_{J} in CM frame the #eta region %d (backward, midrange, forward);x_{J};Entries", i),
                                    22, 0., 1.1);
        hRecoDijetXjCM[i]->Sumw2();
    }


    hRecoDijetPtEtaForward = new TH2D("hRecoDijetPtEtaForward", "Reco dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                      fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                      fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForward->Sumw2();
    hRecoDijetPtEtaForwardCMInLab = new TH2D("hRecoDijetPtEtaForwardCMInLab", "Reco dijet info in CM frame with lab selection (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForwardCMInLab->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForwardCMInLab->Sumw2();
    hRecoDijetPtEtaBackward = new TH2D("hRecoDijetPtEtaBackward", "Reco dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackward->Sumw2();
    hRecoDijetPtEtaBackwardCMInLab = new TH2D("hRecoDijetPtEtaBackwardCMInLab", "Reco dijet info in CM frame with lab selection (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackwardCMInLab->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackwardCMInLab->Sumw2();
    hRecoDijetPtEtaCMForward = new TH2D("hRecoDijetPtEtaCMForward", "Reco dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMForward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMForward->Sumw2();
    hRecoDijetPtEtaForwardLabInCM = new TH2D("hRecoDijetPtEtaForwardLabInCM", "Reco dijet info in lab frame with CM selection (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForwardLabInCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForwardLabInCM->Sumw2();
    hRecoDijetPtEtaCMBackward = new TH2D("hRecoDijetPtEtaCMBackward", "Reco dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMBackward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMBackward->Sumw2();
    hRecoDijetPtEtaBackwardLabInCM = new TH2D("hRecoDijetPtEtaBackwardLabInCM", "Reco dijet info in lab frame with CM selection (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackwardLabInCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackwardLabInCM->Sumw2();
    
    hRecoDijetPtEtaForwardWeighted = new TH2D("hRecoDijetPtEtaForwardWeighted", "Reco dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForwardWeighted->Sumw2();
    hRecoDijetPtEtaBackwardWeighted = new TH2D("hRecoDijetPtEtaBackwardWeighted", "Reco dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackwardWeighted->Sumw2();
    hRecoDijetPtEtaCMForwardWeighted = new TH2D("hRecoDijetPtEtaCMForwardWeighted", "Reco dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);

    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMForwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMForwardWeighted->Sumw2();

    hRecoDijetPtEtaCMBackwardWeighted = new TH2D("hRecoDijetPtEtaCMBackwardWeighted", "Reco dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMBackwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMBackwardWeighted->Sumw2();

    // Reco dijet distributions in different eta regions |eta_CM|<x
    for (int iEta{0}; iEta<6; iEta++) {
        const float etaCut = 1.4 + 0.1*iEta;
        hRecoDijetPtRawEtaForwardArr[iEta] = new TH2D( Form("hRecoDijetPtRawEtaForwardArr_%d", iEta), Form("Reco dijet p_{T}^{raw} vs #eta_{CM}  (forward) |#eta^{dijet}_{CM}|<%.1f;p_{T}^{ave, raw} (GeV);#eta^{dijet}_{CM}", etaCut),
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        hRecoDijetPtRawEtaForwardArr[iEta]->Sumw2();
        hRecoDijetPtRawEtaBackwardArr[iEta] = new TH2D( Form("hRecoDijetPtRawEtaBackwardArr_%d", iEta), Form("Reco dijet p_{T}^{raw} vs #eta_{CM}  (backward) |#eta^{dijet}_{CM}|<%.1f;p_{T}^{ave, raw} (GeV);#eta^{dijet}_{CM}", etaCut),
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        hRecoDijetPtRawEtaBackwardArr[iEta]->Sumw2();
        hRecoDijetPtEtaForwardArr[iEta] = new TH2D( Form("hRecoDijetPtEtaForwardArr_%d", iEta), Form("Reco dijet p_{T} vs #eta in CM frame (forward) |#eta^{dijet}|<%.1f;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}", etaCut),
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        hRecoDijetPtEtaForwardArr[iEta]->Sumw2();
        hRecoDijetPtEtaBackwardArr[iEta] = new TH2D( Form("hRecoDijetPtEtaBackwardArr_%d", iEta), Form("Reco dijet p_{T} vs #eta in CM frame (backward) |#eta^{dijet}|<%.1f;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}", etaCut),
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        hRecoDijetPtEtaBackwardArr[iEta]->Sumw2();
    } // for (int iEta{0}; iEta<6; iEta++)

    //
    // Monte Carlo information
    //
    if (fIsMc) {

        //
        // Gen inclusive jets
        //

        hGenJetCollectionSize = new TH1D("hGenJetCollectionSize","Gen jet collection size;Gen jet collection size;Entries", 100, -0.5, 99.5);
        hGenJetCollectionSize->Sumw2();

        hGenVsRecoJetCollectionSize = new TH2D("hGenVsRecoJetCollectionSize","Reco vs Gen jet collection size;Reco jet collection size;Gen jet collection size", 100, -0.5, 99.5, 100, -0.5, 99.5);
        hGenVsRecoJetCollectionSize->Sumw2();

        hGenLeadJetPtOverPtHatVsLeadJetPt = new TH2D("hGenLeadJetPtOverPtHatVsLeadJetPt", "Lead jet p_{T}^{Gen}/#hat{p}_{T} vs Lead jet p_{T}^{Gen};p_{T}^{Gen} (GeV);p_{T}^{Gen}/#hat{p}_{T}",
                                                            fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenLeadJetPtOverPtHatVsLeadJetPt->Sumw2();
        hGenLeadJetPtOverPtHatVsLeadJetPtWeighted = new TH2D("hGenLeadJetPtOverPtHatVsLeadJetPtWeighted", "Lead jet p_{T}^{Gen}/#hat{p}_{T} vs Lead jet p_{T}^{Gen} weighted;p_{T}^{Gen} (GeV);p_{T}^{Gen}/#hat{p}_{T}",
                                                                    fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenLeadJetPtOverPtHatVsLeadJetPtWeighted->Sumw2();
        hGenDijetPtOverPtHatVsDijetPt = new TH2D("hGenDijetPtOverPtHatVsDijetPt", "Gen dijet p_{T}/#hat{p}_{T} vs gen dijet p_{T};Gen dijet p_{T} (GeV);Gen dijet p_{T}/#hat{p}_{T}",
                                                  fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtOverPtHatVsDijetPt->Sumw2();
        hGenDijetPtOverPtHatVsDijetPtWeighted = new TH2D("hGenDijetPtOverPtHatVsDijetPtWeighted", "Gen dijet p_{T}/#hat{p}_{T} vs gen dijet p_{T} weighted;Gen dijet p_{T} (GeV);Gen dijet p_{T}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtOverPtHatVsDijetPtWeighted->Sumw2();
        hGenDijetPtAveOverPtHatVsDijetPtAve = new TH2D("hGenDijetPtAveOverPtHatVsDijetPtAve", "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs gen dijet p_{T}^{ave};Gen dijet p_{T}^{ave} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtAveOverPtHatVsDijetPtAve->Sumw2();
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted = new TH2D("hGenDijetPtAveOverPtHatVsDijetPtAveWeighted", "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs gen dijet p_{T}^{ave} weighted;Gen dijet p_{T}^{ave} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}",
                                                                fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Sumw2();

        hGenInclusiveJetPt = new TH1D("hGenInclusiveJetPt","Inclusive gen jet;Gen p_{T}^{inclusive} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPt->Sumw2();
        hGenInclusiveJetEta = new TH1D("hGenInclusiveJetEta","Inclusive gen jet;Gen #eta^{inclusive}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenInclusiveJetEta->Sumw2();
        hGenInclusiveJetEtaUnweighted = new TH1D("hGenInclusiveJetEtaUnweighted","Inclusive gen jet unweighted;Gen #eta^{inclusive}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenInclusiveJetEtaUnweighted->Sumw2();
        hGenInclusiveJetPtEta = new TH2D("hGenInclusiveJetPtEta","Gen inclusive jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPtEta->Sumw2();
        hGenInclusiveJetPtEtaCM = new TH2D("hGenInclusiveJetPtEtaCM","Gen inclusive jet acceptance CM frame;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPtEtaCM->Sumw2();
        hGenInclusiveJetPtEtaPtHat = new TH3D("hGenInclusiveJetPtEtaPtHat","Gen inclusive jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        if (fUseVariableBinning) {
            hGenInclusiveJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hGenInclusiveJetPtEtaPtHat->SetBinsLength(-1);
        }
        hGenInclusiveJetPtEtaPtHat->Sumw2();
        hGenLeadJetPtEta = new TH2D("hGenLeadJetPtEta","Gen Lead jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenLeadJetPtEta->Sumw2();
        hGenLeadJetPtEtaCM = new TH2D("hGenLeadJetPtEtaCM","Gen Lead jet acceptance CM frame;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenLeadJetPtEtaCM->Sumw2();
        hGenLeadJetPtEtaPtHat = new TH3D("hGenLeadJetPtEtaPtHat","Gen Lead jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        if (fUseVariableBinning) {
            hGenLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hGenLeadJetPtEtaPtHat->SetBinsLength(-1);
        }
        hGenLeadJetPtEtaPtHat->Sumw2();
        hGenSubLeadJetPtEta = new TH2D("hGenSubLeadJetPtEta","Gen SubLead jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenSubLeadJetPtEta->Sumw2();
        hGenSubLeadJetPtEtaCM = new TH2D("hGenSubLeadJetPtEtaCM","Gen SubLead jet acceptance CM frame;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenSubLeadJetPtEtaCM->Sumw2();
        hGenSubLeadJetPtEtaPtHat = new TH3D("hGenSubLeadJetPtEtaPtHat","Gen SubLead jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        if (fUseVariableBinning) {
            hGenSubLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hGenSubLeadJetPtEtaPtHat->SetBinsLength(-1);
        }
        hGenSubLeadJetPtEtaPtHat->Sumw2();

        //
        // Gen dijets
        //

        // Gen dijet info [10 dimensions]
        // 0 - dijet pt ave, 1 - dijet eta lab, 2 - dijet eta cm, 3 - dijet delta phi,
        // 4 - lead pt, 5 - lead eta lab, 6 - lead eta cm,
        // 7 - sublead pt, 8 - sublead eta lab, 9 - sublead eta cm
        // hGenDijetInfo = new THnSparseD("hGenDijetInfo","Title;Gen p_{T}^{ave} (GeV);Gen #eta^{dijet};Gen #eta^{dijet}_{CM};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV);Gen #eta^{Lead};Gen #eta^{Lead}_{CM};Gen p_{T}^{SubLead} (GeV);Gen #eta^{SubLead};Gen #eta^{SubLead}_{CM}",
        //         10,
        //         bins10D_gen_dijetInfo,
        //         xmin10D_gen_dijetInfo,
        //         xmax10D_gen_dijetInfo);
        // if (fUseVariableBinning) {
        //     hGenDijetInfo->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        //     hGenDijetInfo->GetAxis(2)->Set(dijetEtaBins, dijetEtaVals);
        //     hGenDijetInfo->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
        //     hGenDijetInfo->GetAxis(6)->Set(dijetEtaBins, dijetEtaVals);
        //     hGenDijetInfo->GetAxis(8)->Set(dijetEtaBins, dijetEtaVals);
        //     hGenDijetInfo->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        // }
        // hGenDijetInfo->Sumw2();

        hGenPtLeadPtSublead = new TH2D("hGenPtLeadPtSublead","Lead gen jet pT vs SubLead gen jet pT;Gen p_{T}^{Lead} (GeV);Gen p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenPtLeadPtSublead->Sumw2();
        hGenEtaLeadEtaSublead = new TH2D("hGenEtaLeadEtaSublead","Lead gen jet eta vs SubLead gen jet eta;Gen #eta^{Lead};Gen #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaLeadEtaSublead->Sumw2();
        hGenEtaCMLeadEtaCMSublead = new TH2D("hGenEtaCMLeadEtaCMSublead","Lead gen jet eta in CM vs SubLead gen jet eta in CM;Gen #eta^{Lead}_{CM};Gen #eta^{SubLead}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaCMLeadEtaCMSublead->Sumw2();
        hGenPtLeadPtSubleadMcReweight = new TH2D("hGenPtLeadPtSubleadMcReweight","Lead gen jet pT vs SubLead gen jet pT (MC reweighted to data);Gen p_{T}^{Lead} (GeV);Gen p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenPtLeadPtSubleadMcReweight->Sumw2();
        hGenEtaLeadEtaSubleadMcReweight = new TH2D("hGenEtaLeadEtaSubleadMcReweight","Lead gen jet eta vs SubLead gen jet eta (MC reweighted to data);Gen #eta^{Lead};Gen #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaLeadEtaSubleadMcReweight->Sumw2();
        hGenDijetEta = new TH1D("hGenDijetEta", "Gen dijet #eta;#eta^{dijet}",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEta->Sumw2();

        hGenDijetPtEta = new TH2D("hGenDijetPtEta","Gen dijet info;p_{T}^{ave} (GeV);#eta^{dijet}",
                                      fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                      fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEta->SetBinsLength(-1);
        }
        hGenDijetPtEta->Sumw2();
        hGenDijetPtEtaWeighted = new TH2D("hGenDijetPtEtaWeighted","Gen dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaWeighted->Sumw2();
        hGenDijetPtEtaCMInLab = new TH2D("hGenDijetPtEtaCMInLab","Gen dijet info in CM in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaCMInLab->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMInLab->Sumw2();


        hGenDijetEtaCM = new TH1D("hGenDijetEtaCM", "Gen dijet #eta in CM;#eta^{dijet}_{CM}",
                                  fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEtaCM->Sumw2();
        hGenDijetPtEtaCM = new TH2D("hGenDijetPtEtaCM","Gen dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaCM->Sumw2();
        hGenDijetPtEtaCMWeighted = new TH2D("hGenDijetPtEtaCMWeighted","Gen dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaCMWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMWeighted->Sumw2();
        hGenDijetPtEtaLabInCM = new TH2D("hGenDijetPtEtaLabInCM","Gen dijet info in lab in CM frame;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaLabInCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaLabInCM->Sumw2();

        hGenDijetPtEtaForward = new TH2D("hGenDijetPtEtaForward", "Gen dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForward->SetBinsLength(-1);
        }
        hGenDijetPtEtaForward->Sumw2();
        hGenDijetPtEtaForwardCMInLab = new TH2D("hGenDijetPtEtaForwardCMInLab", "Gen dijet info in CM frame in lab (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForwardCMInLab->SetBinsLength(-1);
        }
        hGenDijetPtEtaForwardCMInLab->Sumw2();
        hGenDijetPtEtaBackward = new TH2D("hGenDijetPtEtaBackward", "Gen dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackward->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackward->Sumw2();
        hGenDijetPtEtaBackwardCMInLab = new TH2D("hGenDijetPtEtaBackwardCMInLab", "Gen dijet info in CM frame in lab (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackwardCMInLab->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackwardCMInLab->Sumw2();
        hGenDijetPtEtaCMForward = new TH2D("hGenDijetPtEtaCMForward", "Gen dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaCMForward->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMForward->Sumw2();
        hGenDijetPtEtaForwardLabInCM = new TH2D("hGenDijetPtEtaForwardLabInCM", "Gen dijet info in lab frame in CM (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForwardLabInCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaForwardLabInCM->Sumw2();
        hGenDijetPtEtaCMBackward = new TH2D("hGenDijetPtEtaCMBackward", "Gen dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaCMBackward->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMBackward->Sumw2();
        hGenDijetPtEtaBackwardLabInCM = new TH2D("hGenDijetPtEtaBackwardLabInCM", "Gen dijet info in lab frame in CM (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackwardLabInCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackwardLabInCM->Sumw2();

        hGenDijetPtEtaForwardWeighted = new TH2D("hGenDijetPtEtaForwardWeighted", "Gen dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaForwardWeighted->Sumw2();
        hGenDijetPtEtaBackwardWeighted = new TH2D("hGenDijetPtEtaBackwardWeighted", "Gen dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackwardWeighted->Sumw2();
        hGenDijetPtEtaCMForwardWeighted = new TH2D("hGenDijetPtEtaCMForwardWeighted", "Gen dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaCMForwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMForwardWeighted->Sumw2();
        hGenDijetPtEtaCMBackwardWeighted = new TH2D("hGenDijetPtEtaCMBackwardWeighted", "Gen dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaCMBackwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMBackwardWeighted->Sumw2();

        // Gen dijet forward-backward jets for different |eta| selections: <1.4, 1.5, 1.6, 1.7, 1.8, 1.9
        for (int iEta{0}; iEta<6; ++iEta) {
            const float etaCut = 1.4 + iEta*0.1;
            hGenDijetPtEtaForwardArr[iEta] = new TH2D(Form("hGenDijetPtEtaForwardArr_%d", iEta), Form("Gen dijet distribuition (forward) in CM frame |#eta^{jet}|<%2.1f;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}", etaCut),
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
            hGenDijetPtEtaForwardArr[iEta]->Sumw2();
            hGenDijetPtEtaBackwardArr[iEta] = new TH2D(Form("hGenDijetPtEtaBackwardArr_%d", iEta), Form("Gen dijet distribuition (backward) in CM frame |#eta^{jet}|<%2.1f;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}", etaCut),
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
            hGenDijetPtEtaBackwardArr[iEta]->Sumw2();
        } // for (int iEta{0}; iEta<6; ++iEta)

        hGenDijetPtAveLeadPtSubLeadPt = new TH3D("hGenDijetPtAveLeadPtSubLeadPt", "Gen dijet pT ave vs Lead pT vs SubLead pT;p_{T}^{ave} (GeV);p_{T}^{Lead} (GeV);p_{T}^{SubLead} (GeV)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hGenDijetPtAveLeadPtSubLeadPt->Sumw2();
        hGenDijetPtAveLeadPtSubLeadPtCM = new TH3D("hGenDijetPtAveLeadPtSubLeadPtCM", "Gen dijet pT ave vs Lead pT vs SubLead pT in CM;p_{T}^{ave} (GeV);p_{T}^{Lead} (GeV);p_{T}^{SubLead} (GeV)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hGenDijetPtAveLeadPtSubLeadPtCM->Sumw2();
        hGenDijetPtAveLeadEtaSubLeadEta = new TH3D("hGenDijetPtAveLeadEtaSubLeadEta", "Gen dijet pT ave vs Lead #eta vs SubLead #eta;p_{T}^{ave} (GeV);#eta^{Lead};#eta^{SubLead}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenDijetPtAveLeadEtaSubLeadEta->Sumw2();
        hGenDijetPtAveLeadEtaSubLeadEtaCM = new TH3D("hGenDijetPtAveLeadEtaSubLeadEtaCM", "Gen dijet pT ave vs Lead #eta vs SubLead #eta in CM;p_{T}^{ave} (GeV);#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenDijetPtAveLeadEtaSubLeadEtaCM->Sumw2();
        hGenDijetEtaLeadEtaSubLeadEta = new TH3D("hGenDijetEtaLeadEtaSubLeadEta", "Gen dijet #eta vs Lead #eta vs SubLead #eta;#eta^{dijet};#eta^{Lead};#eta^{SubLead}",
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenDijetEtaLeadEtaSubLeadEta->Sumw2();
        hGenDijetEtaLeadEtaSubLeadEtaCM = new TH3D("hGenDijetEtaLeadEtaSubLeadEtaCM", "Gen dijet #eta vs Lead #eta vs SubLead #eta in CM;#eta^{dijet}_{CM};#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenDijetEtaLeadEtaSubLeadEtaCM->Sumw2();

        hGenGoodInclusiveJetEtaLabFrame = new TH1D("hGenGoodInclusiveJetEtaLabFrame","Gen good inclusive jet #eta in lab frame;#eta^{Inclusive}",
                                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenGoodInclusiveJetEtaLabFrame->Sumw2();
        hGenGoodInclusiveJetEtaCMFrame = new TH1D("hGenGoodInclusiveJetEtaCMFrame","Gen good inclusive jet #eta in CM frame;#eta^{Inclusive}_{CM}",
                                                  fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenGoodInclusiveJetEtaCMFrame->Sumw2();

        // xPb and xP part
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
        hGenSelectedDijetXPbOverXpEta = new TH2D("hGenSelectedDijetXPbOverXpEta", "Selected Gen Dijet x_{Pb}/x_{p} vs #eta^{dijet};x_{Pb}/x_{p};#eta^{dijet}", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenSelectedDijetXPbOverXpEta->Sumw2();
        hGenSelectedDijetXPbOverXpEtaWeighted = new TH2D("hGenSelectedDijetXPbOverXpEtaWeighted", "Selected Gen Dijet x_{Pb}/x_{p} vs #eta^{dijet} (weighted);x_{Pb}/x_{p};#eta^{dijet}", xPbOverXpBins, xPbOverXpRange[0], xPbOverXpRange[1], fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenSelectedDijetXPbOverXpEtaWeighted->Sumw2();

        // Reco single jets
        hRecoInclusiveUnmatchedJetPtEta = new TH2D("hRecoInclusiveUnmatchedJetPtEta", "Inclusive reco jet unmatched gen pT vs eta;#eta;p_{T} (GeV)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveUnmatchedJetPtEta->Sumw2();

        hRecoInclusiveJetReco2Ref = new THnSparseD("hRecoInclusiveJetReco2Ref","Reconstructed inclusive jets;Reco p_{T, corr}^{Inclusive} (GeV);Reco p_{T, raw}^{Inclusive} (GeV);Ref p_{T}^{Inclusive} (GeV);Reco #eta^{Inclusive};Ref #eta^{Inclusive}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoInclusiveJetReco2Ref->Sumw2();
        hRecoLeadJetReco2Ref = new THnSparseD("hRecoLeadJetReco2Ref","Reconstructed Lead jets;Reco p_{T, corr}^{Lead} (GeV);Reco p_{T, raw}^{Lead} (GeV);Ref p_{T}^{Lead} (GeV);Reco #eta^{Lead};Ref #eta^{Lead}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoLeadJetReco2Ref->Sumw2();
        hRecoSubLeadJetReco2Ref = new THnSparseD("hRecoSubLeadJetReco2Ref","Reconstructed SubLead jets;Reco p_{T, corr}^{SubLead} (GeV);Reco p_{T, raw}^{SubLead} (GeV);Ref p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref #eta^{SubLead}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoSubLeadJetReco2Ref->Sumw2();
        hJESInclusiveJetPtEtaPhi = new THnSparseD("hJESInclusiveJetPtEtaPhi","JES of inclusive jets;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV);#eta^{reco};#phi^{reco} (rad)",
                4,
                bins4D_jet_JESPtEtaPhi,
                xmin4D_jet_JESPtEtaPhi,
                xmax4D_jet_JESPtEtaPhi);
        hJESInclusiveJetPtEtaPhi->Sumw2();

        hRecoLeadJetPtOverPtHatVsLeadJetPt = new TH2D("hRecoLeadJetPtOverPtHatVsLeadJetPt", "Reco Lead jet p_{T}/#hat{p}_{T} vs reco Lead jet p_{T}^{reco};Reco Lead jet p_{T} (GeV);Reco Lead jet p_{T}/#hat{p}_{T}",
                                                            fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoLeadJetPtOverPtHatVsLeadJetPt->Sumw2();
        hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted = new TH2D("hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted", "Reco Lead jet p_{T}/#hat{p}_{T} vs reco Lead jet p_{T}^{reco} weighted;Reco Lead jet p_{T} (GeV);Reco Lead jet p_{T}/#hat{p}_{T}",
                                                                    fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted->Sumw2();
        hRecoDijetPtOverPtHatVsDijetPt = new TH2D("hRecoDijetPtOverPtHatVsDijetPt", "Reco dijet p_{T}/#hat{p}_{T} vs dijet p_{T};Reco dijet p_{T} (GeV);Reco dijet p_{T}/#hat{p}_{T}",
                                                  fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtOverPtHatVsDijetPt->Sumw2();
        hRecoDijetPtOverPtHatVsDijetPtWeighted = new TH2D("hRecoDijetPtOverPtHatVsDijetPtWeighted", "Reco dijet p_{T}/#hat{p}_{T} vs reco dijet p_{T} weighted;Reco dijet p_{T}^{reco} (GeV);Reco dijet p_{T}^{reco}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtOverPtHatVsDijetPtWeighted->Sumw2();
        hRecoDijetPtAveOverPtHatVsDijetPtAve = new TH2D("hRecoDijetPtAveOverPtHatVsDijetPtAve", "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs reco dijet p_{T}^{ave};Reco dijet p_{T}^{ave} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtAveOverPtHatVsDijetPtAve->Sumw2();
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted = new TH2D("hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted", "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs reco dijet p_{T}^{ave} weighted;Reco dijet p_{T}^{ave} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}",
                                                                fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Sumw2();

        // hRecoInclusiveJetJECFactorVsPtEta = new TH3D("hRecoInclusiveJetJECFactorVsPtEta","JEC factor vs p_{T} and #eta;p_{T}^{corr}/p_{T}^{raw};p_{T}^{gen} (GeV);#eta^{gen};JEC factor",
        //                                    20, 0., 2., fPtBins, fPtRange[0], fPtRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        // hRecoInclusiveJetJECFactorVsPtEta->Sumw2(); 
        // hRecoInclusiveJetJEC2FactorVsPtEta = new TH3D("hRecoInclusiveJetJEC2FactorVsPtEta","JEC factor (my) vs p_{T} and #eta;p_{T}^{my, corr}/p_{T}^{raw};p_{T}^{gen} (GeV);#eta^{gen};JEC2 factor",
        //                                    20, 0., 2., fPtBins, fPtRange[0], fPtRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        // hRecoInclusiveJetJEC2FactorVsPtEta->Sumw2();
        // hRecoInclusiveJetPtRawOverPtRefVsPtEta = new TH3D("hRecoInclusiveJetPtRawOverPtRefVsPtEta","p_{T}^{raw}/p_{T}^{ref} vs p_{T} and #eta;p_{T}^{raw}/p_{T}^{ref};p_{T}^{gen} (GeV);#eta^{gen}",
        //                                    20, 0., 2., fPtBins, fPtRange[0], fPtRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        // hRecoInclusiveJetPtRawOverPtRefVsPtEta->Sumw2();
        // hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning = new TH3D("hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning","p_{T}^{raw}/p_{T}^{ref} vs p_{T} and #eta (std binning);p_{T}^{raw}/p_{T}^{ref};p_{T}^{gen} (GeV);#eta^{gen}",
        //                                                             20, 0., 2., 1300, 10., 6510., 104, -5.2, 5.2);
        // if ( fUseVariableBinning ) {
        //     hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->GetZaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        //     hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->SetBinsLength(-1);
        // }
        // hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->Sumw2();
        
        // hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning = new TH3D("hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning","p_{T}^{raw}/p_{T}^{ref} vs p_{T}^{reco} and #eta (std binning);p_{T}^{raw}/p_{T}^{ref};p_{T}^{reco} (GeV);#eta^{reco}",
        //                                                             20, 0., 2., 1300, 10., 6510., 104, -5.2, 5.2);
        // if ( fUseVariableBinning ) {
        //     hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning->GetZaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        //     hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning->SetBinsLength(-1);
        // }
        // hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning->Sumw2();


        hInclusiveJetJESVsPtGen = new TH2D("hInclusiveJetJESVsPtGen","JES vs p_{T}^{gen} for |#eta|<1.4;p_{T}^{gen} (GeV);p_{T}^{reco}/p_{T}^{gen}",
                                           fPtBins, fPtRange[0], fPtRange[1], fJESBins, fJESRange[0], fJESRange[1]);
        hInclusiveJetJESVsPtGen->Sumw2();
        hInclusiveJetJESGenPtGenEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESGenPtGenEtaPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hInclusiveJetJESGenPtGenEtaPtHatWeighted->Sumw2();
        hInclusiveJetJESGenPtGenEtaCMPtHatWeighted = new THnSparseD("hInclusiveJetJESGenPtGenEtaCMPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hInclusiveJetJESGenPtGenEtaCMPtHatWeighted->Sumw2();
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESRecoPtRecoEtaPtHatWeighted","JES vs p_{T}^{reco} vs #eta^{reco} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{reco} (GeV);#eta^{reco};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Sumw2();
        hLeadJetJESGenPtEtaPtHatWeighted = new THnSparseD("hLeadJetJESGenPtEtaPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hLeadJetJESGenPtEtaPtHatWeighted->Sumw2();
        hLeadJetJESGenPtEtaCMPtHatWeighted = new THnSparseD("hLeadJetJESGenPtEtaCMPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hLeadJetJESGenPtEtaCMPtHatWeighted->Sumw2();
        hSubLeadJetJESGenPtEtaPtHatWeighted = new THnSparseD("hSubLeadJetJESGenPtEtaPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hSubLeadJetJESGenPtEtaPtHatWeighted->Sumw2();
        hSubLeadJetJESGenPtEtaCMPtHatWeighted = new THnSparseD("hSubLeadJetJESGenPtEtaCMPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hSubLeadJetJESGenPtEtaCMPtHatWeighted->Sumw2();

        hRecoInclusiveMatchedJetPt = new TH1D("hRecoInclusiveMatchedJetPt","Inclusive reco jet that has matching;p_{T} (GeV)",
                                                    fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveMatchedJetPt->Sumw2();
        hRecoInclusiveMatchedJetPtEta = new TH2D("hRecoInclusiveMatchedJetPtEta","Inclusive reco jet that has matching pT vs eta;#eta;p_{T} (GeV)",
                                                    fEtaBins, fEtaRange[0], fEtaRange[1],
                                                    fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveMatchedJetPtEta->Sumw2();

        hRecoLeadMatchedJetPtEta = new TH2D("hRecoLeadMatchedJetPtEta","Lead jet matched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                              fEtaBins, fEtaRange[0], fEtaRange[1], 
                                              fPtBins, fPtRange[0], fPtRange[1]);
        hRecoLeadMatchedJetPtEta->Sumw2();
        hRecoLeadMatchedJetPtEtaPtHat = new TH3D("hRecoLeadMatchedJetPtEtaPtHat","Lead jet matched p_{T} vs #eta and #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)", 
                                              fEtaBins, fEtaRange[0], fEtaRange[1], 
                                              fPtBins, fPtRange[0], fPtRange[1], 
                                              fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRecoLeadMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRecoLeadMatchedJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRecoLeadMatchedJetPtEtaPtHat->Sumw2();
        hRecoLeadUnmatchedJetPtEta = new TH2D("hRecoLeadUnmatchedJetPtEta","Lead jet unmatched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                              fEtaBins, fEtaRange[0], fEtaRange[1], 
                                              fPtBins, fPtRange[0], fPtRange[1]);
        hRecoLeadUnmatchedJetPtEta->Sumw2();
        hRecoLeadUnmatchedJetPtEtaPtHat = new TH3D("hRecoLeadUnmatchedJetPtEtaPtHat","Lead jet unmatched p_{T} vs #eta and #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)", 
                                                fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                fPtBins, fPtRange[0], fPtRange[1], 
                                                fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRecoLeadUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRecoLeadUnmatchedJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRecoLeadUnmatchedJetPtEtaPtHat->Sumw2();
        hRecoSubLeadMatchedJetPtEta = new TH2D("hRecoSubLeadMatchedJetPtEta","SubLead jet matched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                                 fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                 fPtBins, fPtRange[0], fPtRange[1]);
        hRecoSubLeadMatchedJetPtEta->Sumw2();
        hRecoSubLeadMatchedJetPtEtaPtHat = new TH3D("hRecoSubLeadMatchedJetPtEtaPtHat","SubLead jet matched p_{T} vs #eta and #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)", 
                                                 fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                 fPtBins, fPtRange[0], fPtRange[1], 
                                                 fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRecoSubLeadMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRecoSubLeadMatchedJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRecoSubLeadMatchedJetPtEtaPtHat->Sumw2();
        hRecoSubLeadUnmatchedJetPtEta = new TH2D("hRecoSubLeadUnmatchedJetPtEta","SubLead jet unmatched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                   fPtBins, fPtRange[0], fPtRange[1]);
        hRecoSubLeadUnmatchedJetPtEta->Sumw2();
        hRecoSubLeadUnmatchedJetPtEtaPtHat = new TH3D("hRecoSubLeadUnmatchedJetPtEtaPtHat","SubLead jet unmatched p_{T} vs #eta and #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)", 
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                   fPtBins, fPtRange[0], fPtRange[1], 
                                                   fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRecoSubLeadUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRecoSubLeadUnmatchedJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRecoSubLeadUnmatchedJetPtEtaPtHat->Sumw2();


        // Reco dijet with MC

        // Reco 2 Ref parameters [12 dimensions]
        // 0 - reco dijet ptAve, 1 - dijet eta,
        // 2 - reco lead pt, 3 - lead eta,
        // 4 - reco sublead pt, 5 - sublead eta,
        // 6 - ref dijet ptAve, 7 - dijet eta,
        // 8 - ref lead pt, 9 - lead eta,
        // 10 - ref sublead pt, 11 - sublead eta
        // hReco2RefFull = new THnSparseD("hReco2RefFull",
        //         "Reco to ref correspondence;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Lead} (GeV);Reco #eta^{Lead};Reco p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Lead} (GeV);Ref #eta^{Lead};Ref p_{T}^{SubLead} (GeV);Ref #eta^{SubLead}",
        //         12,
        //         bins12D_ref2reco_info,
        //         xmin12D_ref2reco_info,
        //         xmax12D_ref2reco_info);
        // if ( fUseVariableBinning ) {
        //     hReco2RefFull->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        //     hReco2RefFull->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        //     hReco2RefFull->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
        //     hReco2RefFull->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        //     hReco2RefFull->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        //     hReco2RefFull->GetAxis(11)->Set(dijetEtaBins, dijetEtaVals);
        // } 
        // hReco2RefFull->Sumw2();

        hRecoDijetPtEtaRefDijetPtEta = new THnSparseD("hRecoDijetPtEtaRefDijetPtEta","Reco dijet vs Ref dijet;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}; Ref p_{T}^{ave} (GeV); Ref #eta^{dijet}",
                                                      4,
                                                      bins4D_ref2reco_info,
                                                      xmin4D_ref2reco_info,
                                                      xmax4D_ref2reco_info);
        if ( fUseVariableBinning ) {
            hRecoDijetPtEtaRefDijetPtEta->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
            hRecoDijetPtEtaRefDijetPtEta->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        }
        hRecoDijetPtEtaRefDijetPtEta->Sumw2();
        hRecoDijetPtEtaRefDijetPtEtaWeighted = new THnSparseD("hRecoDijetPtEtaRefDijetPtEtaWeighted","Reco dijet vs Ref dijet;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}; Ref p_{T}^{ave} (GeV); Ref #eta^{dijet}",
                                                      4,
                                                      bins4D_ref2reco_info,
                                                      xmin4D_ref2reco_info,
                                                      xmax4D_ref2reco_info);
        if ( fUseVariableBinning ) {
            hRecoDijetPtEtaRefDijetPtEtaWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
            hRecoDijetPtEtaRefDijetPtEtaWeighted->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        }
        hRecoDijetPtEtaRefDijetPtEtaWeighted->Sumw2();

        //
        // Ref selected dijets
        //

        hRefSelInclusiveJetPt = new TH1D("hRefSelInclusiveJetPt","Ref-selected jet p_{T};Ref p_{T} (GeV);Entries",
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelInclusiveJetPt->Sumw2();
        hRefSelInclusiveJetEta = new TH1D("hRefSelInclusiveJetEta","Ref-selected jet #eta;Ref #eta;Entries",
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefSelInclusiveJetEta->Sumw2();
        hRefSelInclusiveJetEtaUnweighted = new TH1D("hRefSelInclusiveJetEtaUnweighted","Ref-selected jet #eta (unweighted);Ref #eta;Entries",
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefSelInclusiveJetEtaUnweighted->Sumw2();
        hRefSelInclusiveJetPtEta = new TH2D("hRefSelInclusiveJetPtEta","Ref-selected jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelInclusiveJetPtEta->Sumw2();
        hRefSelInclusiveJetPtEtaPtHat = new TH3D("hRefSelInclusiveJetPtEtaPtHat","Ref-selected jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelInclusiveJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefSelInclusiveJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefSelInclusiveJetPtEtaPtHat->Sumw2();
        hRefSelLeadJetPtEta = new TH2D("hRefSelLeadJetPtEta","Ref-selected Lead jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelLeadJetPtEta->Sumw2();
        hRefSelLeadJetPtEtaPtHat = new TH3D("hRefSelLeadJetPtEtaPtHat","Ref-selected Lead jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefSelLeadJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefSelLeadJetPtEtaPtHat->Sumw2();
        hRefSelSubLeadJetPtEta = new TH2D("hRefSelSubLeadJetPtEta","Ref-selected SubLead jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelSubLeadJetPtEta->Sumw2();
        hRefSelSubLeadJetPtEtaPtHat = new TH3D("hRefSelSubLeadJetPtEtaPtHat","Ref-selected SubLead jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelSubLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefSelSubLeadJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefSelSubLeadJetPtEtaPtHat->Sumw2();

        // Reco 2 Ref parameters (RefSelected) [12 dimensions]
        // 0 - reco dijet ptAve, 1 - dijet eta,
        // 2 - reco lead pt, 3 - lead eta,
        // 4 - reco sublead pt, 5 - sublead eta,
        // 6 - ref dijet ptAve, 7 - dijet eta,
        // 8 - ref lead pt, 9 - lead eta,
        // 10 - ref sublead pt, 11 - sublead eta
        // hRefSel2RecoFull = new THnSparseD("hRefSel2RecoFull",
        //         "Reco to ref correspondence (via ref selection) weighted;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Lead} (GeV);Reco #eta^{Lead};Reco p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Lead} (GeV);Ref #eta^{Lead};Ref p_{T}^{SubLead} (GeV);Ref #eta^{SubLead}",
        //         12,
        //         bins12D_ref2reco_info,
        //         xmin12D_ref2reco_info,
        //         xmax12D_ref2reco_info);
        // if ( fUseVariableBinning ) {
        //     hRefSel2RecoFull->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        //     hRefSel2RecoFull->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        //     hRefSel2RecoFull->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
        //     hRefSel2RecoFull->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        //     hRefSel2RecoFull->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        //     hRefSel2RecoFull->GetAxis(11)->Set(dijetEtaBins, dijetEtaVals);
        // } 
        // hRefSel2RecoFull->Sumw2();

        hRefSelDijetEta = new TH1D("hRefSelDijetEta","Ref selected dijets;#eta^{dijet}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEta->Sumw2();
        hRefSelDijetPtEta = new TH2D("hRefSelDijetPtEta","RefSel dijet info;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEta->SetBinsLength(-1);
        }
        hRefSelDijetPtEta->Sumw2();
        hRefSelDijetPtEtaWeighted = new TH2D("hRefSelDijetPtEtaWeighted","RefSel dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaWeighted->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaWeighted->Sumw2();
        hRefSelDijetPtEtaCMInLab = new TH2D("hRefSelDijetPtEtaCMInLab","RefSel dijet info in CM in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaCMInLab->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaCMInLab->Sumw2();
        hRefSelDijetEtaCM = new TH1D("hRefSelDijetEtaCM","Ref selected dijets in CM;#eta^{dijet}_{CM}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEtaCM->Sumw2();
        hRefSelDijetPtEtaCM = new TH2D("hRefSelDijetPtEtaCM","RefSel dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaCM->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaCM->Sumw2();
        hRefSelDijetPtEtaCMWeighted = new TH2D("hRefSelDijetPtEtaCMWeighted","RefSel dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaCMWeighted->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaCMWeighted->Sumw2();
        hRefSelDijetPtEtaLabInCM = new TH2D("hRefSelDijetPtEtaLabInCM","RefSel dijet info in lab in CM frame;p_{T}^{ave} (GeV);#eta^{dijet}_{lab}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaLabInCM->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaLabInCM->Sumw2();

        //
        // Ref jet histograms
        //

        hRefInclusiveJetPt = new TH1D("hRefInclusiveJetPt","Ref jet p_{T};Ref p_{T} (GeV);Entries",
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPt->Sumw2();
        hRefInclusiveJetEta = new TH1D("hRefInclusiveJetEta","Ref jet #eta;Ref #eta;Entries",
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefInclusiveJetEta->Sumw2();
        hRefInclusiveJetEtaUnweighted = new TH1D("hRefInclusiveJetEtaUnweighted","Ref jet #eta (unweighted);Ref #eta;Entries",
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefInclusiveJetEtaUnweighted->Sumw2();
        hRefInclusiveJetPtEta = new TH2D("hRefInclusiveJetPtEta","Ref jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPtEta->Sumw2();
        hRefInclusiveJetPtEtaCM = new TH2D("hRefInclusiveJetPtEtaCM","Ref jet p_{T} vs #eta CM frame;#eta;p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]); 
        hRefInclusiveJetPtEtaCM->Sumw2();
        hRefInclusiveJetPtEtaPtHat = new TH3D("hRefInclusiveJetPtEtaPtHat","Ref jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefInclusiveJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefInclusiveJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefInclusiveJetPtEtaPtHat->Sumw2();
        hRefLeadJetPtEta = new TH2D("hRefLeadJetPtEta","Ref Lead jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefLeadJetPtEta->Sumw2();
        hRefLeadJetPtEtaCM = new TH2D("hRefLeadJetPtEtaCM","Ref Lead jet p_{T} vs #eta CM frame;#eta;p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefLeadJetPtEtaCM->Sumw2();
        hRefLeadJetPtEtaPtHat = new TH3D("hRefLeadJetPtEtaPtHat","Ref Lead jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefLeadJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefLeadJetPtEtaPtHat->Sumw2();
        hRefLeadUnswappedJetPtEta = new TH2D("hRefLeadUnswappedJetPtEta","Ref Lead unswapped jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefLeadUnswappedJetPtEta->Sumw2();
        hRefLeadUnswappedJetPtEtaPtHat = new TH3D("hRefLeadUnswappedJetPtEtaPtHat","Ref Lead unswapped jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefLeadUnswappedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefLeadUnswappedJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefLeadUnswappedJetPtEtaPtHat->Sumw2();
        hRefSubLeadJetPtEta = new TH2D("hRefSubLeadJetPtEta","Ref SubLead jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSubLeadJetPtEta->Sumw2();
        hRefSubLeadJetPtEtaCM = new TH2D("hRefSubLeadJetPtEtaCM","Ref SubLead jet p_{T} vs #eta CM frame;#eta;p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSubLeadJetPtEtaCM->Sumw2();
        hRefSubLeadJetPtEtaPtHat = new TH3D("hRefSubLeadJetPtEtaPtHat","Ref SubLead jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefSubLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefSubLeadJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefSubLeadJetPtEtaPtHat->Sumw2();
        hRefSubLeadUnswappedJetPtEta = new TH2D("hRefSubLeadJetUnswappedPtEta","Ref SubLead unswapped jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSubLeadUnswappedJetPtEta->Sumw2();
        hRefSubLeadUnswappedJetPtEtaPtHat = new TH3D("hRefSubLeadJetUnswappedPtEtaPtHat","Ref SubLead unswapped jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        if ( fUseVariableBinning ) {
            hRefSubLeadUnswappedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
            hRefSubLeadUnswappedJetPtEtaPtHat->SetBinsLength(-1);
        }
        hRefSubLeadUnswappedJetPtEtaPtHat->Sumw2();

        // Ref dijets
        hRefPtLeadPtSublead = new TH2D("hRefPtLeadPtSublead","Ref Lead vs SubLead p_{T};Ref p_{T}^{Lead} (GeV);Ref p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSublead->Sumw2();
        hRefEtaLeadEtaSublead = new TH2D("hRefEtaLeadEtaSublead","Ref Lead vs SubLead #eta;Ref #eta^{Lead};Ref #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaLeadEtaSublead->Sumw2();
        hRefEtaCMLeadEtaCMSublead = new TH2D("hRefEtaCMLeadEtaCMSublead","Ref Lead vs SubLead #eta in CM;Ref #eta^{Lead}_{CM};Ref #eta^{SubLead}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaCMLeadEtaCMSublead->Sumw2();
        hRefPtLeadPtSubleadMcReweight = new TH2D("hRefPtLeadPtSubleadMcReweight","Ref Lead vs SubLead p_{T} (MC reweighted to data);Ref p_{T}^{Lead} (GeV);Ref p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSubleadMcReweight->Sumw2();
        hRefEtaLeadEtaSubleadMcReweight = new TH2D("hRefEtaLeadEtaSubleadMcReweight","Ref Lead vs SubLead #eta (MC reweighted to data);Ref #eta^{Lead};Ref #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaLeadEtaSubleadMcReweight->Sumw2();
        hRefDijetEta = new TH1D("hRefDijetEta","Ref dijet #eta;Ref #eta^{dijet};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEta->Sumw2();


        hRefDijetEtaVsRecoDijetEta = new TH2D("hRefDijetEtaVsRecoDijetEta","Ref dijet #eta vs reco dijet #eta;Reco #eta^{dijet};Ref #eta^{dijet}",
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaVsRecoDijetEta->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T};Reco #eta^{dijet};Ref #eta^{dijet}; Reco p_{T}^{dijet} (GeV)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} weighted;Reco #eta^{dijet};Ref #eta^{dijet}; Reco p_{T}^{dijet} (GeV)",
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Sumw2();
        hRefDijetPtEta = new TH2D("hRefDijetPtEta","Ref dijet info;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEta->SetBinsLength(-1);
        }
        hRefDijetPtEta->Sumw2();
        hRefDijetPtEtaWeighted = new TH2D("hRefDijetPtEtaWeighted","Ref dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaWeighted->Sumw2();
        hRefDijetPtEtaCMInLab = new TH2D("hRefDijetPtEtaCMInLab","Ref dijet info in CM in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaCMInLab->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMInLab->Sumw2();

        hRefDijetEtaCM = new TH1D("hRefDijetEtaCM","Ref dijet #eta in CM;Ref #eta^{dijet}_{CM};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        if (fUseVariableBinning) {
            hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->SetBinsLength(-1);
        }
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} weighted in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV)",
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        if (fUseVariableBinning) {
            hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->SetBinsLength(-1);
        }
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Sumw2();
        hRefDijetPtEtaCM = new TH2D("hRefDijetPtEtaCM","Ref dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaCM->Sumw2();
        hRefDijetPtEtaCMWeighted = new TH2D("hRefDijetPtEtaCMWeighted","Ref dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaCMWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMWeighted->Sumw2();
        hRefDijetPtEtaLabInCM = new TH2D("hRefDijetPtEtaLabInCM","Ref dijet info in lab frame in CM;p_{T}^{ave} (GeV);#eta^{dijet}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaLabInCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaLabInCM->Sumw2();

        hRefDijetPtEtaForward = new TH2D("hRefDijetPtEtaForward", "Ref dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForward->SetBinsLength(-1);
        }
        hRefDijetPtEtaForward->Sumw2();
        hRefDijetPtEtaForwardCMInLab = new TH2D("hRefDijetPtEtaForwardCMInLab", "Ref dijet info in CM frame (forward) in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForwardCMInLab->SetBinsLength(-1);
        }
        hRefDijetPtEtaForwardCMInLab->Sumw2();
        hRefDijetPtEtaBackward = new TH2D("hRefDijetPtEtaBackward", "Ref dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackward->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackward->Sumw2();
        hRefDijetPtEtaBackwardCMInLab = new TH2D("hRefDijetPtEtaBackwardCMInLab", "Ref dijet info in CM frame (backward) in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackwardCMInLab->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackwardCMInLab->Sumw2();
        hRefDijetPtEtaCMForward = new TH2D("hRefDijetPtEtaCMForward", "Ref dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMForward->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMForward->Sumw2();
        hRefDijetPtEtaForwardLabInCM = new TH2D("hRefDijetPtEtaForwardLabInCM", "Ref dijet info in lab frame (forward) in CM;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForwardLabInCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaForwardLabInCM->Sumw2();
        hRefDijetPtEtaCMBackward = new TH2D("hRefDijetPtEtaCMBackward", "Ref dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMBackward->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMBackward->Sumw2();
        hRefDijetPtEtaBackwardLabInCM = new TH2D("hRefDijetPtEtaBackwardLabInCM", "Ref dijet info in lab frame (backward) in CM;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackwardLabInCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackwardLabInCM->Sumw2();
        
        hRefDijetPtEtaForwardWeighted = new TH2D("hRefDijetPtEtaForwardWeighted", "Ref dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaForwardWeighted->Sumw2();
        hRefDijetPtEtaBackwardWeighted = new TH2D("hRefDijetPtEtaBackwardWeighted", "Ref dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                                  fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                  fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackwardWeighted->Sumw2();
        hRefDijetPtEtaCMForwardWeighted = new TH2D("hRefDijetPtEtaCMForwardWeighted", "Ref dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                   fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMForwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMForwardWeighted->Sumw2();

        hRefDijetPtEtaCMBackwardWeighted = new TH2D("hRefDijetPtEtaCMBackwardWeighted", "Ref dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                                    fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                    fDijetEtaFBBins, fDijetEtaFBRange[0], fDijetEtaFBRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMBackwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMBackwardWeighted->Sumw2();
    } // if (fIsMc)
}

//________________
void HistoManagerDiJet::writeOutput() {

    //
    // Event histograms
    //
    hVz->Write();
    hVzWeighted->Write();
    hPtHat->Write();
    hPtHatWeighted->Write();
    hHiBin->Write();
    hHiBinWeighted->Write();

    hVzRecoDijetLab->Write();
    hVzRecoDijetLabWeighted->Write();
    hHiBinRecoDijetLab->Write();
    hHiBinRecoDijetLabWeighted->Write();

    hVzRecoDijetCM->Write();
    hVzRecoDijetCMWeighted->Write();
    hHiBinRecoDijetCM->Write();
    hHiBinRecoDijetCMWeighted->Write();

    if ( fIsMc ) {

        hVzGenDijetLab->Write();
        hVzGenDijetLabWeighted->Write();
        hHiBinGenDijetLab->Write();
        hHiBinGenDijetLabWeighted->Write();

        hVzGenDijetCM->Write();
        hVzGenDijetCMWeighted->Write();
        hHiBinGenDijetCM->Write();
        hHiBinGenDijetCMWeighted->Write();

        hVzRefSelDijetLab->Write();
        hVzRefSelDijetLabWeighted->Write();
        hHiBinRefSelDijetLab->Write();
        hHiBinRefSelDijetLabWeighted->Write();

        hVzRefSelDijetCM->Write();
        hVzRefSelDijetCMWeighted->Write();
        hHiBinRefSelDijetCM->Write();
        hHiBinRefSelDijetCMWeighted->Write();
    }
    //
    // Reco histograms
    //

    for (int i = 0; i < 4; ++i) {
        hRecoInclusiveJetNHF[i]->Write();
        hRecoInclusiveJetNEmF[i]->Write();
        hRecoInclusiveJetNumOfConst[i]->Write();
        hRecoInclusiveJetMUF[i]->Write();
        hRecoInclusiveJetCHF[i]->Write();
        hRecoInclusiveJetChargedMult[i]->Write();
        hRecoInclusiveJetCEmF[i]->Write();
        hRecoInclusiveJetNumOfNeutPart[i]->Write();
    }

    hRecoJetCollectionSize->Write();
    hRecoInclusiveAllJetPtRawEta->Write();
    // hRecoDijetInfo->Write();
    hRecoInclusiveAllJetPt->Write();
    hRecoInclusiveAllJetEta->Write();
    hRecoInclusiveAllJetEtaUnweighted->Write();
    hRecoInclusiveAllJetPtEta->Write();
    hRecoInclusiveAllJetPtEtaCM->Write();
    hRecoInclusiveAllJetPtRawEtaStdBins->Write();
    hRecoInclusiveAllJetPtEtaStdBins->Write();
    hRecoInclusiveAllJetPtEtaPtHat->Write();
    hRecoInclusiveAllJetPtRawEtaPtHat->Write();
    hRecoInclusiveMatchedJetPtEtaPtHat->Write();
    hRecoInclusiveUnmatchedJetPtEtaPtHat->Write();
    for (int i{0}; i<6; ++i) {
        hRecoInclusiveJetEtaRun[i]->Write();
        hRecoLeadJetEtaRun[i]->Write();
        hRecoSubLeadJetEtaRun[i]->Write();
    }
    hRecoPtLeadPtSublead->Write();
    hRecoEtaLeadEtaSublead->Write();
    hRecoEtaCMLeadEtaCMSublead->Write();
    hRecoPtLeadPtSubleadMcReweight->Write();
    hRecoEtaLeadEtaSubleadMcReweight->Write();
    for (int iRun{0}; iRun<6; iRun++) {
        hRecoDijetEtaCMRun[iRun]->Write();
    }

    hRecoDijetPtEtaForward->Write();
    hRecoDijetPtEtaBackward->Write();
    hRecoDijetPtEtaForwardCMInLab->Write();
    hRecoDijetPtEtaBackwardCMInLab->Write();
    hRecoDijetPtEtaCMForward->Write();
    hRecoDijetPtEtaCMBackward->Write();
    hRecoDijetPtEtaForwardLabInCM->Write();
    hRecoDijetPtEtaBackwardLabInCM->Write();
    hRecoDijetPtEtaForwardWeighted->Write();
    hRecoDijetPtEtaBackwardWeighted->Write();
    hRecoDijetPtEtaCMForwardWeighted->Write();
    hRecoDijetPtEtaCMBackwardWeighted->Write();
    for (int i{0}; i<6; ++i) {
        hRecoDijetPtRawEtaForwardArr[i]->Write();
        hRecoDijetPtRawEtaBackwardArr[i]->Write();
        hRecoDijetPtEtaForwardArr[i]->Write();
        hRecoDijetPtEtaBackwardArr[i]->Write();
    }

    hRecoDijetPtEta->Write();
    hRecoDijetPtEtaWeighted->Write();
    hRecoDijetPtEtaCMInLab->Write();
    hRecoDijetPtEtaCM->Write();
    hRecoDijetPtEtaCMWeighted->Write();
    hRecoDijetPtEtaLabInCM->Write();
    if ( fIsMc ) {
        hRecoDijetPtEtaMatched->Write();
        hRecoDijetPtEtaCMMatched->Write();
    }

    hRecoDijetPtAveLeadPtSubLeadPt->Write();
    hRecoDijetPtAveLeadPtSubLeadPtCM->Write();
    hRecoDijetPtAveLeadEtaSubLeadEta->Write();
    hRecoDijetPtAveLeadEtaSubLeadEtaCM->Write();
    hRecoDijetEtaLeadEtaSubLeadEta->Write();
    hRecoDijetEtaLeadEtaSubLeadEtaCM->Write();

    hRecoDijetLeadPtEta->Write();
    hRecoDijetLeadPtEtaStdBins->Write();
    hRecoDijetSubLeadPtEta->Write();
    hRecoDijetSubLeadPtEtaStdBins->Write();
    for (int i{0}; i<3; ++i) {
        hRecoDijetXj[i]->Write();
        hRecoDijetXjCM[i]->Write();
    }

    hRecoLeadAllJetPtEta->Write();
    hRecoLeadAllJetPtEtaCM->Write();
    hRecoLeadAllJetPtRawEtaStdBins->Write();
    hRecoLeadAllJetPtEtaStdBins->Write();
    hRecoLeadAllJetPtEtaPtHat->Write();
    hRecoSubLeadAllJetPtEta->Write();
    hRecoSubLeadAllJetPtEtaCM->Write();
    hRecoSubLeadAllJetPtRawEtaStdBins->Write();
    hRecoSubLeadAllJetPtEtaStdBins->Write();
    hRecoSubLeadAllJetPtEtaPtHat->Write();

    hRecoGoodInclusiveJetEtaLabFrame->Write();
    hRecoGoodInclusiveJetEtaCMFrame->Write();

    hRecoDijetEta->Write();
    hRecoDijetEtaCM->Write();

    if ( fIsMc ) {

        //
        // Gen histograms
        //

        hGenJetCollectionSize->Write();
        hGenVsRecoJetCollectionSize->Write();
        hGenLeadJetPtOverPtHatVsLeadJetPt->Write();
        hGenLeadJetPtOverPtHatVsLeadJetPtWeighted->Write();
        hGenDijetPtOverPtHatVsDijetPt->Write();
        hGenDijetPtOverPtHatVsDijetPtWeighted->Write();
        hGenDijetPtAveOverPtHatVsDijetPtAve->Write();
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Write();
        // hGenDijetInfo->Write();
        hGenInclusiveJetPt->Write();
        hGenInclusiveJetEta->Write();
        hGenInclusiveJetEtaUnweighted->Write();
        hGenInclusiveJetPtEta->Write();
        hGenInclusiveJetPtEtaCM->Write();
        hGenInclusiveJetPtEtaPtHat->Write();
        hGenLeadJetPtEta->Write();
        hGenLeadJetPtEtaCM->Write();
        hGenLeadJetPtEtaPtHat->Write();
        hGenSubLeadJetPtEta->Write();
        hGenSubLeadJetPtEtaCM->Write();
        hGenSubLeadJetPtEtaPtHat->Write();
        hGenPtLeadPtSublead->Write();
        hGenEtaLeadEtaSublead->Write();
        hGenEtaCMLeadEtaCMSublead->Write();
        hGenPtLeadPtSubleadMcReweight->Write();
        hGenEtaLeadEtaSubleadMcReweight->Write();
        hGenDijetEta->Write();
        hGenDijetPtEta->Write();
        hGenDijetPtEtaWeighted->Write();
        hGenDijetPtEtaCMInLab->Write();
        hGenDijetEtaCM->Write();
        hGenDijetPtEtaCM->Write();
        hGenDijetPtEtaCMWeighted->Write();
        hGenDijetPtEtaLabInCM->Write();
        hGenDijetPtEtaForward->Write();
        hGenDijetPtEtaBackward->Write();
        hGenDijetPtEtaForwardCMInLab->Write();
        hGenDijetPtEtaBackwardCMInLab->Write();
        hGenDijetPtEtaCMForward->Write();
        hGenDijetPtEtaCMBackward->Write();
        hGenDijetPtEtaForwardLabInCM->Write();
        hGenDijetPtEtaBackwardLabInCM->Write();
        hGenDijetPtEtaForwardWeighted->Write();
        hGenDijetPtEtaBackwardWeighted->Write();
        hGenDijetPtEtaCMForwardWeighted->Write();
        hGenDijetPtEtaCMBackwardWeighted->Write();
        for (int i{0}; i<6; ++i) {
            hGenDijetPtEtaForwardArr[i]->Write();
            hGenDijetPtEtaBackwardArr[i]->Write();
        }
        hGenDijetPtAveLeadPtSubLeadPt->Write();
        hGenDijetPtAveLeadPtSubLeadPtCM->Write();
        hGenDijetPtAveLeadEtaSubLeadEta->Write();
        hGenDijetPtAveLeadEtaSubLeadEtaCM->Write();
        hGenDijetEtaLeadEtaSubLeadEta->Write();
        hGenDijetEtaLeadEtaSubLeadEtaCM->Write();
        hGenGoodInclusiveJetEtaLabFrame->Write();
        hGenGoodInclusiveJetEtaCMFrame->Write();

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
        hGenSelectedDijetXPbOverXpEta->Write();
        hGenSelectedDijetXPbOverXpEtaWeighted->Write();

        //
        // Ref hisrograms
        //

        hRecoInclusiveJetReco2Ref->Write();
        hRecoLeadJetReco2Ref->Write();
        hRecoSubLeadJetReco2Ref->Write();
        hJESInclusiveJetPtEtaPhi->Write();
        hRecoLeadJetPtOverPtHatVsLeadJetPt->Write();
        hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted->Write();
        hRecoDijetPtOverPtHatVsDijetPt->Write();
        hRecoDijetPtOverPtHatVsDijetPtWeighted->Write();
        hRecoDijetPtAveOverPtHatVsDijetPtAve->Write();
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Write();

        // hRecoInclusiveJetJECFactorVsPtEta->Write();
        // hRecoInclusiveJetJEC2FactorVsPtEta->Write();
        // hRecoInclusiveJetPtRawOverPtRefVsPtEta->Write();
        // hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->Write();
        // hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning->Write();

        hInclusiveJetJESVsPtGen->Write();
        hInclusiveJetJESGenPtGenEtaPtHatWeighted->Write();
        hInclusiveJetJESGenPtGenEtaCMPtHatWeighted->Write();
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Write();
        hLeadJetJESGenPtEtaPtHatWeighted->Write();
        hLeadJetJESGenPtEtaCMPtHatWeighted->Write();
        hSubLeadJetJESGenPtEtaPtHatWeighted->Write();
        hSubLeadJetJESGenPtEtaCMPtHatWeighted->Write();

        hRecoInclusiveMatchedJetPt->Write();
        hRecoInclusiveMatchedJetPtEta->Write();
        hRecoInclusiveUnmatchedJetPtEta->Write();
        hRecoLeadMatchedJetPtEta->Write();
        hRecoLeadMatchedJetPtEtaPtHat->Write();
        hRecoLeadUnmatchedJetPtEta->Write();
        hRecoLeadUnmatchedJetPtEtaPtHat->Write();
        hRecoSubLeadMatchedJetPtEta->Write();
        hRecoSubLeadMatchedJetPtEtaPtHat->Write();
        hRecoSubLeadUnmatchedJetPtEta->Write();
        hRecoSubLeadUnmatchedJetPtEtaPtHat->Write();

        hRefInclusiveJetPt->Write();
        hRefInclusiveJetEta->Write();
        hRefInclusiveJetEtaUnweighted->Write();
        hRefInclusiveJetPtEta->Write();
        hRefInclusiveJetPtEtaCM->Write();
        hRefInclusiveJetPtEtaPtHat->Write();
        hRefLeadJetPtEta->Write();
        hRefLeadJetPtEtaCM->Write();
        hRefLeadJetPtEtaPtHat->Write();
        hRefLeadUnswappedJetPtEta->Write();
        hRefLeadUnswappedJetPtEtaPtHat->Write();
        hRefSubLeadJetPtEta->Write();
        hRefSubLeadJetPtEtaCM->Write();
        hRefSubLeadJetPtEtaPtHat->Write();
        hRefSubLeadUnswappedJetPtEta->Write();
        hRefSubLeadUnswappedJetPtEtaPtHat->Write();

        // hReco2RefFull->Write();

        hRecoDijetPtEtaRefDijetPtEta->Write();
        hRecoDijetPtEtaRefDijetPtEtaWeighted->Write();

        // hRefSel2RecoFull->Write();

        hRefDijetEta->Write();
        hRefDijetEtaVsRecoDijetEta->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Write();
        hRefDijetPtEta->Write();
        hRefDijetPtEtaWeighted->Write();
        hRefDijetPtEtaCMInLab->Write();

        hRefDijetPtEtaForward->Write();
        hRefDijetPtEtaBackward->Write();
        hRefDijetPtEtaForwardCMInLab->Write();
        hRefDijetPtEtaBackwardCMInLab->Write();
        hRefDijetPtEtaCMForward->Write();
        hRefDijetPtEtaCMBackward->Write();
        hRefDijetPtEtaForwardLabInCM->Write();
        hRefDijetPtEtaBackwardLabInCM->Write();
        hRefDijetPtEtaForwardWeighted->Write();
        hRefDijetPtEtaBackwardWeighted->Write();
        hRefDijetPtEtaCMForwardWeighted->Write();
        hRefDijetPtEtaCMBackwardWeighted->Write();

        hRefDijetEtaCM->Write();
        hRefDijetPtEtaCM->Write();
        hRefDijetPtEtaCMWeighted->Write();
        hRefDijetPtEtaLabInCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Write();

        hRefPtLeadPtSublead->Write();
        hRefEtaLeadEtaSublead->Write();
        hRefEtaCMLeadEtaCMSublead->Write();
        hRefPtLeadPtSubleadMcReweight->Write();
        hRefEtaLeadEtaSubleadMcReweight->Write();

        //
        // Ref-selected histograms
        //

        hRefSelInclusiveJetPt->Write();
        hRefSelInclusiveJetEta->Write();
        hRefSelInclusiveJetEtaUnweighted->Write();
        hRefSelInclusiveJetPtEta->Write();
        hRefSelInclusiveJetPtEtaPtHat->Write();
        hRefSelLeadJetPtEtaPtHat->Write();
        hRefSelSubLeadJetPtEtaPtHat->Write();

        hRefSelDijetEta->Write();
        hRefSelDijetPtEta->Write();
        hRefSelDijetPtEtaWeighted->Write();
        hRefSelDijetPtEtaCMInLab->Write();
        hRefSelDijetEtaCM->Write();
        hRefSelDijetPtEtaCM->Write();
        hRefSelDijetPtEtaCMWeighted->Write();
        hRefSelDijetPtEtaLabInCM->Write();
    } // if ( fIsMc )

}
