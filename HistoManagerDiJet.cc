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

    hGenDijetInfo{nullptr},
    hGenInclusiveJetPt{nullptr},
    hGenInclusiveJetEta{nullptr},
    hGenInclusiveJetEtaUnweighted{nullptr},
    hGenInclusiveJetPtEta{nullptr},
    hGenInclusiveJetPtEtaPtHat{nullptr},
    hGenLeadJetPtEta{nullptr},
    hGenLeadJetPtEtaPtHat{nullptr},
    hGenSubLeadJetPtEta{nullptr},
    hGenSubLeadJetPtEtaPtHat{nullptr},
    hGenPtLeadPtSublead{nullptr},
    hGenEtaLeadEtaSublead{nullptr},
    hGenEtaCMLeadEtaCMSublead{nullptr},
    hGenPtLeadPtSubleadMcReweight{nullptr},
    hGenEtaLeadEtaSubleadMcReweight{nullptr},

    hGenDijetEta{nullptr},
    hGenDijetPtEtaPhi{nullptr},
    hGenDijetPtEtaPhiWeighted{nullptr},
    hGenDijetPtEtaPhiCMInLab{nullptr},
    hGenDijetEtaCM{nullptr},
    hGenDijetPtEtaPhiCM{nullptr},
    hGenDijetPtEtaPhiCMWeighted{nullptr},
    hGenDijetPtEtaPhiLabInCM{nullptr},
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

    hGenDijetEta1D{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2D{nullptr},
    hGenDijetEtaForward1D{nullptr},
    hGenDijetEtaBackward1D{nullptr},

    hGenDijetEta1DCM{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DCM{nullptr},
    hGenDijetEtaCMForward1D{nullptr},
    hGenDijetEtaCMBackward1D{nullptr},

    hGenDijetEta1DOldPt{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hGenDijetEtaForward1DOldPt{nullptr},
    hGenDijetEtaBackward1DOldPt{nullptr},

    hGenDijetEta1DOldPtCM{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hGenDijetEtaCMForward1DOldPt{nullptr},
    hGenDijetEtaCMBackward1DOldPt{nullptr},

    hGenDijetEta1DOldPtBinning{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hGenDijetEtaForward1DOldPtBinning{nullptr},
    hGenDijetEtaBackward1DOldPtBinning{nullptr},

    hGenDijetEta1DOldPtBinningCM{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hGenDijetEtaCMForward1DOldPtBinning{nullptr},
    hGenDijetEtaCMBackward1DOldPtBinning{nullptr},

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
    hRecoDijetInfo{nullptr},

    hRecoInclusiveAllJetPt{nullptr},
    hRecoInclusiveAllJetEta{nullptr},
    hRecoInclusiveAllJetEtaUnweighted{nullptr},
    hRecoInclusiveAllJetPtEta{nullptr},
    hRecoInclusiveAllJetPtEtaPtHat{nullptr},
    hRecoInclusiveMatchedJetPtEtaPtHat{nullptr},
    hRecoInclusiveUnmatchedJetPtEtaPtHat{nullptr},

    hRecoPtLeadPtSublead{nullptr},
    hRecoEtaLeadEtaSublead{nullptr},
    hRecoEtaCMLeadEtaCMSublead{nullptr},
    hRecoPtLeadPtSubleadMcReweight{nullptr},
    hRecoEtaLeadEtaSubleadMcReweight{nullptr},

    hRecoDijetPtEta{nullptr},
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
    hRecoDijetPtEtaPhi{nullptr},
    hRecoDijetPtEtaPhiWeighted{nullptr},
    hRecoDijetPtEtaPhiCMInLab{nullptr},
    hRecoDijetPtEtaPhiCM{nullptr},
    hRecoDijetPtEtaPhiCMWeighted{nullptr},
    hRecoDijetPtEtaPhiLabInCM{nullptr},
    hRecoDijetPtEtaPhiMatched{nullptr},
    hRecoDijetPtEtaPhiCMMatched{nullptr},

    hRecoLeadAllJetPtEta{nullptr},
    hRecoLeadAllJetPtEtaPtHat{nullptr},
    hRecoSubLeadAllJetPtEta{nullptr},
    hRecoSubLeadAllJetPtEtaPtHat{nullptr},
    hRecoGoodInclusiveJetEtaLabFrame{nullptr},
    hRecoGoodInclusiveJetEtaCMFrame{nullptr},
    hRecoDijetEta{nullptr},
    hRecoDijetEtaCM{nullptr},

    hRecoDijetEta1D{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2D{nullptr},
    hRecoDijetEtaForward1D{nullptr},
    hRecoDijetEtaBackward1D{nullptr},

    hRecoDijetEta1DCM{nullptr},
    hRecoEtaLeadVsEtaSubLead2DCM{nullptr},
    hRecoDijetEtaCMForward1D{nullptr},
    hRecoDijetEtaCMBackward1D{nullptr},

    hRecoDijetEta1DOldPt{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hRecoDijetEtaForward1DOldPt{nullptr},
    hRecoDijetEtaBackward1DOldPt{nullptr},

    hRecoDijetEta1DOldPtCM{nullptr},
    hRecoEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hRecoDijetEtaCMForward1DOldPt{nullptr},
    hRecoDijetEtaCMBackward1DOldPt{nullptr},

    hRecoDijetEta1DOldPtBinning{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hRecoDijetEtaForward1DOldPtBinning{nullptr},
    hRecoDijetEtaBackward1DOldPtBinning{nullptr},

    hRecoDijetEta1DOldPtBinningCM{nullptr},
    hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hRecoDijetEtaCMForward1DOldPtBinning{nullptr},
    hRecoDijetEtaCMBackward1DOldPtBinning{nullptr},

    //
    // Ref jet histograms
    //

    hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
    hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
    hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
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
    hInclusiveJetJESRecoPtRecoEtaPtHatWeighted{nullptr},

    hLeadJetJESGenPtEtaPtHatWeighted{nullptr},
    hSubLeadJetJESGenPtEtaPtHatWeighted{nullptr},

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
    hRefInclusiveJetPtEtaPtHat{nullptr},
    hRefLeadJetPtEta{nullptr},
    hRefLeadJetPtEtaPtHat{nullptr},
    hRefLeadUnswappedJetPtEta{nullptr},
    hRefLeadUnswappedJetPtEtaPtHat{nullptr},
    hRefSubLeadJetPtEta{nullptr},
    hRefSubLeadJetPtEtaPtHat{nullptr},
    hRefSubLeadUnswappedJetPtEta{nullptr},
    hRefSubLeadUnswappedJetPtEtaPtHat{nullptr},

    hReco2RefFull{nullptr},
    hRecoDijetPtEtaRefDijetPtEta{nullptr},
    hRecoDijetPtEtaRefDijetPtEtaWeighted{nullptr},
    
    hRefSel2RecoFull{nullptr},

    hRefDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted{nullptr},
    hRefDijetPtEtaPhi{nullptr},
    hRefDijetPtEtaPhiWeighted{nullptr},
    hRefDijetPtEtaPhiCMInLab{nullptr},

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
    hRefDijetPtEtaPhiCM{nullptr},
    hRefDijetPtEtaPhiCMWeighted{nullptr},
    hRefDijetPtEtaPhiLabInCM{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted{nullptr},

    hRefPtLeadPtSublead{nullptr},
    hRefEtaLeadEtaSublead{nullptr},
    hRefEtaCMLeadEtaCMSublead{nullptr},
    hRefPtLeadPtSubleadMcReweight{nullptr},
    hRefEtaLeadEtaSubleadMcReweight{nullptr},

    hRefDijetEta1D{nullptr},
    hRefEtaLeadVsEtaSubLead2D{nullptr},
    hRecoVsRefDijetEta2D{nullptr},
    hRecoVsRefLeadJetEta2D{nullptr},
    hRecoVsRefSubLeadJetEta2D{nullptr},
    hRefDijetEtaForward1D{nullptr},
    hRefDijetEtaBackward1D{nullptr},

    hRefDijetEta1DCM{nullptr},
    hRefEtaLeadVsEtaSubLead2DCM{nullptr},
    hRecoVsRefDijetEta2DCM{nullptr},
    hRecoVsRefLeadJetEta2DCM{nullptr},
    hRecoVsRefSubLeadJetEta2DCM{nullptr},
    hRefDijetEtaCMForward1D{nullptr},
    hRefDijetEtaCMBackward1D{nullptr},

    hRefDijetEta1DOldPt{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hRecoVsRefDijetEta2DOldPt{nullptr},
    hRecoVsRefLeadJetEta2DOldPt{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPt{nullptr},
    hRefDijetEtaForward1DOldPt{nullptr},
    hRefDijetEtaBackward1DOldPt{nullptr},

    hRefDijetEta1DOldPtCM{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hRecoVsRefDijetEta2DOldPtCM{nullptr},
    hRecoVsRefLeadJetEta2DOldPtCM{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtCM{nullptr},
    hRefDijetEtaCMForward1DOldPt{nullptr},
    hRefDijetEtaCMBackward1DOldPt{nullptr},

    hRefDijetEta1DOldPtBinning{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hRecoVsRefDijetEta2DOldPtBinning{nullptr},
    hRecoVsRefLeadJetEta2DOldPtBinning{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtBinning{nullptr},
    hRefDijetEtaForward1DOldPtBinning{nullptr},
    hRefDijetEtaBackward1DOldPtBinning{nullptr},

    hRefDijetEta1DOldPtBinningCM{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hRecoVsRefDijetEta2DOldPtBinningCM{nullptr},
    hRecoVsRefLeadJetEta2DOldPtBinningCM{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtBinningCM{nullptr},
    hRefDijetEtaCMForward1DOldPtBinning{nullptr},
    hRefDijetEtaCMBackward1DOldPtBinning{nullptr},

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
    hRefSelDijetPtEtaPhi{nullptr},
    hRefSelDijetPtEtaPhiWeighted{nullptr},
    hRefSelDijetPtEtaPhiCMInLab{nullptr},
    hRefSelDijetEtaCM{nullptr},
    hRefSelDijetPtEtaPhiCM{nullptr},
    hRefSelDijetPtEtaPhiCMWeighted{nullptr},
    hRefSelDijetPtEtaPhiLabInCM{nullptr},

    hRefSelDijetEta1D{nullptr},
    hRefSelRecoDijetEta1D{nullptr},
    hRefSelEtaLeadVsEtaSubLead2D{nullptr},
    hRefSelDijetEtaForward1D{nullptr},
    hRefSelDijetEtaBackward1D{nullptr},

    hRefSelDijetEta1DCM{nullptr},
    hRefSelRecoDijetEta1DCM{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DCM{nullptr},
    hRefSelDijetEtaCMForward1D{nullptr},
    hRefSelDijetEtaCMBackward1D{nullptr},

    hRefSelDijetEta1DOldPt{nullptr},
    hRefSelRecoDijetEta1DOldPt{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hRefSelDijetEtaForward1DOldPt{nullptr},
    hRefSelDijetEtaBackward1DOldPt{nullptr},

    hRefSelDijetEta1DOldPtCM{nullptr},
    hRefSelRecoDijetEta1DOldPtCM{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hRefSelDijetEtaCMForward1DOldPt{nullptr},
    hRefSelDijetEtaCMBackward1DOldPt{nullptr},

    hRefSelDijetEta1DOldPtBinning{nullptr},
    hRefSelRecoDijetEta1DOldPtBinning{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hRefSelDijetEtaForward1DOldPtBinning{nullptr},
    hRefSelDijetEtaBackward1DOldPtBinning{nullptr},

    hRefSelDijetEta1DOldPtBinningCM{nullptr},
    hRefSelRecoDijetEta1DOldPtBinningCM{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hRefSelDijetEtaCMForward1DOldPtBinning{nullptr},
    hRefSelDijetEtaCMBackward1DOldPtBinning{nullptr},

    //
    // Variables
    //
    fIsMc{false}, fUseVariableBinning{true},
    fPtBins{150}, fPtRange{5., 1505.}, 
    fEtaBins{52}, fEtaRange{-5.2, 5.2},
    fPhiBins{16}, fPhiRange{-TMath::Pi(), TMath::Pi()},
    fDijetPtBins{196}, fDijetPtRange{20., 1000.},
    fDijetEtaBins{52}, fDijetEtaRange{-5.2, 5.2},
    fDijetDphiBins{8}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},
    fPtHatBins{100}, fPtHatRange{15., 1015.},
    fFracBins{100}, fFracRange{0., 1.},
    fMultBins{32}, fMultRange{-0.5, 31.5} { 

    double dijetPtVals[17] {  50.,  60.,   70.,  80.,  90.,
                              100., 110.,  120., 130., 140.,
                              150., 160.,  180., 200., 250., 
                              300., 500.};
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    fPtAveBins.assign(dijetPtVals, dijetPtVals + sizeOfPtVals);

    double dijetPtOldVals[7] {25., 55., 75., 95., 115., 150., 400.};
    int sizeOfPtOldVals = sizeof(dijetPtOldVals)/sizeof(dijetPtOldVals[0]);
    fPtAveOldBins.assign(dijetPtOldVals, dijetPtOldVals + sizeOfPtOldVals);
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
        if (hGenDijetInfo) delete hGenDijetInfo;
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
        if (hGenDijetPtEtaPhi) delete hGenDijetPtEtaPhi;
        if (hGenDijetPtEtaPhiWeighted) delete hGenDijetPtEtaPhiWeighted;
        if (hGenDijetPtEtaPhiCMInLab) delete hGenDijetPtEtaPhiCMInLab;
        if (hGenDijetEtaCM) delete hGenDijetEtaCM;
        if (hGenDijetPtEtaPhiCM) delete hGenDijetPtEtaPhiCM;
        if (hGenDijetPtEtaPhiCMWeighted) delete hGenDijetPtEtaPhiCMWeighted;
        if (hGenDijetPtEtaPhiLabInCM) delete hGenDijetPtEtaPhiLabInCM;
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

        for (int i = 0; i < 16; ++i) {
            if (hGenDijetEta1D[i]) { delete hGenDijetEta1D[i]; hGenDijetEta1D[i] = nullptr; }
            if (hGenDijetEtaLeadVsEtaSubLead2D[i]) { delete hGenDijetEtaLeadVsEtaSubLead2D[i]; hGenDijetEtaLeadVsEtaSubLead2D[i] = nullptr; }
            if (hGenDijetEtaForward1D[i]) { delete hGenDijetEtaForward1D[i]; hGenDijetEtaForward1D[i] = nullptr; }
            if (hGenDijetEtaBackward1D[i]) { delete hGenDijetEtaBackward1D[i]; hGenDijetEtaBackward1D[i] = nullptr; }

            if (hGenDijetEta1DCM[i]) { delete hGenDijetEta1DCM[i]; hGenDijetEta1DCM[i] = nullptr; }
            if (hGenDijetEtaLeadVsEtaSubLead2DCM[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DCM[i]; hGenDijetEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
            if (hGenDijetEtaCMForward1D[i]) { delete hGenDijetEtaCMForward1D[i]; hGenDijetEtaCMForward1D[i] = nullptr; }
            if (hGenDijetEtaCMBackward1D[i]) { delete hGenDijetEtaCMBackward1D[i]; hGenDijetEtaCMBackward1D[i] = nullptr; }
        }

        for (int i = 0; i < 6; ++i) {
            if (hGenDijetEta1DOldPt[i]) { delete hGenDijetEta1DOldPt[i]; hGenDijetEta1DOldPt[i] = nullptr; }
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
            if (hGenDijetEtaForward1DOldPt[i]) { delete hGenDijetEtaForward1DOldPt[i]; hGenDijetEtaForward1DOldPt[i] = nullptr; }
            if (hGenDijetEtaBackward1DOldPt[i]) { delete hGenDijetEtaBackward1DOldPt[i]; hGenDijetEtaBackward1DOldPt[i] = nullptr; }

            if (hGenDijetEta1DOldPtCM[i]) { delete hGenDijetEta1DOldPtCM[i]; hGenDijetEta1DOldPtCM[i] = nullptr; }
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
            if (hGenDijetEtaCMForward1DOldPt[i]) { delete hGenDijetEtaCMForward1DOldPt[i]; hGenDijetEtaCMForward1DOldPt[i] = nullptr; }
            if (hGenDijetEtaCMBackward1DOldPt[i]) { delete hGenDijetEtaCMBackward1DOldPt[i]; hGenDijetEtaCMBackward1DOldPt[i] = nullptr; }

            if (hGenDijetEta1DOldPtBinning[i]) { delete hGenDijetEta1DOldPtBinning[i]; hGenDijetEta1DOldPtBinning[i] = nullptr; }
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
            if (hGenDijetEtaForward1DOldPtBinning[i]) { delete hGenDijetEtaForward1DOldPtBinning[i]; hGenDijetEtaForward1DOldPtBinning[i] = nullptr; }
            if (hGenDijetEtaBackward1DOldPtBinning[i]) { delete hGenDijetEtaBackward1DOldPtBinning[i]; hGenDijetEtaBackward1DOldPtBinning[i] = nullptr; }

            if (hGenDijetEta1DOldPtBinningCM[i]) { delete hGenDijetEta1DOldPtBinningCM[i]; hGenDijetEta1DOldPtBinningCM[i] = nullptr; }
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
            if (hGenDijetEtaCMForward1DOldPtBinning[i]) { delete hGenDijetEtaCMForward1DOldPtBinning[i]; hGenDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
            if (hGenDijetEtaCMBackward1DOldPtBinning[i]) { delete hGenDijetEtaCMBackward1DOldPtBinning[i]; hGenDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
        }
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
    if (hRecoDijetInfo) delete hRecoDijetInfo;
    if (hRecoInclusiveAllJetPt) delete hRecoInclusiveAllJetPt;
    if (hRecoInclusiveAllJetEta) delete hRecoInclusiveAllJetEta;
    if (hRecoInclusiveAllJetEtaUnweighted) delete hRecoInclusiveAllJetEtaUnweighted;
    if (hRecoInclusiveAllJetPtEta) delete hRecoInclusiveAllJetPtEta;
    if (hRecoInclusiveAllJetPtEtaPtHat) delete hRecoInclusiveAllJetPtEtaPtHat;
    if (hRecoInclusiveMatchedJetPtEtaPtHat) delete hRecoInclusiveMatchedJetPtEtaPtHat;
    if (hRecoInclusiveUnmatchedJetPtEtaPtHat) delete hRecoInclusiveUnmatchedJetPtEtaPtHat;
    if (hRecoPtLeadPtSublead) delete hRecoPtLeadPtSublead;
    if (hRecoEtaLeadEtaSublead) delete hRecoEtaLeadEtaSublead;
    if (hRecoEtaCMLeadEtaCMSublead) delete hRecoEtaCMLeadEtaCMSublead;
    if (hRecoPtLeadPtSubleadMcReweight) delete hRecoPtLeadPtSubleadMcReweight;
    if (hRecoEtaLeadEtaSubleadMcReweight) delete hRecoEtaLeadEtaSubleadMcReweight;
    if (hRecoDijetPtEta) delete hRecoDijetPtEta;
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
    if (hRecoDijetPtEtaPhi) delete hRecoDijetPtEtaPhi;
    if (hRecoDijetPtEtaPhiWeighted) delete hRecoDijetPtEtaPhiWeighted;
    if (hRecoDijetPtEtaPhiCMInLab) delete hRecoDijetPtEtaPhiCMInLab;
    if (hRecoDijetPtEtaPhiCM) delete hRecoDijetPtEtaPhiCM;
    if (hRecoDijetPtEtaPhiCMWeighted) delete hRecoDijetPtEtaPhiCMWeighted;
    if (hRecoDijetPtEtaPhiLabInCM) delete hRecoDijetPtEtaPhiLabInCM;
    if (hRecoDijetPtEtaPhiMatched) delete hRecoDijetPtEtaPhiMatched;
    if (hRecoDijetPtEtaPhiCMMatched) delete hRecoDijetPtEtaPhiCMMatched;
    if (hRecoLeadAllJetPtEta) delete hRecoLeadAllJetPtEta;
    if (hRecoLeadAllJetPtEtaPtHat) delete hRecoLeadAllJetPtEtaPtHat;
    if (hRecoSubLeadAllJetPtEta) delete hRecoSubLeadAllJetPtEta;
    if (hRecoSubLeadAllJetPtEtaPtHat) delete hRecoSubLeadAllJetPtEtaPtHat;
    if (hRecoGoodInclusiveJetEtaLabFrame) delete hRecoGoodInclusiveJetEtaLabFrame;
    if (hRecoGoodInclusiveJetEtaCMFrame) delete hRecoGoodInclusiveJetEtaCMFrame;
    if (hRecoDijetEta) delete hRecoDijetEta;
    if (hRecoDijetEtaCM) delete hRecoDijetEtaCM;

    // New ptAve and eta binning
    for (int i = 0; i < 16; ++i) {
        if (hRecoDijetEta1D[i]) { delete hRecoDijetEta1D[i]; hRecoDijetEta1D[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2D[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2D[i]; hRecoDijetEtaLeadVsEtaSubLead2D[i] = nullptr; }
        if (hRecoDijetEtaForward1D[i]) { delete hRecoDijetEtaForward1D[i]; hRecoDijetEtaForward1D[i] = nullptr; }
        if (hRecoDijetEtaBackward1D[i]) { delete hRecoDijetEtaBackward1D[i]; hRecoDijetEtaBackward1D[i] = nullptr; }

        if (hRecoDijetEta1DCM[i]) { delete hRecoDijetEta1DCM[i]; hRecoDijetEta1DCM[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DCM[i]) { delete hRecoEtaLeadVsEtaSubLead2DCM[i]; hRecoEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
        if (hRecoDijetEtaCMForward1D[i]) { delete hRecoDijetEtaCMForward1D[i]; hRecoDijetEtaCMForward1D[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1D[i]) { delete hRecoDijetEtaCMBackward1D[i]; hRecoDijetEtaCMBackward1D[i] = nullptr; }
    }

    // Old ptAve and new eta binning
    for (int i = 0; i < 6; ++i) {
        if (hRecoDijetEta1DOldPt[i]) { delete hRecoDijetEta1DOldPt[i]; hRecoDijetEta1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]; hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
        if (hRecoDijetEtaForward1DOldPt[i]) { delete hRecoDijetEtaForward1DOldPt[i]; hRecoDijetEtaForward1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaBackward1DOldPt[i]) { delete hRecoDijetEtaBackward1DOldPt[i]; hRecoDijetEtaBackward1DOldPt[i] = nullptr; }

        if (hRecoDijetEta1DOldPtCM[i]) { delete hRecoDijetEta1DOldPtCM[i]; hRecoDijetEta1DOldPtCM[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]; hRecoEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DOldPt[i]) { delete hRecoDijetEtaCMForward1DOldPt[i]; hRecoDijetEtaCMForward1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DOldPt[i]) { delete hRecoDijetEtaCMBackward1DOldPt[i]; hRecoDijetEtaCMBackward1DOldPt[i] = nullptr; }
    }

    // Old ptAve and old eta binning
    for (int i = 0; i < 6; ++i) {
        if (hRecoDijetEta1DOldPtBinning[i]) { delete hRecoDijetEta1DOldPtBinning[i]; hRecoDijetEta1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]; hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaForward1DOldPtBinning[i]) { delete hRecoDijetEtaForward1DOldPtBinning[i]; hRecoDijetEtaForward1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaBackward1DOldPtBinning[i]) { delete hRecoDijetEtaBackward1DOldPtBinning[i]; hRecoDijetEtaBackward1DOldPtBinning[i] = nullptr; }

        if (hRecoDijetEta1DOldPtBinningCM[i]) { delete hRecoDijetEta1DOldPtBinningCM[i]; hRecoDijetEta1DOldPtBinningCM[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DOldPtBinning[i]) { delete hRecoDijetEtaCMForward1DOldPtBinning[i]; hRecoDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DOldPtBinning[i]) { delete hRecoDijetEtaCMBackward1DOldPtBinning[i]; hRecoDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
    }

    if ( fIsMc ) {

        //
        // Ref jet histograms
        //
        if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen;
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
        if ( hInclusiveJetJESRecoPtRecoEtaPtHatWeighted ) delete hInclusiveJetJESRecoPtRecoEtaPtHatWeighted;
        if ( hLeadJetJESGenPtEtaPtHatWeighted ) delete hLeadJetJESGenPtEtaPtHatWeighted;
        if ( hSubLeadJetJESGenPtEtaPtHatWeighted ) delete hSubLeadJetJESGenPtEtaPtHatWeighted;

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

        if (hReco2RefFull) delete hReco2RefFull;
        if (hRecoDijetPtEtaRefDijetPtEta) delete hRecoDijetPtEtaRefDijetPtEta;
        if (hRecoDijetPtEtaRefDijetPtEtaWeighted) delete hRecoDijetPtEtaRefDijetPtEtaWeighted;
        if (hRefSel2RecoFull) delete hRefSel2RecoFull;

        if (hRefDijetEta) delete hRefDijetEta;
        if (hRefDijetEtaVsRecoDijetEta) delete hRefDijetEtaVsRecoDijetEta;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted;
        if (hRefDijetPtEtaPhi) delete hRefDijetPtEtaPhi;
        if (hRefDijetPtEtaPhiWeighted) delete hRefDijetPtEtaPhiWeighted;
        if (hRefDijetPtEtaPhiCMInLab) delete hRefDijetPtEtaPhiCMInLab;

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
        if (hRefDijetPtEtaPhiCM) delete hRefDijetPtEtaPhiCM;
        if (hRefDijetPtEtaPhiCMWeighted) delete hRefDijetPtEtaPhiCMWeighted;
        if (hRefDijetPtEtaPhiLabInCM) delete hRefDijetPtEtaPhiLabInCM;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted;

        if (hRefPtLeadPtSublead) delete hRefPtLeadPtSublead;
        if (hRefEtaLeadEtaSublead) delete hRefEtaLeadEtaSublead;
        if (hRefEtaCMLeadEtaCMSublead) delete hRefEtaCMLeadEtaCMSublead;
        if (hRefPtLeadPtSubleadMcReweight) delete hRefPtLeadPtSubleadMcReweight;
        if (hRefEtaLeadEtaSubleadMcReweight) delete hRefEtaLeadEtaSubleadMcReweight;

        // New ptAve and eta binning
        for (int i = 0; i < 16; ++i) {
            if (hRefDijetEta1D[i]) { delete hRefDijetEta1D[i]; hRefDijetEta1D[i] = nullptr; }
            if (hRefEtaLeadVsEtaSubLead2D[i]) { delete hRefEtaLeadVsEtaSubLead2D[i]; hRefEtaLeadVsEtaSubLead2D[i] = nullptr; }
            if (hRecoVsRefDijetEta2D[i]) { delete hRecoVsRefDijetEta2D[i]; hRecoVsRefDijetEta2D[i] = nullptr; }
            if (hRecoVsRefLeadJetEta2D[i]) { delete hRecoVsRefLeadJetEta2D[i]; hRecoVsRefLeadJetEta2D[i] = nullptr; }
            if (hRecoVsRefSubLeadJetEta2D[i]) { delete hRecoVsRefSubLeadJetEta2D[i]; hRecoVsRefSubLeadJetEta2D[i] = nullptr; }
            if (hRefDijetEtaForward1D[i]) { delete hRefDijetEtaForward1D[i]; hRefDijetEtaForward1D[i] = nullptr; }
            if (hRefDijetEtaBackward1D[i]) { delete hRefDijetEtaBackward1D[i]; hRefDijetEtaBackward1D[i] = nullptr; }

            if (hRefDijetEta1DCM[i]) { delete hRefDijetEta1DCM[i]; hRefDijetEta1DCM[i] = nullptr; }
            if (hRefEtaLeadVsEtaSubLead2DCM[i]) { delete hRefEtaLeadVsEtaSubLead2DCM[i]; hRefEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
            if (hRecoVsRefDijetEta2DCM[i]) { delete hRecoVsRefDijetEta2DCM[i]; hRecoVsRefDijetEta2DCM[i] = nullptr; }
            if (hRecoVsRefLeadJetEta2DCM[i]) { delete hRecoVsRefLeadJetEta2DCM[i]; hRecoVsRefLeadJetEta2DCM[i] = nullptr; }
            if (hRecoVsRefSubLeadJetEta2DCM[i]) { delete hRecoVsRefSubLeadJetEta2DCM[i]; hRecoVsRefSubLeadJetEta2DCM[i] = nullptr; }
            if (hRefDijetEtaCMForward1D[i]) { delete hRefDijetEtaCMForward1D[i]; hRefDijetEtaCMForward1D[i] = nullptr; }
            if (hRefDijetEtaCMBackward1D[i]) { delete hRefDijetEtaCMBackward1D[i]; hRefDijetEtaCMBackward1D[i] = nullptr; }
        }

        // Old ptAve and new eta binning
        for (int i = 0; i < 6; ++i) {
            if (hRefDijetEta1DOldPt[i]) { delete hRefDijetEta1DOldPt[i]; hRefDijetEta1DOldPt[i] = nullptr; }
            if (hRefEtaLeadVsEtaSubLead2DOldPt[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPt[i]; hRefEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
            if (hRecoVsRefDijetEta2DOldPt[i]) { delete hRecoVsRefDijetEta2DOldPt[i]; hRecoVsRefDijetEta2DOldPt[i] = nullptr; }
            if (hRecoVsRefLeadJetEta2DOldPt[i]) { delete hRecoVsRefLeadJetEta2DOldPt[i]; hRecoVsRefLeadJetEta2DOldPt[i] = nullptr; }
            if (hRecoVsRefSubLeadJetEta2DOldPt[i]) { delete hRecoVsRefSubLeadJetEta2DOldPt[i]; hRecoVsRefSubLeadJetEta2DOldPt[i] = nullptr; }
            if (hRefDijetEtaForward1DOldPt[i]) { delete hRefDijetEtaForward1DOldPt[i]; hRefDijetEtaForward1DOldPt[i] = nullptr; }
            if (hRefDijetEtaBackward1DOldPt[i]) { delete hRefDijetEtaBackward1DOldPt[i]; hRefDijetEtaBackward1DOldPt[i] = nullptr; }

            if (hRefDijetEta1DOldPtCM[i]) { delete hRefDijetEta1DOldPtCM[i]; hRefDijetEta1DOldPtCM[i] = nullptr; }
            if (hRefEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPtCM[i]; hRefEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
            if (hRecoVsRefDijetEta2DOldPtCM[i]) { delete hRecoVsRefDijetEta2DOldPtCM[i]; hRecoVsRefDijetEta2DOldPtCM[i] = nullptr; }
            if (hRecoVsRefLeadJetEta2DOldPtCM[i]) { delete hRecoVsRefLeadJetEta2DOldPtCM[i]; hRecoVsRefLeadJetEta2DOldPtCM[i] = nullptr; }
            if (hRecoVsRefSubLeadJetEta2DOldPtCM[i]) { delete hRecoVsRefSubLeadJetEta2DOldPtCM[i]; hRecoVsRefSubLeadJetEta2DOldPtCM[i] = nullptr; }
            if (hRefDijetEtaCMForward1DOldPt[i]) { delete hRefDijetEtaCMForward1DOldPt[i]; hRefDijetEtaCMForward1DOldPt[i] = nullptr; }
            if (hRefDijetEtaCMBackward1DOldPt[i]) { delete hRefDijetEtaCMBackward1DOldPt[i]; hRefDijetEtaCMBackward1DOldPt[i] = nullptr; }
        }

        // Old ptAve and old eta binning
        for (int i = 0; i < 6; ++i) {
            if (hRefDijetEta1DOldPtBinning[i]) { delete hRefDijetEta1DOldPtBinning[i]; hRefDijetEta1DOldPtBinning[i] = nullptr; }
            if (hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]; hRefEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
            if (hRecoVsRefDijetEta2DOldPtBinning[i]) { delete hRecoVsRefDijetEta2DOldPtBinning[i]; hRecoVsRefDijetEta2DOldPtBinning[i] = nullptr; }
            if (hRecoVsRefLeadJetEta2DOldPtBinning[i]) { delete hRecoVsRefLeadJetEta2DOldPtBinning[i]; hRecoVsRefLeadJetEta2DOldPtBinning[i] = nullptr; }
            if (hRecoVsRefSubLeadJetEta2DOldPtBinning[i]) { delete hRecoVsRefSubLeadJetEta2DOldPtBinning[i]; hRecoVsRefSubLeadJetEta2DOldPtBinning[i] = nullptr; }
            if (hRefDijetEtaForward1DOldPtBinning[i]) { delete hRefDijetEtaForward1DOldPtBinning[i]; hRefDijetEtaForward1DOldPtBinning[i] = nullptr; }
            if (hRefDijetEtaBackward1DOldPtBinning[i]) { delete hRefDijetEtaBackward1DOldPtBinning[i]; hRefDijetEtaBackward1DOldPtBinning[i] = nullptr; }

            if (hRefDijetEta1DOldPtBinningCM[i]) { delete hRefDijetEta1DOldPtBinningCM[i]; hRefDijetEta1DOldPtBinningCM[i] = nullptr; }
            if (hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
            if (hRecoVsRefDijetEta2DOldPtBinningCM[i]) { delete hRecoVsRefDijetEta2DOldPtBinningCM[i]; hRecoVsRefDijetEta2DOldPtBinningCM[i] = nullptr; }
            if (hRecoVsRefLeadJetEta2DOldPtBinningCM[i]) { delete hRecoVsRefLeadJetEta2DOldPtBinningCM[i]; hRecoVsRefLeadJetEta2DOldPtBinningCM[i] = nullptr; }
            if (hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]) { delete hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]; hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i] = nullptr; }
            if (hRefDijetEtaCMForward1DOldPtBinning[i]) { delete hRefDijetEtaCMForward1DOldPtBinning[i]; hRefDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
            if (hRefDijetEtaCMBackward1DOldPtBinning[i]) { delete hRefDijetEtaCMBackward1DOldPtBinning[i]; hRefDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
        }

        // Ref-selected jet histograms

        if (hRefSelInclusiveJetPt) delete hRefSelInclusiveJetPt;
        if (hRefSelInclusiveJetEta) delete hRefSelInclusiveJetEta;
        if (hRefSelInclusiveJetEtaUnweighted) delete hRefSelInclusiveJetEtaUnweighted;
        if (hRefSelInclusiveJetPtEta) delete hRefSelInclusiveJetPtEta;
        if (hRefSelInclusiveJetPtEtaPtHat) delete hRefSelInclusiveJetPtEtaPtHat;
        if (hRefSelLeadJetPtEtaPtHat) delete hRefSelLeadJetPtEtaPtHat;
        if (hRefSelSubLeadJetPtEtaPtHat) delete hRefSelSubLeadJetPtEtaPtHat;

        if (hRefSelDijetEta) delete hRefSelDijetEta;
        if (hRefSelDijetPtEtaPhi) delete hRefSelDijetPtEtaPhi;
        if (hRefSelDijetPtEtaPhiWeighted) delete hRefSelDijetPtEtaPhiWeighted;
        if (hRefDijetPtEtaPhiCMInLab) delete hRefDijetPtEtaPhiCMInLab;
        if (hRefSelDijetEtaCM) delete hRefSelDijetEtaCM;
        if (hRefSelDijetPtEtaPhiCM) delete hRefSelDijetPtEtaPhiCM;
        if (hRefSelDijetPtEtaPhiCMWeighted) delete hRefSelDijetPtEtaPhiCMWeighted;
        if (hRefSelDijetPtEtaPhiLabInCM) delete hRefSelDijetPtEtaPhiLabInCM;

        // New ptAve and eta binning
        for (int i = 0; i < 16; ++i) {
            if (hRefSelDijetEta1D[i]) { delete hRefSelDijetEta1D[i]; hRefSelDijetEta1D[i] = nullptr; }
            if (hRefSelRecoDijetEta1D[i]) { delete hRefSelRecoDijetEta1D[i]; hRefSelRecoDijetEta1D[i] = nullptr; }
            if (hRefSelEtaLeadVsEtaSubLead2D[i]) { delete hRefSelEtaLeadVsEtaSubLead2D[i]; hRefSelEtaLeadVsEtaSubLead2D[i] = nullptr; }
            if (hRefSelDijetEtaForward1D[i]) { delete hRefSelDijetEtaForward1D[i]; hRefSelDijetEtaForward1D[i] = nullptr; }
            if (hRefSelDijetEtaBackward1D[i]) { delete hRefSelDijetEtaBackward1D[i]; hRefSelDijetEtaBackward1D[i] = nullptr; }

            if (hRefSelDijetEta1DCM[i]) { delete hRefSelDijetEta1DCM[i]; hRefSelDijetEta1DCM[i] = nullptr; }
            if (hRefSelRecoDijetEta1DCM[i]) { delete hRefSelRecoDijetEta1DCM[i]; hRefSelRecoDijetEta1DCM[i] = nullptr; }
            if (hRefSelEtaLeadVsEtaSubLead2DCM[i]) { delete hRefSelEtaLeadVsEtaSubLead2DCM[i]; hRefSelEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
            if (hRefSelDijetEtaCMForward1D[i]) { delete hRefSelDijetEtaCMForward1D[i]; hRefSelDijetEtaCMForward1D[i] = nullptr; }
            if (hRefSelDijetEtaCMBackward1D[i]) { delete hRefSelDijetEtaCMBackward1D[i]; hRefSelDijetEtaCMBackward1D[i] = nullptr; }
        }

        // Old ptAve and new eta binning
        for (int i = 0; i < 6; ++i) {
            if (hRefSelDijetEta1DOldPt[i]) { delete hRefSelDijetEta1DOldPt[i]; hRefSelDijetEta1DOldPt[i] = nullptr; }
            if (hRefSelRecoDijetEta1DOldPt[i]) { delete hRefSelRecoDijetEta1DOldPt[i]; hRefSelRecoDijetEta1DOldPt[i] = nullptr; }
            if (hRefSelEtaLeadVsEtaSubLead2DOldPt[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPt[i]; hRefSelEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
            if (hRefSelDijetEtaForward1DOldPt[i]) { delete hRefSelDijetEtaForward1DOldPt[i]; hRefSelDijetEtaForward1DOldPt[i] = nullptr; }
            if (hRefSelDijetEtaBackward1DOldPt[i]) { delete hRefSelDijetEtaBackward1DOldPt[i]; hRefSelDijetEtaBackward1DOldPt[i] = nullptr; }

            if (hRefSelDijetEta1DOldPtCM[i]) { delete hRefSelDijetEta1DOldPtCM[i]; hRefSelDijetEta1DOldPtCM[i] = nullptr; }
            if (hRefSelRecoDijetEta1DOldPtCM[i]) { delete hRefSelRecoDijetEta1DOldPtCM[i]; hRefSelRecoDijetEta1DOldPtCM[i] = nullptr; }
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
            if (hRefSelDijetEtaCMForward1DOldPt[i]) { delete hRefSelDijetEtaCMForward1DOldPt[i]; hRefSelDijetEtaCMForward1DOldPt[i] = nullptr; }
            if (hRefSelDijetEtaCMBackward1DOldPt[i]) { delete hRefSelDijetEtaCMBackward1DOldPt[i]; hRefSelDijetEtaCMBackward1DOldPt[i] = nullptr; }
        }

        // Old ptAve and old eta binning
        for (int i = 0; i < 6; ++i) {
            if (hRefSelDijetEta1DOldPtBinning[i]) { delete hRefSelDijetEta1DOldPtBinning[i]; hRefSelDijetEta1DOldPtBinning[i] = nullptr; }
            if (hRefSelRecoDijetEta1DOldPtBinning[i]) { delete hRefSelRecoDijetEta1DOldPtBinning[i]; hRefSelRecoDijetEta1DOldPtBinning[i] = nullptr; }
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
            if (hRefSelDijetEtaForward1DOldPtBinning[i]) { delete hRefSelDijetEtaForward1DOldPtBinning[i]; hRefSelDijetEtaForward1DOldPtBinning[i] = nullptr; }
            if (hRefSelDijetEtaBackward1DOldPtBinning[i]) { delete hRefSelDijetEtaBackward1DOldPtBinning[i]; hRefSelDijetEtaBackward1DOldPtBinning[i] = nullptr; }

            if (hRefSelDijetEta1DOldPtBinningCM[i]) { delete hRefSelDijetEta1DOldPtBinningCM[i]; hRefSelDijetEta1DOldPtBinningCM[i] = nullptr; }
            if (hRefSelRecoDijetEta1DOldPtBinningCM[i]) { delete hRefSelRecoDijetEta1DOldPtBinningCM[i]; hRefSelRecoDijetEta1DOldPtBinningCM[i] = nullptr; }
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
            if (hRefSelDijetEtaCMForward1DOldPtBinning[i]) { delete hRefSelDijetEtaCMForward1DOldPtBinning[i]; hRefSelDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
            if (hRefSelDijetEtaCMBackward1DOldPtBinning[i]) { delete hRefSelDijetEtaCMBackward1DOldPtBinning[i]; hRefSelDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
        }
    } // if (fIsMc)
}

//________________
void HistoManagerDiJet::init() {

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

    const int dijetEtaFBBins{13};
    double dijetEtaFBVals[dijetEtaFBBins+1] { 0.0,  0.2,  0.4,  0.6,  0.8,  
                                              1.0,  1.2,  1.4,  1.6,  1.8,  
                                              2.0,  2.2,  2.4,  3.0 };


    // const int dijetPtBins{16};
    // double dijetPtVals[dijetPtBins+1] {  50.,  60.,   70.,  80.,  90.,
    //                                      100., 110.,  120., 130., 140.,
    //                                      150., 160.,  180., 200., 250., 
    //                                      300., 500.};

    // Old binning convention
    // const int dijetPtOldBins{6};
    // double dijetPtOldVals[dijetPtOldBins+1] {25., 55., 75., 95., 115., 150., 400.}; // 6 bins
    
    int    prescale = 2;

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

    //
    // Gen
    //
    int    bins9D_gen_dijetInfo[9]
    { fDijetPtBins, fDijetEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    double xmin9D_gen_dijetInfo[9]
    { fDijetPtRange[0], fDijetEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0]  };
    double xmax9D_gen_dijetInfo[9]
    { fDijetPtRange[1], fDijetEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

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
    int    bins9D_dijet_info[9]
    { fDijetPtBins, fDijetEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    double xmin9D_dijet_info[9]
    { fDijetPtRange[0], fDijetEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0] };
    double xmax9D_dijet_info[9]
    { fDijetPtRange[1], fDijetEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

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

    hRecoLeadAllJetPtEtaPtHat = new TH3D("hRecoLeadAllJetPtEtaPtHat","Lead jet (matched+unmatched) p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                         prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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

    hRecoSubLeadAllJetPtEtaPtHat = new TH3D("hRecoSubLeadAllJetPtEtaPtHat","SubLead jet (matched+unmatched) p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
    hRecoInclusiveAllJetPtEtaPtHat = new TH3D("hRecoInclusiveAllJetPtEtaPtHat","Reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                           prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                           fPtBins, fPtRange[0], fPtRange[1],
                                           fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveAllJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveAllJetPtEtaPtHat->Sumw2();
    hRecoInclusiveMatchedJetPtEtaPtHat = new TH3D("hRecoInclusiveMatchedJetPtEtaPtHat","Matched reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                         prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveMatchedJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveMatchedJetPtEtaPtHat->Sumw2();
    hRecoInclusiveUnmatchedJetPtEtaPtHat = new TH3D("hRecoInclusiveUnmatchedJetPtEtaPtHat","Unmatched reco jet p_{T} vs #eta vs #hat{p}_{T}};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                         prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    if ( fUseVariableBinning) {
        hRecoInclusiveUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveUnmatchedJetPtEtaPtHat->SetBinsLength(-1);
    }
    hRecoInclusiveUnmatchedJetPtEtaPtHat->Sumw2();

    hRecoGoodInclusiveJetEtaLabFrame = new TH1D("hRecoGoodInclusiveJetEtaLabFrame","Reco good jet #eta in lab frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaLabFrame->Sumw2();
    hRecoGoodInclusiveJetEtaCMFrame = new TH1D("hRecoGoodInclusiveJetEtaCMFrame","Reco good jet #eta in CM frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaCMFrame->Sumw2();

    // Dijet histograms

    // Reco dijet info [9 dimensions]
    // 0 - dijet pt, 1 - dijet eta, 2 - dijet phi,
    // 3 - lead pt, 4 - lead eta, 5 - lead phi,
    // 6 - sublead pt, 7 - sublead eta, 8 - sublead phi    
    hRecoDijetInfo = new THnSparseD("hRecoDijetInfo",
            "Reconstructed dijet and jet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi^{dijet} (rad);p_{T}^{Lead} (GeV);#eta^{Lead};#phi^{Lead} (rad);p_{T}^{SubLead} (GeV);#eta^{SubLead};#phi^{SubLead} (rad)",
            9,
            bins9D_dijet_info,
            xmin9D_dijet_info,
            xmax9D_dijet_info);
    if ( fUseVariableBinning ) {
        hRecoDijetInfo->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetInfo->GetAxis(4)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetInfo->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
    }
    hRecoDijetInfo->Sumw2();
    

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
    hRecoDijetEta = new TH1D("hRecoDijetEta","Reco dijet #eta;Reco #eta^{dijet};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEta->Sumw2();

    // New ptAve and eta binning
    for (unsigned int i{0}; i<fPtAveBins.size()-1; ++i) {
        double ptAveLow = fPtAveBins.at(i);
        double ptAveHi = fPtAveBins.at(i+1);
        
        // Lab frame
        hRecoDijetEta1D[i] = new TH1D(Form("hRecoDijetEta1D_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                      prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRecoDijetEta1D[i]->SetBinsLength(-1);
        }
        hRecoDijetEta1D[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2D_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}", i, ptAveLow, ptAveHi),
                                                      fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2D[i]->Sumw2();
        hRecoDijetEtaForward1D[i] = new TH1D(Form("hRecoDijetEtaForward1D_%d",i),Form("Reco #eta^{dijet}_{forward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                             fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRecoDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRecoDijetEtaForward1D[i]->SetBinsLength(-1);
        }
        hRecoDijetEtaForward1D[i]->Sumw2();
        hRecoDijetEtaBackward1D[i] = new TH1D(Form("hRecoDijetEtaBackward1D_%d",i),Form("Reco #eta^{dijet}_{backward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                              fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRecoDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRecoDijetEtaBackward1D[i]->SetBinsLength(-1);
        }
        hRecoDijetEtaBackward1D[i]->Sumw2();

        // CM frame
        hRecoDijetEta1DCM[i] = new TH1D(Form("hRecoDijetEta1DCM_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                        prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRecoDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRecoDijetEta1DCM[i]->SetBinsLength(-1);
        }
        hRecoDijetEta1DCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DCM_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}", i, ptAveLow, ptAveHi),
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
        hRecoDijetEtaCMForward1D[i] = new TH1D(Form("hRecoDijetEtaCMForward1D_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                               fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRecoDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRecoDijetEtaCMForward1D[i]->SetBinsLength(-1);
        }
        hRecoDijetEtaCMForward1D[i]->Sumw2();
        hRecoDijetEtaCMBackward1D[i] = new TH1D(Form("hRecoDijetEtaCMBackward1D_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRecoDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRecoDijetEtaCMBackward1D[i]->SetBinsLength(-1);
        }   
        hRecoDijetEtaCMBackward1D[i]->Sumw2();
    } // for (int i{0}; i<fPtAveBins.size()-2; ++i)


    for (unsigned int i{0}; i<fPtAveOldBins.size()-1; ++i) {

        double ptAveLow = fPtAveOldBins.at(i);
        double ptAveHi = fPtAveOldBins.at(i+1);

        // Lab frame
        hRecoDijetEta1DOldPt[i] = new TH1D(Form("hRecoDijetEta1DOldPt_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                        dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPt[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Reco dijet #eta lead vs sublead in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
        hRecoDijetEtaForward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaForward1DOldPt_%d",i),Form("Reco #eta^{dijet} forward in %d in the lab frame for %3.0f<p_{T}^{ave} (GeV)<%3.0f;Reco #eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DOldPt[i]->Sumw2();
        hRecoDijetEtaBackward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaBackward1DOldPt_%d",i),Form("Reco #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DOldPt[i]->Sumw2();

        hRecoDijetEta1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinning_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                  dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Reco dijet #eta lead vs sublead in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                              fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaForward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet} forward  in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;Reco #eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaBackward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DOldPtBinning[i]->Sumw2();

        // CM frame
        hRecoDijetEta1DOldPtCM[i] = new TH1D(Form("hRecoDijetEta1DOldPtCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                            dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPtCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
        hRecoDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaCMForward1DOldPt_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DOldPt[i]->Sumw2();
        hRecoDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DOldPt_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DOldPt[i]->Sumw2();

        hRecoDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinningCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                  dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinningCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                              fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
        hRecoDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaCMForward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
    }

    hRecoDijetPtEta = new TH2D("hRecoDijetPtEta", "Reco dijet #eta vs p_{T};p_{T}^{ave} (GeV);#eta^{dijet}", 
                               fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                               fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEta->Sumw2();
    hRecoDijetPtEtaPhi = new TH3D("hRecoDijetPtEtaPhi","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhi->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhi->Sumw2();
    hRecoDijetPtEtaPhiWeighted = new TH3D("hRecoDijetPtEtaPhiWeighted","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiWeighted->Sumw2();
    hRecoDijetPtEtaPhiCMInLab = new TH3D("hRecoDijetPtEtaPhiCMInLab","Reco dijet info in CM with lab frame selection;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiCMInLab->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiCMInLab->Sumw2();
    hRecoDijetEtaCM = new TH1D("hRecoDijetEtaCM","Reco dijet #eta in CM;Reco #eta^{dijet}_{CM};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEtaCM->Sumw2();
    hRecoDijetPtEtaPhiCM = new TH3D("hRecoDijetPtEtaPhiCM","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiCM->Sumw2();
    hRecoDijetPtEtaPhiCMWeighted = new TH3D("hRecoDijetPtEtaPhiCMWeighted","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiCMWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiCMWeighted->Sumw2();
    hRecoDijetPtEtaPhiLabInCM = new TH3D("hRecoDijetPtEtaPhiLabInCM","Reco dijet info in lab with CM selection;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiLabInCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiLabInCM->Sumw2();
    hRecoDijetPtEtaPhiMatched = new TH3D("hRecoDijetPtEtaPhiMatched","Matched reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiMatched->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiMatched->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiMatched->Sumw2();
    hRecoDijetPtEtaPhiCMMatched = new TH3D("hRecoDijetPtEtaPhiCMMatched","Matched reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaPhiCMMatched->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaPhiCMMatched->SetBinsLength(-1);
    }
    hRecoDijetPtEtaPhiCMMatched->Sumw2();

    hRecoDijetPtEtaForward = new TH2D("hRecoDijetPtEtaForward", "Reco dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForward->Sumw2();
    hRecoDijetPtEtaForwardCMInLab = new TH2D("hRecoDijetPtEtaForwardCMInLab", "Reco dijet info in CM frame with lab selection (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForwardCMInLab->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForwardCMInLab->Sumw2();
    hRecoDijetPtEtaBackward = new TH2D("hRecoDijetPtEtaBackward", "Reco dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackward->Sumw2();
    hRecoDijetPtEtaBackwardCMInLab = new TH2D("hRecoDijetPtEtaBackwardCMInLab", "Reco dijet info in CM frame with lab selection (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackwardCMInLab->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackwardCMInLab->Sumw2();
    hRecoDijetPtEtaCMForward = new TH2D("hRecoDijetPtEtaCMForward", "Reco dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMForward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMForward->Sumw2();
    hRecoDijetPtEtaForwardLabInCM = new TH2D("hRecoDijetPtEtaForwardLabInCM", "Reco dijet info in lab frame with CM selection (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForwardLabInCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForwardLabInCM->Sumw2();
    hRecoDijetPtEtaCMBackward = new TH2D("hRecoDijetPtEtaCMBackward", "Reco dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMBackward->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMBackward->Sumw2();
    hRecoDijetPtEtaBackwardLabInCM = new TH2D("hRecoDijetPtEtaBackwardLabInCM", "Reco dijet info in lab frame with CM selection (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackwardLabInCM->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackwardLabInCM->Sumw2();
    
    hRecoDijetPtEtaForwardWeighted = new TH2D("hRecoDijetPtEtaForwardWeighted", "Reco dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaForwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaForwardWeighted->Sumw2();
    hRecoDijetPtEtaBackwardWeighted = new TH2D("hRecoDijetPtEtaBackwardWeighted", "Reco dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaBackwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaBackwardWeighted->Sumw2();
    hRecoDijetPtEtaCMForwardWeighted = new TH2D("hRecoDijetPtEtaCMForwardWeighted", "Reco dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);

    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMForwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMForwardWeighted->Sumw2();

    hRecoDijetPtEtaCMBackwardWeighted = new TH2D("hRecoDijetPtEtaCMBackwardWeighted", "Reco dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if ( fUseVariableBinning ) {
        hRecoDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetPtEtaCMBackwardWeighted->SetBinsLength(-1);
    }
    hRecoDijetPtEtaCMBackwardWeighted->Sumw2();

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
        hGenInclusiveJetPtEtaPtHat = new TH3D("hGenInclusiveJetPtEtaPtHat","Gen inclusive jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
        hGenLeadJetPtEtaPtHat = new TH3D("hGenLeadJetPtEtaPtHat","Gen Lead jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
        hGenSubLeadJetPtEtaPtHat = new TH3D("hGenSubLeadJetPtEtaPtHat","Gen SubLead jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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

        // Gen dijet info [9 dimensions]
        // 0 - dijet pt, 1 - dijet eta, 2 - dijet phi,
        // 3 - lead pt, 4 - lead eta, 5 - lead phi,
        // 6 - sublead pt, 7 - sublead eta, 8 - sublead phi
        hGenDijetInfo = new THnSparseD("hGenDijetInfo","Title;Gen p_{T}^{ave} (GeV);Gen #eta^{dijet};Gen #phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
                9, 
                bins9D_gen_dijetInfo,
                xmin9D_gen_dijetInfo,
                xmax9D_gen_dijetInfo);
        if (fUseVariableBinning) {
            hGenDijetInfo->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetInfo->GetAxis(4)->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetInfo->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        }
        hGenDijetInfo->Sumw2();

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

        for (unsigned int i=0; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hGenDijetEta1D[i] = new TH1D(Form("hGenDijetEta1D_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                         prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) {
            hGenDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1D[i]->SetBinsLength(-1);
            }
            hGenDijetEta1D[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2D_%d",i), Form("Gen #eta^{dijet} Lead vs SubLead #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hGenDijetEtaForward1D[i] = new TH1D(Form("hGenDijetEtaForward1D_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) {
            hGenDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1D[i]->SetBinsLength(-1);
            }
            hGenDijetEtaForward1D[i]->Sumw2();
            hGenDijetEtaBackward1D[i] = new TH1D(Form("hGenDijetEtaBackward1D_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                             fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) {
            hGenDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1D[i]->SetBinsLength(-1);
            }
            hGenDijetEtaBackward1D[i]->Sumw2();
            hGenDijetEta1DCM[i] = new TH1D(Form("hGenDijetEta1DCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                           prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) {
            hGenDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DCM[i]->SetBinsLength(-1);
            }
            hGenDijetEta1DCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DCM_%d",i), Form("Gen #eta^{dijet} Lead vs SubLead #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hGenDijetEtaCMForward1D[i] = new TH1D(Form("hGenDijetEtaCMForward1D_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                              fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) {
            hGenDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1D[i]->SetBinsLength(-1);
            }
            hGenDijetEtaCMForward1D[i]->Sumw2();
            hGenDijetEtaCMBackward1D[i] = new TH1D(Form("hGenDijetEtaCMBackward1D_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) {
            hGenDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1D[i]->SetBinsLength(-1);
            }
            hGenDijetEtaCMBackward1D[i]->Sumw2();
        }

        for (unsigned int i=0; i<fPtAveOldBins.size()-1; i++) {

            double ptAveLow = fPtAveOldBins.at(i);
            double ptAveHi = fPtAveOldBins.at(i+1);
            // New eta binning
            hGenDijetEta1DOldPt[i] = new TH1D(Form("hGenDijetEta1DOldPt_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                              dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPt[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPt_%d",i), Form("Gen #eta^{dijet} Lead vs SubLead #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                      fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hGenDijetEtaForward1DOldPt[i] = new TH1D(Form("hGenDijetEtaForward1DOldPt_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                              dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DOldPt[i]->Sumw2();
            hGenDijetEtaBackward1DOldPt[i] = new TH1D(Form("hGenDijetEtaBackward1DOldPt_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                              dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DOldPt[i]->Sumw2();

        
            hGenDijetEta1DOldPtCM[i] = new TH1D(Form("hGenDijetEta1DOldPtCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                              dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPtCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtCM_%d",i), Form("Gen #eta^{dijet} Lead vs SubLead #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                      fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hGenDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hGenDijetEtaCMForward1DOldPt_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                   dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DOldPt[i]->Sumw2();
            hGenDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hGenDijetEtaCMBackward1DOldPt_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DOldPt[i]->Sumw2();

            // Old eta binning
            hGenDijetEta1DOldPtBinning[i] = new TH1D(Form("hGenDijetEta1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                 dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning_%d",i), Form("Gen #eta^{dijet} Lead vs SubLead #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                         fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hGenDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaForward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                     dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaBackward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                     dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DOldPtBinning[i]->Sumw2();

            hGenDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hGenDijetEta1DOldPtBinningCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                 dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinningCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i), Form("Gen #eta^{dijet} Lead vs SubLead #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                         fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hGenDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaCMForward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                     dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaCMBackward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                     dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
        }

        hGenDijetPtEtaPhi = new TH3D("hGenDijetPtEtaPhi","Gen dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                      fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                      fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                      fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaPhi->SetBinsLength(-1);
        }
        hGenDijetPtEtaPhi->Sumw2();
        hGenDijetPtEtaPhiWeighted = new TH3D("hGenDijetPtEtaPhiWeighted","Gen dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaPhiWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaPhiWeighted->Sumw2();
        hGenDijetPtEtaPhiCMInLab = new TH3D("hGenDijetPtEtaPhiCMInLab","Gen dijet info in CM in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaPhiCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaPhiCMInLab->SetBinsLength(-1);
        }
        hGenDijetPtEtaPhiCMInLab->Sumw2();


        hGenDijetEtaCM = new TH1D("hGenDijetEtaCM", "Gen dijet #eta in CM;#eta^{dijet}_{CM}",
                                  fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEtaCM->Sumw2();
        hGenDijetPtEtaPhiCM = new TH3D("hGenDijetPtEtaPhiCM","Gen dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaPhiCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaPhiCM->Sumw2();
        hGenDijetPtEtaPhiCMWeighted = new TH3D("hGenDijetPtEtaPhiCMWeighted","Gen dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaPhiCMWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaPhiCMWeighted->Sumw2();
        hGenDijetPtEtaPhiLabInCM = new TH3D("hGenDijetPtEtaPhiLabInCM","Gen dijet info in lab in CM frame;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hGenDijetPtEtaPhiLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetPtEtaPhiLabInCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaPhiLabInCM->Sumw2();

        hGenDijetPtEtaForward = new TH2D("hGenDijetPtEtaForward", "Gen dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaForward->Sumw2();
        hGenDijetPtEtaForwardCMInLab = new TH2D("hGenDijetPtEtaForwardCMInLab", "Gen dijet info in CM frame in lab (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForwardCMInLab->SetBinsLength(-1);
        }
        hGenDijetPtEtaForwardCMInLab->Sumw2();
        hGenDijetPtEtaBackward = new TH2D("hGenDijetPtEtaBackward", "Gen dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaBackward->Sumw2();
        hGenDijetPtEtaBackwardCMInLab = new TH2D("hGenDijetPtEtaBackwardCMInLab", "Gen dijet info in CM frame in lab (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackwardCMInLab->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackwardCMInLab->Sumw2();
        hGenDijetPtEtaCMForward = new TH2D("hGenDijetPtEtaCMForward", "Gen dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaCMForward->Sumw2();
        hGenDijetPtEtaForwardLabInCM = new TH2D("hGenDijetPtEtaForwardLabInCM", "Gen dijet info in lab frame in CM (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForwardLabInCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaForwardLabInCM->Sumw2();
        hGenDijetPtEtaCMBackward = new TH2D("hGenDijetPtEtaCMBackward", "Gen dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetPtEtaCMBackward->Sumw2();
        hGenDijetPtEtaBackwardLabInCM = new TH2D("hGenDijetPtEtaBackwardLabInCM", "Gen dijet info in lab frame in CM (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackwardLabInCM->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackwardLabInCM->Sumw2();

        hGenDijetPtEtaForwardWeighted = new TH2D("hGenDijetPtEtaForwardWeighted", "Gen dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaForwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaForwardWeighted->Sumw2();
        hGenDijetPtEtaBackwardWeighted = new TH2D("hGenDijetPtEtaBackwardWeighted", "Gen dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaBackwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaBackwardWeighted->Sumw2();
        hGenDijetPtEtaCMForwardWeighted = new TH2D("hGenDijetPtEtaCMForwardWeighted", "Gen dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaCMForwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMForwardWeighted->Sumw2();
        hGenDijetPtEtaCMBackwardWeighted = new TH2D("hGenDijetPtEtaCMBackwardWeighted", "Gen dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hGenDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetPtEtaCMBackwardWeighted->SetBinsLength(-1);
        }
        hGenDijetPtEtaCMBackwardWeighted->Sumw2();



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

        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed inclusive jets;Reco p_{T, corr}^{Inclusive} (GeV);Reco p_{T, raw}^{Inclusive} (GeV);Ref p_{T}^{Inclusive} (GeV);Reco #eta^{Inclusive};Ref #eta^{Inclusive}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
        hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed Lead jets;Reco p_{T, corr}^{Lead} (GeV);Reco p_{T, raw}^{Lead} (GeV);Ref p_{T}^{Lead} (GeV);Reco #eta^{Lead};Ref #eta^{Lead}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
        hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed SubLead jets;Reco p_{T, corr}^{SubLead} (GeV);Reco p_{T, raw}^{SubLead} (GeV);Ref p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref #eta^{SubLead}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
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
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESRecoPtRecoEtaPtHatWeighted","JES vs p_{T}^{reco} vs #eta^{reco} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{reco} (GeV);#eta^{reco};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Sumw2();
        hLeadJetJESGenPtEtaPtHatWeighted = new THnSparseD("hLeadJetJESGenPtEtaPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hLeadJetJESGenPtEtaPtHatWeighted->Sumw2();
        hSubLeadJetJESGenPtEtaPtHatWeighted = new THnSparseD("hSubLeadJetJESGenPtEtaPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hSubLeadJetJESGenPtEtaPtHatWeighted->Sumw2();

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
                                              prescale * fEtaBins, fEtaRange[0], fEtaRange[1], 
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
                                                prescale *  fEtaBins, fEtaRange[0], fEtaRange[1], 
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
                                                 prescale * fEtaBins, fEtaRange[0], fEtaRange[1], 
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
                                                   prescale * fEtaBins, fEtaRange[0], fEtaRange[1], 
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
        hReco2RefFull = new THnSparseD("hReco2RefFull",
                "Reco to ref correspondence;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Lead} (GeV);Reco #eta^{Lead};Reco p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Lead} (GeV);Ref #eta^{Lead};Ref p_{T}^{SubLead} (GeV);Ref #eta^{SubLead}",
                12,
                bins12D_ref2reco_info,
                xmin12D_ref2reco_info,
                xmax12D_ref2reco_info);
        if ( fUseVariableBinning ) {
            hReco2RefFull->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
            hReco2RefFull->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
            hReco2RefFull->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
            hReco2RefFull->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
            hReco2RefFull->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
            hReco2RefFull->GetAxis(11)->Set(dijetEtaBins, dijetEtaVals);
        } 
        hReco2RefFull->Sumw2();

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
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
        hRefSel2RecoFull = new THnSparseD("hRefSel2RecoFull",
                "Reco to ref correspondence (via ref selection) weighted;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Lead} (GeV);Reco #eta^{Lead};Reco p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Lead} (GeV);Ref #eta^{Lead};Ref p_{T}^{SubLead} (GeV);Ref #eta^{SubLead}",
                12,
                bins12D_ref2reco_info,
                xmin12D_ref2reco_info,
                xmax12D_ref2reco_info);
        if ( fUseVariableBinning ) {
            hRefSel2RecoFull->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
            hRefSel2RecoFull->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
            hRefSel2RecoFull->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
            hRefSel2RecoFull->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
            hRefSel2RecoFull->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
            hRefSel2RecoFull->GetAxis(11)->Set(dijetEtaBins, dijetEtaVals);
        } 
        hRefSel2RecoFull->Sumw2();
        hRefSelDijetEta = new TH1D("hRefSelDijetEta","Ref selected dijets;#eta^{dijet}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEta->Sumw2();
        hRefSelDijetPtEtaPhi = new TH3D("hRefSelDijetPtEtaPhi","RefSel dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaPhi->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaPhi->Sumw2();
        hRefSelDijetPtEtaPhiWeighted = new TH3D("hRefSelDijetPtEtaPhiWeighted","RefSel dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaPhiWeighted->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaPhiWeighted->Sumw2();
        hRefSelDijetPtEtaPhiCMInLab = new TH3D("hRefSelDijetPtEtaPhiCMInLab","RefSel dijet info in CM in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaPhiCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaPhiCMInLab->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaPhiCMInLab->Sumw2();
        hRefSelDijetEtaCM = new TH1D("hRefSelDijetEtaCM","Ref selected dijets in CM;#eta^{dijet}_{CM}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEtaCM->Sumw2();
        hRefSelDijetPtEtaPhiCM = new TH3D("hRefSelDijetPtEtaPhiCM","RefSel dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaPhiCM->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaPhiCM->Sumw2();
        hRefSelDijetPtEtaPhiCMWeighted = new TH3D("hRefSelDijetPtEtaPhiCMWeighted","RefSel dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaPhiCMWeighted->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaPhiCMWeighted->Sumw2();
        hRefSelDijetPtEtaPhiLabInCM = new TH3D("hRefSelDijetPtEtaPhiLabInCM","RefSel dijet info in lab in CM frame;p_{T}^{ave} (GeV);#eta^{dijet}_{lab};#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fPhiBins, fPhiRange[0], fPhiRange[1] );
        if ( fUseVariableBinning ) {
            hRefSelDijetPtEtaPhiLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetPtEtaPhiLabInCM->SetBinsLength(-1);
        }
        hRefSelDijetPtEtaPhiLabInCM->Sumw2();

        // New pT and eta binning
        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hRefSelDijetEta1D[i] = new TH1D(Form("hRefSelDijetEta1D_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1D[i]->SetBinsLength(-1);
            }
            hRefSelDijetEta1D[i]->Sumw2();
            hRefSelRecoDijetEta1D[i] = new TH1D(Form("hRefSelRecoDijetEta1D_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1D[i]->SetBinsLength(-1);
            }
            hRefSelRecoDijetEta1D[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2D_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hRefSelDijetEtaForward1D[i] = new TH1D(Form("hRefSelDijetEtaForward1D_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                               fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1D[i]->SetBinsLength(-1);
            }
            hRefSelDijetEtaForward1D[i]->Sumw2();
            hRefSelDijetEtaBackward1D[i] = new TH1D(Form("hRefSelDijetEtaBackward1D_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1D[i]->SetBinsLength(-1);
            }
            hRefSelDijetEtaBackward1D[i]->Sumw2();


            hRefSelDijetEta1DCM[i] = new TH1D(Form("hRefSelDijetEta1DCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                              prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DCM[i]->SetBinsLength(-1);
            }
            hRefSelDijetEta1DCM[i]->Sumw2();
            hRefSelRecoDijetEta1DCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                              prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelRecoDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DCM[i]->SetBinsLength(-1);
            }
            hRefSelRecoDijetEta1DCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DCM_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hRefSelDijetEtaCMForward1D[i] = new TH1D(Form("hRefSelDijetEtaCMForward1D_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                 fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1D[i]->SetBinsLength(-1);
            }
            hRefSelDijetEtaCMForward1D[i]->Sumw2();
            hRefSelDijetEtaCMBackward1D[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1D_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                   fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefSelDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1D[i]->SetBinsLength(-1);
            }
            hRefSelDijetEtaCMBackward1D[i]->Sumw2();
        } // for (int i{0}; i<fPtAveBins.size()-1; i++)

        // Old pT binning
        for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {
            double ptAveLow = fPtAveOldBins.at(i);
            double ptAveHi = fPtAveOldBins.at(i+1);

            // New eta binning
            hRefSelDijetEta1DOldPt[i] = new TH1D(Form("hRefSelDijetEta1DOldPt_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                 dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DOldPt[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPt[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPt_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                 dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DOldPt[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hRefSelDijetEtaForward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaForward1DOldPt_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                       dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DOldPt[i]->Sumw2();
            hRefSelDijetEtaBackward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaBackward1DOldPt_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                       dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DOldPt[i]->Sumw2();

            hRefSelDijetEta1DOldPtCM[i] = new TH1D(Form("hRefSelDijetEta1DOldPtCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                 dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DOldPtCM[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                 dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DOldPtCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hRefSelDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DOldPt_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                       dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DOldPt[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DOldPt_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                       dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DOldPt[i]->Sumw2();
            
            // Old eta binning
            hRefSelDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEta1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                    dijetEtaOldBins, dijetEtaOldVals);
            hRefSelDijetEta1DOldPtBinning[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtBinning_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                    dijetEtaOldBins, dijetEtaOldVals);
            hRefSelRecoDijetEta1DOldPtBinning[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                           fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaForward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaBackward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DOldPtBinning[i]->Sumw2();
            
            // Old eta binning
            hRefSelDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRefSelDijetEta1DOldPtBinningCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaOldBins, dijetEtaOldVals);
            hRefSelDijetEta1DOldPtBinningCM[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtBinningCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaOldBins, dijetEtaOldVals);
            hRefSelRecoDijetEta1DOldPtBinningCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hRefSelDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                               dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                               dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();

        } // for (int i{0}; i<fPtAveOldBins.size()-1; i++)


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
        hRefInclusiveJetPtEtaPtHat = new TH3D("hRefInclusiveJetPtEtaPtHat","Ref jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
        hRefLeadJetPtEtaPtHat = new TH3D("hRefLeadJetPtEtaPtHat","Ref Lead jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
        hRefSubLeadJetPtEtaPtHat = new TH3D("hRefSubLeadJetPtEtaPtHat","Ref SubLead jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Reco p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
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

        // New pT and eta binning
        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hRefDijetEta1D[i] = new TH1D(Form("hRefDijetEta1D_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                         prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1D[i]->SetBinsLength(-1);
            }
            hRefDijetEta1D[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2D_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hRecoVsRefDijetEta2D[i] = new TH2D(Form("hRecoVsRefDijetEta2D_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2D[i]->Sumw2();
            hRecoVsRefLeadJetEta2D[i] = new TH2D(Form("hRecoVsRefLeadJetEta2D_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2D[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2D[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2D_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2D[i]->Sumw2();
            hRefDijetEtaForward1D[i] = new TH1D(Form("hRefDijetEtaForward1D_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1D[i]->SetBinsLength(-1);
            }
            hRefDijetEtaForward1D[i]->Sumw2();
            hRefDijetEtaBackward1D[i] = new TH1D(Form("hRefDijetEtaBackward1D_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                             fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1D[i]->SetBinsLength(-1);
            }
            hRefDijetEtaBackward1D[i]->Sumw2();

            // Ref dijets in CM frame
            hRefDijetEta1DCM[i] = new TH1D(Form("hRefDijetEta1DCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                           prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DCM[i]->SetBinsLength(-1);
            }
            hRefDijetEta1DCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DCM_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hRecoVsRefDijetEta2DCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DCM_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DCM_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DCM[i]->Sumw2();
            hRefDijetEtaCMForward1D[i] = new TH1D(Form("hRefDijetEtaCMForward1D_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                              fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1D[i]->SetBinsLength(-1);
            }
            hRefDijetEtaCMForward1D[i]->Sumw2();
            hRefDijetEtaCMBackward1D[i] = new TH1D(Form("hRefDijetEtaCMBackward1D_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                               fDijetEtaBins, 0., fDijetEtaRange[1]);
            if ( fUseVariableBinning ) {
            hRefDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1D[i]->SetBinsLength(-1);
            }
            hRefDijetEtaCMBackward1D[i]->Sumw2();
        }

        // Old pT binning
        for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {

            double ptAveLow = fPtAveOldBins.at(i);
            double ptAveHi = fPtAveOldBins.at(i+1);
            // New eta binning
            hRefDijetEta1DOldPt[i] = new TH1D(Form("hRefDijetEta1DOldPt_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                         dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPt[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPt[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPt[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPt[i]->Sumw2();
            hRefDijetEtaForward1DOldPt[i] = new TH1D(Form("hRefDijetEtaForward1DOldPt_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DOldPt[i]->Sumw2();
            hRefDijetEtaBackward1DOldPt[i] = new TH1D(Form("hRefDijetEtaBackward1DOldPt_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DOldPt[i]->Sumw2();

            hRefDijetEta1DOldPtCM[i] = new TH1D(Form("hRefDijetEta1DOldPtCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                         dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPtCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtCM[i]->Sumw2();
            hRefDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hRefDijetEtaCMForward1DOldPt_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DOldPt[i]->Sumw2();
            hRefDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hRefDijetEtaCMBackward1DOldPt_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DOldPt[i]->Sumw2();

            // Old eta binning
            hRefDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefDijetEta1DOldPtBinning_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinning[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);                                  
            hRecoVsRefLeadJetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i]->Sumw2();
            hRefDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaForward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DOldPtBinning[i]->Sumw2();
            hRefDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaBackward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DOldPtBinning[i]->Sumw2();

            hRefDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRefDijetEta1DOldPtBinningCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinningCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]->Sumw2();
            hRefDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaCMForward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
            hRefDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaCMBackward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                            dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
        }

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
        hRefDijetPtEtaPhi = new TH3D("hRefDijetPtEtaPhi","Ref dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) {
            hRefDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaPhi->SetBinsLength(-1);
        }
        hRefDijetPtEtaPhi->Sumw2();
        hRefDijetPtEtaPhiWeighted = new TH3D("hRefDijetPtEtaPhiWeighted","Ref dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) {
            hRefDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaPhiWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaPhiWeighted->Sumw2();
        hRefDijetPtEtaPhiCMInLab = new TH3D("hRefDijetPtEtaPhiCMInLab","Ref dijet info in CM in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) {
            hRefDijetPtEtaPhiCMInLab->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaPhiCMInLab->SetBinsLength(-1);
        }
        hRefDijetPtEtaPhiCMInLab->Sumw2();

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
        hRefDijetPtEtaPhiCM = new TH3D("hRefDijetPtEtaPhiCM","Ref dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) {
            hRefDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaPhiCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaPhiCM->Sumw2();
        hRefDijetPtEtaPhiCMWeighted = new TH3D("hRefDijetPtEtaPhiCMWeighted","Ref dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) {
            hRefDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaPhiCMWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaPhiCMWeighted->Sumw2();
        hRefDijetPtEtaPhiLabInCM = new TH3D("hRefDijetPtEtaPhiLabInCM","Ref dijet info in lab frame in CM;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) {
            hRefDijetPtEtaPhiLabInCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetPtEtaPhiLabInCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaPhiLabInCM->Sumw2();

        hRefDijetPtEtaForward = new TH2D("hRefDijetPtEtaForward", "Ref dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForward->SetBinsLength(-1);
        }
        hRefDijetPtEtaForward->Sumw2();
        hRefDijetPtEtaForwardCMInLab = new TH2D("hRefDijetPtEtaForwardCMInLab", "Ref dijet info in CM frame (forward) in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForwardCMInLab->SetBinsLength(-1);
        }
        hRefDijetPtEtaForwardCMInLab->Sumw2();
        hRefDijetPtEtaBackward = new TH2D("hRefDijetPtEtaBackward", "Ref dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackward->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackward->Sumw2();
        hRefDijetPtEtaBackwardCMInLab = new TH2D("hRefDijetPtEtaBackwardCMInLab", "Ref dijet info in CM frame (backward) in lab frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackwardCMInLab->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackwardCMInLab->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackwardCMInLab->Sumw2();
        hRefDijetPtEtaCMForward = new TH2D("hRefDijetPtEtaCMForward", "Ref dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMForward->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMForward->Sumw2();
        hRefDijetPtEtaForwardLabInCM = new TH2D("hRefDijetPtEtaForwardLabInCM", "Ref dijet info in lab frame (forward) in CM;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForwardLabInCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaForwardLabInCM->Sumw2();
        hRefDijetPtEtaCMBackward = new TH2D("hRefDijetPtEtaCMBackward", "Ref dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMBackward->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMBackward->Sumw2();
        hRefDijetPtEtaBackwardLabInCM = new TH2D("hRefDijetPtEtaBackwardLabInCM", "Ref dijet info in lab frame (backward) in CM;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackwardLabInCM->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackwardLabInCM->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackwardLabInCM->Sumw2();
        
        hRefDijetPtEtaForwardWeighted = new TH2D("hRefDijetPtEtaForwardWeighted", "Ref dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaForwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaForwardWeighted->Sumw2();
        hRefDijetPtEtaBackwardWeighted = new TH2D("hRefDijetPtEtaBackwardWeighted", "Ref dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaBackwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaBackwardWeighted->Sumw2();
        hRefDijetPtEtaCMForwardWeighted = new TH2D("hRefDijetPtEtaCMForwardWeighted", "Ref dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) {
            hRefDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetPtEtaCMForwardWeighted->SetBinsLength(-1);
        }
        hRefDijetPtEtaCMForwardWeighted->Sumw2();

        hRefDijetPtEtaCMBackwardWeighted = new TH2D("hRefDijetPtEtaCMBackwardWeighted", "Ref dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
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
    hRecoDijetInfo->Write();
    hRecoInclusiveAllJetPt->Write();
    hRecoInclusiveAllJetEta->Write();
    hRecoInclusiveAllJetEtaUnweighted->Write();
    hRecoInclusiveAllJetPtEta->Write();
    hRecoInclusiveAllJetPtEtaPtHat->Write();
    hRecoInclusiveMatchedJetPtEtaPtHat->Write();
    hRecoInclusiveUnmatchedJetPtEtaPtHat->Write();
    hRecoPtLeadPtSublead->Write();
    hRecoEtaLeadEtaSublead->Write();
    hRecoEtaCMLeadEtaCMSublead->Write();
    hRecoPtLeadPtSubleadMcReweight->Write();
    hRecoEtaLeadEtaSubleadMcReweight->Write();
    hRecoDijetPtEta->Write();
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

    hRecoDijetPtEtaPhi->Write();
    hRecoDijetPtEtaPhiWeighted->Write();
    hRecoDijetPtEtaPhiCMInLab->Write();
    hRecoDijetPtEtaPhiCM->Write();
    hRecoDijetPtEtaPhiCMWeighted->Write();
    hRecoDijetPtEtaPhiLabInCM->Write();
    if ( fIsMc ) {
        hRecoDijetPtEtaPhiMatched->Write();
        hRecoDijetPtEtaPhiCMMatched->Write();
    }

    hRecoLeadAllJetPtEta->Write();
    hRecoLeadAllJetPtEtaPtHat->Write();
    hRecoSubLeadAllJetPtEta->Write();
    hRecoSubLeadAllJetPtEtaPtHat->Write();

    hRecoGoodInclusiveJetEtaLabFrame->Write();
    hRecoGoodInclusiveJetEtaCMFrame->Write();

    hRecoDijetEta->Write();
    hRecoDijetEtaCM->Write();

    for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
        hRecoDijetEta1D[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2D[i]->Write();
        hRecoDijetEtaForward1D[i]->Write();
        hRecoDijetEtaBackward1D[i]->Write();

        hRecoDijetEta1DCM[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DCM[i]->Write();
        hRecoDijetEtaCMForward1D[i]->Write();
        hRecoDijetEtaCMBackward1D[i]->Write();
    }

    for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
        hRecoDijetEta1DOldPt[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]->Write();
        hRecoDijetEtaForward1DOldPt[i]->Write();
        hRecoDijetEtaBackward1DOldPt[i]->Write();
        hRecoDijetEta1DOldPtCM[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
        hRecoDijetEtaCMForward1DOldPt[i]->Write();
        hRecoDijetEtaCMBackward1DOldPt[i]->Write();
        hRecoDijetEta1DOldPtBinning[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
        hRecoDijetEtaForward1DOldPtBinning[i]->Write();
        hRecoDijetEtaBackward1DOldPtBinning[i]->Write();
        hRecoDijetEta1DOldPtBinningCM[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
        hRecoDijetEtaCMForward1DOldPtBinning[i]->Write();
        hRecoDijetEtaCMBackward1DOldPtBinning[i]->Write();
    }


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
        hGenDijetInfo->Write();
        hGenInclusiveJetPt->Write();
        hGenInclusiveJetEta->Write();
        hGenInclusiveJetEtaUnweighted->Write();
        hGenInclusiveJetPtEta->Write();
        hGenInclusiveJetPtEtaPtHat->Write();
        hGenLeadJetPtEta->Write();
        hGenLeadJetPtEtaPtHat->Write();
        hGenSubLeadJetPtEta->Write();
        hGenSubLeadJetPtEtaPtHat->Write();
        hGenPtLeadPtSublead->Write();
        hGenEtaLeadEtaSublead->Write();
        hGenEtaCMLeadEtaCMSublead->Write();
        hGenPtLeadPtSubleadMcReweight->Write();
        hGenEtaLeadEtaSubleadMcReweight->Write();
        hGenDijetEta->Write();
        hGenDijetPtEtaPhi->Write();
        hGenDijetPtEtaPhiWeighted->Write();
        hGenDijetPtEtaPhiCMInLab->Write();
        hGenDijetEtaCM->Write();
        hGenDijetPtEtaPhiCM->Write();
        hGenDijetPtEtaPhiCMWeighted->Write();
        hGenDijetPtEtaPhiLabInCM->Write();
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

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            hGenDijetEta1D[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2D[i]->Write();
            hGenDijetEtaForward1D[i]->Write();
            hGenDijetEtaBackward1D[i]->Write();
            hGenDijetEta1DCM[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DCM[i]->Write();
            hGenDijetEtaCMForward1D[i]->Write();
            hGenDijetEtaCMBackward1D[i]->Write();
        }

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            hGenDijetEta1DOldPt[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]->Write();
            hGenDijetEtaForward1DOldPt[i]->Write();
            hGenDijetEtaBackward1DOldPt[i]->Write();
            hGenDijetEta1DOldPtCM[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
            hGenDijetEtaCMForward1DOldPt[i]->Write();
            hGenDijetEtaCMBackward1DOldPt[i]->Write();
            hGenDijetEta1DOldPtBinning[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
            hGenDijetEtaForward1DOldPtBinning[i]->Write();
            hGenDijetEtaBackward1DOldPtBinning[i]->Write();
            hGenDijetEta1DOldPtBinningCM[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
            hGenDijetEtaCMForward1DOldPtBinning[i]->Write();
            hGenDijetEtaCMBackward1DOldPtBinning[i]->Write();
        }

        //
        // Ref hisrograms
        //

        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoLeadJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoSubLeadJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
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
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Write();
        hLeadJetJESGenPtEtaPtHatWeighted->Write();
        hSubLeadJetJESGenPtEtaPtHatWeighted->Write();

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
        hRefInclusiveJetPtEtaPtHat->Write();
        hRefLeadJetPtEta->Write();
        hRefLeadJetPtEtaPtHat->Write();
        hRefLeadUnswappedJetPtEta->Write();
        hRefLeadUnswappedJetPtEtaPtHat->Write();
        hRefSubLeadJetPtEta->Write();
        hRefSubLeadJetPtEtaPtHat->Write();
        hRefSubLeadUnswappedJetPtEta->Write();
        hRefSubLeadUnswappedJetPtEtaPtHat->Write();

        hReco2RefFull->Write();

        hRecoDijetPtEtaRefDijetPtEta->Write();
        hRecoDijetPtEtaRefDijetPtEtaWeighted->Write();

        hRefSel2RecoFull->Write();

        hRefDijetEta->Write();
        hRefDijetEtaVsRecoDijetEta->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Write();
        hRefDijetPtEtaPhi->Write();
        hRefDijetPtEtaPhiWeighted->Write();
        hRefDijetPtEtaPhiCMInLab->Write();

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
        hRefDijetPtEtaPhiCM->Write();
        hRefDijetPtEtaPhiCMWeighted->Write();
        hRefDijetPtEtaPhiLabInCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Write();

        hRefPtLeadPtSublead->Write();
        hRefEtaLeadEtaSublead->Write();
        hRefEtaCMLeadEtaCMSublead->Write();
        hRefPtLeadPtSubleadMcReweight->Write();
        hRefEtaLeadEtaSubleadMcReweight->Write();

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            hRefDijetEta1D[i]->Write();
            hRefEtaLeadVsEtaSubLead2D[i]->Write();
            hRecoVsRefDijetEta2D[i]->Write();
            hRecoVsRefLeadJetEta2D[i]->Write();
            hRecoVsRefSubLeadJetEta2D[i]->Write();
            hRefDijetEtaForward1D[i]->Write();
            hRefDijetEtaBackward1D[i]->Write();

            hRefDijetEta1DCM[i]->Write();
            hRefEtaLeadVsEtaSubLead2DCM[i]->Write();
            hRecoVsRefDijetEta2DCM[i]->Write();
            hRecoVsRefLeadJetEta2DCM[i]->Write();
            hRecoVsRefSubLeadJetEta2DCM[i]->Write();
            hRefDijetEtaCMForward1D[i]->Write();
            hRefDijetEtaCMBackward1D[i]->Write();
        }

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            hRefDijetEta1DOldPt[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPt[i]->Write();
            hRecoVsRefDijetEta2DOldPt[i]->Write();
            hRecoVsRefLeadJetEta2DOldPt[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPt[i]->Write();
            hRefDijetEtaForward1DOldPt[i]->Write();
            hRefDijetEtaBackward1DOldPt[i]->Write();

            hRefDijetEta1DOldPtCM[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
            hRecoVsRefDijetEta2DOldPtCM[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtCM[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtCM[i]->Write();
            hRefDijetEtaCMForward1DOldPt[i]->Write();
            hRefDijetEtaCMBackward1DOldPt[i]->Write();

            hRefDijetEta1DOldPtBinning[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
            hRecoVsRefDijetEta2DOldPtBinning[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtBinning[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i]->Write();
            hRefDijetEtaForward1DOldPtBinning[i]->Write();
            hRefDijetEtaBackward1DOldPtBinning[i]->Write();

            hRefDijetEta1DOldPtBinningCM[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
            hRecoVsRefDijetEta2DOldPtBinningCM[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]->Write();
            hRefDijetEtaCMForward1DOldPtBinning[i]->Write();
            hRefDijetEtaCMBackward1DOldPtBinning[i]->Write();
        }


        //
        // Ref-seletected histograms
        //

        hRefSelInclusiveJetPt->Write();
        hRefSelInclusiveJetEta->Write();
        hRefSelInclusiveJetEtaUnweighted->Write();
        hRefSelInclusiveJetPtEta->Write();
        hRefSelInclusiveJetPtEtaPtHat->Write();
        hRefSelLeadJetPtEtaPtHat->Write();
        hRefSelSubLeadJetPtEtaPtHat->Write();

        hRefSelDijetEta->Write();
        hRefSelDijetPtEtaPhi->Write();
        hRefSelDijetPtEtaPhiWeighted->Write();
        hRefSelDijetPtEtaPhiCMInLab->Write();
        hRefSelDijetEtaCM->Write();
        hRefSelDijetPtEtaPhiCM->Write();
        hRefSelDijetPtEtaPhiCMWeighted->Write();
        hRefSelDijetPtEtaPhiLabInCM->Write();

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            hRefSelDijetEta1D[i]->Write();
            hRefSelRecoDijetEta1D[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2D[i]->Write();
            hRefSelDijetEtaForward1D[i]->Write();
            hRefSelDijetEtaBackward1D[i]->Write();
            hRefSelDijetEta1DCM[i]->Write();
            hRefSelRecoDijetEta1DCM[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DCM[i]->Write();
            hRefSelDijetEtaCMForward1D[i]->Write();
            hRefSelDijetEtaCMBackward1D[i]->Write();
        }

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            hRefSelDijetEta1DOldPt[i]->Write();
            hRefSelRecoDijetEta1DOldPt[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i]->Write();
            hRefSelDijetEtaForward1DOldPt[i]->Write();
            hRefSelDijetEtaBackward1DOldPt[i]->Write();
            hRefSelDijetEta1DOldPtCM[i]->Write();
            hRefSelRecoDijetEta1DOldPtCM[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
            hRefSelDijetEtaCMForward1DOldPt[i]->Write();
            hRefSelDijetEtaCMBackward1DOldPt[i]->Write();
            hRefSelDijetEta1DOldPtBinning[i]->Write();
            hRefSelRecoDijetEta1DOldPtBinning[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
            hRefSelDijetEtaForward1DOldPtBinning[i]->Write();
            hRefSelDijetEtaBackward1DOldPtBinning[i]->Write();
            hRefSelDijetEta1DOldPtBinningCM[i]->Write();
            hRefSelRecoDijetEta1DOldPtBinningCM[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
            hRefSelDijetEtaCMForward1DOldPtBinning[i]->Write();
            hRefSelDijetEtaCMBackward1DOldPtBinning[i]->Write();
        }

    } // if ( fIsMc )

}
