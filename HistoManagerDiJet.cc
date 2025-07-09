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
HistoManagerDiJet::HistoManagerDiJet() : BaseHistoManager(),

    // Event histograms
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
    // Gen jet histograms
    //

    hGenJetCollectionSize{nullptr},    
    hGenVsRecoJetCollectionSize{nullptr},
    hGenLeadingJetPtOverPtHatVsLeadingJetPt{nullptr},
    hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted{nullptr},
    hGenDijetPtOverPtHatVsDijetPt{nullptr},
    hGenDijetPtOverPtHatVsDijetPtWeighted{nullptr},
    hGenDijetPtAveOverPtHatVsDijetPtAve{nullptr},
    hGenDijetPtAveOverPtHatVsDijetPtAveWeighted{nullptr},

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

    hGenDijetLeadSubLead{nullptr},
    hGenDijetLeadSubLeadCM{nullptr},

    hGenPtLeadPtSublead{nullptr},
    hGenEtaLeadEtaSublead{nullptr},
    hGenDijetEta{nullptr},
    hGenDijetPtEtaPhi{nullptr},
    hGenDijetPtEtaPhiWeighted{nullptr},

    hGenPtLeadPtSubleadCM{nullptr},
    hGenEtaCMLeadEtaCMSublead{nullptr},
    hGenDijetEtaCM{nullptr},
    hGenDijetPtEtaPhiCM{nullptr},
    hGenDijetPtEtaPhiCMWeighted{nullptr},

    hGenDijetPtEtaForward{nullptr},
    hGenDijetPtEtaBackward{nullptr},
    hGenDijetPtEtaCMForward{nullptr},
    hGenDijetPtEtaCMBackward{nullptr},
    hGenDijetPtEtaForwardWeighted{nullptr},
    hGenDijetPtEtaBackwardWeighted{nullptr},
    hGenDijetPtEtaCMForwardWeighted{nullptr},
    hGenDijetPtEtaCMBackwardWeighted{nullptr},

    hGenDijetEta1D{nullptr},
    hGenDijetEta1DWeighted{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2D{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DWeighted{nullptr},
    hGenDijetEtaForward1D{nullptr},
    hGenDijetEtaForward1DWeighted{nullptr},
    hGenDijetEtaBackward1D{nullptr},
    hGenDijetEtaBackward1DWeighted{nullptr},

    hGenDijetEta1DCM{nullptr},
    hGenDijetEta1DCMWeighted{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DCM{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DCMWeighted{nullptr},
    hGenDijetEtaCMForward1D{nullptr},
    hGenDijetEtaCMForward1DWeighted{nullptr},
    hGenDijetEtaCMBackward1D{nullptr},
    hGenDijetEtaCMBackward1DWeighted{nullptr},

    hGenDijetEta1DOldPt{nullptr},
    hGenDijetEta1DOldPtWeighted{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted{nullptr},
    hGenDijetEtaForward1DOldPt{nullptr},
    hGenDijetEtaForward1DOldPtWeighted{nullptr},
    hGenDijetEtaBackward1DOldPt{nullptr},
    hGenDijetEtaBackward1DOldPtWeighted{nullptr},

    hGenDijetEta1DOldPtCM{nullptr},
    hGenDijetEta1DOldPtCMWeighted{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted{nullptr},
    hGenDijetEtaCMForward1DOldPt{nullptr},
    hGenDijetEtaCMForward1DOldPtWeighted{nullptr},
    hGenDijetEtaCMBackward1DOldPt{nullptr},
    hGenDijetEtaCMBackward1DOldPtWeighted{nullptr},

    hGenDijetEta1DOldPtBinning{nullptr},
    hGenDijetEta1DOldPtBinningWeighted{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted{nullptr},
    hGenDijetEtaForward1DOldPtBinning{nullptr},
    hGenDijetEtaForward1DOldPtBinningWeighted{nullptr},
    hGenDijetEtaBackward1DOldPtBinning{nullptr},
    hGenDijetEtaBackward1DOldPtBinningWeighted{nullptr},

    hGenDijetEta1DOldPtBinningCM{nullptr},
    hGenDijetEta1DOldPtBinningCMWeighted{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted{nullptr},
    hGenDijetEtaCMForward1DOldPtBinning{nullptr},
    hGenDijetEtaCMForward1DOldPtBinningWeighted{nullptr},
    hGenDijetEtaCMBackward1DOldPtBinning{nullptr},
    hGenDijetEtaCMBackward1DOldPtBinningWeighted{nullptr},

    //
    // Reco jet histograms
    //

    hRecoJetCollectionSize{nullptr},

    hRecoDijetLeadSubLead{nullptr},
    hRecoDijetLeadSubLeadCM{nullptr},

    hRecoPtLeadPtSublead{nullptr},
    hRecoEtaLeadEtaSublead{nullptr},
    hRecoPtLeadPtSubleadMcReweight{nullptr},
    hRecoEtaLeadEtaSubleadMcReweight{nullptr},
    hRecoDijetEta{nullptr},
    hRecoDijetPtEta{nullptr},
    hRecoDijetPtEtaForward{nullptr},
    hRecoDijetPtEtaBackward{nullptr},
    hRecoDijetPtEtaForwardWeighted{nullptr},
    hRecoDijetPtEtaBackwardWeighted{nullptr},
    hRecoDijetPtEtaPhi{nullptr},
    hRecoDijetPtEtaPhiWeighted{nullptr},

    hRecoPtLeadPtSubleadCM{nullptr},
    hRecoEtaCMLeadEtaCMSublead{nullptr},
    hRecoPtLeadPtSubleadCMMcReweight{nullptr},
    hRecoEtaCMLeadEtaCMSubleadMcReweight{nullptr},
    hRecoDijetEtaCM{nullptr},
    hRecoDijetPtEtaCM{nullptr},
    hRecoDijetPtEtaCMForward{nullptr},
    hRecoDijetPtEtaCMBackward{nullptr},
    hRecoDijetPtEtaCMForwardWeighted{nullptr},
    hRecoDijetPtEtaCMBackwardWeighted{nullptr},
    hRecoDijetPtEtaPhiCM{nullptr},
    hRecoDijetPtEtaPhiCMWeighted{nullptr},

    hRecoDijetEta1D{nullptr},
    hRecoDijetEta1DWeighted{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2D{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DWeighted{nullptr},
    hRecoDijetEtaForward1D{nullptr},
    hRecoDijetEtaForward1DWeighted{nullptr},
    hRecoDijetEtaBackward1D{nullptr},
    hRecoDijetEtaBackward1DWeighted{nullptr},

    hRecoDijetEta1DCM{nullptr},
    hRecoDijetEta1DCMWeighted{nullptr},
    hRecoEtaLeadVsEtaSubLead2DCM{nullptr},
    hRecoEtaLeadVsEtaSubLead2DCMWeighted{nullptr},
    hRecoDijetEtaCMForward1D{nullptr},
    hRecoDijetEtaCMForward1DWeighted{nullptr},
    hRecoDijetEtaCMBackward1D{nullptr},
    hRecoDijetEtaCMBackward1DWeighted{nullptr},

    hRecoDijetEta1DOldPt{nullptr},
    hRecoDijetEta1DOldPtWeighted{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted{nullptr},
    hRecoDijetEtaForward1DOldPt{nullptr},
    hRecoDijetEtaForward1DOldPtWeighted{nullptr},
    hRecoDijetEtaBackward1DOldPt{nullptr},
    hRecoDijetEtaBackward1DOldPtWeighted{nullptr},

    hRecoDijetEta1DOldPtCM{nullptr},
    hRecoDijetEta1DOldPtCMWeighted{nullptr},
    hRecoEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted{nullptr},
    hRecoDijetEtaCMForward1DOldPt{nullptr},
    hRecoDijetEtaCMForward1DOldPtWeighted{nullptr},
    hRecoDijetEtaCMBackward1DOldPt{nullptr},
    hRecoDijetEtaCMBackward1DOldPtWeighted{nullptr},

    hRecoDijetEta1DOldPtBinning{nullptr},
    hRecoDijetEta1DOldPtBinningWeighted{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted{nullptr},
    hRecoDijetEtaForward1DOldPtBinning{nullptr},
    hRecoDijetEtaForward1DOldPtBinningWeighted{nullptr},
    hRecoDijetEtaBackward1DOldPtBinning{nullptr},
    hRecoDijetEtaBackward1DOldPtBinningWeighted{nullptr},

    hRecoDijetEta1DOldPtBinningCM{nullptr},
    hRecoDijetEta1DOldPtBinningCMWeighted{nullptr},
    hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted{nullptr},
    hRecoDijetEtaCMForward1DOldPtBinning{nullptr},
    hRecoDijetEtaCMForward1DOldPtBinningWeighted{nullptr},
    hRecoDijetEtaCMBackward1DOldPtBinning{nullptr},
    hRecoDijetEtaCMBackward1DOldPtBinningWeighted{nullptr},

    //
    // Ref jet histograms
    //

    hRecoLeadingJetPtOverPtHatVsLeadingJetPt{nullptr},
    hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted{nullptr},
    hRecoDijetPtOverPtHatVsDijetPt{nullptr},
    hRecoDijetPtOverPtHatVsDijetPtWeighted{nullptr},
    hRecoDijetPtAveOverPtHatVsDijetPtAve{nullptr},
    hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted{nullptr},

    hReco2RefDijetLeadSubLead{nullptr},
    hReco2RefDijetLeadSubLeadCM{nullptr},
    hRecoDijetPtEtaRefDijetPtEta{nullptr},
    hRecoDijetPtEtaRefDijetPtEtaCM{nullptr},

    hRefPtLeadPtSublead{nullptr},
    hRefEtaLeadEtaSublead{nullptr},
    hRefPtLeadPtSubleadMcReweight{nullptr},
    hRefEtaLeadEtaSubleadMcReweight{nullptr},

    hRefPtLeadPtSubleadCM{nullptr},
    hRefEtaCMLeadEtaCMSublead{nullptr},
    hRefPtLeadPtSubleadCMMcReweight{nullptr},
    hRefEtaCMLeadEtaCMSubleadMcReweight{nullptr},

    hRefDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEta{nullptr},
    hRefDijetPtVsRecoDijetPt{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt{nullptr},
    hRefDijetPtEtaPhi{nullptr},
    hRefDijetPtEtaPhiWeighted{nullptr},
    hRefDijetPtEtaForward{nullptr},
    hRefDijetPtEtaBackward{nullptr},
    hRefDijetPtEtaForwardWeighted{nullptr},
    hRefDijetPtEtaBackwardWeighted{nullptr},

    hRefDijetEtaCM{nullptr},
    hRefDijetEtaVsRecoDijetEtaCM{nullptr},
    hRefDijetPtVsRecoDijetPtCM{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM{nullptr},
    hRefDijetPtEtaPhiCM{nullptr},
    hRefDijetPtEtaPhiCMWeighted{nullptr},
    hRefDijetPtEtaCMForward{nullptr},
    hRefDijetPtEtaCMBackward{nullptr},
    hRefDijetPtEtaCMForwardWeighted{nullptr},
    hRefDijetPtEtaCMBackwardWeighted{nullptr},

    hRefDijetEta1D{nullptr},
    hRefDijetEta1DWeighted{nullptr},
    hRefEtaLeadVsEtaSubLead2D{nullptr},
    hRefEtaLeadVsEtaSubLead2DWeighted{nullptr},
    hRecoVsRefDijetEta2D{nullptr},
    hRecoVsRefDijetEta2DWeighted{nullptr},
    hRecoVsRefLeadJetEta2D{nullptr},
    hRecoVsRefLeadJetEta2DWeighted{nullptr},
    hRecoVsRefSubLeadJetEta2D{nullptr},
    hRecoVsRefSubLeadJetEta2DWeighted{nullptr},
    hRefDijetEtaForward1D{nullptr},
    hRefDijetEtaForward1DWeighted{nullptr},
    hRefDijetEtaBackward1D{nullptr},
    hRefDijetEtaBackward1DWeighted{nullptr},

    hRefDijetEta1DCM{nullptr},
    hRefDijetEta1DCMWeighted{nullptr},
    hRefEtaLeadVsEtaSubLead2DCM{nullptr},
    hRefEtaLeadVsEtaSubLead2DCMWeighted{nullptr},
    hRecoVsRefDijetEta2DCM{nullptr},
    hRecoVsRefDijetEta2DCMWeighted{nullptr},
    hRecoVsRefLeadJetEta2DCM{nullptr},
    hRecoVsRefLeadJetEta2DCMWeighted{nullptr},
    hRecoVsRefSubLeadJetEta2DCM{nullptr},
    hRecoVsRefSubLeadJetEta2DCMWeighted{nullptr},
    hRefDijetEtaCMForward1D{nullptr},
    hRefDijetEtaCMForward1DWeighted{nullptr},
    hRefDijetEtaCMBackward1D{nullptr},
    hRefDijetEtaCMBackward1DWeighted{nullptr},

    hRefDijetEta1DOldPt{nullptr},
    hRefDijetEta1DOldPtWeighted{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtWeighted{nullptr},
    hRecoVsRefDijetEta2DOldPt{nullptr},
    hRecoVsRefDijetEta2DOldPtWeighted{nullptr},
    hRecoVsRefLeadJetEta2DOldPt{nullptr},
    hRecoVsRefLeadJetEta2DOldPtWeighted{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPt{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtWeighted{nullptr},
    hRefDijetEtaForward1DOldPt{nullptr},
    hRefDijetEtaForward1DOldPtWeighted{nullptr},
    hRefDijetEtaBackward1DOldPt{nullptr},
    hRefDijetEtaBackward1DOldPtWeighted{nullptr},

    hRefDijetEta1DOldPtCM{nullptr},
    hRefDijetEta1DOldPtCMWeighted{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted{nullptr},
    hRecoVsRefDijetEta2DOldPtCM{nullptr},
    hRecoVsRefDijetEta2DOldPtCMWeighted{nullptr},
    hRecoVsRefLeadJetEta2DOldPtCM{nullptr},
    hRecoVsRefLeadJetEta2DOldPtCMWeighted{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtCM{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtCMWeighted{nullptr},
    hRefDijetEtaCMForward1DOldPt{nullptr},
    hRefDijetEtaCMForward1DOldPtWeighted{nullptr},
    hRefDijetEtaCMBackward1DOldPt{nullptr},
    hRefDijetEtaCMBackward1DOldPtWeighted{nullptr},

    hRefDijetEta1DOldPtBinning{nullptr},
    hRefDijetEta1DOldPtBinningWeighted{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted{nullptr},
    hRecoVsRefDijetEta2DOldPtBinning{nullptr},
    hRecoVsRefDijetEta2DOldPtBinningWeighted{nullptr},
    hRecoVsRefLeadJetEta2DOldPtBinning{nullptr},
    hRecoVsRefLeadJetEta2DOldPtBinningWeighted{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtBinning{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted{nullptr},
    hRefDijetEtaForward1DOldPtBinning{nullptr},
    hRefDijetEtaForward1DOldPtBinningWeighted{nullptr},
    hRefDijetEtaBackward1DOldPtBinning{nullptr},
    hRefDijetEtaBackward1DOldPtBinningWeighted{nullptr},

    hRefDijetEta1DOldPtBinningCM{nullptr},
    hRefDijetEta1DOldPtBinningCMWeighted{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted{nullptr},
    hRecoVsRefDijetEta2DOldPtBinningCM{nullptr},
    hRecoVsRefDijetEta2DOldPtBinningCMWeighted{nullptr},
    hRecoVsRefLeadJetEta2DOldPtBinningCM{nullptr},
    hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtBinningCM{nullptr},
    hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted{nullptr},
    hRefDijetEtaCMForward1DOldPtBinning{nullptr},
    hRefDijetEtaCMForward1DOldPtBinningWeighted{nullptr},
    hRefDijetEtaCMBackward1DOldPtBinning{nullptr},
    hRefDijetEtaCMBackward1DOldPtBinningWeighted{nullptr},

    //
    // Ref-selected jet histograms
    //

    hRefSelDijetEta{nullptr},
    hRefSelDijetPtEtaPhi{nullptr},
    hRefSelDijetPtEtaPhiWeighted{nullptr},
    hRefSelDijetEtaCM{nullptr},
    hRefSelDijetPtEtaPhiCM{nullptr},
    hRefSelDijetPtEtaPhiCMWeighted{nullptr},

    hRefSelDijetEta1D{nullptr},
    hRefSelDijetEta1DWeighted{nullptr},
    hRefSelRecoDijetEta1D{nullptr},
    hRefSelRecoDijetEta1DWeighted{nullptr},
    hRefSelEtaLeadVsEtaSubLead2D{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DWeighted{nullptr},
    hRefSelDijetEtaForward1D{nullptr},
    hRefSelDijetEtaForward1DWeighted{nullptr},
    hRefSelDijetEtaBackward1D{nullptr},
    hRefSelDijetEtaBackward1DWeighted{nullptr},

    hRefSelDijetEta1DCM{nullptr},
    hRefSelDijetEta1DCMWeighted{nullptr},
    hRefSelRecoDijetEta1DCM{nullptr},
    hRefSelRecoDijetEta1DCMWeighted{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DCM{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DCMWeighted{nullptr},
    hRefSelDijetEtaCMForward1D{nullptr},
    hRefSelDijetEtaCMForward1DWeighted{nullptr},
    hRefSelDijetEtaCMBackward1D{nullptr},
    hRefSelDijetEtaCMBackward1DWeighted{nullptr},

    hRefSelDijetEta1DOldPt{nullptr},
    hRefSelDijetEta1DOldPtWeighted{nullptr},
    hRefSelRecoDijetEta1DOldPt{nullptr},
    hRefSelRecoDijetEta1DOldPtWeighted{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPt{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted{nullptr},
    hRefSelDijetEtaForward1DOldPt{nullptr},
    hRefSelDijetEtaForward1DOldPtWeighted{nullptr},
    hRefSelDijetEtaBackward1DOldPt{nullptr},
    hRefSelDijetEtaBackward1DOldPtWeighted{nullptr},

    hRefSelDijetEta1DOldPtCM{nullptr},
    hRefSelDijetEta1DOldPtCMWeighted{nullptr},
    hRefSelRecoDijetEta1DOldPtCM{nullptr},
    hRefSelRecoDijetEta1DOldPtCMWeighted{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtCM{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted{nullptr},
    hRefSelDijetEtaCMForward1DOldPt{nullptr},
    hRefSelDijetEtaCMForward1DOldPtWeighted{nullptr},
    hRefSelDijetEtaCMBackward1DOldPt{nullptr},
    hRefSelDijetEtaCMBackward1DOldPtWeighted{nullptr},

    hRefSelDijetEta1DOldPtBinning{nullptr},
    hRefSelDijetEta1DOldPtBinningWeighted{nullptr},
    hRefSelRecoDijetEta1DOldPtBinning{nullptr},
    hRefSelRecoDijetEta1DOldPtBinningWeighted{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtBinning{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted{nullptr},
    hRefSelDijetEtaForward1DOldPtBinning{nullptr},
    hRefSelDijetEtaForward1DOldPtBinningWeighted{nullptr},
    hRefSelDijetEtaBackward1DOldPtBinning{nullptr},
    hRefSelDijetEtaBackward1DOldPtBinningWeighted{nullptr},

    hRefSelDijetEta1DOldPtBinningCM{nullptr},
    hRefSelDijetEta1DOldPtBinningCMWeighted{nullptr},
    hRefSelRecoDijetEta1DOldPtBinningCM{nullptr},
    hRefSelRecoDijetEta1DOldPtBinningCMWeighted{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM{nullptr},
    hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted{nullptr},
    hRefSelDijetEtaCMForward1DOldPtBinning{nullptr},
    hRefSelDijetEtaCMForward1DOldPtBinningWeighted{nullptr},
    hRefSelDijetEtaCMBackward1DOldPtBinning{nullptr},
    hRefSelDijetEtaCMBackward1DOldPtBinningWeighted{nullptr},

    //
    // Variables
    //
    fIsMc{false}, 
    fPtBins{150}, fPtRange{5., 1505.}, 
    fEtaBins{36}, fEtaRange{-3.6, 3.6},
    fPhiBins{8}, fPhiRange{-TMath::Pi(), TMath::Pi()},
    fDijetPtBins{196}, fDijetPtRange{20., 1000.},
    fDijetEtaBins{36}, fDijetEtaRange{-3.6, 3.6},
    fDijetDphiBins{8}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},
    fPtHatBins{100}, fPtHatRange{15., 1015.},
    fUseVariableBinning{true}
{ 

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

    //
    // Event histograms
    //

    if (hVz) { delete hVz; hVz = nullptr; }
    if (hVzWeighted) { delete hVzWeighted; hVzWeighted = nullptr; }
    if (hPtHat) { delete hPtHat; hPtHat = nullptr; }
    if (hPtHatWeighted) { delete hPtHatWeighted; hPtHatWeighted = nullptr; }
    if (hHiBin) { delete hHiBin; hHiBin = nullptr; }
    if (hHiBinWeighted) { delete hHiBinWeighted; hHiBinWeighted = nullptr; }

    if (hVzGenDijetLab) { delete hVzGenDijetLab; hVzGenDijetLab = nullptr; }
    if (hVzGenDijetLabWeighted) { delete hVzGenDijetLabWeighted; hVzGenDijetLabWeighted = nullptr; }
    if (hHiBinGenDijetLab) { delete hHiBinGenDijetLab; hHiBinGenDijetLab = nullptr; }
    if (hHiBinGenDijetLabWeighted) { delete hHiBinGenDijetLabWeighted; hHiBinGenDijetLabWeighted = nullptr; }

    if (hVzGenDijetCM) { delete hVzGenDijetCM; hVzGenDijetCM = nullptr; }
    if (hVzGenDijetCMWeighted) { delete hVzGenDijetCMWeighted; hVzGenDijetCMWeighted = nullptr; }
    if (hHiBinGenDijetCM) { delete hHiBinGenDijetCM; hHiBinGenDijetCM = nullptr; }
    if (hHiBinGenDijetCMWeighted) { delete hHiBinGenDijetCMWeighted; hHiBinGenDijetCMWeighted = nullptr; }

    if (hVzRecoDijetLab) { delete hVzRecoDijetLab; hVzRecoDijetLab = nullptr; }
    if (hVzRecoDijetLabWeighted) { delete hVzRecoDijetLabWeighted; hVzRecoDijetLabWeighted = nullptr; }
    if (hHiBinRecoDijetLab) { delete hHiBinRecoDijetLab; hHiBinRecoDijetLab = nullptr; }
    if (hHiBinRecoDijetLabWeighted) { delete hHiBinRecoDijetLabWeighted; hHiBinRecoDijetLabWeighted = nullptr; }

    if (hVzRecoDijetCM) { delete hVzRecoDijetCM; hVzRecoDijetCM = nullptr; }
    if (hVzRecoDijetCMWeighted) { delete hVzRecoDijetCMWeighted; hVzRecoDijetCMWeighted = nullptr; }
    if (hHiBinRecoDijetCM) { delete hHiBinRecoDijetCM; hHiBinRecoDijetCM = nullptr; }
    if (hHiBinRecoDijetCMWeighted) { delete hHiBinRecoDijetCMWeighted; hHiBinRecoDijetCMWeighted = nullptr; }

    if (hVzRefSelDijetLab) { delete hVzRefSelDijetLab; hVzRefSelDijetLab = nullptr; }
    if (hVzRefSelDijetLabWeighted) { delete hVzRefSelDijetLabWeighted; hVzRefSelDijetLabWeighted = nullptr; }
    if (hHiBinRefSelDijetLab) { delete hHiBinRefSelDijetLab; hHiBinRefSelDijetLab = nullptr; }
    if (hHiBinRefSelDijetLabWeighted) { delete hHiBinRefSelDijetLabWeighted; hHiBinRefSelDijetLabWeighted = nullptr; }

    if (hVzRefSelDijetCM) { delete hVzRefSelDijetCM; hVzRefSelDijetCM = nullptr; }
    if (hVzRefSelDijetCMWeighted) { delete hVzRefSelDijetCMWeighted; hVzRefSelDijetCMWeighted = nullptr; }
    if (hHiBinRefSelDijetCM) { delete hHiBinRefSelDijetCM; hHiBinRefSelDijetCM = nullptr; }
    if (hHiBinRefSelDijetCMWeighted) { delete hHiBinRefSelDijetCMWeighted; hHiBinRefSelDijetCMWeighted = nullptr; }

    //
    // Gen jet histograms
    //

    if (hGenJetCollectionSize) { delete hGenJetCollectionSize; hGenJetCollectionSize = nullptr; }
    if (hGenVsRecoJetCollectionSize) { delete hGenVsRecoJetCollectionSize; hGenVsRecoJetCollectionSize = nullptr; }

    if (hGenLeadingJetPtOverPtHatVsLeadingJetPt) { delete hGenLeadingJetPtOverPtHatVsLeadingJetPt; hGenLeadingJetPtOverPtHatVsLeadingJetPt = nullptr; }
    if (hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted) { delete hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted; hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted = nullptr; }
    if (hGenDijetPtOverPtHatVsDijetPt) { delete hGenDijetPtOverPtHatVsDijetPt; hGenDijetPtOverPtHatVsDijetPt = nullptr; }
    if (hGenDijetPtOverPtHatVsDijetPtWeighted) { delete hGenDijetPtOverPtHatVsDijetPtWeighted; hGenDijetPtOverPtHatVsDijetPtWeighted = nullptr; }
    if (hGenDijetPtAveOverPtHatVsDijetPtAve) { delete hGenDijetPtAveOverPtHatVsDijetPtAve; hGenDijetPtAveOverPtHatVsDijetPtAve = nullptr; }
    if (hGenDijetPtAveOverPtHatVsDijetPtAveWeighted) { delete hGenDijetPtAveOverPtHatVsDijetPtAveWeighted; hGenDijetPtAveOverPtHatVsDijetPtAveWeighted = nullptr; }

    if (hGenInclusiveDijetDetaCM) { delete hGenInclusiveDijetDetaCM; hGenInclusiveDijetDetaCM = nullptr; }
    if (hGenInclusiveDijetDetaCMWeighted) { delete hGenInclusiveDijetDetaCMWeighted; hGenInclusiveDijetDetaCMWeighted = nullptr; }
    if (hGenInclusiveDijetDetaCMPt) { delete hGenInclusiveDijetDetaCMPt; hGenInclusiveDijetDetaCMPt = nullptr; }
    if (hGenInclusiveDijetDetaCMPtWeighted) { delete hGenInclusiveDijetDetaCMPtWeighted; hGenInclusiveDijetDetaCMPtWeighted = nullptr; }
    if (hGenInclusiveDijetEtaDetaCMPt) { delete hGenInclusiveDijetEtaDetaCMPt; hGenInclusiveDijetEtaDetaCMPt = nullptr; }
    if (hGenInclusiveDijetEtaDetaCMPtWeighted) { delete hGenInclusiveDijetEtaDetaCMPtWeighted; hGenInclusiveDijetEtaDetaCMPtWeighted = nullptr; }
    if (hGenInclusiveDijetXPb) { delete hGenInclusiveDijetXPb; hGenInclusiveDijetXPb = nullptr; }
    if (hGenInclusiveDijetXPbWeighted) { delete hGenInclusiveDijetXPbWeighted; hGenInclusiveDijetXPbWeighted = nullptr; }
    if (hGenInclusiveDijetXp) { delete hGenInclusiveDijetXp; hGenInclusiveDijetXp = nullptr; }
    if (hGenInclusiveDijetXpWeighted) { delete hGenInclusiveDijetXpWeighted; hGenInclusiveDijetXpWeighted = nullptr; }
    if (hGenInclusiveDijetXPbOverXp) { delete hGenInclusiveDijetXPbOverXp; hGenInclusiveDijetXPbOverXp = nullptr; }
    if (hGenInclusiveDijetXPbOverXpWeighted) { delete hGenInclusiveDijetXPbOverXpWeighted; hGenInclusiveDijetXPbOverXpWeighted = nullptr; }
    if (hGenInclusiveDijetXPbOverXpEta) { delete hGenInclusiveDijetXPbOverXpEta; hGenInclusiveDijetXPbOverXpEta = nullptr; }
    if (hGenInclusiveDijetXPbOverXpEtaWeighted) { delete hGenInclusiveDijetXPbOverXpEtaWeighted; hGenInclusiveDijetXPbOverXpEtaWeighted = nullptr; }

    if (hGenSelectedDijetDetaCM) { delete hGenSelectedDijetDetaCM; hGenSelectedDijetDetaCM = nullptr; }
    if (hGenSelectedDijetDetaCMWeighted) { delete hGenSelectedDijetDetaCMWeighted; hGenSelectedDijetDetaCMWeighted = nullptr; }
    if (hGenSelectedDijetDetaCMPt) { delete hGenSelectedDijetDetaCMPt; hGenSelectedDijetDetaCMPt = nullptr; }
    if (hGenSelectedDijetDetaCMPtWeighted) { delete hGenSelectedDijetDetaCMPtWeighted; hGenSelectedDijetDetaCMPtWeighted = nullptr; }
    if (hGenSelectedDijetEtaDetaCMPt) { delete hGenSelectedDijetEtaDetaCMPt; hGenSelectedDijetEtaDetaCMPt = nullptr; }
    if (hGenSelectedDijetEtaDetaCMPtWeighted) { delete hGenSelectedDijetEtaDetaCMPtWeighted; hGenSelectedDijetEtaDetaCMPtWeighted = nullptr; }
    if (hGenSelectedDijetXPb) { delete hGenSelectedDijetXPb; hGenSelectedDijetXPb = nullptr; }
    if (hGenSelectedDijetXPbWeighted) { delete hGenSelectedDijetXPbWeighted; hGenSelectedDijetXPbWeighted = nullptr; }
    if (hGenSelectedDijetXp) { delete hGenSelectedDijetXp; hGenSelectedDijetXp = nullptr; }
    if (hGenSelectedDijetXpWeighted) { delete hGenSelectedDijetXpWeighted; hGenSelectedDijetXpWeighted = nullptr; }
    if (hGenSelectedDijetXPbOverXp) { delete hGenSelectedDijetXPbOverXp; hGenSelectedDijetXPbOverXp = nullptr; }
    if (hGenSelectedDijetXPbOverXpWeighted) { delete hGenSelectedDijetXPbOverXpWeighted; hGenSelectedDijetXPbOverXpWeighted = nullptr; }
    if (hGenSelectedDijetXPbOverXpEta) { delete hGenSelectedDijetXPbOverXpEta; hGenSelectedDijetXPbOverXpEta = nullptr; }
    if (hGenSelectedDijetXPbOverXpEtaWeighted) { delete hGenSelectedDijetXPbOverXpEtaWeighted; hGenSelectedDijetXPbOverXpEtaWeighted = nullptr; }

    if (hGenDijetLeadSubLead) { delete hGenDijetLeadSubLead; hGenDijetLeadSubLead = nullptr; }
    if (hGenDijetLeadSubLeadCM) { delete hGenDijetLeadSubLeadCM; hGenDijetLeadSubLeadCM = nullptr; }

    if (hGenPtLeadPtSublead) { delete hGenPtLeadPtSublead; hGenPtLeadPtSublead = nullptr; }
    if (hGenEtaLeadEtaSublead) { delete hGenEtaLeadEtaSublead; hGenEtaLeadEtaSublead = nullptr; }
    if (hGenDijetEta) { delete hGenDijetEta; hGenDijetEta = nullptr; }
    if (hGenDijetPtEtaPhi) { delete hGenDijetPtEtaPhi; hGenDijetPtEtaPhi = nullptr; }
    if (hGenDijetPtEtaPhiWeighted) { delete hGenDijetPtEtaPhiWeighted; hGenDijetPtEtaPhiWeighted = nullptr; }

    if (hGenPtLeadPtSubleadCM) { delete hGenPtLeadPtSubleadCM; hGenPtLeadPtSubleadCM = nullptr; }
    if (hGenEtaCMLeadEtaCMSublead) { delete hGenEtaCMLeadEtaCMSublead; hGenEtaCMLeadEtaCMSublead = nullptr; }
    if (hGenDijetEtaCM) { delete hGenDijetEtaCM; hGenDijetEtaCM = nullptr; }
    if (hGenDijetPtEtaPhiCM) { delete hGenDijetPtEtaPhiCM; hGenDijetPtEtaPhiCM = nullptr; }
    if (hGenDijetPtEtaPhiCMWeighted) { delete hGenDijetPtEtaPhiCMWeighted; hGenDijetPtEtaPhiCMWeighted = nullptr; }

    if (hGenDijetPtEtaForward) { delete hGenDijetPtEtaForward; hGenDijetPtEtaForward = nullptr; }
    if (hGenDijetPtEtaBackward) { delete hGenDijetPtEtaBackward; hGenDijetPtEtaBackward = nullptr; }
    if (hGenDijetPtEtaCMForward) { delete hGenDijetPtEtaCMForward; hGenDijetPtEtaCMForward = nullptr; }
    if (hGenDijetPtEtaCMBackward) { delete hGenDijetPtEtaCMBackward; hGenDijetPtEtaCMBackward = nullptr; }
    if (hGenDijetPtEtaForwardWeighted) { delete hGenDijetPtEtaForwardWeighted; hGenDijetPtEtaForwardWeighted = nullptr; }
    if (hGenDijetPtEtaBackwardWeighted) { delete hGenDijetPtEtaBackwardWeighted; hGenDijetPtEtaBackwardWeighted = nullptr; }
    if (hGenDijetPtEtaCMForwardWeighted) { delete hGenDijetPtEtaCMForwardWeighted; hGenDijetPtEtaCMForwardWeighted = nullptr; }
    if (hGenDijetPtEtaCMBackwardWeighted) { delete hGenDijetPtEtaCMBackwardWeighted; hGenDijetPtEtaCMBackwardWeighted = nullptr; }


    for (int i = 0; i < 16; ++i) {
        if (hGenDijetEta1D[i]) { delete hGenDijetEta1D[i]; hGenDijetEta1D[i] = nullptr; }
        if (hGenDijetEta1DWeighted[i]) { delete hGenDijetEta1DWeighted[i]; hGenDijetEta1DWeighted[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2D[i]) { delete hGenDijetEtaLeadVsEtaSubLead2D[i]; hGenDijetEtaLeadVsEtaSubLead2D[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DWeighted[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DWeighted[i]; hGenDijetEtaLeadVsEtaSubLead2DWeighted[i] = nullptr; }
        if (hGenDijetEtaForward1D[i]) { delete hGenDijetEtaForward1D[i]; hGenDijetEtaForward1D[i] = nullptr; }
        if (hGenDijetEtaForward1DWeighted[i]) { delete hGenDijetEtaForward1DWeighted[i]; hGenDijetEtaForward1DWeighted[i] = nullptr; }
        if (hGenDijetEtaBackward1D[i]) { delete hGenDijetEtaBackward1D[i]; hGenDijetEtaBackward1D[i] = nullptr; }
        if (hGenDijetEtaBackward1DWeighted[i]) { delete hGenDijetEtaBackward1DWeighted[i]; hGenDijetEtaBackward1DWeighted[i] = nullptr; }

        if (hGenDijetEta1DCM[i]) { delete hGenDijetEta1DCM[i]; hGenDijetEta1DCM[i] = nullptr; }
        if (hGenDijetEta1DCMWeighted[i]) { delete hGenDijetEta1DCMWeighted[i]; hGenDijetEta1DCMWeighted[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DCM[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DCM[i]; hGenDijetEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i]; hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr; }
        if (hGenDijetEtaCMForward1D[i]) { delete hGenDijetEtaCMForward1D[i]; hGenDijetEtaCMForward1D[i] = nullptr; }
        if (hGenDijetEtaCMForward1DWeighted[i]) { delete hGenDijetEtaCMForward1DWeighted[i]; hGenDijetEtaCMForward1DWeighted[i] = nullptr; }
        if (hGenDijetEtaCMBackward1D[i]) { delete hGenDijetEtaCMBackward1D[i]; hGenDijetEtaCMBackward1D[i] = nullptr; }
        if (hGenDijetEtaCMBackward1DWeighted[i]) { delete hGenDijetEtaCMBackward1DWeighted[i]; hGenDijetEtaCMBackward1DWeighted[i] = nullptr; }
    }

    for (int i = 0; i < 6; ++i) {
        if (hGenDijetEta1DOldPt[i]) { delete hGenDijetEta1DOldPt[i]; hGenDijetEta1DOldPt[i] = nullptr; }
        if (hGenDijetEta1DOldPtWeighted[i]) { delete hGenDijetEta1DOldPtWeighted[i]; hGenDijetEta1DOldPtWeighted[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr; }
        if (hGenDijetEtaForward1DOldPt[i]) { delete hGenDijetEtaForward1DOldPt[i]; hGenDijetEtaForward1DOldPt[i] = nullptr; }
        if (hGenDijetEtaForward1DOldPtWeighted[i]) { delete hGenDijetEtaForward1DOldPtWeighted[i]; hGenDijetEtaForward1DOldPtWeighted[i] = nullptr; }
        if (hGenDijetEtaBackward1DOldPt[i]) { delete hGenDijetEtaBackward1DOldPt[i]; hGenDijetEtaBackward1DOldPt[i] = nullptr; }
        if (hGenDijetEtaBackward1DOldPtWeighted[i]) { delete hGenDijetEtaBackward1DOldPtWeighted[i]; hGenDijetEtaBackward1DOldPtWeighted[i] = nullptr; }

        if (hGenDijetEta1DOldPtCM[i]) { delete hGenDijetEta1DOldPtCM[i]; hGenDijetEta1DOldPtCM[i] = nullptr; }
        if (hGenDijetEta1DOldPtCMWeighted[i]) { delete hGenDijetEta1DOldPtCMWeighted[i]; hGenDijetEta1DOldPtCMWeighted[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr; }
        if (hGenDijetEtaCMForward1DOldPt[i]) { delete hGenDijetEtaCMForward1DOldPt[i]; hGenDijetEtaCMForward1DOldPt[i] = nullptr; }
        if (hGenDijetEtaCMForward1DOldPtWeighted[i]) { delete hGenDijetEtaCMForward1DOldPtWeighted[i]; hGenDijetEtaCMForward1DOldPtWeighted[i] = nullptr; }
        if (hGenDijetEtaCMBackward1DOldPt[i]) { delete hGenDijetEtaCMBackward1DOldPt[i]; hGenDijetEtaCMBackward1DOldPt[i] = nullptr; }
        if (hGenDijetEtaCMBackward1DOldPtWeighted[i]) { delete hGenDijetEtaCMBackward1DOldPtWeighted[i]; hGenDijetEtaCMBackward1DOldPtWeighted[i] = nullptr; }

        if (hGenDijetEta1DOldPtBinning[i]) { delete hGenDijetEta1DOldPtBinning[i]; hGenDijetEta1DOldPtBinning[i] = nullptr; }
        if (hGenDijetEta1DOldPtBinningWeighted[i]) { delete hGenDijetEta1DOldPtBinningWeighted[i]; hGenDijetEta1DOldPtBinningWeighted[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr; }
        if (hGenDijetEtaForward1DOldPtBinning[i]) { delete hGenDijetEtaForward1DOldPtBinning[i]; hGenDijetEtaForward1DOldPtBinning[i] = nullptr; }
        if (hGenDijetEtaForward1DOldPtBinningWeighted[i]) { delete hGenDijetEtaForward1DOldPtBinningWeighted[i]; hGenDijetEtaForward1DOldPtBinningWeighted[i] = nullptr; }
        if (hGenDijetEtaBackward1DOldPtBinning[i]) { delete hGenDijetEtaBackward1DOldPtBinning[i]; hGenDijetEtaBackward1DOldPtBinning[i] = nullptr; }
        if (hGenDijetEtaBackward1DOldPtBinningWeighted[i]) { delete hGenDijetEtaBackward1DOldPtBinningWeighted[i]; hGenDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr; }

        if (hGenDijetEta1DOldPtBinningCM[i]) { delete hGenDijetEta1DOldPtBinningCM[i]; hGenDijetEta1DOldPtBinningCM[i] = nullptr; }
        if (hGenDijetEta1DOldPtBinningCMWeighted[i]) { delete hGenDijetEta1DOldPtBinningCMWeighted[i]; hGenDijetEta1DOldPtBinningCMWeighted[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
        if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) { delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]; hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr; }
        if (hGenDijetEtaCMForward1DOldPtBinning[i]) { delete hGenDijetEtaCMForward1DOldPtBinning[i]; hGenDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
        if (hGenDijetEtaCMForward1DOldPtBinningWeighted[i]) { delete hGenDijetEtaCMForward1DOldPtBinningWeighted[i]; hGenDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr; }
        if (hGenDijetEtaCMBackward1DOldPtBinning[i]) { delete hGenDijetEtaCMBackward1DOldPtBinning[i]; hGenDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
        if (hGenDijetEtaCMBackward1DOldPtBinningWeighted[i]) { delete hGenDijetEtaCMBackward1DOldPtBinningWeighted[i]; hGenDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr; }
    }

    //
    // Reco jet histograms
    //

    if (hRecoJetCollectionSize) { delete hRecoJetCollectionSize; hRecoJetCollectionSize = nullptr; }

    if (hRecoDijetLeadSubLead) { delete hRecoDijetLeadSubLead; hRecoDijetLeadSubLead = nullptr; }
    if (hRecoDijetLeadSubLeadCM) { delete hRecoDijetLeadSubLeadCM; hRecoDijetLeadSubLeadCM = nullptr; }

    if (hRecoPtLeadPtSublead) { delete hRecoPtLeadPtSublead; hRecoPtLeadPtSublead = nullptr; }
    if (hRecoEtaLeadEtaSublead) { delete hRecoEtaLeadEtaSublead; hRecoEtaLeadEtaSublead = nullptr; }
    if (hRecoPtLeadPtSubleadMcReweight) { delete hRecoPtLeadPtSubleadMcReweight; hRecoPtLeadPtSubleadMcReweight = nullptr; }
    if (hRecoEtaLeadEtaSubleadMcReweight) { delete hRecoEtaLeadEtaSubleadMcReweight; hRecoEtaLeadEtaSubleadMcReweight = nullptr; }
    if (hRecoDijetEta) { delete hRecoDijetEta; hRecoDijetEta = nullptr; }
    if (hRecoDijetPtEta) { delete hRecoDijetPtEta; hRecoDijetPtEta = nullptr; }
    if (hRecoDijetPtEtaForward) { delete hRecoDijetPtEtaForward; hRecoDijetPtEtaForward = nullptr; }
    if (hRecoDijetPtEtaBackward) { delete hRecoDijetPtEtaBackward; hRecoDijetPtEtaBackward = nullptr; }
    if (hRecoDijetPtEtaForwardWeighted) { delete hRecoDijetPtEtaForwardWeighted; hRecoDijetPtEtaForwardWeighted = nullptr; }
    if (hRecoDijetPtEtaBackwardWeighted) { delete hRecoDijetPtEtaBackwardWeighted; hRecoDijetPtEtaBackwardWeighted = nullptr; }
    if (hRecoDijetPtEtaPhi) { delete hRecoDijetPtEtaPhi; hRecoDijetPtEtaPhi = nullptr; }
    if (hRecoDijetPtEtaPhiWeighted) { delete hRecoDijetPtEtaPhiWeighted; hRecoDijetPtEtaPhiWeighted = nullptr; }

    if (hRecoPtLeadPtSubleadCM) { delete hRecoPtLeadPtSubleadCM; hRecoPtLeadPtSubleadCM = nullptr; }
    if (hRecoEtaCMLeadEtaCMSublead) { delete hRecoEtaCMLeadEtaCMSublead; hRecoEtaCMLeadEtaCMSublead = nullptr; }
    if (hRecoPtLeadPtSubleadCMMcReweight) { delete hRecoPtLeadPtSubleadCMMcReweight; hRecoPtLeadPtSubleadCMMcReweight = nullptr; }
    if (hRecoEtaCMLeadEtaCMSubleadMcReweight) { delete hRecoEtaCMLeadEtaCMSubleadMcReweight; hRecoEtaCMLeadEtaCMSubleadMcReweight = nullptr; }
    if (hRecoDijetEtaCM) { delete hRecoDijetEtaCM; hRecoDijetEtaCM = nullptr; }
    if (hRecoDijetPtEtaCM) { delete hRecoDijetPtEtaCM; hRecoDijetPtEtaCM = nullptr; }
    if (hRecoDijetPtEtaCMForward) { delete hRecoDijetPtEtaCMForward; hRecoDijetPtEtaCMForward = nullptr; }
    if (hRecoDijetPtEtaCMBackward) { delete hRecoDijetPtEtaCMBackward; hRecoDijetPtEtaCMBackward = nullptr; }
    if (hRecoDijetPtEtaCMForwardWeighted) { delete hRecoDijetPtEtaCMForwardWeighted; hRecoDijetPtEtaCMForwardWeighted = nullptr; }
    if (hRecoDijetPtEtaCMBackwardWeighted) { delete hRecoDijetPtEtaCMBackwardWeighted; hRecoDijetPtEtaCMBackwardWeighted = nullptr; }
    if (hRecoDijetPtEtaPhiCM) { delete hRecoDijetPtEtaPhiCM; hRecoDijetPtEtaPhiCM = nullptr; }
    if (hRecoDijetPtEtaPhiCMWeighted) { delete hRecoDijetPtEtaPhiCMWeighted; hRecoDijetPtEtaPhiCMWeighted = nullptr; }

    for (int i = 0; i < 16; ++i) {
        if (hRecoDijetEta1D[i]) { delete hRecoDijetEta1D[i]; hRecoDijetEta1D[i] = nullptr; }
        if (hRecoDijetEta1DWeighted[i]) { delete hRecoDijetEta1DWeighted[i]; hRecoDijetEta1DWeighted[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2D[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2D[i]; hRecoDijetEtaLeadVsEtaSubLead2D[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i]; hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i] = nullptr; }
        if (hRecoDijetEtaForward1D[i]) { delete hRecoDijetEtaForward1D[i]; hRecoDijetEtaForward1D[i] = nullptr; }
        if (hRecoDijetEtaForward1DWeighted[i]) { delete hRecoDijetEtaForward1DWeighted[i]; hRecoDijetEtaForward1DWeighted[i] = nullptr; }
        if (hRecoDijetEtaBackward1D[i]) { delete hRecoDijetEtaBackward1D[i]; hRecoDijetEtaBackward1D[i] = nullptr; }
        if (hRecoDijetEtaBackward1DWeighted[i]) { delete hRecoDijetEtaBackward1DWeighted[i]; hRecoDijetEtaBackward1DWeighted[i] = nullptr; }

        if (hRecoDijetEta1DCM[i]) { delete hRecoDijetEta1DCM[i]; hRecoDijetEta1DCM[i] = nullptr; }
        if (hRecoDijetEta1DCMWeighted[i]) { delete hRecoDijetEta1DCMWeighted[i]; hRecoDijetEta1DCMWeighted[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DCM[i]) { delete hRecoEtaLeadVsEtaSubLead2DCM[i]; hRecoEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DCMWeighted[i]) { delete hRecoEtaLeadVsEtaSubLead2DCMWeighted[i]; hRecoEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr; }
        if (hRecoDijetEtaCMForward1D[i]) { delete hRecoDijetEtaCMForward1D[i]; hRecoDijetEtaCMForward1D[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DWeighted[i]) { delete hRecoDijetEtaCMForward1DWeighted[i]; hRecoDijetEtaCMForward1DWeighted[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1D[i]) { delete hRecoDijetEtaCMBackward1D[i]; hRecoDijetEtaCMBackward1D[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DWeighted[i]) { delete hRecoDijetEtaCMBackward1DWeighted[i]; hRecoDijetEtaCMBackward1DWeighted[i] = nullptr; }
    }

    for (int i = 0; i < 6; ++i) {
        if (hRecoDijetEta1DOldPt[i]) { delete hRecoDijetEta1DOldPt[i]; hRecoDijetEta1DOldPt[i] = nullptr; }
        if (hRecoDijetEta1DOldPtWeighted[i]) { delete hRecoDijetEta1DOldPtWeighted[i]; hRecoDijetEta1DOldPtWeighted[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]; hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]; hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr; }
        if (hRecoDijetEtaForward1DOldPt[i]) { delete hRecoDijetEtaForward1DOldPt[i]; hRecoDijetEtaForward1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaForward1DOldPtWeighted[i]) { delete hRecoDijetEtaForward1DOldPtWeighted[i]; hRecoDijetEtaForward1DOldPtWeighted[i] = nullptr; }
        if (hRecoDijetEtaBackward1DOldPt[i]) { delete hRecoDijetEtaBackward1DOldPt[i]; hRecoDijetEtaBackward1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaBackward1DOldPtWeighted[i]) { delete hRecoDijetEtaBackward1DOldPtWeighted[i]; hRecoDijetEtaBackward1DOldPtWeighted[i] = nullptr; }

        if (hRecoDijetEta1DOldPtCM[i]) { delete hRecoDijetEta1DOldPtCM[i]; hRecoDijetEta1DOldPtCM[i] = nullptr; }
        if (hRecoDijetEta1DOldPtCMWeighted[i]) { delete hRecoDijetEta1DOldPtCMWeighted[i]; hRecoDijetEta1DOldPtCMWeighted[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]; hRecoEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) { delete hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]; hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DOldPt[i]) { delete hRecoDijetEtaCMForward1DOldPt[i]; hRecoDijetEtaCMForward1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DOldPtWeighted[i]) { delete hRecoDijetEtaCMForward1DOldPtWeighted[i]; hRecoDijetEtaCMForward1DOldPtWeighted[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DOldPt[i]) { delete hRecoDijetEtaCMBackward1DOldPt[i]; hRecoDijetEtaCMBackward1DOldPt[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DOldPtWeighted[i]) { delete hRecoDijetEtaCMBackward1DOldPtWeighted[i]; hRecoDijetEtaCMBackward1DOldPtWeighted[i] = nullptr; }

        if (hRecoDijetEta1DOldPtBinning[i]) { delete hRecoDijetEta1DOldPtBinning[i]; hRecoDijetEta1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEta1DOldPtBinningWeighted[i]) { delete hRecoDijetEta1DOldPtBinningWeighted[i]; hRecoDijetEta1DOldPtBinningWeighted[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]; hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) { delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]; hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr; }
        if (hRecoDijetEtaForward1DOldPtBinning[i]) { delete hRecoDijetEtaForward1DOldPtBinning[i]; hRecoDijetEtaForward1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaForward1DOldPtBinningWeighted[i]) { delete hRecoDijetEtaForward1DOldPtBinningWeighted[i]; hRecoDijetEtaForward1DOldPtBinningWeighted[i] = nullptr; }
        if (hRecoDijetEtaBackward1DOldPtBinning[i]) { delete hRecoDijetEtaBackward1DOldPtBinning[i]; hRecoDijetEtaBackward1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaBackward1DOldPtBinningWeighted[i]) { delete hRecoDijetEtaBackward1DOldPtBinningWeighted[i]; hRecoDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr; }

        if (hRecoDijetEta1DOldPtBinningCM[i]) { delete hRecoDijetEta1DOldPtBinningCM[i]; hRecoDijetEta1DOldPtBinningCM[i] = nullptr; }
        if (hRecoDijetEta1DOldPtBinningCMWeighted[i]) { delete hRecoDijetEta1DOldPtBinningCMWeighted[i]; hRecoDijetEta1DOldPtBinningCMWeighted[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
        if (hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) { delete hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]; hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DOldPtBinning[i]) { delete hRecoDijetEtaCMForward1DOldPtBinning[i]; hRecoDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaCMForward1DOldPtBinningWeighted[i]) { delete hRecoDijetEtaCMForward1DOldPtBinningWeighted[i]; hRecoDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DOldPtBinning[i]) { delete hRecoDijetEtaCMBackward1DOldPtBinning[i]; hRecoDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
        if (hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i]) { delete hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i]; hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr; }
    }

    //
    // Ref jet histograms
    //

    if (hRecoLeadingJetPtOverPtHatVsLeadingJetPt) { delete hRecoLeadingJetPtOverPtHatVsLeadingJetPt; hRecoLeadingJetPtOverPtHatVsLeadingJetPt = nullptr; }
    if (hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted) { delete hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted; hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted = nullptr; }
    if (hRecoDijetPtOverPtHatVsDijetPt) { delete hRecoDijetPtOverPtHatVsDijetPt; hRecoDijetPtOverPtHatVsDijetPt = nullptr; }
    if (hRecoDijetPtOverPtHatVsDijetPtWeighted) { delete hRecoDijetPtOverPtHatVsDijetPtWeighted; hRecoDijetPtOverPtHatVsDijetPtWeighted = nullptr; }
    if (hRecoDijetPtAveOverPtHatVsDijetPtAve) { delete hRecoDijetPtAveOverPtHatVsDijetPtAve; hRecoDijetPtAveOverPtHatVsDijetPtAve = nullptr; }
    if (hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted) { delete hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted; hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted = nullptr; }

    if (hReco2RefDijetLeadSubLead) { delete hReco2RefDijetLeadSubLead; hReco2RefDijetLeadSubLead = nullptr; }
    if (hReco2RefDijetLeadSubLeadCM) { delete hReco2RefDijetLeadSubLeadCM; hReco2RefDijetLeadSubLeadCM = nullptr; }
    if (hRecoDijetPtEtaRefDijetPtEta) { delete hRecoDijetPtEtaRefDijetPtEta; hRecoDijetPtEtaRefDijetPtEta = nullptr; }
    if (hRecoDijetPtEtaRefDijetPtEtaCM) { delete hRecoDijetPtEtaRefDijetPtEtaCM; hRecoDijetPtEtaRefDijetPtEtaCM = nullptr; }

    if (hRefPtLeadPtSublead) { delete hRefPtLeadPtSublead; hRefPtLeadPtSublead = nullptr; }
    if (hRefEtaLeadEtaSublead) { delete hRefEtaLeadEtaSublead; hRefEtaLeadEtaSublead = nullptr; }
    if (hRefPtLeadPtSubleadMcReweight) { delete hRefPtLeadPtSubleadMcReweight; hRefPtLeadPtSubleadMcReweight = nullptr; }
    if (hRefEtaLeadEtaSubleadMcReweight) { delete hRefEtaLeadEtaSubleadMcReweight; hRefEtaLeadEtaSubleadMcReweight = nullptr; }

    if (hRefPtLeadPtSubleadCM) { delete hRefPtLeadPtSubleadCM; hRefPtLeadPtSubleadCM = nullptr; }
    if (hRefEtaCMLeadEtaCMSublead) { delete hRefEtaCMLeadEtaCMSublead; hRefEtaCMLeadEtaCMSublead = nullptr; }
    if (hRefPtLeadPtSubleadCMMcReweight) { delete hRefPtLeadPtSubleadCMMcReweight; hRefPtLeadPtSubleadCMMcReweight = nullptr; }
    if (hRefEtaCMLeadEtaCMSubleadMcReweight) { delete hRefEtaCMLeadEtaCMSubleadMcReweight; hRefEtaCMLeadEtaCMSubleadMcReweight = nullptr; }

    if (hRefDijetEta) { delete hRefDijetEta; hRefDijetEta = nullptr; }
    if (hRefDijetEtaVsRecoDijetEta) { delete hRefDijetEtaVsRecoDijetEta; hRefDijetEtaVsRecoDijetEta = nullptr; }
    if (hRefDijetPtVsRecoDijetPt) { delete hRefDijetPtVsRecoDijetPt; hRefDijetPtVsRecoDijetPt = nullptr; }
    if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt) { delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt; hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt = nullptr; }
    if (hRefDijetPtEtaPhi) { delete hRefDijetPtEtaPhi; hRefDijetPtEtaPhi = nullptr; }
    if (hRefDijetPtEtaPhiWeighted) { delete hRefDijetPtEtaPhiWeighted; hRefDijetPtEtaPhiWeighted = nullptr; }
    if (hRefDijetPtEtaForward) { delete hRefDijetPtEtaForward; hRefDijetPtEtaForward = nullptr; }
    if (hRefDijetPtEtaBackward) { delete hRefDijetPtEtaBackward; hRefDijetPtEtaBackward = nullptr; }
    if (hRefDijetPtEtaForwardWeighted) { delete hRefDijetPtEtaForwardWeighted; hRefDijetPtEtaForwardWeighted = nullptr; }
    if (hRefDijetPtEtaBackwardWeighted) { delete hRefDijetPtEtaBackwardWeighted; hRefDijetPtEtaBackwardWeighted = nullptr; }

    if (hRefDijetEtaCM) { delete hRefDijetEtaCM; hRefDijetEtaCM = nullptr; }
    if (hRefDijetEtaVsRecoDijetEtaCM) { delete hRefDijetEtaVsRecoDijetEtaCM; hRefDijetEtaVsRecoDijetEtaCM = nullptr; }
    if (hRefDijetPtVsRecoDijetPtCM) { delete hRefDijetPtVsRecoDijetPtCM; hRefDijetPtVsRecoDijetPtCM = nullptr; }
    if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM) { delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM; hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM = nullptr; }
    if (hRefDijetPtEtaPhiCM) { delete hRefDijetPtEtaPhiCM; hRefDijetPtEtaPhiCM = nullptr; }
    if (hRefDijetPtEtaPhiCMWeighted) { delete hRefDijetPtEtaPhiCMWeighted; hRefDijetPtEtaPhiCMWeighted = nullptr; }
    if (hRefDijetPtEtaCMForward) { delete hRefDijetPtEtaCMForward; hRefDijetPtEtaCMForward = nullptr; }
    if (hRefDijetPtEtaCMBackward) { delete hRefDijetPtEtaCMBackward; hRefDijetPtEtaCMBackward = nullptr; }
    if (hRefDijetPtEtaCMForwardWeighted) { delete hRefDijetPtEtaCMForwardWeighted; hRefDijetPtEtaCMForwardWeighted = nullptr; }
    if (hRefDijetPtEtaCMBackwardWeighted) { delete hRefDijetPtEtaCMBackwardWeighted; hRefDijetPtEtaCMBackwardWeighted = nullptr; }

    for (int i = 0; i < 16; ++i) {
        if (hRefDijetEta1D[i]) { delete hRefDijetEta1D[i]; hRefDijetEta1D[i] = nullptr; }
        if (hRefDijetEta1DWeighted[i]) { delete hRefDijetEta1DWeighted[i]; hRefDijetEta1DWeighted[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2D[i]) { delete hRefEtaLeadVsEtaSubLead2D[i]; hRefEtaLeadVsEtaSubLead2D[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DWeighted[i]) { delete hRefEtaLeadVsEtaSubLead2DWeighted[i]; hRefEtaLeadVsEtaSubLead2DWeighted[i] = nullptr; }
        if (hRecoVsRefDijetEta2D[i]) { delete hRecoVsRefDijetEta2D[i]; hRecoVsRefDijetEta2D[i] = nullptr; }
        if (hRecoVsRefDijetEta2DWeighted[i]) { delete hRecoVsRefDijetEta2DWeighted[i]; hRecoVsRefDijetEta2DWeighted[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2D[i]) { delete hRecoVsRefLeadJetEta2D[i]; hRecoVsRefLeadJetEta2D[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DWeighted[i]) { delete hRecoVsRefLeadJetEta2DWeighted[i]; hRecoVsRefLeadJetEta2DWeighted[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2D[i]) { delete hRecoVsRefSubLeadJetEta2D[i]; hRecoVsRefSubLeadJetEta2D[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DWeighted[i]) { delete hRecoVsRefSubLeadJetEta2DWeighted[i]; hRecoVsRefSubLeadJetEta2DWeighted[i] = nullptr; }
        if (hRefDijetEtaForward1D[i]) { delete hRefDijetEtaForward1D[i]; hRefDijetEtaForward1D[i] = nullptr; }
        if (hRefDijetEtaForward1DWeighted[i]) { delete hRefDijetEtaForward1DWeighted[i]; hRefDijetEtaForward1DWeighted[i] = nullptr; }
        if (hRefDijetEtaBackward1D[i]) { delete hRefDijetEtaBackward1D[i]; hRefDijetEtaBackward1D[i] = nullptr; }
        if (hRefDijetEtaBackward1DWeighted[i]) { delete hRefDijetEtaBackward1DWeighted[i]; hRefDijetEtaBackward1DWeighted[i] = nullptr; }

        if (hRefDijetEta1DCM[i]) { delete hRefDijetEta1DCM[i]; hRefDijetEta1DCM[i] = nullptr; }
        if (hRefDijetEta1DCMWeighted[i]) { delete hRefDijetEta1DCMWeighted[i]; hRefDijetEta1DCMWeighted[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DCM[i]) { delete hRefEtaLeadVsEtaSubLead2DCM[i]; hRefEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DCMWeighted[i]) { delete hRefEtaLeadVsEtaSubLead2DCMWeighted[i]; hRefEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr; }
        if (hRecoVsRefDijetEta2DCM[i]) { delete hRecoVsRefDijetEta2DCM[i]; hRecoVsRefDijetEta2DCM[i] = nullptr; }
        if (hRecoVsRefDijetEta2DCMWeighted[i]) { delete hRecoVsRefDijetEta2DCMWeighted[i]; hRecoVsRefDijetEta2DCMWeighted[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DCM[i]) { delete hRecoVsRefLeadJetEta2DCM[i]; hRecoVsRefLeadJetEta2DCM[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DCMWeighted[i]) { delete hRecoVsRefLeadJetEta2DCMWeighted[i]; hRecoVsRefLeadJetEta2DCMWeighted[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DCM[i]) { delete hRecoVsRefSubLeadJetEta2DCM[i]; hRecoVsRefSubLeadJetEta2DCM[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DCMWeighted[i]) { delete hRecoVsRefSubLeadJetEta2DCMWeighted[i]; hRecoVsRefSubLeadJetEta2DCMWeighted[i] = nullptr; }
        if (hRefDijetEtaCMForward1D[i]) { delete hRefDijetEtaCMForward1D[i]; hRefDijetEtaCMForward1D[i] = nullptr; }
        if (hRefDijetEtaCMForward1DWeighted[i]) { delete hRefDijetEtaCMForward1DWeighted[i]; hRefDijetEtaCMForward1DWeighted[i] = nullptr; }
        if (hRefDijetEtaCMBackward1D[i]) { delete hRefDijetEtaCMBackward1D[i]; hRefDijetEtaCMBackward1D[i] = nullptr; }
        if (hRefDijetEtaCMBackward1DWeighted[i]) { delete hRefDijetEtaCMBackward1DWeighted[i]; hRefDijetEtaCMBackward1DWeighted[i] = nullptr; }
    }

    for (int i = 0; i < 6; ++i) {
        if (hRefDijetEta1DOldPt[i]) { delete hRefDijetEta1DOldPt[i]; hRefDijetEta1DOldPt[i] = nullptr; }
        if (hRefDijetEta1DOldPtWeighted[i]) { delete hRefDijetEta1DOldPtWeighted[i]; hRefDijetEta1DOldPtWeighted[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DOldPt[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPt[i]; hRefEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i]; hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr; }
        if (hRecoVsRefDijetEta2DOldPt[i]) { delete hRecoVsRefDijetEta2DOldPt[i]; hRecoVsRefDijetEta2DOldPt[i] = nullptr; }
        if (hRecoVsRefDijetEta2DOldPtWeighted[i]) { delete hRecoVsRefDijetEta2DOldPtWeighted[i]; hRecoVsRefDijetEta2DOldPtWeighted[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DOldPt[i]) { delete hRecoVsRefLeadJetEta2DOldPt[i]; hRecoVsRefLeadJetEta2DOldPt[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DOldPtWeighted[i]) { delete hRecoVsRefLeadJetEta2DOldPtWeighted[i]; hRecoVsRefLeadJetEta2DOldPtWeighted[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DOldPt[i]) { delete hRecoVsRefSubLeadJetEta2DOldPt[i]; hRecoVsRefSubLeadJetEta2DOldPt[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DOldPtWeighted[i]) { delete hRecoVsRefSubLeadJetEta2DOldPtWeighted[i]; hRecoVsRefSubLeadJetEta2DOldPtWeighted[i] = nullptr; }
        if (hRefDijetEtaForward1DOldPt[i]) { delete hRefDijetEtaForward1DOldPt[i]; hRefDijetEtaForward1DOldPt[i] = nullptr; }
        if (hRefDijetEtaForward1DOldPtWeighted[i]) { delete hRefDijetEtaForward1DOldPtWeighted[i]; hRefDijetEtaForward1DOldPtWeighted[i] = nullptr; }
        if (hRefDijetEtaBackward1DOldPt[i]) { delete hRefDijetEtaBackward1DOldPt[i]; hRefDijetEtaBackward1DOldPt[i] = nullptr; }
        if (hRefDijetEtaBackward1DOldPtWeighted[i]) { delete hRefDijetEtaBackward1DOldPtWeighted[i]; hRefDijetEtaBackward1DOldPtWeighted[i] = nullptr; }

        if (hRefDijetEta1DOldPtCM[i]) { delete hRefDijetEta1DOldPtCM[i]; hRefDijetEta1DOldPtCM[i] = nullptr; }
        if (hRefDijetEta1DOldPtCMWeighted[i]) { delete hRefDijetEta1DOldPtCMWeighted[i]; hRefDijetEta1DOldPtCMWeighted[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPtCM[i]; hRefEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
        if (hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) { delete hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]; hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr; }
        if (hRecoVsRefDijetEta2DOldPtCM[i]) { delete hRecoVsRefDijetEta2DOldPtCM[i]; hRecoVsRefDijetEta2DOldPtCM[i] = nullptr; }
        if (hRecoVsRefDijetEta2DOldPtCMWeighted[i]) { delete hRecoVsRefDijetEta2DOldPtCMWeighted[i]; hRecoVsRefDijetEta2DOldPtCMWeighted[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DOldPtCM[i]) { delete hRecoVsRefLeadJetEta2DOldPtCM[i]; hRecoVsRefLeadJetEta2DOldPtCM[i] = nullptr; }
        if (hRecoVsRefLeadJetEta2DOldPtCMWeighted[i]) { delete hRecoVsRefLeadJetEta2DOldPtCMWeighted[i]; hRecoVsRefLeadJetEta2DOldPtCMWeighted[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DOldPtCM[i]) { delete hRecoVsRefSubLeadJetEta2DOldPtCM[i]; hRecoVsRefSubLeadJetEta2DOldPtCM[i] = nullptr; }
        if (hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i]) { delete hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i]; hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i] = nullptr; }
        if (hRefDijetEtaCMForward1DOldPt[i]) { delete hRefDijetEtaCMForward1DOldPt[i]; hRefDijetEtaCMForward1DOldPt[i] = nullptr; }
        if (hRefDijetEtaCMForward1DOldPtWeighted[i]) { delete hRefDijetEtaCMForward1DOldPtWeighted[i]; hRefDijetEtaCMForward1DOldPtWeighted[i] = nullptr; }
        if (hRefDijetEtaCMBackward1DOldPt[i]) { delete hRefDijetEtaCMBackward1DOldPt[i]; hRefDijetEtaCMBackward1DOldPt[i] = nullptr; }
        if (hRefDijetEtaCMBackward1DOldPtWeighted[i]) { delete hRefDijetEtaCMBackward1DOldPtWeighted[i]; hRefDijetEtaCMBackward1DOldPtWeighted[i] = nullptr; }
    }

    //
    // Ref-selected jet histograms
    //

    if (hRefSelDijetEta) { delete hRefSelDijetEta; hRefSelDijetEta = nullptr; }
    if (hRefSelDijetPtEtaPhi) { delete hRefSelDijetPtEtaPhi; hRefSelDijetPtEtaPhi = nullptr; }
    if (hRefSelDijetPtEtaPhiWeighted) { delete hRefSelDijetPtEtaPhiWeighted; hRefSelDijetPtEtaPhiWeighted = nullptr; }

    if (hRefSelDijetEtaCM) { delete hRefSelDijetEtaCM; hRefSelDijetEtaCM = nullptr; }
    if (hRefSelDijetPtEtaPhiCM) { delete hRefSelDijetPtEtaPhiCM; hRefSelDijetPtEtaPhiCM = nullptr; }
    if (hRefSelDijetPtEtaPhiCMWeighted) { delete hRefSelDijetPtEtaPhiCMWeighted; hRefSelDijetPtEtaPhiCMWeighted = nullptr; }

    for (int i = 0; i < 16; ++i) {
        if (hRefSelDijetEta1D[i]) { delete hRefSelDijetEta1D[i]; hRefSelDijetEta1D[i] = nullptr; }
        if (hRefSelDijetEta1DWeighted[i]) { delete hRefSelDijetEta1DWeighted[i]; hRefSelDijetEta1DWeighted[i] = nullptr; }
        if (hRefSelRecoDijetEta1D[i]) { delete hRefSelRecoDijetEta1D[i]; hRefSelRecoDijetEta1D[i] = nullptr; }
        if (hRefSelRecoDijetEta1DWeighted[i]) { delete hRefSelRecoDijetEta1DWeighted[i]; hRefSelRecoDijetEta1DWeighted[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2D[i]) { delete hRefSelEtaLeadVsEtaSubLead2D[i]; hRefSelEtaLeadVsEtaSubLead2D[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DWeighted[i]) { delete hRefSelEtaLeadVsEtaSubLead2DWeighted[i]; hRefSelEtaLeadVsEtaSubLead2DWeighted[i] = nullptr; }
        if (hRefSelDijetEtaForward1D[i]) { delete hRefSelDijetEtaForward1D[i]; hRefSelDijetEtaForward1D[i] = nullptr; }
        if (hRefSelDijetEtaForward1DWeighted[i]) { delete hRefSelDijetEtaForward1DWeighted[i]; hRefSelDijetEtaForward1DWeighted[i] = nullptr; }
        if (hRefSelDijetEtaBackward1D[i]) { delete hRefSelDijetEtaBackward1D[i]; hRefSelDijetEtaBackward1D[i] = nullptr; }
        if (hRefSelDijetEtaBackward1DWeighted[i]) { delete hRefSelDijetEtaBackward1DWeighted[i]; hRefSelDijetEtaBackward1DWeighted[i] = nullptr; }

        if (hRefSelDijetEta1DCM[i]) { delete hRefSelDijetEta1DCM[i]; hRefSelDijetEta1DCM[i] = nullptr; }
        if (hRefSelDijetEta1DCMWeighted[i]) { delete hRefSelDijetEta1DCMWeighted[i]; hRefSelDijetEta1DCMWeighted[i] = nullptr; }
        if (hRefSelRecoDijetEta1DCM[i]) { delete hRefSelRecoDijetEta1DCM[i]; hRefSelRecoDijetEta1DCM[i] = nullptr; }
        if (hRefSelRecoDijetEta1DCMWeighted[i]) { delete hRefSelRecoDijetEta1DCMWeighted[i]; hRefSelRecoDijetEta1DCMWeighted[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DCM[i]) { delete hRefSelEtaLeadVsEtaSubLead2DCM[i]; hRefSelEtaLeadVsEtaSubLead2DCM[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i]) { delete hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i]; hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr; }
        if (hRefSelDijetEtaCMForward1D[i]) { delete hRefSelDijetEtaCMForward1D[i]; hRefSelDijetEtaCMForward1D[i] = nullptr; }
        if (hRefSelDijetEtaCMForward1DWeighted[i]) { delete hRefSelDijetEtaCMForward1DWeighted[i]; hRefSelDijetEtaCMForward1DWeighted[i] = nullptr; }
        if (hRefSelDijetEtaCMBackward1D[i]) { delete hRefSelDijetEtaCMBackward1D[i]; hRefSelDijetEtaCMBackward1D[i] = nullptr; }
        if (hRefSelDijetEtaCMBackward1DWeighted[i]) { delete hRefSelDijetEtaCMBackward1DWeighted[i]; hRefSelDijetEtaCMBackward1DWeighted[i] = nullptr; }
    }

    for (int i = 0; i < 6; ++i) {
        if (hRefSelDijetEta1DOldPt[i]) { delete hRefSelDijetEta1DOldPt[i]; hRefSelDijetEta1DOldPt[i] = nullptr; }
        if (hRefSelDijetEta1DOldPtWeighted[i]) { delete hRefSelDijetEta1DOldPtWeighted[i]; hRefSelDijetEta1DOldPtWeighted[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPt[i]) { delete hRefSelRecoDijetEta1DOldPt[i]; hRefSelRecoDijetEta1DOldPt[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtWeighted[i]) { delete hRefSelRecoDijetEta1DOldPtWeighted[i]; hRefSelRecoDijetEta1DOldPtWeighted[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPt[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPt[i]; hRefSelEtaLeadVsEtaSubLead2DOldPt[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr; }
        if (hRefSelDijetEtaForward1DOldPt[i]) { delete hRefSelDijetEtaForward1DOldPt[i]; hRefSelDijetEtaForward1DOldPt[i] = nullptr; }
        if (hRefSelDijetEtaForward1DOldPtWeighted[i]) { delete hRefSelDijetEtaForward1DOldPtWeighted[i]; hRefSelDijetEtaForward1DOldPtWeighted[i] = nullptr; }
        if (hRefSelDijetEtaBackward1DOldPt[i]) { delete hRefSelDijetEtaBackward1DOldPt[i]; hRefSelDijetEtaBackward1DOldPt[i] = nullptr; }
        if (hRefSelDijetEtaBackward1DOldPtWeighted[i]) { delete hRefSelDijetEtaBackward1DOldPtWeighted[i]; hRefSelDijetEtaBackward1DOldPtWeighted[i] = nullptr; }

        if (hRefSelDijetEta1DOldPtCM[i]) { delete hRefSelDijetEta1DOldPtCM[i]; hRefSelDijetEta1DOldPtCM[i] = nullptr; }
        if (hRefSelDijetEta1DOldPtCMWeighted[i]) { delete hRefSelDijetEta1DOldPtCMWeighted[i]; hRefSelDijetEta1DOldPtCMWeighted[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtCM[i]) { delete hRefSelRecoDijetEta1DOldPtCM[i]; hRefSelRecoDijetEta1DOldPtCM[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtCMWeighted[i]) { delete hRefSelRecoDijetEta1DOldPtCMWeighted[i]; hRefSelRecoDijetEta1DOldPtCMWeighted[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr; }
        if (hRefSelDijetEtaCMForward1DOldPt[i]) { delete hRefSelDijetEtaCMForward1DOldPt[i]; hRefSelDijetEtaCMForward1DOldPt[i] = nullptr; }
        if (hRefSelDijetEtaCMForward1DOldPtWeighted[i]) { delete hRefSelDijetEtaCMForward1DOldPtWeighted[i]; hRefSelDijetEtaCMForward1DOldPtWeighted[i] = nullptr; }
        if (hRefSelDijetEtaCMBackward1DOldPt[i]) { delete hRefSelDijetEtaCMBackward1DOldPt[i]; hRefSelDijetEtaCMBackward1DOldPt[i] = nullptr; }
        if (hRefSelDijetEtaCMBackward1DOldPtWeighted[i]) { delete hRefSelDijetEtaCMBackward1DOldPtWeighted[i]; hRefSelDijetEtaCMBackward1DOldPtWeighted[i] = nullptr; }

        if (hRefSelDijetEta1DOldPtBinning[i]) { delete hRefSelDijetEta1DOldPtBinning[i]; hRefSelDijetEta1DOldPtBinning[i] = nullptr; }
        if (hRefSelDijetEta1DOldPtBinningWeighted[i]) { delete hRefSelDijetEta1DOldPtBinningWeighted[i]; hRefSelDijetEta1DOldPtBinningWeighted[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtBinning[i]) { delete hRefSelRecoDijetEta1DOldPtBinning[i]; hRefSelRecoDijetEta1DOldPtBinning[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtBinningWeighted[i]) { delete hRefSelRecoDijetEta1DOldPtBinningWeighted[i]; hRefSelRecoDijetEta1DOldPtBinningWeighted[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr; }
        if (hRefSelDijetEtaForward1DOldPtBinning[i]) { delete hRefSelDijetEtaForward1DOldPtBinning[i]; hRefSelDijetEtaForward1DOldPtBinning[i] = nullptr; }
        if (hRefSelDijetEtaForward1DOldPtBinningWeighted[i]) { delete hRefSelDijetEtaForward1DOldPtBinningWeighted[i]; hRefSelDijetEtaForward1DOldPtBinningWeighted[i] = nullptr; }
        if (hRefSelDijetEtaBackward1DOldPtBinning[i]) { delete hRefSelDijetEtaBackward1DOldPtBinning[i]; hRefSelDijetEtaBackward1DOldPtBinning[i] = nullptr; }
        if (hRefSelDijetEtaBackward1DOldPtBinningWeighted[i]) { delete hRefSelDijetEtaBackward1DOldPtBinningWeighted[i]; hRefSelDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr; }

        if (hRefSelDijetEta1DOldPtBinningCM[i]) { delete hRefSelDijetEta1DOldPtBinningCM[i]; hRefSelDijetEta1DOldPtBinningCM[i] = nullptr; }
        if (hRefSelDijetEta1DOldPtBinningCMWeighted[i]) { delete hRefSelDijetEta1DOldPtBinningCMWeighted[i]; hRefSelDijetEta1DOldPtBinningCMWeighted[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtBinningCM[i]) { delete hRefSelRecoDijetEta1DOldPtBinningCM[i]; hRefSelRecoDijetEta1DOldPtBinningCM[i] = nullptr; }
        if (hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i]) { delete hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i]; hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr; }
        if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) { delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]; hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr; }
        if (hRefSelDijetEtaCMForward1DOldPtBinning[i]) { delete hRefSelDijetEtaCMForward1DOldPtBinning[i]; hRefSelDijetEtaCMForward1DOldPtBinning[i] = nullptr; }
        if (hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i]) { delete hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i]; hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr; }
        if (hRefSelDijetEtaCMBackward1DOldPtBinning[i]) { delete hRefSelDijetEtaCMBackward1DOldPtBinning[i]; hRefSelDijetEtaCMBackward1DOldPtBinning[i] = nullptr; }
        if (hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i]) { delete hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i]; hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr; }
    }

}

//________________
void HistoManagerDiJet::init() {

    const int dijetEtaBins{28};
    double dijetEtaVals[dijetEtaBins+1] { -4.0, -3.0, -2.4, -2.2, 
                                          -2.0, -1.8, -1.6, -1.4, -1.2, 
                                          -1.0, -0.8, -0.6, -0.4, -0.2,  
                                           0.0,  0.2,  0.4,  0.6,  0.8,  
                                           1.0,  1.2,  1.4,  1.6,  1.8,  
                                           2.0,  2.2,  2.4,  3.0,  4.0 };

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
    // int    fJESBins{100}; 
    // double fJESRange[2] {0., 2.};
    int    ptHatBins{100};
    double ptHatRange[2] {0., 1000.};
    int    dEtaBins{50}; 
    double dEtaRange[2] {-0.05, 0.05};

    // [11 dimensions]
    // 0 - Dijet ptAve, 1 - dijet eta, 2 - dijet delta eta, 3 - dijet dphi, 4 - dijet phi,
    // 5 - Lead pt, 6 - lead eta, 7 - lead phi,
    // 8 - Sublead pt, 9 - sublead eta, 10 - sublead phi
    int    bins11D_dijet_params[11]
    { fDijetPtBins, fDijetEtaBins, dEtaBins, fDijetDphiBins, fPhiBins, fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins };
    double xmin11D_dijet_params[11]
    { fDijetPtRange[0], fDijetEtaRange[0], dEtaRange[0], fDijetDphiRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0] };
    double xmax11D_dijet_params[11]
    { fDijetPtRange[1], fDijetEtaRange[1], dEtaRange[1], fDijetDphiRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1] };

    // [12 dimensions]
    // 0 - Reco dijet ptAve, 1 - dijet eta, 2 - lead pt, 3 - lead eta,
    // 4 - sublead pt, 5 - sublead eta,
    // 6 - Ref dijet ptAve, 7 - dijet eta, 8 - Ref lead pt, 9 - lead eta,
    // 10 - sublead pt, 11 - sublead eta
    int    bins12D_dijet_reco2ref[12]
    { fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins, fDijetPtBins, fDijetEtaBins, fPtBins, fEtaBins, fPtBins, fEtaBins };
    double xmin12D_dijet_reco2ref[12]
    { fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0], fPtRange[0], fEtaRange[0], fPtRange[0], fEtaRange[0] };
    double xmax12D_dijet_reco2ref[12]
    { fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1], fPtRange[1], fEtaRange[1], fPtRange[1], fEtaRange[1] };

    // [4 dimensions]
    // 0 - reco dijet ptAve, 1 - reco dijet eta,
    // 2 - ref dijet ptAve, 3 - ref dijet eta   
    int    bins4D_dijet_reco2ref[4] { fDijetPtBins, fDijetEtaBins, fDijetPtBins, fDijetEtaBins };
    double xmin4D_dijet_reco2ref[4] { fDijetPtRange[0], fDijetEtaRange[0], fDijetPtRange[0], fDijetEtaRange[0] };
    double xmax4D_dijet_reco2ref[4] { fDijetPtRange[1], fDijetEtaRange[1], fDijetPtRange[1], fDijetEtaRange[1] };


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

    hRecoJetCollectionSize = new TH1D("hRecoJetCollectionSize","Reco jet collection size;Number of jets;Entries", 100, -0.5, 99.5);
    hRecoJetCollectionSize->Sumw2();

    // Weighted dijet distribution (dijet, lead, sublead) in the laboratory frame [11 dimensions]
    // 0 - Dijet ptAve, 1 - dijet eta, 2 - dijet delta eta, 3 - dijet dphi, 4 - dijet phi,
    // 5 - Lead pt, 6 - lead eta, 7 - lead phi,
    // 8 - Sublead pt, 9 - sublead eta, 10 - sublead phi
    hRecoDijetLeadSubLead = new THnSparseD("hRecoDijetLeadSubLead",
            "Reconstructed dijet and jet info in the lab frame;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#eta^{dijet};#Delta#phi^{dijet} (rad);#phi^{dijet} (rad);p_{T}^{Lead} (GeV);#eta^{Lead};#phi^{Lead} (rad);p_{T}^{SubLead} (GeV);#eta^{SubLead};#phi^{SubLead} (rad)",
            11,
            bins11D_dijet_params,
            xmin11D_dijet_params,
            xmax11D_dijet_params);
    if (fUseVariableBinning) hRecoDijetLeadSubLead->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    if (fUseVariableBinning) hRecoDijetLeadSubLead->GetAxis(6)->Set(dijetEtaBins, dijetEtaVals);
    if (fUseVariableBinning) hRecoDijetLeadSubLead->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetLeadSubLead->Sumw2();
    // Weighted dijet distribution (dijet, lead, sublead) in the center-of-mass frame [11 dimensions]
    // 0 - Dijet ptAve, 1 - dijet eta, 2 - dijet delta eta, 3 - dijet dphi, 4 - dijet phi,
    // 5 - Lead pt, 6 - lead eta, 7 - lead phi,
    // 8 - Sublead pt, 9 - sublead eta, 10 - sublead phi
    hRecoDijetLeadSubLeadCM = new THnSparseD("hRecoDijetLeadSubLeadCM",
            "Reconstructed dijet and jet info in the CM frame;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#eta^{dijet}_{CM};#Delta#phi^{dijet} (rad);#phi^{dijet} (rad);p_{T}^{Lead} (GeV);#eta^{Lead}_{CM};#phi^{Lead} (rad);p_{T}^{SubLead} (GeV);#eta^{SubLead}_{CM};#phi^{SubLead} (rad)",
            11,
            bins11D_dijet_params,
            xmin11D_dijet_params,
            xmax11D_dijet_params);
    if (fUseVariableBinning) hRecoDijetLeadSubLeadCM->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    if (fUseVariableBinning) hRecoDijetLeadSubLeadCM->GetAxis(6)->Set(dijetEtaBins, dijetEtaVals);
    if (fUseVariableBinning) hRecoDijetLeadSubLeadCM->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetLeadSubLeadCM->Sumw2();


    // Lab frame
    hRecoPtLeadPtSublead = new TH2D("hRecoPtLeadPtSublead","Reco leading vs subleading p_{T};Reco p_{T}^{Lead} (GeV);Reco p_{T}^{SubLead} (GeV)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSublead->Sumw2();
    hRecoEtaLeadEtaSublead = new TH2D("hRecoEtaLeadEtaSublead","Reco leading vs subleading #eta;Reco #eta^{Lead};Reco #eta^{SubLead}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSublead->Sumw2();
    hRecoPtLeadPtSubleadMcReweight = new TH2D("hRecoPtLeadPtSubleadMcReweight","Reco leading vs subleading p_{T} (MC reweighted to data);Reco p_{T}^{Lead} (GeV);Reco p_{T}^{SubLead} (GeV)",
                                     fPtBins, fPtRange[0], fPtRange[1],
                                     fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSubleadMcReweight->Sumw2();
    hRecoEtaLeadEtaSubleadMcReweight = new TH2D("hRecoEtaLeadEtaSubleadMcReweight","Reco leading vs subleading #eta (MC reweighted to data);Reco #eta^{Lead};Reco #eta^{SubLead}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaLeadEtaSubleadMcReweight->Sumw2();
    hRecoDijetEta = new TH1D("hRecoDijetEta","Reco dijet #eta;Reco #eta^{dijet};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEta->Sumw2();
    hRecoDijetPtEta = new TH2D("hRecoDijetPtEta", "Reco dijet #eta vs p_{T};p_{T}^{ave} (GeV);#eta^{dijet}", 
                               fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                               fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEta->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEta->Sumw2();
    hRecoDijetPtEtaForward = new TH2D("hRecoDijetPtEtaForward", "Reco dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaForward->Sumw2();
    hRecoDijetPtEtaBackward = new TH2D("hRecoDijetPtEtaBackward", "Reco dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaBackward->Sumw2();
    hRecoDijetPtEtaForwardWeighted = new TH2D("hRecoDijetPtEtaForwardWeighted", "Reco dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaForwardWeighted->Sumw2();
    hRecoDijetPtEtaBackwardWeighted = new TH2D("hRecoDijetPtEtaBackwardWeighted", "Reco dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaBackwardWeighted->Sumw2();
    hRecoDijetPtEtaPhi = new TH3D("hRecoDijetPtEtaPhi","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fPhiBins, fPhiRange[0], fPhiRange[1] );
    if (fUseVariableBinning) hRecoDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaPhi->Sumw2();
    hRecoDijetPtEtaPhiWeighted = new TH3D("hRecoDijetPtEtaPhiWeighted","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fPhiBins, fPhiRange[0], fPhiRange[1] );
    if (fUseVariableBinning) hRecoDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaPhiWeighted->Sumw2();

    // Center-of-mass frame
    hRecoPtLeadPtSubleadCM = new TH2D("hRecoPtLeadPtSubleadCM","Reco leading vs subleading p_{T} in CM;Reco p_{T}^{Lead}_{CM} (GeV);Reco p_{T}^{SubLead}_{CM} (GeV)",
                                       fPtBins, fPtRange[0], fPtRange[1],
                                       fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSubleadCM->Sumw2();
    hRecoEtaCMLeadEtaCMSublead = new TH2D("hRecoEtaCMLeadEtaCMSublead","Reco leading vs subleading #eta in CM;Reco #eta^{Lead}_{CM};Reco #eta^{SubLead}_{CM}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaCMLeadEtaCMSublead->Sumw2();
    hRecoPtLeadPtSubleadCMMcReweight = new TH2D("hRecoPtLeadPtSubleadCMMcReweight","Reco leading vs subleading p_{T} in CM (MC reweighted to data);Reco p_{T}^{Lead}_{CM} (GeV);Reco p_{T}^{SubLead}_{CM} (GeV)",
                                       fPtBins, fPtRange[0], fPtRange[1],
                                       fPtBins, fPtRange[0], fPtRange[1]);
    hRecoPtLeadPtSubleadCMMcReweight->Sumw2();
    hRecoEtaCMLeadEtaCMSubleadMcReweight = new TH2D("hRecoEtaCMLeadEtaCMSubleadMcReweight","Reco leading vs subleading #eta in CM (MC reweighted to data);Reco #eta^{Lead}_{CM};Reco #eta^{SubLead}_{CM}",
                                       fEtaBins, fEtaRange[0], fEtaRange[1],
                                       fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoEtaCMLeadEtaCMSubleadMcReweight->Sumw2();

    hRecoDijetEtaCM = new TH1D("hRecoDijetEtaCM","Reco dijet #eta in CM;Reco #eta^{dijet}_{CM};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEtaCM->Sumw2();
    hRecoDijetPtEtaCM = new TH2D("hRecoDijetPtEtaCM", "Reco dijet #eta vs p_{T} in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}", 
                               fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                               fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaCM->Sumw2();
    hRecoDijetPtEtaCMForward = new TH2D("hRecoDijetPtEtaCMForward", "Reco dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMForward->Sumw2();
    hRecoDijetPtEtaCMBackward = new TH2D("hRecoDijetPtEtaCMBackward", "Reco dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMBackward->Sumw2();
    hRecoDijetPtEtaCMForwardWeighted = new TH2D("hRecoDijetPtEtaCMForwardWeighted", "Reco dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMForwardWeighted->Sumw2();
    hRecoDijetPtEtaCMBackwardWeighted = new TH2D("hRecoDijetPtEtaCMBackwardWeighted", "Reco dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
    if (fUseVariableBinning) hRecoDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
    hRecoDijetPtEtaCMBackwardWeighted->Sumw2();

    hRecoDijetPtEtaPhiCM = new TH3D("hRecoDijetPtEtaPhiCM","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fPhiBins, fPhiRange[0], fPhiRange[1] );
    if (fUseVariableBinning) hRecoDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaPhiCM->Sumw2();
    hRecoDijetPtEtaPhiCMWeighted = new TH3D("hRecoDijetPtEtaPhiCMWeighted","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                       fPhiBins, fPhiRange[0], fPhiRange[1] );
    if (fUseVariableBinning) hRecoDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaPhiCMWeighted->Sumw2();


    // New ptAve and eta binning
    for (unsigned int i{0}; i<fPtAveBins.size()-1; ++i) {
        double ptAveLow = fPtAveBins.at(i);
        double ptAveHi = fPtAveBins.at(i+1);
        
        // Lab frame
        hRecoDijetEta1D[i] = new TH1D(Form("hRecoDijetEta1D_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                      prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1D[i]->Sumw2();
        
        hRecoDijetEta1DWeighted[i] = new TH1D(Form("hRecoDijetEta1DWeighted_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                              prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DWeighted[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2D_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}", i, ptAveLow, ptAveHi),
                                                      fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2D[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DWeighted_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}", i, ptAveLow, ptAveHi),
                                                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaForward1D[i] = new TH1D(Form("hRecoDijetEtaForward1D_%d",i),Form("Reco #eta^{dijet}_{forward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                             fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1D[i]->Sumw2();
        hRecoDijetEtaForward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaForward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{forward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                     fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DWeighted[i]->Sumw2();
        hRecoDijetEtaBackward1D[i] = new TH1D(Form("hRecoDijetEtaBackward1D_%d",i),Form("Reco #eta^{dijet}_{backward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                              fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1D[i]->Sumw2();
        hRecoDijetEtaBackward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaBackward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{backward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                      fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DWeighted[i]->Sumw2();

        // CM frame
        hRecoDijetEta1DCM[i] = new TH1D(Form("hRecoDijetEta1DCM_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                        prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DCM[i]->Sumw2();
        hRecoDijetEta1DCMWeighted[i] = new TH1D(Form("hRecoDijetEta1DCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DCMWeighted[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DCM_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}", i, ptAveLow, ptAveHi),
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}", i, ptAveLow, ptAveHi),
                                                           fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
        hRecoDijetEtaCMForward1D[i] = new TH1D(Form("hRecoDijetEtaCMForward1D_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                               fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1D[i]->Sumw2();
        hRecoDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaCMForward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                       fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DWeighted[i]->Sumw2();
        hRecoDijetEtaCMBackward1D[i] = new TH1D(Form("hRecoDijetEtaCMBackward1D_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1D[i]->Sumw2();
        hRecoDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM", i, ptAveLow, ptAveHi),
                                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRecoDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DWeighted[i]->Sumw2();

    } // for (int i{0}; i<fPtAveBins.size()-2; ++i)


    for (unsigned int i{0}; i<fPtAveOldBins.size()-1; ++i) {

        double ptAveLow = fPtAveOldBins.at(i);
        double ptAveHi = fPtAveOldBins.at(i+1);

        // Lab frame
        hRecoDijetEta1DOldPt[i] = new TH1D(Form("hRecoDijetEta1DOldPt_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                        dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPt[i]->Sumw2();
        hRecoDijetEta1DOldPtWeighted[i] = new TH1D(Form("hRecoDijetEta1DOldPtWeighted_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPtWeighted[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Reco dijet #eta lead vs sublead in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted_%d",i),Form("Reco dijet #eta lead vs sublead in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Sumw2();
        hRecoDijetEtaForward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaForward1DOldPt_%d",i),Form("Reco #eta^{dijet} forward in %d in the lab frame for %3.0f<p_{T}^{ave} (GeV)<%3.0f;Reco #eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DOldPt[i]->Sumw2();
        hRecoDijetEtaForward1DOldPtWeighted[i] = new TH1D(Form("hRecoDijetEtaForward1DOldPtWeighted_%d",i),Form("Reco #eta^{dijet} forward in %d in the lab frame for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DOldPtWeighted[i]->Sumw2();
        hRecoDijetEtaBackward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaBackward1DOldPt_%d",i),Form("Reco #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DOldPt[i]->Sumw2();
        hRecoDijetEtaBackward1DOldPtWeighted[i] = new TH1D(Form("hRecoDijetEtaBackward1DOldPtWeighted_%d",i),Form("Reco #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DOldPtWeighted[i]->Sumw2();


        hRecoDijetEta1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinning_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                  dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinning[i]->Sumw2();
        hRecoDijetEta1DOldPtBinningWeighted[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinningWeighted_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                          dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinningWeighted[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Reco dijet #eta lead vs sublead in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                              fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted_%d",i),Form("Reco dijet #eta lead vs sublead in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                  fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]); 
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Sumw2();
        hRecoDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaForward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet} forward  in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;Reco #eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaForward1DOldPtBinningWeighted[i] = new TH1D(Form("hRecoDijetEtaForward1DOldPtBinningWeighted_%d",i),Form("Reco #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                          dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DOldPtBinningWeighted[i]->Sumw2();
        hRecoDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaBackward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hRecoDijetEtaBackward1DOldPtBinningWeighted_%d",i),Form("Reco #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                          dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DOldPtBinningWeighted[i]->Sumw2();

        // CM frame
        hRecoDijetEta1DOldPtCM[i] = new TH1D(Form("hRecoDijetEta1DOldPtCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                            dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPtCM[i]->Sumw2();
        hRecoDijetEta1DOldPtCMWeighted[i] = new TH1D(Form("hRecoDijetEta1DOldPtCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DOldPtCMWeighted[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Sumw2();
        hRecoDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaCMForward1DOldPt_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DOldPt[i]->Sumw2();
        hRecoDijetEtaCMForward1DOldPtWeighted[i] = new TH1D(Form("hRecoDijetEtaCMForward1DOldPtWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                    dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DOldPtWeighted[i]->Sumw2();
        hRecoDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DOldPt_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DOldPt[i]->Sumw2();
        hRecoDijetEtaCMBackward1DOldPtWeighted[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DOldPtWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                    dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DOldPtWeighted[i]->Sumw2();

        
        hRecoDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinningCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                  dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinningCM[i]->Sumw2();
        hRecoDijetEta1DOldPtBinningCMWeighted[i] = new TH1D(Form("hRecoDijetEta1DOldPtBinningCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                          dijetEtaOldBins, dijetEtaOldVals);
        hRecoDijetEta1DOldPtBinningCMWeighted[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                              fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Sumw2();
        hRecoDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaCMForward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaCMForward1DOldPtBinningWeighted[i] = new TH1D(Form("hRecoDijetEtaCMForward1DOldPtBinningWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                          dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DOldPtBinningWeighted[i]->Sumw2();
        hRecoDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DOldPtBinning_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                      dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
        hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DOldPtBinningWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                          dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i]->Sumw2();
    }

    
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

        hGenLeadingJetPtOverPtHatVsLeadingJetPt = new TH2D("hGenLeadingJetPtOverPtHatVsLeadingJetPt", "Leading jet p_{T}^{Gen}/#hat{p}_{T} vs leading jet p_{T}^{Gen};p_{T}^{Gen} (GeV);p_{T}^{Gen}/#hat{p}_{T}",
                                                            fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenLeadingJetPtOverPtHatVsLeadingJetPt->Sumw2();
        hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted = new TH2D("hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted", "Leading jet p_{T}^{Gen}/#hat{p}_{T} vs leading jet p_{T}^{Gen} weighted;p_{T}^{Gen} (GeV);p_{T}^{Gen}/#hat{p}_{T}",
                                                                    fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Sumw2();
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
        
        // Inclusive dijets
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

        // Inclusive dijets that passed the selection
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

        // Weighted dijet distribution (dijet, lead, sublead) in the laboratory frame [11 dimensions]
        // 0 - Dijet ptAve, 1 - dijet eta, 2 - dijet delta eta, 3 - dijet dphi, 4 - dijet phi,
        // 5 - Lead pt, 6 - lead eta, 7 - lead phi,
        // 8 - Sublead pt, 9 - sublead eta, 10 - sublead phi
        hGenDijetLeadSubLead = new THnSparseD("hGenDijetLeadSubLead",
            "Gen dijet and jet info in the lab frame;Gen p_{T}^{ave} (GeV);Gen #eta^{dijet};Gen #Delta#eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen #phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{SubLead} (GeV);Gen #eta^{SubLead};Gen #phi^{SubLead} (rad)",
            11,
            bins11D_dijet_params,
            xmin11D_dijet_params,
            xmax11D_dijet_params);
        if (fUseVariableBinning) hGenDijetLeadSubLead->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hGenDijetLeadSubLead->GetAxis(6)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hGenDijetLeadSubLead->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetLeadSubLead->Sumw2();
        // Weighted dijet distribution (dijet, lead, sublead) in the center-of-mass frame [11 dimensions]
        // 0 - Dijet ptAve, 1 - dijet eta, 2 - dijet delta eta, 3 - dijet dphi, 4 - dijet phi,
        // 5 - Lead pt, 6 - lead eta, 7 - lead phi,
        // 8 - Sublead pt, 9 - sublead eta, 10 - sublead phi        
        hGenDijetLeadSubLeadCM = new THnSparseD("hGenDijetLeadSubLeadCM",
            "Gen dijet and jet info in the CM frame;Gen p_{T}^{ave} (GeV);Gen #eta^{dijet}_{CM};Gen #Delta#eta^{dijet}_{CM};Gen #Delta#phi^{dijet} (rad);Gen #phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV);Gen #eta^{Lead}_{CM};Gen #phi^{Lead} (rad);Gen p_{T}^{SubLead} (GeV);Gen #eta^{SubLead}_{CM};Gen #phi^{SubLead} (rad)",
            11,
            bins11D_dijet_params,
            xmin11D_dijet_params,
            xmax11D_dijet_params);
        if (fUseVariableBinning) hGenDijetLeadSubLeadCM->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hGenDijetLeadSubLeadCM->GetAxis(6)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hGenDijetLeadSubLeadCM->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetLeadSubLeadCM->Sumw2();

        // Lab frame
        hGenPtLeadPtSublead = new TH2D("hGenPtLeadPtSublead","Leading gen jet pT vs subleading gen jet pT;Gen p_{T}^{Lead} (GeV);Gen p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenPtLeadPtSublead->Sumw2();
        hGenEtaLeadEtaSublead = new TH2D("hGenEtaLeadEtaSublead","Leading gen jet eta vs subleading gen jet eta;Gen #eta^{Lead};Gen #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaLeadEtaSublead->Sumw2();
        hGenDijetEta = new TH1D("hGenDijetEta", "Gen dijet #eta;#eta^{dijet}",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEta->Sumw2();
        hGenDijetPtEtaPhi = new TH3D("hGenDijetPtEtaPhi","Gen dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi^{dijet} (rad)",
                                      fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                      fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                      fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hGenDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaPhi->Sumw2();
        hGenDijetPtEtaPhiWeighted = new TH3D("hGenDijetPtEtaPhiWeighted","Gen dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hGenDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaPhiWeighted->Sumw2();
        hGenDijetPtEtaForward = new TH2D("hGenDijetPtEtaForward", "Gen dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaForward->Sumw2();
        hGenDijetPtEtaBackward = new TH2D("hGenDijetPtEtaBackward", "Gen dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaBackward->Sumw2();
        hGenDijetPtEtaForwardWeighted = new TH2D("hGenDijetPtEtaForwardWeighted", "Gen dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaForwardWeighted->Sumw2();
        hGenDijetPtEtaBackwardWeighted = new TH2D("hGenDijetPtEtaBackwardWeighted", "Gen dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaBackwardWeighted->Sumw2();


        // Center-of-mass frame
        hGenPtLeadPtSubleadCM = new TH2D("hGenPtLeadPtSubleadCM","Leading gen jet pT in CM vs subleading gen jet pT in CM;Gen p_{T}^{Lead}_{CM} (GeV);Gen p_{T}^{SubLead}_{CM} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenPtLeadPtSubleadCM->Sumw2();
        hGenEtaCMLeadEtaCMSublead = new TH2D("hGenEtaCMLeadEtaCMSublead","Leading gen jet eta in CM vs subleading gen jet eta in CM;Gen #eta^{Lead}_{CM};Gen #eta^{SubLead}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenEtaCMLeadEtaCMSublead->Sumw2();
        hGenDijetEtaCM = new TH1D("hGenDijetEtaCM", "Gen dijet #eta in CM;#eta^{dijet}_{CM}",
                                  fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEtaCM->Sumw2();
        hGenDijetPtEtaPhiCM = new TH3D("hGenDijetPtEtaPhiCM","Gen dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hGenDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaPhiCM->Sumw2();
        hGenDijetPtEtaPhiCMWeighted = new TH3D("hGenDijetPtEtaPhiCMWeighted","Gen dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hGenDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaPhiCMWeighted->Sumw2();
        hGenDijetPtEtaCMForward = new TH2D("hGenDijetPtEtaCMForward", "Gen dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMForward->Sumw2();
        hGenDijetPtEtaCMBackward = new TH2D("hGenDijetPtEtaCMBackward", "Gen dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMBackward->Sumw2();
        hGenDijetPtEtaCMForwardWeighted = new TH2D("hGenDijetPtEtaCMForwardWeighted", "Gen dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMForwardWeighted->Sumw2();
        hGenDijetPtEtaCMBackwardWeighted = new TH2D("hGenDijetPtEtaCMBackwardWeighted", "Gen dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hGenDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hGenDijetPtEtaCMBackwardWeighted->Sumw2();


        for (int i=0; i< (int)(fPtAveBins.size()-1); i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hGenDijetEta1D[i] = new TH1D(Form("hGenDijetEta1D_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1D[i]->Sumw2();
            hGenDijetEta1DWeighted[i] = new TH1D(Form("hGenDijetEta1DWeighted_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DWeighted[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2D_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                            fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DWeighted[i]->Sumw2();
            hGenDijetEtaForward1D[i] = new TH1D(Form("hGenDijetEtaForward1D_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1D[i]->Sumw2();
            hGenDijetEtaForward1DWeighted[i] = new TH1D(Form("hGenDijetEtaForward1DWeighted_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DWeighted[i]->Sumw2();
            hGenDijetEtaBackward1D[i] = new TH1D(Form("hGenDijetEtaBackward1D_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1D[i]->Sumw2();
            hGenDijetEtaBackward1DWeighted[i] = new TH1D(Form("hGenDijetEtaBackward1DWeighted_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DWeighted[i]->Sumw2();

            hGenDijetEta1DCM[i] = new TH1D(Form("hGenDijetEta1DCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                           prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DCM[i]->Sumw2();
            hGenDijetEta1DCMWeighted[i] = new TH1D(Form("hGenDijetEta1DCMWeighted_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DCMWeighted[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DCM_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DCMWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                            fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
            hGenDijetEtaCMForward1D[i] = new TH1D(Form("hGenDijetEtaCMForward1D_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                  fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1D[i]->Sumw2();
            hGenDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hGenDijetEtaCMForward1DWeighted_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DWeighted[i]->Sumw2();
            hGenDijetEtaCMBackward1D[i] = new TH1D(Form("hGenDijetEtaCMBackward1D_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1D[i]->Sumw2();
            hGenDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hGenDijetEtaCMBackward1DWeighted_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hGenDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DWeighted[i]->Sumw2();
        }

        for (int i=0; i<(int)(fPtAveOldBins.size()-1); i++) {

            double ptAveLow = fPtAveOldBins.at(i);
            double ptAveHi = fPtAveOldBins.at(i+1);
            // New eta binning
            hGenDijetEta1DOldPt[i] = new TH1D(Form("hGenDijetEta1DOldPt_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                              dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPt[i]->Sumw2();
            hGenDijetEta1DOldPtWeighted[i] = new TH1D(Form("hGenDijetEta1DOldPtWeighted_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                      dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPtWeighted[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPt_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                              fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                  fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Sumw2();
            hGenDijetEtaForward1DOldPt[i] = new TH1D(Form("hGenDijetEtaForward1DOldPt_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                  dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DOldPt[i]->Sumw2();
            hGenDijetEtaForward1DOldPtWeighted[i] = new TH1D(Form("hGenDijetEtaForward1DOldPtWeighted_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                          dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DOldPtWeighted[i]->Sumw2();
            hGenDijetEtaBackward1DOldPt[i] = new TH1D(Form("hGenDijetEtaBackward1DOldPt_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                  dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DOldPt[i]->Sumw2();
            hGenDijetEtaBackward1DOldPtWeighted[i] = new TH1D(Form("hGenDijetEtaBackward1DOldPtWeighted_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                          dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DOldPtWeighted[i]->Sumw2();

        
            hGenDijetEta1DOldPtCM[i] = new TH1D(Form("hGenDijetEta1DOldPtCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                              dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPtCM[i]->Sumw2();
            hGenDijetEta1DOldPtCMWeighted[i] = new TH1D(Form("hGenDijetEta1DOldPtCMWeighted_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                      dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DOldPtCMWeighted[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtCM_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                              fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                                  fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Sumw2();
            hGenDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hGenDijetEtaCMForward1DOldPt_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                       dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DOldPt[i]->Sumw2();
            hGenDijetEtaCMForward1DOldPtWeighted[i] = new TH1D(Form("hGenDijetEtaCMForward1DOldPtWeighted_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                               dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DOldPtWeighted[i]->Sumw2();
            hGenDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hGenDijetEtaCMBackward1DOldPt_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DOldPt[i]->Sumw2();
            hGenDijetEtaCMBackward1DOldPtWeighted[i] = new TH1D(Form("hGenDijetEtaCMBackward1DOldPtWeighted_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                                dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DOldPtWeighted[i]->Sumw2();

            // Old eta binning
            hGenDijetEta1DOldPtBinning[i] = new TH1D(Form("hGenDijetEta1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                     dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinning[i]->Sumw2();
            hGenDijetEta1DOldPtBinningWeighted[i] = new TH1D(Form("hGenDijetEta1DOldPtBinningWeighted_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                             dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinningWeighted[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                     fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                             fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Sumw2();
            hGenDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaForward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                             dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaForward1DOldPtBinningWeighted[i] = new TH1D(Form("hGenDijetEtaForward1DOldPtBinningWeighted_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                                 dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaForward1DOldPtBinningWeighted[i]->Sumw2();
            hGenDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaBackward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                             dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hGenDijetEtaBackward1DOldPtBinningWeighted_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                                 dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DOldPtBinningWeighted[i]->Sumw2();

            hGenDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hGenDijetEta1DOldPtBinningCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                     dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinningCM[i]->Sumw2();
            hGenDijetEta1DOldPtBinningCMWeighted[i] = new TH1D(Form("hGenDijetEta1DOldPtBinningCMWeighted_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             dijetEtaOldBins, dijetEtaOldVals);
            hGenDijetEta1DOldPtBinningCMWeighted[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                                     fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                                             fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Sumw2();
            hGenDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaCMForward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaCMForward1DOldPtBinningWeighted[i] = new TH1D(Form("hGenDijetEtaCMForward1DOldPtBinningWeighted_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                                 dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DOldPtBinningWeighted[i]->Sumw2();
            hGenDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hGenDijetEtaCMBackward1DOldPtBinning_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
            hGenDijetEtaCMBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hGenDijetEtaCMBackward1DOldPtBinningWeighted_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                                 dijetEtaFBBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DOldPtBinningWeighted[i]->Sumw2();
        }

        //
        // Ref jet histograms
        //

        hRecoLeadingJetPtOverPtHatVsLeadingJetPt = new TH2D("hRecoLeadingJetPtOverPtHatVsLeadingJetPt", "Reco leading jet p_{T}/#hat{p}_{T} vs reco leading jet p_{T}^{reco};Reco leading jet p_{T} (GeV);Reco leading jet p_{T}/#hat{p}_{T}",
                                                            fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Sumw2();
        hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted = new TH2D("hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted", "Reco leading jet p_{T}/#hat{p}_{T} vs reco leading jet p_{T}^{reco} weighted;Reco leading jet p_{T} (GeV);Reco leading jet p_{T}/#hat{p}_{T}",
                                                                    fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Sumw2();
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

        // Reco to ref dijet correspondence in the laboratory frame [12 dimensions]
        // 0 - Reco dijet ptAve, 1 - dijet eta, 2 - lead pt, 3 - lead eta,
        // 4 - sublead pt, 5 - sublead eta,
        // 6 - Ref dijet ptAve, 7 - dijet eta, 8 - Ref lead pt, 9 - lead eta,
        // 10 - sublead pt, 11 - sublead eta
        hReco2RefDijetLeadSubLead = new THnSparseD("hReco2RefDijetLeadSubLead",
            "Reco to ref correspondence in the lab frame;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet};Reco p_{T}^{Lead} (GeV);Reco #eta^{Lead};Reco p_{T}^{SubLead} (GeV);Reco #eta^{SubLead};Ref p_{T}^{ave} (GeV);Ref #eta^{dijet};Ref p_{T}^{Lead} (GeV);Ref #eta^{Lead};Ref p_{T}^{SubLead} (GeV);Ref #eta^{SubLead}",
            12,
            bins12D_dijet_reco2ref,
            xmin12D_dijet_reco2ref,
            xmax12D_dijet_reco2ref);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLead->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLead->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLead->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLead->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLead->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLead->GetAxis(11)->Set(dijetEtaBins, dijetEtaVals);
        hReco2RefDijetLeadSubLead->Sumw2();

        // Reco to ref dijet correspondence in the center-of-mass frame [12 dimensions]
        // 0 - Reco dijet ptAve, 1 - dijet eta, 2 - lead pt, 3 - lead eta,
        // 4 - sublead pt, 5 - sublead eta,
        // 6 - Ref dijet ptAve, 7 - dijet eta, 8 - Ref lead pt, 9 - lead eta,
        // 10 - sublead pt, 11 - sublead eta
        hReco2RefDijetLeadSubLeadCM = new THnSparseD("hReco2RefDijetLeadSubLeadCM",
            "Reco to ref correspondence in the CM frame;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}_{CM};Reco p_{T}^{Lead} (GeV);Reco #eta^{Lead}_{CM};Reco p_{T}^{SubLead} (GeV);Reco #eta^{SubLead}_{CM};Ref p_{T}^{ave} (GeV);Ref #eta^{dijet}_{CM};Ref p_{T}^{Lead} (GeV);Ref #eta^{Lead}_{CM};Ref p_{T}^{SubLead} (GeV);Ref #eta^{SubLead}_{CM}",
            12,
            bins12D_dijet_reco2ref,
            xmin12D_dijet_reco2ref,
            xmax12D_dijet_reco2ref);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLeadCM->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLeadCM->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLeadCM->GetAxis(5)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLeadCM->GetAxis(7)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLeadCM->GetAxis(9)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hReco2RefDijetLeadSubLeadCM->GetAxis(11)->Set(dijetEtaBins, dijetEtaVals);
        hReco2RefDijetLeadSubLeadCM->Sumw2();

        // Simple reco to ref dijet correspondence in the laboratory frame [4 dimensions]
        // 0 - reco dijet ptAve, 1 - reco dijet eta, 
        // 2 - ref dijet ptAve, 3 - ref dijet eta
        hRecoDijetPtEtaRefDijetPtEta = new THnSparseD("hRecoDijetPtEtaRefDijetPtEta",
            "Reco dijet vs Ref dijet in the lab frame;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}; Ref p_{T}^{ave} (GeV); Ref #eta^{dijet}",
            4,
            bins4D_dijet_reco2ref,
            xmin4D_dijet_reco2ref,
            xmax4D_dijet_reco2ref);
        if (fUseVariableBinning) hRecoDijetPtEtaRefDijetPtEta->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hRecoDijetPtEtaRefDijetPtEta->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaRefDijetPtEta->Sumw2();


        // Simple reco to ref dijet correspondence in the center-of-mass frame [4 dimensions]
        // 0 - reco dijet ptAve, 1 - reco dijet eta,
        // 2 - ref dijet ptAve, 3 - ref dijet eta
        hRecoDijetPtEtaRefDijetPtEtaCM = new THnSparseD("hRecoDijetPtEtaRefDijetPtEtaCM",
            "Reco dijet vs Ref dijet in the CM frame;Reco p_{T}^{ave} (GeV);Reco #eta^{dijet}_{CM}; Ref p_{T}^{ave} (GeV); Ref #eta^{dijet}_{CM}",
            4,
            bins4D_dijet_reco2ref,
            xmin4D_dijet_reco2ref,
            xmax4D_dijet_reco2ref);
        if (fUseVariableBinning) hRecoDijetPtEtaRefDijetPtEtaCM->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        if (fUseVariableBinning) hRecoDijetPtEtaRefDijetPtEtaCM->GetAxis(3)->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetPtEtaRefDijetPtEtaCM->Sumw2();


        hRefPtLeadPtSublead = new TH2D("hRefPtLeadPtSublead","Ref leading vs subleading p_{T};Ref p_{T}^{Lead} (GeV);Ref p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSublead->Sumw2();
        hRefEtaLeadEtaSublead = new TH2D("hRefEtaLeadEtaSublead","Ref leading vs subleading #eta;Ref #eta^{Lead};Ref #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaLeadEtaSublead->Sumw2();
        hRefPtLeadPtSubleadMcReweight = new TH2D("hRefPtLeadPtSubleadMcReweight","Ref leading vs subleading p_{T} (MC reweighted to data);Ref p_{T}^{Lead} (GeV);Ref p_{T}^{SubLead} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSubleadMcReweight->Sumw2();
        hRefEtaLeadEtaSubleadMcReweight = new TH2D("hRefEtaLeadEtaSubleadMcReweight","Ref leading vs subleading #eta (MC reweighted to data);Ref #eta^{Lead};Ref #eta^{SubLead}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaLeadEtaSubleadMcReweight->Sumw2();


        hRefPtLeadPtSubleadCM = new TH2D("hRefPtLeadPtSubleadCM","Ref leading vs subleading p_{T} in CM;Ref p_{T}^{Lead}_{CM} (GeV);Ref p_{T}^{SubLead}_{CM} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSubleadCM->Sumw2();
        hRefEtaCMLeadEtaCMSublead = new TH2D("hRefEtaCMLeadEtaCMSublead","Ref leading vs subleading #eta in CM;Ref #eta^{Lead}_{CM};Ref #eta^{SubLead}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaCMLeadEtaCMSublead->Sumw2();
        hRefPtLeadPtSubleadCMMcReweight = new TH2D("hRefPtLeadPtSubleadCMMcReweight","Ref leading vs subleading p_{T} in CM (MC reweighted to data);Ref p_{T}^{Lead}_{CM} (GeV);Ref p_{T}^{SubLead}_{CM} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefPtLeadPtSubleadCMMcReweight->Sumw2();
        hRefEtaCMLeadEtaCMSubleadMcReweight = new TH2D("hRefEtaCMLeadEtaCMSubleadMcReweight","Ref leading vs subleading #eta in CM (MC reweighted to data);Ref #eta^{Lead}_{CM};Ref #eta^{SubLead}_{CM}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefEtaCMLeadEtaCMSubleadMcReweight->Sumw2();


        hRefDijetEta = new TH1D("hRefDijetEta","Ref dijet #eta;Ref #eta^{dijet};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEta->Sumw2();
        hRefDijetEtaVsRecoDijetEta = new TH2D("hRefDijetEtaVsRecoDijetEta","Ref dijet #eta vs reco dijet #eta;Reco #eta^{dijet};Ref #eta^{dijet}",
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaVsRecoDijetEta->Sumw2();
        hRefDijetPtVsRecoDijetPt = new TH2D("hRefDijetPtVsRecoDijetPt","Ref dijet p_{T} vs reco dijet p_{T};Reco p_{T}^{dijet} (GeV);Ref p_{T}^{dijet} (GeV)",
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRefDijetPtVsRecoDijetPt->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T};Reco #eta^{dijet};Ref #eta^{dijet}; Reco p_{T}^{dijet} (GeV)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Sumw2();
        hRefDijetPtEtaPhi = new TH3D("hRefDijetPtEtaPhi","Ref dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefDijetPtEtaPhi->GetYaxis()->Set(dijetPtBins, dijetPtVals);
        hRefDijetPtEtaPhi->Sumw2();
        hRefDijetPtEtaPhiWeighted = new TH3D("hRefDijetPtEtaPhiWeighted","Ref dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetPtBins, dijetPtVals);
        hRefDijetPtEtaPhiWeighted->Sumw2();
        hRefDijetPtEtaForward = new TH2D("hRefDijetPtEtaForward", "Ref dijet info in lab frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaForward->Sumw2();
        hRefDijetPtEtaBackward = new TH2D("hRefDijetPtEtaBackward", "Ref dijet info in lab frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaBackward->Sumw2();
        hRefDijetPtEtaForwardWeighted = new TH2D("hRefDijetPtEtaForwardWeighted", "Ref dijet info in lab frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaForwardWeighted->Sumw2();
        hRefDijetPtEtaBackwardWeighted = new TH2D("hRefDijetPtEtaBackwardWeighted", "Ref dijet info in lab frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaBackwardWeighted->Sumw2();


        hRefDijetEtaCM = new TH1D("hRefDijetEtaCM","Ref dijet #eta in CM;Ref #eta^{dijet}_{CM};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaCM = new TH2D("hRefDijetEtaVsRecoDijetEtaCM","Ref dijet #eta vs reco dijet #eta in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaVsRecoDijetEtaCM->Sumw2();
        hRefDijetPtVsRecoDijetPtCM = new TH2D("hRefDijetPtVsRecoDijetPtCM","Ref dijet p_{T} vs reco dijet p_{T} in CM;Reco p_{T}^{dijet}_{CM} (GeV);Ref p_{T}^{dijet}_{CM} (GeV)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1]);
        hRefDijetPtVsRecoDijetPtCM->Sumw2();

        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Sumw2();
        hRefDijetPtEtaPhiCM = new TH3D("hRefDijetPtEtaPhiCM","Ref dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefDijetPtEtaPhiCM->GetYaxis()->Set(dijetPtBins, dijetPtVals);
        hRefDijetPtEtaPhiCM->Sumw2();
        hRefDijetPtEtaPhiCMWeighted = new TH3D("hRefDijetPtEtaPhiCMWeighted","Ref dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetPtBins, dijetPtVals);
        hRefDijetPtEtaPhiCMWeighted->Sumw2();
        hRefDijetPtEtaCMForward = new TH2D("hRefDijetPtEtaCMForward", "Ref dijet info in CM frame (forward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaCMForward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMForward->Sumw2();
        hRefDijetPtEtaCMBackward = new TH2D("hRefDijetPtEtaCMBackward", "Ref dijet info in CM frame (backward);p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaCMBackward->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMBackward->Sumw2();
        hRefDijetPtEtaCMForwardWeighted = new TH2D("hRefDijetPtEtaCMForwardWeighted", "Ref dijet info in CM frame (forward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaCMForwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMForwardWeighted->Sumw2();

        hRefDijetPtEtaCMBackwardWeighted = new TH2D("hRefDijetPtEtaCMBackwardWeighted", "Ref dijet info in CM frame (backward) weighted;p_{T}^{ave} (GeV);#eta^{dijet}_{CM}",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, 0., fDijetEtaRange[1]);
        if (fUseVariableBinning) hRefDijetPtEtaCMBackwardWeighted->GetYaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRefDijetPtEtaCMBackwardWeighted->Sumw2();       

        // New pT and eta binning
        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hRefDijetEta1D[i] = new TH1D(Form("hRefDijetEta1D_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1D[i]->Sumw2();
            hRefDijetEta1DWeighted[i] = new TH1D(Form("hRefDijetEta1DWeighted_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DWeighted[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2D_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DWeighted_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2D[i] = new TH2D(Form("hRecoVsRefDijetEta2D_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2D[i]->Sumw2();
            hRecoVsRefDijetEta2DWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DWeighted_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2D[i] = new TH2D(Form("hRecoVsRefLeadJetEta2D_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2D[i]->Sumw2();
            hRecoVsRefLeadJetEta2DWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DWeighted_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2D[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2D_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2D[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DWeighted_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DWeighted[i]->Sumw2();
            hRefDijetEtaForward1D[i] = new TH1D(Form("hRefDijetEtaForward1D_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1D[i]->Sumw2();
            hRefDijetEtaForward1DWeighted[i] = new TH1D(Form("hRefDijetEtaForward1DWeighted_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DWeighted[i]->Sumw2();
            hRefDijetEtaBackward1D[i] = new TH1D(Form("hRefDijetEtaBackward1D_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1D[i]->Sumw2();
            hRefDijetEtaBackward1DWeighted[i] = new TH1D(Form("hRefDijetEtaBackward1DWeighted_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DWeighted[i]->Sumw2();

            // Ref dijets in CM frame
            hRefDijetEta1DCM[i] = new TH1D(Form("hRefDijetEta1DCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                           prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DCM[i]->Sumw2();
            hRefDijetEta1DCMWeighted[i] = new TH1D(Form("hRefDijetEta1DCMWeighted_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DCMWeighted[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DCM_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DCMWeighted_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DCM[i]->Sumw2();
            hRecoVsRefDijetEta2DCMWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DCMWeighted_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DCMWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DCM_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DCMWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DCMWeighted_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DCMWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DCM_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DCMWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DCMWeighted_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DCMWeighted[i]->Sumw2();
            hRefDijetEtaCMForward1D[i] = new TH1D(Form("hRefDijetEtaCMForward1D_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                  fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1D[i]->Sumw2();
            hRefDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hRefDijetEtaCMForward1DWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DWeighted[i]->Sumw2();
            hRefDijetEtaCMBackward1D[i] = new TH1D(Form("hRefDijetEtaCMBackward1D_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1D[i]->Sumw2();
            hRefDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hRefDijetEtaCMBackward1DWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DWeighted[i]->Sumw2();
        }

        // Old pT binning
        for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {

            double ptAveLow = fPtAveOldBins.at(i);
            double ptAveHi = fPtAveOldBins.at(i+1);
            // New eta binning
            hRefDijetEta1DOldPt[i] = new TH1D(Form("hRefDijetEta1DOldPt_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPt[i]->Sumw2();
            hRefDijetEta1DOldPtWeighted[i] = new TH1D(Form("hRefDijetEta1DOldPtWeighted_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                      dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPtWeighted[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtWeighted_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPt[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtWeighted_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPt[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtWeighted_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPt[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtWeighted_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtWeighted[i]->Sumw2();
            hRefDijetEtaForward1DOldPt[i] = new TH1D(Form("hRefDijetEtaForward1DOldPt_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DOldPt[i]->Sumw2();
            hRefDijetEtaForward1DOldPtWeighted[i] = new TH1D(Form("hRefDijetEtaForward1DOldPtWeighted_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DOldPtWeighted[i]->Sumw2();
            hRefDijetEtaBackward1DOldPt[i] = new TH1D(Form("hRefDijetEtaBackward1DOldPt_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DOldPt[i]->Sumw2();
            hRefDijetEtaBackward1DOldPtWeighted[i] = new TH1D(Form("hRefDijetEtaBackward1DOldPtWeighted_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DOldPtWeighted[i]->Sumw2();


            hRefDijetEta1DOldPtCM[i] = new TH1D(Form("hRefDijetEta1DOldPtCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPtCM[i]->Sumw2();
            hRefDijetEta1DOldPtCMWeighted[i] = new TH1D(Form("hRefDijetEta1DOldPtCMWeighted_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                      dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DOldPtCMWeighted[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtCMWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtCMWeighted_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtCMWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtCMWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtCMWeighted_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtCMWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtCMWeighted_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i]->Sumw2();
            hRefDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hRefDijetEtaCMForward1DOldPt_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DOldPt[i]->Sumw2();
            hRefDijetEtaCMForward1DOldPtWeighted[i] = new TH1D(Form("hRefDijetEtaCMForward1DOldPtWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DOldPtWeighted[i]->Sumw2();
            hRefDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hRefDijetEtaCMBackward1DOldPt_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DOldPt[i]->Sumw2();
            hRefDijetEtaCMBackward1DOldPtWeighted[i] = new TH1D(Form("hRefDijetEtaCMBackward1DOldPtWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DOldPtWeighted[i]->Sumw2();

            // Old eta binning
            hRefDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefDijetEta1DOldPtBinning_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinning[i]->Sumw2();
            hRefDijetEta1DOldPtBinningWeighted[i] = new TH1D(Form("hRefDijetEta1DOldPtBinningWeighted_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinningWeighted[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted_%d",i),Form("Ref #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningWeighted_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);                                  
            hRecoVsRefLeadJetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningWeighted_%d",i),Form("Reco vs Ref #eta^{Lead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Lead};Ref #eta^{Lead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted_%d",i),Form("Reco vs Ref #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{SubLead};Ref #eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i]->Sumw2();
            hRefDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaForward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DOldPtBinning[i]->Sumw2();
            hRefDijetEtaForward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefDijetEtaForward1DOldPtBinningWeighted_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DOldPtBinningWeighted[i]->Sumw2();
            hRefDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaBackward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DOldPtBinning[i]->Sumw2();
            hRefDijetEtaBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefDijetEtaBackward1DOldPtBinningWeighted_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DOldPtBinningWeighted[i]->Sumw2();


            hRefDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRefDijetEta1DOldPtBinningCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinningCM[i]->Sumw2();
            hRefDijetEta1DOldPtBinningCMWeighted[i] = new TH1D(Form("hRefDijetEta1DOldPtBinningCMWeighted_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                               dijetEtaOldBins, dijetEtaOldVals);
            hRefDijetEta1DOldPtBinningCMWeighted[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted_%d",i),Form("Ref #eta^{Lead}_{CM} vs #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningCMWeighted_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted_%d",i),Form("Reco vs Ref #eta^{Lead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Lead}_{CM};Ref #eta^{Lead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted_%d",i),Form("Reco vs Ref #eta^{SubLead}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{SubLead}_{CM};Ref #eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i]->Sumw2();
            hRefDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaCMForward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
            hRefDijetEtaCMForward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefDijetEtaCMForward1DOldPtBinningWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DOldPtBinningWeighted[i]->Sumw2();
            hRefDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hRefDijetEtaCMBackward1DOldPtBinning_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
            hRefDijetEtaCMBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefDijetEtaCMBackward1DOldPtBinningWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DOldPtBinningWeighted[i]->Sumw2();
        }

        //
        // Ref-selected jet histograms
        //


        hRefSelDijetEta = new TH1D("hRefSelDijetEta","Ref selected dijets;#eta^{dijet}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEta->Sumw2();
        hRefSelDijetPtEtaPhi = new TH3D("hRefSelDijetPtEtaPhi","Ref selected dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefSelDijetPtEtaPhi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaPhi->Sumw2();
        hRefSelDijetPtEtaPhiWeighted = new TH3D("hRefSelDijetPtEtaPhiWeighted","Ref selected dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefSelDijetPtEtaPhiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaPhiWeighted->Sumw2();

        hRefSelDijetEtaCM = new TH1D("hRefSelDijetEtaCM","Ref selected dijets in CM;#eta^{dijet}_{CM}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEtaCM->Sumw2();
        hRefSelDijetPtEtaPhiCM = new TH3D("hRefSelDijetPtEtaPhiCM","Ref selected dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefSelDijetPtEtaPhiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaPhiCM->Sumw2();
        hRefSelDijetPtEtaPhiCMWeighted = new TH3D("hRefSelDijetPtEtaPhiCMWeighted","Ref selected dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fPhiBins, fPhiRange[0], fPhiRange[1] );
        if (fUseVariableBinning) hRefSelDijetPtEtaPhiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRefSelDijetPtEtaPhiCMWeighted->Sumw2();

        // New pT and eta binning
        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hRefSelDijetEta1D[i] = new TH1D(Form("hRefSelDijetEta1D_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                            prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1D[i]->Sumw2();
            hRefSelDijetEta1DWeighted[i] = new TH1D(Form("hRefSelDijetEta1DWeighted_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                  prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DWeighted[i]->Sumw2();
            hRefSelRecoDijetEta1D[i] = new TH1D(Form("hRefSelRecoDijetEta1D_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1D[i]->Sumw2();
            hRefSelRecoDijetEta1DWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                  prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelRecoDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DWeighted[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2D_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DWeighted_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DWeighted[i]->Sumw2();
            hRefSelDijetEtaForward1D[i] = new TH1D(Form("hRefSelDijetEtaForward1D_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1D[i]->Sumw2();
            hRefSelDijetEtaForward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaForward1DWeighted_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DWeighted[i]->Sumw2();
            hRefSelDijetEtaBackward1D[i] = new TH1D(Form("hRefSelDijetEtaBackward1D_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1D[i]->Sumw2();
            hRefSelDijetEtaBackward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaBackward1DWeighted_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                            fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DWeighted[i]->Sumw2();


            hRefSelDijetEta1DCM[i] = new TH1D(Form("hRefSelDijetEta1DCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                              prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DCM[i]->Sumw2();
            hRefSelDijetEta1DCMWeighted[i] = new TH1D(Form("hRefSelDijetEta1DCMWeighted_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                      prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DCMWeighted[i]->Sumw2();
            hRefSelRecoDijetEta1DCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                  prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelRecoDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DCM[i]->Sumw2();
            hRefSelRecoDijetEta1DCMWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DCMWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          prescale * fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelRecoDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DCMWeighted[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DCM_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DCMWeighted_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
            hRefSelDijetEtaCMForward1D[i] = new TH1D(Form("hRefSelDijetEtaCMForward1D_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                     fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1D[i]->Sumw2();
            hRefSelDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DWeighted[i]->Sumw2();
            hRefSelDijetEtaCMBackward1D[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1D_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                           fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1D[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             fDijetEtaBins, 0., fDijetEtaRange[1]);
            if (fUseVariableBinning) hRefSelDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DWeighted[i]->Sumw2();
        } // for (int i{0}; i<fPtAveBins.size()-1; i++)

        // Old pT binning
        for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {
            double ptAveLow = fPtAveOldBins.at(i);
            double ptAveHi = fPtAveOldBins.at(i+1);

            // New eta binning
            hRefSelDijetEta1DOldPt[i] = new TH1D(Form("hRefSelDijetEta1DOldPt_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DOldPt[i]->Sumw2();
            hRefSelDijetEta1DOldPtWeighted[i] = new TH1D(Form("hRefSelDijetEta1DOldPtWeighted_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DOldPtWeighted[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPt[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPt_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DOldPt[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DOldPtWeighted[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Sumw2();
            hRefSelDijetEtaForward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaForward1DOldPt_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DOldPt[i]->Sumw2();
            hRefSelDijetEtaForward1DOldPtWeighted[i] = new TH1D(Form("hRefSelDijetEtaForward1DOldPtWeighted_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DOldPtWeighted[i]->Sumw2();
            hRefSelDijetEtaBackward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaBackward1DOldPt_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DOldPt[i]->Sumw2();
            hRefSelDijetEtaBackward1DOldPtWeighted[i] = new TH1D(Form("hRefSelDijetEtaBackward1DOldPtWeighted_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DOldPtWeighted[i]->Sumw2();

            hRefSelDijetEta1DOldPtCM[i] = new TH1D(Form("hRefSelDijetEta1DOldPtCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DOldPtCM[i]->Sumw2();
            hRefSelDijetEta1DOldPtCMWeighted[i] = new TH1D(Form("hRefSelDijetEta1DOldPtCMWeighted_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DOldPtCMWeighted[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DOldPtCM[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtCMWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtCMWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                         dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DOldPtCMWeighted[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Sumw2();
            hRefSelDijetEtaCMForward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DOldPt_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DOldPt[i]->Sumw2();
            hRefSelDijetEtaCMForward1DOldPtWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DOldPtWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DOldPtWeighted[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DOldPt[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DOldPt_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DOldPt[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DOldPtWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DOldPtWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DOldPtWeighted[i]->Sumw2();
            
            // Old eta binning
            hRefSelDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEta1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelDijetEta1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEta1DOldPtBinningWeighted[i] = new TH1D(Form("hRefSelDijetEta1DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelDijetEta1DOldPtBinningWeighted[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtBinning[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtBinning_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelRecoDijetEta1DOldPtBinning[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtBinningWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtBinningWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelRecoDijetEta1DOldPtBinningWeighted[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Sumw2();
            hRefSelDijetEtaForward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaForward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                           dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaForward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefSelDijetEtaForward1DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                           dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaForward1DOldPtBinningWeighted[i]->Sumw2();
            hRefSelDijetEtaBackward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaBackward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                           dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefSelDijetEtaBackward1DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                           dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaBackward1DOldPtBinningWeighted[i]->Sumw2();


            hRefSelDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRefSelDijetEta1DOldPtBinningCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelDijetEta1DOldPtBinningCM[i]->Sumw2();
            hRefSelDijetEta1DOldPtBinningCMWeighted[i] = new TH1D(Form("hRefSelDijetEta1DOldPtBinningCMWeighted_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelDijetEta1DOldPtBinningCMWeighted[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtBinningCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtBinningCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelRecoDijetEta1DOldPtBinningCM[i]->Sumw2();
            hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DOldPtBinningCMWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals);
            hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted_%d",i),Form("Ref selected #eta^{Lead} vs #eta^{SubLead} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Sumw2();
            hRefSelDijetEtaCMForward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DOldPtBinning[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DOldPtBinning_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DOldPtBinning[i]->Sumw2();
            hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   dijetEtaFBBins, dijetEtaFBVals);
            hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i]->Sumw2();

        } // for (int i{0}; i<fPtAveOldBins.size()-1; i++)
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



    hRecoJetCollectionSize->Write();
    hRecoDijetLeadSubLead->Write();
    hRecoDijetLeadSubLeadCM->Write();

    hRecoPtLeadPtSublead->Write();
    hRecoEtaLeadEtaSublead->Write();
    hRecoPtLeadPtSubleadMcReweight->Write();
    hRecoEtaLeadEtaSubleadMcReweight->Write();
    hRecoDijetEta->Write();
    hRecoDijetPtEta->Write();
    hRecoDijetPtEtaForward->Write();
    hRecoDijetPtEtaBackward->Write();
    hRecoDijetPtEtaForwardWeighted->Write();
    hRecoDijetPtEtaBackwardWeighted->Write();
    hRecoDijetPtEtaPhi->Write();
    hRecoDijetPtEtaPhiWeighted->Write();

    hRecoPtLeadPtSubleadCM->Write();
    hRecoEtaCMLeadEtaCMSublead->Write();
    hRecoPtLeadPtSubleadCMMcReweight->Write();
    hRecoEtaCMLeadEtaCMSubleadMcReweight->Write();
    hRecoDijetEtaCM->Write();
    hRecoDijetPtEtaCM->Write();
    hRecoDijetPtEtaCMForward->Write();
    hRecoDijetPtEtaCMBackward->Write();
    hRecoDijetPtEtaCMForwardWeighted->Write();
    hRecoDijetPtEtaCMBackwardWeighted->Write();
    hRecoDijetPtEtaPhiCM->Write();
    hRecoDijetPtEtaPhiCMWeighted->Write();


    for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
        hRecoDijetEta1D[i]->Write();
        hRecoDijetEta1DWeighted[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2D[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i]->Write();
        hRecoDijetEtaForward1D[i]->Write();
        hRecoDijetEtaForward1DWeighted[i]->Write();
        hRecoDijetEtaBackward1D[i]->Write();
        hRecoDijetEtaBackward1DWeighted[i]->Write();

        hRecoDijetEta1DCM[i]->Write();
        hRecoDijetEta1DCMWeighted[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DCM[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DCMWeighted[i]->Write();
        hRecoDijetEtaCMForward1D[i]->Write();
        hRecoDijetEtaCMForward1DWeighted[i]->Write();
        hRecoDijetEtaCMBackward1D[i]->Write();
        hRecoDijetEtaCMBackward1DWeighted[i]->Write();
    }

    for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
        hRecoDijetEta1DOldPt[i]->Write();
        hRecoDijetEta1DOldPtWeighted[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Write();
        hRecoDijetEtaForward1DOldPt[i]->Write();
        hRecoDijetEtaForward1DOldPtWeighted[i]->Write();
        hRecoDijetEtaBackward1DOldPt[i]->Write();
        hRecoDijetEtaBackward1DOldPtWeighted[i]->Write();
        hRecoDijetEta1DOldPtCM[i]->Write();
        hRecoDijetEta1DOldPtCMWeighted[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Write();
        hRecoDijetEtaCMForward1DOldPt[i]->Write();
        hRecoDijetEtaCMForward1DOldPtWeighted[i]->Write();
        hRecoDijetEtaCMBackward1DOldPt[i]->Write();
        hRecoDijetEtaCMBackward1DOldPtWeighted[i]->Write();
        hRecoDijetEta1DOldPtBinning[i]->Write();
        hRecoDijetEta1DOldPtBinningWeighted[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Write();
        hRecoDijetEtaForward1DOldPtBinning[i]->Write();
        hRecoDijetEtaForward1DOldPtBinningWeighted[i]->Write();
        hRecoDijetEtaBackward1DOldPtBinning[i]->Write();
        hRecoDijetEtaBackward1DOldPtBinningWeighted[i]->Write();
        hRecoDijetEta1DOldPtBinningCM[i]->Write();
        hRecoDijetEta1DOldPtBinningCMWeighted[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Write();
        hRecoDijetEtaCMForward1DOldPtBinning[i]->Write();
        hRecoDijetEtaCMForward1DOldPtBinningWeighted[i]->Write();
        hRecoDijetEtaCMBackward1DOldPtBinning[i]->Write();
        hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i]->Write();
    }


    if ( fIsMc ) {

        //
        // Gen histograms
        //

        hGenJetCollectionSize->Write();
        hGenVsRecoJetCollectionSize->Write();

        hGenLeadingJetPtOverPtHatVsLeadingJetPt->Write();
        hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Write();
        hGenDijetPtOverPtHatVsDijetPt->Write();
        hGenDijetPtOverPtHatVsDijetPtWeighted->Write();
        hGenDijetPtAveOverPtHatVsDijetPtAve->Write();
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Write();

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

        hGenDijetLeadSubLead->Write();
        hGenDijetLeadSubLeadCM->Write();

        hGenPtLeadPtSublead->Write();
        hGenEtaLeadEtaSublead->Write();
        hGenDijetEta->Write();
        hGenDijetPtEtaPhi->Write();
        hGenDijetPtEtaPhiWeighted->Write();

        hGenPtLeadPtSubleadCM->Write();
        hGenEtaCMLeadEtaCMSublead->Write();
        hGenDijetEtaCM->Write();
        hGenDijetPtEtaPhiCM->Write();
        hGenDijetPtEtaPhiCMWeighted->Write();

        hGenDijetPtEtaForward->Write();
        hGenDijetPtEtaBackward->Write();
        hGenDijetPtEtaCMForward->Write();
        hGenDijetPtEtaCMBackward->Write();
        hGenDijetPtEtaForwardWeighted->Write();
        hGenDijetPtEtaBackwardWeighted->Write();
        hGenDijetPtEtaCMForwardWeighted->Write();
        hGenDijetPtEtaCMBackwardWeighted->Write();

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            hGenDijetEta1D[i]->Write();
            hGenDijetEta1DWeighted[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2D[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DWeighted[i]->Write();
            hGenDijetEtaForward1D[i]->Write();
            hGenDijetEtaForward1DWeighted[i]->Write();
            hGenDijetEtaBackward1D[i]->Write();
            hGenDijetEtaBackward1DWeighted[i]->Write();
            hGenDijetEta1DCM[i]->Write();
            hGenDijetEta1DCMWeighted[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DCM[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i]->Write();
            hGenDijetEtaCMForward1D[i]->Write();
            hGenDijetEtaCMForward1DWeighted[i]->Write();
            hGenDijetEtaCMBackward1D[i]->Write();
            hGenDijetEtaCMBackward1DWeighted[i]->Write();
        }

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            hGenDijetEta1DOldPt[i]->Write();
            hGenDijetEta1DOldPtWeighted[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Write();
            hGenDijetEtaForward1DOldPt[i]->Write();
            hGenDijetEtaForward1DOldPtWeighted[i]->Write();
            hGenDijetEtaBackward1DOldPt[i]->Write();
            hGenDijetEtaBackward1DOldPtWeighted[i]->Write();
            hGenDijetEta1DOldPtCM[i]->Write();
            hGenDijetEta1DOldPtCMWeighted[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Write();
            hGenDijetEtaCMForward1DOldPt[i]->Write();
            hGenDijetEtaCMForward1DOldPtWeighted[i]->Write();
            hGenDijetEtaCMBackward1DOldPt[i]->Write();
            hGenDijetEtaCMBackward1DOldPtWeighted[i]->Write();
            hGenDijetEta1DOldPtBinning[i]->Write();
            hGenDijetEta1DOldPtBinningWeighted[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Write();
            hGenDijetEtaForward1DOldPtBinning[i]->Write();
            hGenDijetEtaForward1DOldPtBinningWeighted[i]->Write();
            hGenDijetEtaBackward1DOldPtBinning[i]->Write();
            hGenDijetEtaBackward1DOldPtBinningWeighted[i]->Write();
            hGenDijetEta1DOldPtBinningCM[i]->Write();
            hGenDijetEta1DOldPtBinningCMWeighted[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
            hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Write();
            hGenDijetEtaCMForward1DOldPtBinning[i]->Write();
            hGenDijetEtaCMForward1DOldPtBinningWeighted[i]->Write();
            hGenDijetEtaCMBackward1DOldPtBinning[i]->Write();
            hGenDijetEtaCMBackward1DOldPtBinningWeighted[i]->Write();
        }

        //
        // Ref hisrograms
        //

        hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Write();
        hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Write();
        hRecoDijetPtOverPtHatVsDijetPt->Write();
        hRecoDijetPtOverPtHatVsDijetPtWeighted->Write();
        hRecoDijetPtAveOverPtHatVsDijetPtAve->Write();
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Write();

        hReco2RefDijetLeadSubLead->Write();
        hReco2RefDijetLeadSubLeadCM->Write();
        hRecoDijetPtEtaRefDijetPtEta->Write();
        hRecoDijetPtEtaRefDijetPtEtaCM->Write();

        hRefPtLeadPtSublead->Write();
        hRefEtaLeadEtaSublead->Write();
        hRefPtLeadPtSubleadMcReweight->Write();
        hRefEtaLeadEtaSubleadMcReweight->Write();

        hRefPtLeadPtSubleadCM->Write();
        hRefEtaCMLeadEtaCMSublead->Write();
        hRefPtLeadPtSubleadCMMcReweight->Write();
        hRefEtaCMLeadEtaCMSubleadMcReweight->Write();

        hRefDijetEta->Write();
        hRefDijetEtaVsRecoDijetEta->Write();
        hRefDijetPtVsRecoDijetPt->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Write();
        hRefDijetPtEtaPhi->Write();
        hRefDijetPtEtaPhiWeighted->Write();
        hRefDijetPtEtaForward->Write();
        hRefDijetPtEtaBackward->Write();
        hRefDijetPtEtaForwardWeighted->Write();
        hRefDijetPtEtaBackwardWeighted->Write();

        hRefDijetEtaCM->Write();
        hRefDijetEtaVsRecoDijetEtaCM->Write();
        hRefDijetPtVsRecoDijetPtCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Write();
        hRefDijetPtEtaPhiCM->Write();
        hRefDijetPtEtaPhiCMWeighted->Write();
        hRefDijetPtEtaCMForward->Write();
        hRefDijetPtEtaCMBackward->Write();
        hRefDijetPtEtaCMForwardWeighted->Write();
        hRefDijetPtEtaCMBackwardWeighted->Write();

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            hRefDijetEta1D[i]->Write();
            hRefDijetEta1DWeighted[i]->Write();
            hRefEtaLeadVsEtaSubLead2D[i]->Write();
            hRefEtaLeadVsEtaSubLead2DWeighted[i]->Write();
            hRecoVsRefDijetEta2D[i]->Write();
            hRecoVsRefDijetEta2DWeighted[i]->Write();
            hRecoVsRefLeadJetEta2D[i]->Write();
            hRecoVsRefLeadJetEta2DWeighted[i]->Write();
            hRecoVsRefSubLeadJetEta2D[i]->Write();
            hRecoVsRefSubLeadJetEta2DWeighted[i]->Write();
            hRefDijetEtaForward1D[i]->Write();
            hRefDijetEtaForward1DWeighted[i]->Write();
            hRefDijetEtaBackward1D[i]->Write();
            hRefDijetEtaBackward1DWeighted[i]->Write();

            hRefDijetEta1DCM[i]->Write();
            hRefDijetEta1DCMWeighted[i]->Write();
            hRefEtaLeadVsEtaSubLead2DCM[i]->Write();
            hRefEtaLeadVsEtaSubLead2DCMWeighted[i]->Write();
            hRecoVsRefDijetEta2DCM[i]->Write();
            hRecoVsRefDijetEta2DCMWeighted[i]->Write();
            hRecoVsRefLeadJetEta2DCM[i]->Write();
            hRecoVsRefLeadJetEta2DCMWeighted[i]->Write();
            hRecoVsRefSubLeadJetEta2DCM[i]->Write();
            hRecoVsRefSubLeadJetEta2DCMWeighted[i]->Write();
            hRefDijetEtaCMForward1D[i]->Write();
            hRefDijetEtaCMForward1DWeighted[i]->Write();
            hRefDijetEtaCMBackward1D[i]->Write();
            hRefDijetEtaCMBackward1DWeighted[i]->Write();
        }

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            hRefDijetEta1DOldPt[i]->Write();
            hRefDijetEta1DOldPtWeighted[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPt[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Write();
            hRecoVsRefDijetEta2DOldPt[i]->Write();
            hRecoVsRefDijetEta2DOldPtWeighted[i]->Write();
            hRecoVsRefLeadJetEta2DOldPt[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtWeighted[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPt[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtWeighted[i]->Write();
            hRefDijetEtaForward1DOldPt[i]->Write();
            hRefDijetEtaForward1DOldPtWeighted[i]->Write();
            hRefDijetEtaBackward1DOldPt[i]->Write();
            hRefDijetEtaBackward1DOldPtWeighted[i]->Write();

            hRefDijetEta1DOldPtCM[i]->Write();
            hRefDijetEta1DOldPtCMWeighted[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Write();
            hRecoVsRefDijetEta2DOldPtCM[i]->Write();
            hRecoVsRefDijetEta2DOldPtCMWeighted[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtCM[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtCMWeighted[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtCM[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i]->Write();
            hRefDijetEtaCMForward1DOldPt[i]->Write();
            hRefDijetEtaCMForward1DOldPtWeighted[i]->Write();
            hRefDijetEtaCMBackward1DOldPt[i]->Write();
            hRefDijetEtaCMBackward1DOldPtWeighted[i]->Write();

            hRefDijetEta1DOldPtBinning[i]->Write();
            hRefDijetEta1DOldPtBinningWeighted[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Write();
            hRecoVsRefDijetEta2DOldPtBinning[i]->Write();
            hRecoVsRefDijetEta2DOldPtBinningWeighted[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtBinning[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i]->Write();
            hRefDijetEtaForward1DOldPtBinning[i]->Write();
            hRefDijetEtaForward1DOldPtBinningWeighted[i]->Write();
            hRefDijetEtaBackward1DOldPtBinning[i]->Write();
            hRefDijetEtaBackward1DOldPtBinningWeighted[i]->Write();

            hRefDijetEta1DOldPtBinningCM[i]->Write();
            hRefDijetEta1DOldPtBinningCMWeighted[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Write();
            hRecoVsRefDijetEta2DOldPtBinningCM[i]->Write();
            hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i]->Write();
            hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]->Write();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i]->Write();
            hRefDijetEtaCMForward1DOldPtBinning[i]->Write();
            hRefDijetEtaCMForward1DOldPtBinningWeighted[i]->Write();
            hRefDijetEtaCMBackward1DOldPtBinning[i]->Write();
            hRefDijetEtaCMBackward1DOldPtBinningWeighted[i]->Write();
        }

        //
        // Ref-selected histograms
        //

        hRefSelDijetEta->Write();
        hRefSelDijetPtEtaPhi->Write();
        hRefSelDijetPtEtaPhiWeighted->Write();

        hRefSelDijetEtaCM->Write();
        hRefSelDijetPtEtaPhiCM->Write();
        hRefSelDijetPtEtaPhiCMWeighted->Write();

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            hRefSelDijetEta1D[i]->Write();
            hRefSelDijetEta1DWeighted[i]->Write();
            hRefSelRecoDijetEta1D[i]->Write();
            hRefSelRecoDijetEta1DWeighted[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2D[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DWeighted[i]->Write();
            hRefSelDijetEtaForward1D[i]->Write();
            hRefSelDijetEtaForward1DWeighted[i]->Write();
            hRefSelDijetEtaBackward1D[i]->Write();
            hRefSelDijetEtaBackward1DWeighted[i]->Write();
            hRefSelDijetEta1DCM[i]->Write();
            hRefSelDijetEta1DCMWeighted[i]->Write();
            hRefSelRecoDijetEta1DCM[i]->Write();
            hRefSelRecoDijetEta1DCMWeighted[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DCM[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i]->Write();
            hRefSelDijetEtaCMForward1D[i]->Write();
            hRefSelDijetEtaCMForward1DWeighted[i]->Write();
            hRefSelDijetEtaCMBackward1D[i]->Write();
            hRefSelDijetEtaCMBackward1DWeighted[i]->Write();
        }

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            hRefSelDijetEta1DOldPt[i]->Write();
            hRefSelDijetEta1DOldPtWeighted[i]->Write();
            hRefSelRecoDijetEta1DOldPt[i]->Write();
            hRefSelRecoDijetEta1DOldPtWeighted[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Write();
            hRefSelDijetEtaForward1DOldPt[i]->Write();
            hRefSelDijetEtaForward1DOldPtWeighted[i]->Write();
            hRefSelDijetEtaBackward1DOldPt[i]->Write();
            hRefSelDijetEtaBackward1DOldPtWeighted[i]->Write();
            hRefSelDijetEta1DOldPtCM[i]->Write();
            hRefSelDijetEta1DOldPtCMWeighted[i]->Write();
            hRefSelRecoDijetEta1DOldPtCM[i]->Write();
            hRefSelRecoDijetEta1DOldPtCMWeighted[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Write();
            hRefSelDijetEtaCMForward1DOldPt[i]->Write();
            hRefSelDijetEtaCMForward1DOldPtWeighted[i]->Write();
            hRefSelDijetEtaCMBackward1DOldPt[i]->Write();
            hRefSelDijetEtaCMBackward1DOldPtWeighted[i]->Write();
            hRefSelDijetEta1DOldPtBinning[i]->Write();
            hRefSelDijetEta1DOldPtBinningWeighted[i]->Write();
            hRefSelRecoDijetEta1DOldPtBinning[i]->Write();
            hRefSelRecoDijetEta1DOldPtBinningWeighted[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Write();
            hRefSelDijetEtaForward1DOldPtBinning[i]->Write();
            hRefSelDijetEtaForward1DOldPtBinningWeighted[i]->Write();
            hRefSelDijetEtaBackward1DOldPtBinning[i]->Write();
            hRefSelDijetEtaBackward1DOldPtBinningWeighted[i]->Write();
            hRefSelDijetEta1DOldPtBinningCM[i]->Write();
            hRefSelDijetEta1DOldPtBinningCMWeighted[i]->Write();
            hRefSelRecoDijetEta1DOldPtBinningCM[i]->Write();
            hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Write();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Write();
            hRefSelDijetEtaCMForward1DOldPtBinning[i]->Write();
            hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i]->Write();
            hRefSelDijetEtaCMBackward1DOldPtBinning[i]->Write();
            hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i]->Write();
        }
    } // if ( fIsMc )
}
