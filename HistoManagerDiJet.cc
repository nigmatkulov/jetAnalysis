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
    //
    // Variables
    //
    fIsMc{false}, 
    fPtBins{150}, fPtRange{5., 1505.}, 
    fEtaBins{52}, fEtaRange{-5.2, 5.2},
    fPhiBins{32}, fPhiRange{-TMath::Pi(), TMath::Pi()},
    fDijetPtBins{196}, fDijetPtRange{20., 1000.},
    fDijetEtaBins{48}, fDijetEtaRange{-4.8, 4.8},
    fDijetDphiBins{16}, fDijetDphiRange{-TMath::Pi(), TMath::Pi()},
    fPtHatBins{100}, fPtHatRange{15., 1015.},
    fFracBins{100}, fFracRange{0., 1.},
    fMultBins{32}, fMultRange{-0.5, 31.5},

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
    hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi{nullptr},
    hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted{nullptr},
    hGenInclusiveJetPt{nullptr},
    hGenInclusiveJetPtEta{nullptr},
    hGenInclusiveJetPtEtaPtHat{nullptr},
    hGenPtLeadPtSublead{nullptr},
    hGenEtaLeadEtaSublead{nullptr},
    hGenEtaCMLeadEtaCMSublead{nullptr},
    hGenPtLeadPtSubleadMcReweight{nullptr},
    hGenEtaLeadEtaSubleadMcReweight{nullptr},

    hGenDijetEta{nullptr},
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

    hRecoInclusiveAllJetPtVsEta{nullptr},
    hRecoInclusiveJetPtRawVsEta{nullptr},
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi{nullptr},
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted{nullptr},

    hRecoInclusiveJetPt{nullptr},
    hRecoInclusiveJetPtEta{nullptr},
    hRecoInclusiveJetPtEtaPtHat{nullptr},
    hRecoPtLeadPtSublead{nullptr},
    hRecoEtaLeadEtaSublead{nullptr},
    hRecoEtaCMLeadEtaCMSublead{nullptr},
    hRecoPtLeadPtSubleadMcReweight{nullptr},
    hRecoEtaLeadEtaSubleadMcReweight{nullptr},

    hRecoDijetPtEta{nullptr},
    hRecoDijetPtEtaForward{nullptr},
    hRecoDijetPtEtaBackward{nullptr},
    hRecoDijetPtEtaCMForward{nullptr},
    hRecoDijetPtEtaCMBackward{nullptr},
    hRecoDijetPtEtaForwardWeighted{nullptr},
    hRecoDijetPtEtaBackwardWeighted{nullptr},
    hRecoDijetPtEtaCMForwardWeighted{nullptr},
    hRecoDijetPtEtaCMBackwardWeighted{nullptr},
    hRecoDijetPtEtaDphi{nullptr},
    hRecoDijetPtEtaDphiWeighted{nullptr},
    hRecoDijetPtEtaDphiCM{nullptr},
    hRecoDijetPtEtaDphiCMWeighted{nullptr},

    hRecoLeadJetAllPtVsEta{nullptr},
    hRecoSubLeadJetAllPtVsEta{nullptr},
    hRecoGoodInclusiveJetEtaLabFrame{nullptr},
    hRecoGoodInclusiveJetEtaCMFrame{nullptr},
    hRecoDijetEta{nullptr},
    hRecoDijetEtaCM{nullptr},

    //
    // Ref jet histograms
    //

    hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
    hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
    hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
    hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},
    hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen{nullptr},
    hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted{nullptr},

    hJESInclusiveJetPtEtaPhi{nullptr},
    hJESInclusiveJetPtEtaPhiWeighted{nullptr},
    hInclusiveJetJESVsPtGen{nullptr},
    hInclusiveJetJESGenPtGenEtaPtHatWeighted{nullptr},
    hInclusiveJetJESRecoPtRecoEtaPtHatWeighted{nullptr},

    hRecoInclusiveJetJECFactorVsPtEta{nullptr},
    hRecoInclusiveJetJEC2FactorVsPtEta{nullptr},
    hRecoInclusiveJetPtRawOverPtRefVsPtEta{nullptr},
    hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning{nullptr},

    hRecoInclusiveMatchedJetPt{nullptr},
    hRecoInclusiveMatchedJetPtVsEta{nullptr},
    hRecoInclusiveUnmatchedJetPtVsEta{nullptr},
    hRecoLeadJetMatchedPtVsEta{nullptr},
    hRecoLeadJetUnmatchedPtVsEta{nullptr},
    hRecoSubLeadJetMatchedPtVsEta{nullptr},
    hRecoSubLeadJetUnmatchedPtVsEta{nullptr},

    hRefInclusiveJetPt{nullptr},
    hRefInclusiveJetPtEta{nullptr},
    hRefInclusiveJetPtEtaPtHat{nullptr},

    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta{nullptr},
    hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},
    hRecoDijetPtEtaRefDijetPtEta{nullptr},
    hRecoDijetPtEtaRefDijetPtEtaWeighted{nullptr},
    hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted{nullptr},

    hRefSelInclusiveJetPt{nullptr},
    hRefSelInclusiveJetPtEta{nullptr},
    hRefSelInclusiveJetPtEtaPtHat{nullptr},

    hRefDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEta{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt{nullptr},
    hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted{nullptr},
    hRefDijetPtEtaDphi{nullptr},
    hRefDijetPtEtaDphiWeighted{nullptr},
    hRefDijetPtEtaForward{nullptr},
    hRefDijetPtEtaBackward{nullptr},
    hRefDijetPtEtaCMForward{nullptr},
    hRefDijetPtEtaCMBackward{nullptr},
    hRefDijetPtEtaForwardWeighted{nullptr},
    hRefDijetPtEtaBackwardWeighted{nullptr},
    hRefDijetPtEtaCMForwardWeighted{nullptr},
    hRefDijetPtEtaCMBackwardWeighted{nullptr},

    hRefDijetEtaCM{nullptr},
    hRefDijetPtEtaDphiCM{nullptr},
    hRefDijetPtEtaDphiCMWeighted{nullptr},
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

    hRefSelDijetEta{nullptr},
    hRefSelDijetPtEtaDphi{nullptr},
    hRefSelDijetPtEtaDphiWeighted{nullptr},
    hRefSelDijetEtaCM{nullptr},
    hRefSelDijetPtEtaDphiCM{nullptr},
    hRefSelDijetPtEtaDphiCMWeighted{nullptr}
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

    // New ptAve binning
    for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {

        // Gen jets
        hGenDijetEta1D[i] = nullptr;
        hGenDijetEta1DWeighted[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2D[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DWeighted[i] = nullptr;
        hGenDijetEtaForward1D[i] = nullptr;
        hGenDijetEtaForward1DWeighted[i] = nullptr;
        hGenDijetEtaBackward1D[i] = nullptr;
        hGenDijetEtaBackward1DWeighted[i] = nullptr;
        hGenDijetEta1DCM[i] = nullptr;
        hGenDijetEta1DCMWeighted[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DCM[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr;
        hGenDijetEtaCMForward1D[i] = nullptr;
        hGenDijetEtaCMForward1DWeighted[i] = nullptr;
        hGenDijetEtaCMBackward1D[i] = nullptr;
        hGenDijetEtaCMBackward1DWeighted[i] = nullptr;

        // Reco jets
        hRecoDijetEta1D[i] = nullptr;
        hRecoDijetEta1DWeighted[i] = nullptr;
        hRecoDijetEtaLeadVsEtaSubLead2D[i] = nullptr;
        hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i] = nullptr;
        hRecoDijetEtaForward1D[i] = nullptr;
        hRecoDijetEtaForward1DWeighted[i] = nullptr;
        hRecoDijetEtaBackward1D[i] = nullptr;
        hRecoDijetEtaBackward1DWeighted[i] = nullptr;
        hRecoDijetEta1DCM[i] = nullptr;
        hRecoDijetEta1DCMWeighted[i] = nullptr;
        hRecoEtaLeadVsEtaSubLead2DCM[i] = nullptr;
        hRecoEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr;
        hRecoDijetEtaCMForward1D[i] = nullptr;
        hRecoDijetEtaCMForward1DWeighted[i] = nullptr;
        hRecoDijetEtaCMBackward1D[i] = nullptr;
        hRecoDijetEtaCMBackward1DWeighted[i] = nullptr;

        // Ref jets
        hRefDijetEta1D[i] = nullptr;
        hRefDijetEta1DWeighted[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2D[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DWeighted[i] = nullptr;
        hRecoVsRefDijetEta2D[i] = nullptr;
        hRecoVsRefDijetEta2DWeighted[i] = nullptr;
        hRecoVsRefLeadJetEta2D[i] = nullptr;
        hRecoVsRefLeadJetEta2DWeighted[i] = nullptr;
        hRecoVsRefSubLeadJetEta2D[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DWeighted[i] = nullptr;
        hRefDijetEtaForward1D[i] = nullptr;
        hRefDijetEtaForward1DWeighted[i] = nullptr;
        hRefDijetEtaBackward1D[i] = nullptr;
        hRefDijetEtaBackward1DWeighted[i] = nullptr;

        hRefDijetEta1DCM[i] = nullptr;
        hRefDijetEta1DCMWeighted[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DCM[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr;
        hRecoVsRefDijetEta2DCM[i] = nullptr;
        hRecoVsRefDijetEta2DCMWeighted[i] = nullptr;
        hRecoVsRefLeadJetEta2DCM[i] = nullptr;
        hRecoVsRefLeadJetEta2DCMWeighted[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DCM[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DCMWeighted[i] = nullptr;
        hRefDijetEtaCMForward1D[i] = nullptr;
        hRefDijetEtaCMForward1DWeighted[i] = nullptr;
        hRefDijetEtaCMBackward1D[i] = nullptr;
        hRefDijetEtaCMBackward1DWeighted[i] = nullptr;

        // Ref-selected jets
        hRefSelDijetEta1D[i] = nullptr;
        hRefSelDijetEta1DWeighted[i] = nullptr;
        hRefSelRecoDijetEta1D[i] = nullptr;
        hRefSelRecoDijetEta1DWeighted[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2D[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DWeighted[i] = nullptr;
        hRefSelDijetEtaForward1D[i] = nullptr;
        hRefSelDijetEtaForward1DWeighted[i] = nullptr;
        hRefSelDijetEtaBackward1D[i] = nullptr;
        hRefSelDijetEtaBackward1DWeighted[i] = nullptr;

        hRefSelDijetEta1DCM[i] = nullptr;
        hRefSelDijetEta1DCMWeighted[i] = nullptr;
        hRefSelRecoDijetEta1DCM[i] = nullptr;
        hRefSelRecoDijetEta1DCMWeighted[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DCM[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i] = nullptr;
        hRefSelDijetEtaCMForward1D[i] = nullptr;
        hRefSelDijetEtaCMForward1DWeighted[i] = nullptr;
        hRefSelDijetEtaCMBackward1D[i] = nullptr;
        hRefSelDijetEtaCMBackward1DWeighted[i] = nullptr;

    } // for (int i{0}; i<fPtAveBins.size()-1; i++)

    // Old ptAve binning
    for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {

        // Gen jets
        hGenDijetEta1DOldPt[i] = nullptr;
        hGenDijetEta1DOldPtWeighted[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPt[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr;
        hGenDijetEtaForward1DOldPt[i] = nullptr;
        hGenDijetEtaForward1DOldPtWeighted[i] = nullptr;
        hGenDijetEtaBackward1DOldPt[i] = nullptr;
        hGenDijetEtaBackward1DOldPtWeighted[i] = nullptr;

        hGenDijetEta1DOldPtCM[i] = nullptr;
        hGenDijetEta1DOldPtCMWeighted[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr;
        hGenDijetEtaCMForward1DOldPt[i] = nullptr;
        hGenDijetEtaCMForward1DOldPtWeighted[i] = nullptr;
        hGenDijetEtaCMBackward1DOldPt[i] = nullptr;
        hGenDijetEtaCMBackward1DOldPtWeighted[i] = nullptr;
        hGenDijetEta1DOldPtBinning[i] = nullptr;
        hGenDijetEta1DOldPtBinningWeighted[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr;
        hGenDijetEtaForward1DOldPtBinning[i] = nullptr;
        hGenDijetEtaForward1DOldPtBinningWeighted[i] = nullptr;
        hGenDijetEtaBackward1DOldPtBinning[i] = nullptr;
        hGenDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr;
        hGenDijetEta1DOldPtBinningCM[i] = nullptr;
        hGenDijetEta1DOldPtBinningCMWeighted[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr;
        hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr;
        hGenDijetEtaCMForward1DOldPtBinning[i] = nullptr;
        hGenDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr;
        hGenDijetEtaCMBackward1DOldPtBinning[i] = nullptr;
        hGenDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr;

        // Reco jets
        hRecoDijetEta1DOldPt[i] = nullptr;
        hRecoDijetEta1DOldPtWeighted[i] = nullptr;
        hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i] = nullptr;
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr;
        hRecoDijetEtaForward1DOldPt[i] = nullptr;
        hRecoDijetEtaForward1DOldPtWeighted[i] = nullptr;
        hRecoDijetEtaBackward1DOldPt[i] = nullptr;
        hRecoDijetEtaBackward1DOldPtWeighted[i] = nullptr;
        hRecoDijetEta1DOldPtCM[i] = nullptr;
        hRecoDijetEta1DOldPtCMWeighted[i] = nullptr;
        hRecoEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr;
        hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr;
        hRecoDijetEtaCMForward1DOldPt[i] = nullptr;
        hRecoDijetEtaCMForward1DOldPtWeighted[i] = nullptr;
        hRecoDijetEtaCMBackward1DOldPt[i] = nullptr;
        hRecoDijetEtaCMBackward1DOldPtWeighted[i] = nullptr;
        hRecoDijetEta1DOldPtBinning[i] = nullptr;
        hRecoDijetEta1DOldPtBinningWeighted[i] = nullptr;
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr;
        hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr;
        hRecoDijetEtaForward1DOldPtBinning[i] = nullptr;
        hRecoDijetEtaForward1DOldPtBinningWeighted[i] = nullptr;
        hRecoDijetEtaBackward1DOldPtBinning[i] = nullptr;
        hRecoDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr;
        hRecoDijetEta1DOldPtBinningCM[i] = nullptr;
        hRecoDijetEta1DOldPtBinningCMWeighted[i] = nullptr;
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr;
        hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr;
        hRecoDijetEtaCMForward1DOldPtBinning[i] = nullptr;
        hRecoDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr;
        hRecoDijetEtaCMBackward1DOldPtBinning[i] = nullptr;
        hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr;

        // Ref jets
        hRefDijetEta1DOldPt[i] = nullptr;
        hRefDijetEta1DOldPtWeighted[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPt[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr;
        hRecoVsRefDijetEta2DOldPt[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtWeighted[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPt[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtWeighted[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPt[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtWeighted[i] = nullptr;
        hRefDijetEtaForward1DOldPt[i] = nullptr;
        hRefDijetEtaForward1DOldPtWeighted[i] = nullptr;
        hRefDijetEtaBackward1DOldPt[i] = nullptr;
        hRefDijetEtaBackward1DOldPtWeighted[i] = nullptr;
        hRefDijetEta1DOldPtCM[i] = nullptr;
        hRefDijetEta1DOldPtCMWeighted[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtCM[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtCMWeighted[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtCM[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtCMWeighted[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtCM[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i] = nullptr;
        hRefDijetEtaCMForward1DOldPt[i] = nullptr;
        hRefDijetEtaCMForward1DOldPtWeighted[i] = nullptr;
        hRefDijetEtaCMBackward1DOldPt[i] = nullptr;
        hRefDijetEtaCMBackward1DOldPtWeighted[i] = nullptr;
        hRefDijetEta1DOldPtBinning[i] = nullptr;
        hRefDijetEta1DOldPtBinningWeighted[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtBinning[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtBinningWeighted[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtBinning[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtBinning[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i] = nullptr;
        hRefDijetEtaForward1DOldPtBinning[i] = nullptr;
        hRefDijetEtaForward1DOldPtBinningWeighted[i] = nullptr;
        hRefDijetEtaBackward1DOldPtBinning[i] = nullptr;
        hRefDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr;
        hRefDijetEta1DOldPtBinningCM[i] = nullptr;
        hRefDijetEta1DOldPtBinningCMWeighted[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr;
        hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtBinningCM[i] = nullptr;
        hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtBinningCM[i] = nullptr;
        hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i] = nullptr;
        hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i] = nullptr;
        hRefDijetEtaCMForward1DOldPtBinning[i] = nullptr;
        hRefDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr;
        hRefDijetEtaCMBackward1DOldPtBinning[i] = nullptr;
        hRefDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr;

        // Ref-selected jets
        hRefSelDijetEta1DOldPt[i] = nullptr;
        hRefSelDijetEta1DOldPtWeighted[i] = nullptr;
        hRefSelRecoDijetEta1DOldPt[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtWeighted[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPt[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i] = nullptr;
        hRefSelDijetEtaForward1DOldPt[i] = nullptr;
        hRefSelDijetEtaForward1DOldPtWeighted[i] = nullptr;
        hRefSelDijetEtaBackward1DOldPt[i] = nullptr;
        hRefSelDijetEtaBackward1DOldPtWeighted[i] = nullptr;
        hRefSelDijetEta1DOldPtCM[i] = nullptr;
        hRefSelDijetEta1DOldPtCMWeighted[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtCM[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtCMWeighted[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = nullptr;
        hRefSelDijetEtaCMForward1DOldPt[i] = nullptr;
        hRefSelDijetEtaCMForward1DOldPtWeighted[i] = nullptr;
        hRefSelDijetEtaCMBackward1DOldPt[i] = nullptr;
        hRefSelDijetEtaCMBackward1DOldPtWeighted[i] = nullptr;
        hRefSelDijetEta1DOldPtBinning[i] = nullptr;
        hRefSelDijetEta1DOldPtBinningWeighted[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtBinning[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtBinningWeighted[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = nullptr;
        hRefSelDijetEtaForward1DOldPtBinning[i] = nullptr;
        hRefSelDijetEtaForward1DOldPtBinningWeighted[i] = nullptr;
        hRefSelDijetEtaBackward1DOldPtBinning[i] = nullptr;
        hRefSelDijetEtaBackward1DOldPtBinningWeighted[i] = nullptr;
        hRefSelDijetEta1DOldPtBinningCM[i] = nullptr;
        hRefSelDijetEta1DOldPtBinningCMWeighted[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtBinningCM[i] = nullptr;
        hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = nullptr;
        hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = nullptr;
        hRefSelDijetEtaCMForward1DOldPtBinning[i] = nullptr;
        hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i] = nullptr;
        hRefSelDijetEtaCMBackward1DOldPtBinning[i] = nullptr;
        hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i] = nullptr;
    } // for (int i{0}; i<fPtAveOldBins.size()-2; i++)
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
        if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
        if (hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
        if (hGenInclusiveJetPt) delete hGenInclusiveJetPt;
        if (hGenInclusiveJetPtEta) delete hGenInclusiveJetPtEta;
        if (hGenInclusiveJetPtEtaPtHat) delete hGenInclusiveJetPtEtaPtHat;
        if (hGenPtLeadPtSublead) delete hGenPtLeadPtSublead;
        if (hGenEtaLeadEtaSublead) delete hGenEtaLeadEtaSublead;
        if (hGenEtaCMLeadEtaCMSublead) delete hGenEtaCMLeadEtaCMSublead;
        if (hGenPtLeadPtSubleadMcReweight) delete hGenPtLeadPtSubleadMcReweight;
        if (hGenEtaLeadEtaSubleadMcReweight) delete hGenEtaLeadEtaSubleadMcReweight;
        if (hGenDijetEta) delete hGenDijetEta;
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

        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            if (hGenDijetEta1D[i]) delete hGenDijetEta1D[i];
            if (hGenDijetEta1DWeighted[i]) delete hGenDijetEta1DWeighted[i];
            if (hGenDijetEtaLeadVsEtaSubLead2D[i]) delete hGenDijetEtaLeadVsEtaSubLead2D[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DWeighted[i]) delete hGenDijetEtaLeadVsEtaSubLead2DWeighted[i];
            if (hGenDijetEtaForward1D[i]) delete hGenDijetEtaForward1D[i];
            if (hGenDijetEtaForward1DWeighted[i]) delete hGenDijetEtaForward1DWeighted[i];
            if (hGenDijetEtaBackward1D[i]) delete hGenDijetEtaBackward1D[i];
            if (hGenDijetEtaBackward1DWeighted[i]) delete hGenDijetEtaBackward1DWeighted[i];
            if (hGenDijetEta1DCM[i]) delete hGenDijetEta1DCM[i];
            if (hGenDijetEta1DCMWeighted[i]) delete hGenDijetEta1DCMWeighted[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DCM[i]) delete hGenDijetEtaLeadVsEtaSubLead2DCM[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i]) delete hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i];
            if (hGenDijetEtaCMForward1D[i]) delete hGenDijetEtaCMForward1D[i];
            if (hGenDijetEtaCMForward1DWeighted[i]) delete hGenDijetEtaCMForward1DWeighted[i];
            if (hGenDijetEtaCMBackward1D[i]) delete hGenDijetEtaCMBackward1D[i];
            if (hGenDijetEtaCMBackward1DWeighted[i]) delete hGenDijetEtaCMBackward1DWeighted[i];
        }

        for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {
            if (hGenDijetEta1DOldPt[i]) delete hGenDijetEta1DOldPt[i];
            if (hGenDijetEta1DOldPtWeighted[i]) delete hGenDijetEta1DOldPtWeighted[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPt[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPt[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i];
            if (hGenDijetEtaForward1DOldPt[i]) delete hGenDijetEtaForward1DOldPt[i];
            if (hGenDijetEtaForward1DOldPtWeighted[i]) delete hGenDijetEtaForward1DOldPtWeighted[i];
            if (hGenDijetEtaBackward1DOldPt[i]) delete hGenDijetEtaBackward1DOldPt[i];
            if (hGenDijetEtaBackward1DOldPtWeighted[i]) delete hGenDijetEtaBackward1DOldPtWeighted[i];
            if (hGenDijetEta1DOldPtCM[i]) delete hGenDijetEta1DOldPtCM[i];
            if (hGenDijetEta1DOldPtCMWeighted[i]) delete hGenDijetEta1DOldPtCMWeighted[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[i];
            if (hGenDijetEtaCMForward1DOldPt[i]) delete hGenDijetEtaCMForward1DOldPt[i];
            if (hGenDijetEtaCMForward1DOldPtWeighted[i]) delete hGenDijetEtaCMForward1DOldPtWeighted[i];
            if (hGenDijetEtaCMBackward1DOldPt[i]) delete hGenDijetEtaCMBackward1DOldPt[i];
            if (hGenDijetEtaCMBackward1DOldPtWeighted[i]) delete hGenDijetEtaCMBackward1DOldPtWeighted[i];
            if (hGenDijetEta1DOldPtBinning[i]) delete hGenDijetEta1DOldPtBinning[i];
            if (hGenDijetEta1DOldPtBinningWeighted[i]) delete hGenDijetEta1DOldPtBinningWeighted[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i];
            if (hGenDijetEtaForward1DOldPtBinning[i]) delete hGenDijetEtaForward1DOldPtBinning[i];
            if (hGenDijetEtaForward1DOldPtBinningWeighted[i]) delete hGenDijetEtaForward1DOldPtBinningWeighted[i];
            if (hGenDijetEtaBackward1DOldPtBinning[i]) delete hGenDijetEtaBackward1DOldPtBinning[i];
            if (hGenDijetEtaBackward1DOldPtBinningWeighted[i]) delete hGenDijetEtaBackward1DOldPtBinningWeighted[i];
            if (hGenDijetEta1DOldPtBinningCM[i]) delete hGenDijetEta1DOldPtBinningCM[i];
            if (hGenDijetEta1DOldPtBinningCMWeighted[i]) delete hGenDijetEta1DOldPtBinningCMWeighted[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[i];
            if (hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) delete hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i];
            if (hGenDijetEtaCMForward1DOldPtBinning[i]) delete hGenDijetEtaCMForward1DOldPtBinning[i];
            if (hGenDijetEtaCMForward1DOldPtBinningWeighted[i]) delete hGenDijetEtaCMForward1DOldPtBinningWeighted[i];
            if (hGenDijetEtaCMBackward1DOldPtBinning[i]) delete hGenDijetEtaCMBackward1DOldPtBinning[i];
            if (hGenDijetEtaCMBackward1DOldPtBinningWeighted[i]) delete hGenDijetEtaCMBackward1DOldPtBinningWeighted[i];
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

    if (hRecoInclusiveAllJetPtVsEta) delete hRecoInclusiveAllJetPtVsEta;
    if (hRecoInclusiveJetPtRawVsEta) delete hRecoInclusiveJetPtRawVsEta;
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    if (hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted) delete hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    if (hRecoInclusiveJetPt) delete hRecoInclusiveJetPt;
    if (hRecoInclusiveJetPtEta) delete hRecoInclusiveJetPtEta;
    if (hRecoInclusiveJetPtEtaPtHat) delete hRecoInclusiveJetPtEtaPtHat;
    if (hRecoPtLeadPtSublead) delete hRecoPtLeadPtSublead;
    if (hRecoEtaLeadEtaSublead) delete hRecoEtaLeadEtaSublead;
    if (hRecoEtaCMLeadEtaCMSublead) delete hRecoEtaCMLeadEtaCMSublead;
    if (hRecoPtLeadPtSubleadMcReweight) delete hRecoPtLeadPtSubleadMcReweight;
    if (hRecoEtaLeadEtaSubleadMcReweight) delete hRecoEtaLeadEtaSubleadMcReweight;
    if (hRecoDijetPtEta) delete hRecoDijetPtEta;
    if (hRecoDijetPtEtaForward) delete hRecoDijetPtEtaForward;
    if (hRecoDijetPtEtaBackward) delete hRecoDijetPtEtaBackward;
    if (hRecoDijetPtEtaCMForward) delete hRecoDijetPtEtaCMForward;
    if (hRecoDijetPtEtaCMBackward) delete hRecoDijetPtEtaCMBackward;
    if (hRecoDijetPtEtaForwardWeighted) delete hRecoDijetPtEtaForwardWeighted;
    if (hRecoDijetPtEtaBackwardWeighted) delete hRecoDijetPtEtaBackwardWeighted;
    if (hRecoDijetPtEtaCMForwardWeighted) delete hRecoDijetPtEtaCMForwardWeighted;
    if (hRecoDijetPtEtaCMBackwardWeighted) delete hRecoDijetPtEtaCMBackwardWeighted;
    if (hRecoDijetPtEtaDphi) delete hRecoDijetPtEtaDphi;
    if (hRecoDijetPtEtaDphiWeighted) delete hRecoDijetPtEtaDphiWeighted;
    if (hRecoDijetPtEtaDphiCM) delete hRecoDijetPtEtaDphiCM;
    if (hRecoDijetPtEtaDphiCMWeighted) delete hRecoDijetPtEtaDphiCMWeighted;
    if (hRecoLeadJetAllPtVsEta) delete hRecoLeadJetAllPtVsEta;
    if (hRecoSubLeadJetAllPtVsEta) delete hRecoSubLeadJetAllPtVsEta;
    if (hRecoGoodInclusiveJetEtaLabFrame) delete hRecoGoodInclusiveJetEtaLabFrame;
    if (hRecoGoodInclusiveJetEtaCMFrame) delete hRecoGoodInclusiveJetEtaCMFrame;
    if (hRecoDijetEta) delete hRecoDijetEta;
    if (hRecoDijetEtaCM) delete hRecoDijetEtaCM;


    for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
        if (hRecoDijetEta1D[i]) delete hRecoDijetEta1D[i];
        if (hRecoDijetEta1DWeighted[i]) delete hRecoDijetEta1DWeighted[i];
        if (hRecoDijetEtaLeadVsEtaSubLead2D[i]) delete hRecoDijetEtaLeadVsEtaSubLead2D[i];
        if (hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i]) delete hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i];
        if (hRecoDijetEtaForward1D[i]) delete hRecoDijetEtaForward1D[i];
        if (hRecoDijetEtaForward1DWeighted[i]) delete hRecoDijetEtaForward1DWeighted[i];
        if (hRecoDijetEtaBackward1D[i]) delete hRecoDijetEtaBackward1D[i];
        if (hRecoDijetEtaBackward1DWeighted[i]) delete hRecoDijetEtaBackward1DWeighted[i];
        if (hRecoDijetEta1DCM[i]) delete hRecoDijetEta1DCM[i];
        if (hRecoDijetEta1DCMWeighted[i]) delete hRecoDijetEta1DCMWeighted[i];
        if (hRecoEtaLeadVsEtaSubLead2DCM[i]) delete hRecoEtaLeadVsEtaSubLead2DCM[i];
        if (hRecoEtaLeadVsEtaSubLead2DCMWeighted[i]) delete hRecoEtaLeadVsEtaSubLead2DCMWeighted[i];
        if (hRecoDijetEtaCMForward1D[i]) delete hRecoDijetEtaCMForward1D[i];
        if (hRecoDijetEtaCMForward1DWeighted[i]) delete hRecoDijetEtaCMForward1DWeighted[i];
        if (hRecoDijetEtaCMBackward1D[i]) delete hRecoDijetEtaCMBackward1D[i];
        if (hRecoDijetEtaCMBackward1DWeighted[i]) delete hRecoDijetEtaCMBackward1DWeighted[i];
    }

    for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {
        if (hRecoDijetEta1DOldPt[i]) delete hRecoDijetEta1DOldPt[i];
        if (hRecoDijetEta1DOldPtWeighted[i]) delete hRecoDijetEta1DOldPtWeighted[i];
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i]) delete hRecoDijetEtaLeadVsEtaSubLead2DOldPt[i];
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i]) delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[i];
        if (hRecoDijetEtaForward1DOldPt[i]) delete hRecoDijetEtaForward1DOldPt[i];
        if (hRecoDijetEtaForward1DOldPtWeighted[i]) delete hRecoDijetEtaForward1DOldPtWeighted[i];
        if (hRecoDijetEtaBackward1DOldPt[i]) delete hRecoDijetEtaBackward1DOldPt[i];
        if (hRecoDijetEtaBackward1DOldPtWeighted[i]) delete hRecoDijetEtaBackward1DOldPtWeighted[i];

        if (hRecoDijetEta1DOldPtCM[i]) delete hRecoDijetEta1DOldPtCM[i];
        if (hRecoDijetEta1DOldPtCMWeighted[i]) delete hRecoDijetEta1DOldPtCMWeighted[i];
        if (hRecoEtaLeadVsEtaSubLead2DOldPtCM[i]) delete hRecoEtaLeadVsEtaSubLead2DOldPtCM[i];
        if (hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) delete hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[i];
        if (hRecoDijetEtaCMForward1DOldPt[i]) delete hRecoDijetEtaCMForward1DOldPt[i];
        if (hRecoDijetEtaCMForward1DOldPtWeighted[i]) delete hRecoDijetEtaCMForward1DOldPtWeighted[i];
        if (hRecoDijetEtaCMBackward1DOldPt[i]) delete hRecoDijetEtaCMBackward1DOldPt[i];
        if (hRecoDijetEtaCMBackward1DOldPtWeighted[i]) delete hRecoDijetEtaCMBackward1DOldPtWeighted[i];

        if (hRecoDijetEta1DOldPtBinning[i]) delete hRecoDijetEta1DOldPtBinning[i];
        if (hRecoDijetEta1DOldPtBinningWeighted[i]) delete hRecoDijetEta1DOldPtBinningWeighted[i];
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i]) delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[i];
        if (hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) delete hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i];
        if (hRecoDijetEtaForward1DOldPtBinning[i]) delete hRecoDijetEtaForward1DOldPtBinning[i];
        if (hRecoDijetEtaForward1DOldPtBinningWeighted[i]) delete hRecoDijetEtaForward1DOldPtBinningWeighted[i];
        if (hRecoDijetEtaBackward1DOldPtBinning[i]) delete hRecoDijetEtaBackward1DOldPtBinning[i];
        if (hRecoDijetEtaBackward1DOldPtBinningWeighted[i]) delete hRecoDijetEtaBackward1DOldPtBinningWeighted[i];

        if (hRecoDijetEta1DOldPtBinningCM[i]) delete hRecoDijetEta1DOldPtBinningCM[i];
        if (hRecoDijetEta1DOldPtBinningCMWeighted[i]) delete hRecoDijetEta1DOldPtBinningCMWeighted[i];
        if (hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) delete hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[i];
        if (hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) delete hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i];
        if (hRecoDijetEtaCMForward1DOldPtBinning[i]) delete hRecoDijetEtaCMForward1DOldPtBinning[i];
        if (hRecoDijetEtaCMForward1DOldPtBinningWeighted[i]) delete hRecoDijetEtaCMForward1DOldPtBinningWeighted[i];
        if (hRecoDijetEtaCMBackward1DOldPtBinning[i]) delete hRecoDijetEtaCMBackward1DOldPtBinning[i];
        if (hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i]) delete hRecoDijetEtaCMBackward1DOldPtBinningWeighted[i];
    }

    if ( fIsMc ) {
        // Ref jet histograms
        if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
        if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
        if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
        if (hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted) delete hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
        if (hJESInclusiveJetPtEtaPhi) delete hJESInclusiveJetPtEtaPhi;
        if (hJESInclusiveJetPtEtaPhiWeighted) delete hJESInclusiveJetPtEtaPhiWeighted;

        if (hRecoInclusiveJetJECFactorVsPtEta) delete hRecoInclusiveJetJECFactorVsPtEta;
        if (hRecoInclusiveJetJEC2FactorVsPtEta) delete hRecoInclusiveJetJEC2FactorVsPtEta;
        if (hRecoInclusiveJetPtRawOverPtRefVsPtEta) delete hRecoInclusiveJetPtRawOverPtRefVsPtEta;
        if (hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning) delete hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning;

        if ( hInclusiveJetJESVsPtGen ) delete hInclusiveJetJESVsPtGen;
        if ( hInclusiveJetJESGenPtGenEtaPtHatWeighted ) delete hInclusiveJetJESGenPtGenEtaPtHatWeighted;
        if ( hInclusiveJetJESRecoPtRecoEtaPtHatWeighted ) delete hInclusiveJetJESRecoPtRecoEtaPtHatWeighted;

        if (hRecoInclusiveMatchedJetPt) delete hRecoInclusiveMatchedJetPt;
        if (hRecoInclusiveMatchedJetPtVsEta) delete hRecoInclusiveMatchedJetPtVsEta;
        if (hRecoInclusiveUnmatchedJetPtVsEta) delete hRecoInclusiveUnmatchedJetPtVsEta;
        if (hRecoLeadJetMatchedPtVsEta) delete hRecoLeadJetMatchedPtVsEta;
        if (hRecoLeadJetUnmatchedPtVsEta) delete hRecoLeadJetUnmatchedPtVsEta;
        if (hRecoSubLeadJetMatchedPtVsEta) delete hRecoSubLeadJetMatchedPtVsEta;
        if (hRecoSubLeadJetUnmatchedPtVsEta) delete hRecoSubLeadJetUnmatchedPtVsEta;

        if (hRefInclusiveJetPt) delete hRefInclusiveJetPt;
        if (hRefInclusiveJetPtEta) delete hRefInclusiveJetPtEta;
        if (hRefInclusiveJetPtEtaPtHat) delete hRefInclusiveJetPtEtaPtHat;

        if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta;
        if (hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;
        if (hRecoDijetPtEtaRefDijetPtEta) delete hRecoDijetPtEtaRefDijetPtEta;
        if (hRecoDijetPtEtaRefDijetPtEtaWeighted) delete hRecoDijetPtEtaRefDijetPtEtaWeighted;
        if (hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted) delete hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;

        if (hRefDijetEta) delete hRefDijetEta;
        if (hRefDijetEtaVsRecoDijetEta) delete hRefDijetEtaVsRecoDijetEta;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted;
        if (hRefDijetPtEtaDphi) delete hRefDijetPtEtaDphi;
        if (hRefDijetPtEtaDphiWeighted) delete hRefDijetPtEtaDphiWeighted;

        if (hRefDijetPtEtaForward) delete hRefDijetPtEtaForward;
        if (hRefDijetPtEtaBackward) delete hRefDijetPtEtaBackward;
        if (hRefDijetPtEtaCMForward) delete hRefDijetPtEtaCMForward;
        if (hRefDijetPtEtaCMBackward) delete hRefDijetPtEtaCMBackward;
        if (hRefDijetPtEtaForwardWeighted) delete hRefDijetPtEtaForwardWeighted;
        if (hRefDijetPtEtaBackwardWeighted) delete hRefDijetPtEtaBackwardWeighted;
        if (hRefDijetPtEtaCMForwardWeighted) delete hRefDijetPtEtaCMForwardWeighted;
        if (hRefDijetPtEtaCMBackwardWeighted) delete hRefDijetPtEtaCMBackwardWeighted;

        if (hRefDijetEtaCM) delete hRefDijetEtaCM;
        if (hRefDijetPtEtaDphiCM) delete hRefDijetPtEtaDphiCM;
        if (hRefDijetPtEtaDphiCMWeighted) delete hRefDijetPtEtaDphiCMWeighted;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM;
        if (hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted) delete hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted;

        if (hRefPtLeadPtSublead) delete hRefPtLeadPtSublead;
        if (hRefEtaLeadEtaSublead) delete hRefEtaLeadEtaSublead;
        if (hRefEtaCMLeadEtaCMSublead) delete hRefEtaCMLeadEtaCMSublead;
        if (hRefPtLeadPtSubleadMcReweight) delete hRefPtLeadPtSubleadMcReweight;
        if (hRefEtaLeadEtaSubleadMcReweight) delete hRefEtaLeadEtaSubleadMcReweight;

        for (unsigned int i = 0; i < fPtAveBins.size()-1; ++i) {
            if (hRefDijetEta1D[i]) delete hRefDijetEta1D[i];
            if (hRefDijetEta1DWeighted[i]) delete hRefDijetEta1DWeighted[i];
            if (hRefEtaLeadVsEtaSubLead2D[i]) delete hRefEtaLeadVsEtaSubLead2D[i];
            if (hRefEtaLeadVsEtaSubLead2DWeighted[i]) delete hRefEtaLeadVsEtaSubLead2DWeighted[i];
            if (hRecoVsRefDijetEta2D[i]) delete hRecoVsRefDijetEta2D[i];
            if (hRecoVsRefDijetEta2DWeighted[i]) delete hRecoVsRefDijetEta2DWeighted[i];
            if (hRecoVsRefLeadJetEta2D[i]) delete hRecoVsRefLeadJetEta2D[i];
            if (hRecoVsRefLeadJetEta2DWeighted[i]) delete hRecoVsRefLeadJetEta2DWeighted[i];
            if (hRecoVsRefSubLeadJetEta2D[i]) delete hRecoVsRefSubLeadJetEta2D[i];
            if (hRecoVsRefSubLeadJetEta2DWeighted[i]) delete hRecoVsRefSubLeadJetEta2DWeighted[i];
            if (hRefDijetEtaForward1D[i]) delete hRefDijetEtaForward1D[i];
            if (hRefDijetEtaForward1DWeighted[i]) delete hRefDijetEtaForward1DWeighted[i];
            if (hRefDijetEtaBackward1D[i]) delete hRefDijetEtaBackward1D[i];
            if (hRefDijetEtaBackward1DWeighted[i]) delete hRefDijetEtaBackward1DWeighted[i];
            
            if (hRefDijetEta1DCM[i]) delete hRefDijetEta1DCM[i];
            if (hRefDijetEta1DCMWeighted[i]) delete hRefDijetEta1DCMWeighted[i];
            if (hRefEtaLeadVsEtaSubLead2DCM[i]) delete hRefEtaLeadVsEtaSubLead2DCM[i];
            if (hRefEtaLeadVsEtaSubLead2DCMWeighted[i]) delete hRefEtaLeadVsEtaSubLead2DCMWeighted[i];
            if (hRecoVsRefDijetEta2DCM[i]) delete hRecoVsRefDijetEta2DCM[i];
            if (hRecoVsRefDijetEta2DCMWeighted[i]) delete hRecoVsRefDijetEta2DCMWeighted[i];
            if (hRecoVsRefLeadJetEta2DCM[i]) delete hRecoVsRefLeadJetEta2DCM[i];
            if (hRecoVsRefLeadJetEta2DCMWeighted[i]) delete hRecoVsRefLeadJetEta2DCMWeighted[i];
            if (hRecoVsRefSubLeadJetEta2DCM[i]) delete hRecoVsRefSubLeadJetEta2DCM[i];
            if (hRecoVsRefSubLeadJetEta2DCMWeighted[i]) delete hRecoVsRefSubLeadJetEta2DCMWeighted[i];
            if (hRefDijetEtaCMForward1D[i]) delete hRefDijetEtaCMForward1D[i];
            if (hRefDijetEtaCMForward1DWeighted[i]) delete hRefDijetEtaCMForward1DWeighted[i];
            if (hRefDijetEtaCMBackward1D[i]) delete hRefDijetEtaCMBackward1D[i];
            if (hRefDijetEtaCMBackward1DWeighted[i]) delete hRefDijetEtaCMBackward1DWeighted[i];
        } // for (int i = 0; i < 16; ++i)

        for (unsigned int i = 0; i < fPtAveOldBins.size()-1; ++i) {
            if (hRefDijetEta1DOldPt[i]) delete hRefDijetEta1DOldPt[i];
            if (hRefDijetEta1DOldPtWeighted[i]) delete hRefDijetEta1DOldPtWeighted[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPt[i]) delete hRefEtaLeadVsEtaSubLead2DOldPt[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i];
            if (hRecoVsRefDijetEta2DOldPt[i]) delete hRecoVsRefDijetEta2DOldPt[i];
            if (hRecoVsRefDijetEta2DOldPtWeighted[i]) delete hRecoVsRefDijetEta2DOldPtWeighted[i];
            if (hRecoVsRefLeadJetEta2DOldPt[i]) delete hRecoVsRefLeadJetEta2DOldPt[i];
            if (hRecoVsRefLeadJetEta2DOldPtWeighted[i]) delete hRecoVsRefLeadJetEta2DOldPtWeighted[i];
            if (hRecoVsRefSubLeadJetEta2DOldPt[i]) delete hRecoVsRefSubLeadJetEta2DOldPt[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtWeighted[i]) delete hRecoVsRefSubLeadJetEta2DOldPtWeighted[i];
            if (hRefDijetEtaForward1DOldPt[i]) delete hRefDijetEtaForward1DOldPt[i];
            if (hRefDijetEtaForward1DOldPtWeighted[i]) delete hRefDijetEtaForward1DOldPtWeighted[i];
            if (hRefDijetEtaBackward1DOldPt[i]) delete hRefDijetEtaBackward1DOldPt[i];
            if (hRefDijetEtaBackward1DOldPtWeighted[i]) delete hRefDijetEtaBackward1DOldPtWeighted[i];

            if (hRefDijetEta1DOldPtCM[i]) delete hRefDijetEta1DOldPtCM[i];
            if (hRefDijetEta1DOldPtCMWeighted[i]) delete hRefDijetEta1DOldPtCMWeighted[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtCM[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtCM[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i];
            if (hRecoVsRefDijetEta2DOldPtCM[i]) delete hRecoVsRefDijetEta2DOldPtCM[i];
            if (hRecoVsRefDijetEta2DOldPtCMWeighted[i]) delete hRecoVsRefDijetEta2DOldPtCMWeighted[i];
            if (hRecoVsRefLeadJetEta2DOldPtCM[i]) delete hRecoVsRefLeadJetEta2DOldPtCM[i];
            if (hRecoVsRefLeadJetEta2DOldPtCMWeighted[i]) delete hRecoVsRefLeadJetEta2DOldPtCMWeighted[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtCM[i]) delete hRecoVsRefSubLeadJetEta2DOldPtCM[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i]) delete hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i];
            if (hRefDijetEtaCMForward1DOldPt[i]) delete hRefDijetEtaCMForward1DOldPt[i];
            if (hRefDijetEtaCMForward1DOldPtWeighted[i]) delete hRefDijetEtaCMForward1DOldPtWeighted[i];
            if (hRefDijetEtaCMBackward1DOldPt[i]) delete hRefDijetEtaCMBackward1DOldPt[i];
            if (hRefDijetEtaCMBackward1DOldPtWeighted[i]) delete hRefDijetEtaCMBackward1DOldPtWeighted[i];

            if (hRefDijetEta1DOldPtBinning[i]) delete hRefDijetEta1DOldPtBinning[i];
            if (hRefDijetEta1DOldPtBinningWeighted[i]) delete hRefDijetEta1DOldPtBinningWeighted[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtBinning[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i];
            if (hRecoVsRefDijetEta2DOldPtBinning[i]) delete hRecoVsRefDijetEta2DOldPtBinning[i];
            if (hRecoVsRefDijetEta2DOldPtBinningWeighted[i]) delete hRecoVsRefDijetEta2DOldPtBinningWeighted[i];
            if (hRecoVsRefLeadJetEta2DOldPtBinning[i]) delete hRecoVsRefLeadJetEta2DOldPtBinning[i];
            if (hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i]) delete hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtBinning[i]) delete hRecoVsRefSubLeadJetEta2DOldPtBinning[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i]) delete hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i];
            if (hRefDijetEtaForward1DOldPtBinning[i]) delete hRefDijetEtaForward1DOldPtBinning[i];
            if (hRefDijetEtaForward1DOldPtBinningWeighted[i]) delete hRefDijetEtaForward1DOldPtBinningWeighted[i];
            if (hRefDijetEtaBackward1DOldPtBinning[i]) delete hRefDijetEtaBackward1DOldPtBinning[i];
            if (hRefDijetEtaBackward1DOldPtBinningWeighted[i]) delete hRefDijetEtaBackward1DOldPtBinningWeighted[i];

            if (hRefDijetEta1DOldPtBinningCM[i]) delete hRefDijetEta1DOldPtBinningCM[i];
            if (hRefDijetEta1DOldPtBinningCMWeighted[i]) delete hRefDijetEta1DOldPtBinningCMWeighted[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i];
            if (hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) delete hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i];
            if (hRecoVsRefDijetEta2DOldPtBinningCM[i]) delete hRecoVsRefDijetEta2DOldPtBinningCM[i];
            if (hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i]) delete hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i];
            if (hRecoVsRefLeadJetEta2DOldPtBinningCM[i]) delete hRecoVsRefLeadJetEta2DOldPtBinningCM[i];
            if (hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i]) delete hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]) delete hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i];
            if (hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i]) delete hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i];
            if (hRefDijetEtaCMForward1DOldPtBinning[i]) delete hRefDijetEtaCMForward1DOldPtBinning[i];
            if (hRefDijetEtaCMForward1DOldPtBinningWeighted[i]) delete hRefDijetEtaCMForward1DOldPtBinningWeighted[i];
            if (hRefDijetEtaCMBackward1DOldPtBinning[i]) delete hRefDijetEtaCMBackward1DOldPtBinning[i];
            if (hRefDijetEtaCMBackward1DOldPtBinningWeighted[i]) delete hRefDijetEtaCMBackward1DOldPtBinningWeighted[i];
        } // for (int i = 0; i < 5; ++i)

        // Ref-selected jet histograms

        if (hRefSelInclusiveJetPt) delete hRefSelInclusiveJetPt;
        if (hRefSelInclusiveJetPtEta) delete hRefSelInclusiveJetPt;
        if (hRefSelInclusiveJetPtEtaPtHat) delete hRefSelInclusiveJetPt;

        if (hRefSelDijetEta) delete hRefSelDijetEta;
        if (hRefSelDijetPtEtaDphi) delete hRefSelDijetPtEtaDphi;
        if (hRefSelDijetPtEtaDphiWeighted) delete hRefSelDijetPtEtaDphiWeighted;
        if (hRefSelDijetEtaCM) delete hRefSelDijetEtaCM;
        if (hRefSelDijetPtEtaDphiCM) delete hRefSelDijetPtEtaDphiCM;
        if (hRefSelDijetPtEtaDphiCMWeighted) delete hRefSelDijetPtEtaDphiCMWeighted;

        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            if (hRefSelDijetEta1D[i]) delete hRefSelDijetEta1D[i];
            if (hRefSelDijetEta1DWeighted[i]) delete hRefSelDijetEta1DWeighted[i];
            if (hRefSelRecoDijetEta1D[i]) delete hRefSelRecoDijetEta1D[i];
            if (hRefSelRecoDijetEta1DWeighted[i]) delete hRefSelRecoDijetEta1DWeighted[i];
            if (hRefSelEtaLeadVsEtaSubLead2D[i]) delete hRefSelEtaLeadVsEtaSubLead2D[i];
            if (hRefSelEtaLeadVsEtaSubLead2DWeighted[i]) delete hRefSelEtaLeadVsEtaSubLead2DWeighted[i];
            if (hRefSelDijetEtaForward1D[i]) delete hRefSelDijetEtaForward1D[i];
            if (hRefSelDijetEtaForward1DWeighted[i]) delete hRefSelDijetEtaForward1DWeighted[i];
            if (hRefSelDijetEtaBackward1D[i]) delete hRefSelDijetEtaBackward1D[i];
            if (hRefSelDijetEtaBackward1DWeighted[i]) delete hRefSelDijetEtaBackward1DWeighted[i];

            if (hRefSelDijetEta1DCM[i]) delete hRefSelDijetEta1DCM[i];
            if (hRefSelDijetEta1DCMWeighted[i]) delete hRefSelDijetEta1DCMWeighted[i];
            if (hRefSelRecoDijetEta1DCM[i]) delete hRefSelRecoDijetEta1DCM[i];
            if (hRefSelRecoDijetEta1DCMWeighted[i]) delete hRefSelRecoDijetEta1DCMWeighted[i];
            if (hRefSelEtaLeadVsEtaSubLead2DCM[i]) delete hRefSelEtaLeadVsEtaSubLead2DCM[i];
            if (hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i]) delete hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i];
            if (hRefSelDijetEtaCMForward1D[i]) delete hRefSelDijetEtaCMForward1D[i];
            if (hRefSelDijetEtaCMForward1DWeighted[i]) delete hRefSelDijetEtaCMForward1DWeighted[i];
            if (hRefSelDijetEtaCMBackward1D[i]) delete hRefSelDijetEtaCMBackward1D[i];
            if (hRefSelDijetEtaCMBackward1DWeighted[i]) delete hRefSelDijetEtaCMBackward1DWeighted[i];
        }

        for (unsigned int i{0}; i<fPtAveOldBins.size()-1; i++) {
            if (hRefSelDijetEta1DOldPt[i]) delete hRefSelDijetEta1DOldPt[i];
            if (hRefSelDijetEta1DOldPtWeighted[i]) delete hRefSelDijetEta1DOldPtWeighted[i];
            if (hRefSelRecoDijetEta1DOldPt[i]) delete hRefSelRecoDijetEta1DOldPt[i];
            if (hRefSelRecoDijetEta1DOldPtWeighted[i]) delete hRefSelRecoDijetEta1DOldPtWeighted[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPt[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPt[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i];
            if (hRefSelDijetEtaForward1DOldPt[i]) delete hRefSelDijetEtaForward1DOldPt[i];
            if (hRefSelDijetEtaForward1DOldPtWeighted[i]) delete hRefSelDijetEtaForward1DOldPtWeighted[i];
            if (hRefSelDijetEtaBackward1DOldPt[i]) delete hRefSelDijetEtaBackward1DOldPt[i];
            if (hRefSelDijetEtaBackward1DOldPtWeighted[i]) delete hRefSelDijetEtaBackward1DOldPtWeighted[i];

            if (hRefSelDijetEta1DOldPtCM[i]) delete hRefSelDijetEta1DOldPtCM[i];
            if (hRefSelDijetEta1DOldPtCMWeighted[i]) delete hRefSelDijetEta1DOldPtCMWeighted[i];
            if (hRefSelRecoDijetEta1DOldPtCM[i]) delete hRefSelRecoDijetEta1DOldPtCM[i];
            if (hRefSelRecoDijetEta1DOldPtCMWeighted[i]) delete hRefSelRecoDijetEta1DOldPtCMWeighted[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i];
            if (hRefSelDijetEtaCMForward1DOldPt[i]) delete hRefSelDijetEtaCMForward1DOldPt[i];
            if (hRefSelDijetEtaCMForward1DOldPtWeighted[i]) delete hRefSelDijetEtaCMForward1DOldPtWeighted[i];
            if (hRefSelDijetEtaCMBackward1DOldPt[i]) delete hRefSelDijetEtaCMBackward1DOldPt[i];
            if (hRefSelDijetEtaCMBackward1DOldPtWeighted[i]) delete hRefSelDijetEtaCMBackward1DOldPtWeighted[i];

            if (hRefSelDijetEta1DOldPtBinning[i]) delete hRefSelDijetEta1DOldPtBinning[i];
            if (hRefSelDijetEta1DOldPtBinningWeighted[i]) delete hRefSelDijetEta1DOldPtBinningWeighted[i];
            if (hRefSelRecoDijetEta1DOldPtBinning[i]) delete hRefSelRecoDijetEta1DOldPtBinning[i];
            if (hRefSelRecoDijetEta1DOldPtBinningWeighted[i]) delete hRefSelRecoDijetEta1DOldPtBinningWeighted[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i];
            if (hRefSelDijetEtaForward1DOldPtBinning[i]) delete hRefSelDijetEtaForward1DOldPtBinning[i];
            if (hRefSelDijetEtaForward1DOldPtBinningWeighted[i]) delete hRefSelDijetEtaForward1DOldPtBinningWeighted[i];
            if (hRefSelDijetEtaBackward1DOldPtBinning[i]) delete hRefSelDijetEtaBackward1DOldPtBinning[i];
            if (hRefSelDijetEtaBackward1DOldPtBinningWeighted[i]) delete hRefSelDijetEtaBackward1DOldPtBinningWeighted[i];

            if (hRefSelDijetEta1DOldPtBinningCM[i]) delete hRefSelDijetEta1DOldPtBinningCM[i];
            if (hRefSelDijetEta1DOldPtBinningCMWeighted[i]) delete hRefSelDijetEta1DOldPtBinningCMWeighted[i];
            if (hRefSelRecoDijetEta1DOldPtBinningCM[i]) delete hRefSelRecoDijetEta1DOldPtBinningCM[i];
            if (hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i]) delete hRefSelRecoDijetEta1DOldPtBinningCMWeighted[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i];
            if (hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]) delete hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i];
            if (hRefSelDijetEtaCMForward1DOldPtBinning[i]) delete hRefSelDijetEtaCMForward1DOldPtBinning[i];
            if (hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i]) delete hRefSelDijetEtaCMForward1DOldPtBinningWeighted[i];
            if (hRefSelDijetEtaCMBackward1DOldPtBinning[i]) delete hRefSelDijetEtaCMBackward1DOldPtBinning[i];
            if (hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i]) delete hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[i];
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


    const int dijetPtBins{16};
    double dijetPtVals[dijetPtBins+1] {  50.,  60.,   70.,  80.,  90.,
                                         100., 110.,  120., 130., 140.,
                                         150., 160.,  180., 200., 250., 
                                         300., 500.};

    // Old binning convention
    const int dijetPtOldBins{6};
    double dijetPtOldVals[dijetPtBins+1] {25., 55., 75., 95., 115., 150., 400.}; // 6 bins
    
    int    prescale = 2;

    int xBins = 10000;
    double xRange[2] = {0.0001, 1.};
    int xPbOverXpBins = 10000;
    double xPbOverXpRange[2] = {0.0001, 1000.};

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

    int    bins4D_jet_JESPtEtaPtHat[4] { fJESBins, fPtBins, fEtaBins, fPtHatBins };
    double xmin4D_jet_JESPtEtaPtHat[4] { fJESRange[0], fPtRange[0], fEtaRange[0], fPtHatRange[0] };
    double xmax4D_jet_JESPtEtaPtHat[4] { fJESRange[1], fPtRange[1], fEtaRange[1], fPtHatRange[1] };

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

    hRecoLeadJetAllPtVsEta = new TH2D("hRecoLeadJetAllPtVsEta","Leading jet all p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoLeadJetAllPtVsEta->Sumw2();

    hRecoInclusiveJetPtRawVsEta = new TH2D("hRecoInclusiveJetPtRawVsEta","Reco jet raw p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveJetPtRawVsEta->Sumw2();

    hRecoSubLeadJetAllPtVsEta = new TH2D("hRecoSubLeadJetAllPtVsEta","Subleading jet all p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1], 
                                            fPtBins, fPtRange[0], fPtRange[1]);
    hRecoSubLeadJetAllPtVsEta->Sumw2();


    hRecoInclusiveJetPt = new TH1D("hRecoInclusiveJetPt","Reco jet p_{T};Reco p_{T} (GeV);Entries",
                                   fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveJetPt->Sumw2();
    hRecoInclusiveJetPtEta = new TH2D("hRecoInclusiveJetPtEta","Reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                      fEtaBins, fEtaRange[0], fEtaRange[1],
                                      fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveJetPtEta->Sumw2();
    hRecoInclusiveJetPtEtaPtHat = new TH3D("hRecoInclusiveJetPtEtaPtHat","Reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                           fEtaBins, fEtaRange[0], fEtaRange[1],
                                           fPtBins, fPtRange[0], fPtRange[1],
                                           fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hRecoInclusiveJetPtEtaPtHat->Sumw2();
    hRecoInclusiveAllJetPtVsEta = new TH2D("hRecoInclusiveAllJetPtVsEta", "Inclusive reco jet pT vs eta;#eta;p_{T} (GeV)",
                                           fEtaBins, fEtaRange[0], fEtaRange[1],
                                           fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtVsEta->Sumw2();

    hRecoGoodInclusiveJetEtaLabFrame = new TH1D("hRecoGoodInclusiveJetEtaLabFrame","Reco good jet #eta in lab frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaLabFrame->Sumw2();
    hRecoGoodInclusiveJetEtaCMFrame = new TH1D("hRecoGoodInclusiveJetEtaCMFrame","Reco good jet #eta in CM frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaCMFrame->Sumw2();

    // Dijet histograms
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi = new THnSparseD("hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi",
            "Reconstructed dijet and jet info;p_{T}^{dijet} (GeV);#eta^{dijet};#Delta#phi^{dijet} (rad);p_{T}^{Leading} (GeV);#eta^{Leading};#phi^{Leading} (rad);p_{T}^{Subleading} (GeV);#eta^{Subleading};#phi^{Subleading} (rad)",
            9,
            bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmax9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Sumw2();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted = new THnSparseD("hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted",
            "Reconstructed dijet and jet info weighted;p_{T}^{dijet} (GeV);#eta^{dijet};#Delta#phi^{dijet} (rad);p_{T}^{Leading} (GeV);#eta^{Leading};#phi^{Leading} (rad);p_{T}^{Subleading} (GeV);#eta^{Subleading};#phi^{Subleading} (rad)",
            9,
            bins9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmin9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi,
            xmax9D_dijet_PtEtaDphiPtEtaPhiPtEtaPhi);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Sumw2();
    

    hRecoPtLeadPtSublead = new TH2D("hRecoPtLeadPtSublead","Reco leading vs subleading p_{T};Reco p_{T}^{Leading} (GeV);Reco p_{T}^{Subleading} (GeV)",
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
    hRecoPtLeadPtSubleadMcReweight = new TH2D("hRecoPtLeadPtSubleadMcReweight","Reco leading vs subleading p_{T} (MC reweighted to data);Reco p_{T}^{Leading} (GeV);Reco p_{T}^{Subleading} (GeV)",
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

    // New ptAve and eta binning
    for (unsigned int i{0}; i<fPtAveBins.size()-1; ++i) {
        double ptAveLow = fPtAveBins.at(i);
        double ptAveHi = fPtAveBins.at(i+1);
        
        // Lab frame
        hRecoDijetEta1D[i] = new TH1D(Form("hRecoDijetEta1D_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                      2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEta1D[i]->Sumw2();
        // hRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DWeighted[i] = new TH1D(Form("hRecoDijetEta1DWeighted_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                              2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEta1DWeighted[i]->Sumw2();
        // hRecoDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2D_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}", i, ptAveLow, ptAveHi),
                                                      fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaLeadVsEtaSubLead2D[i]->Sumw2();
        hRecoDijetEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hRecoDijetEtaLeadVsEtaSubLead2DWeighted_%d",i),Form("Reco #eta^{dijet} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}", i, ptAveLow, ptAveHi),
                                                                fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEtaForward1D[i] = new TH1D(Form("hRecoDijetEtaForward1D_%d",i),Form("Reco #eta^{dijet}_{forward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                             fEtaBins / 2, 0., fEtaRange[1]);
        hRecoDijetEtaForward1D[i]->Sumw2();
        //hRecoDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaForward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaForward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{forward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                     fEtaBins / 2, 0., fEtaRange[1]);
        hRecoDijetEtaForward1DWeighted[i]->Sumw2();
        // hRecoDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1D[i] = new TH1D(Form("hRecoDijetEtaBackward1D_%d",i),Form("Reco #eta^{dijet}_{backward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                              fEtaBins / 2, 0., fEtaRange[1]);
        hRecoDijetEtaBackward1D[i]->Sumw2();
        //hRecoDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaBackward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaBackward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{backward} in the lab frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}", i, ptAveLow, ptAveHi),
                                                      fEtaBins / 2, 0., fEtaRange[1]);
        hRecoDijetEtaBackward1DWeighted[i]->Sumw2();
        //hRecoDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);

        // CM frame
        hRecoDijetEta1DCM[i] = new TH1D(Form("hRecoDijetEta1DCM_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                        2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEta1DCM[i]->Sumw2();
        //hRecoDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoDijetEta1DCMWeighted[i] = new TH1D(Form("hRecoDijetEta1DCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoDijetEta1DCMWeighted[i]->Sumw2();
        //hRecoDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
        hRecoEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DCM_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}", i, ptAveLow, ptAveHi),
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
        hRecoEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hRecoEtaLeadVsEtaSubLead2DCMWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}", i, ptAveLow, ptAveHi),
                                                           fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
        hRecoDijetEtaCMForward1D[i] = new TH1D(Form("hRecoDijetEtaCMForward1D_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                               fEtaBins, 0., fEtaRange[1]);
        hRecoDijetEtaCMForward1D[i]->Sumw2();
        //hRecoDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaCMForward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} forward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                       fEtaBins, 0., fEtaRange[1]);
        hRecoDijetEtaCMForward1DWeighted[i]->Sumw2();
        //hRecoDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1D[i] = new TH1D(Form("hRecoDijetEtaCMBackward1D_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}", i, ptAveLow, ptAveHi),
                                                fEtaBins, 0., fEtaRange[1]);
        hRecoDijetEtaCMBackward1D[i]->Sumw2();
        //hRecoDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
        hRecoDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hRecoDijetEtaCMBackward1DWeighted_%d",i),Form("Reco #eta^{dijet}_{CM} backward in the CM frame in %d for %3.0f<p_{T}^{ave} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM", i, ptAveLow, ptAveHi),
                                                        fEtaBins, 0., fEtaRange[1]);
        hRecoDijetEtaCMBackward1DWeighted[i]->Sumw2();
        //hRecoDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);

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

    hRecoDijetPtEta = new TH2D("hRecoDijetPtEta", "Reco dijet #eta vs p_{T};p_{T}^{ave} (GeV);#eta^{dijet}", 
                               fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                               fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetPtEta->Sumw2();
    hRecoDijetPtEtaDphi = new TH3D("hRecoDijetPtEtaDphi","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphi->Sumw2();
    hRecoDijetPtEtaDphiWeighted = new TH3D("hRecoDijetPtEtaDphiWeighted","Reco dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                           fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                           fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                           fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphiWeighted->Sumw2();
    hRecoDijetEtaCM = new TH1D("hRecoDijetEtaCM","Reco dijet #eta in CM;Reco #eta^{dijet}_{CM};Entries",
                             fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
    hRecoDijetEtaCM->Sumw2();
    hRecoDijetPtEtaDphiCM = new TH3D("hRecoDijetPtEtaDphiCM","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                   fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
    hRecoDijetPtEtaDphiCM->Sumw2();
    hRecoDijetPtEtaDphiCMWeighted = new TH3D("hRecoDijetPtEtaDphiCMWeighted","Reco dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
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

        hGenInclusiveJetPt = new TH1D("hGenInclusiveJetPt","Inclusive gen jet;Gen p_{T}^{inclusive} (GeV)",
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPt->Sumw2();
        hGenInclusiveJetPtEta = new TH2D("hGenInclusiveJetPtEta","Gen inclusive jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPtEta->Sumw2();
        hGenInclusiveJetPtEtaPtHat = new TH3D("hGenInclusiveJetPtEtaPtHat","Gen inclusive jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        hGenInclusiveJetPtEtaPtHat->Sumw2();

        //
        // Gen dijets
        //
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi = new THnSparseD("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi","Title;Gen p_{T}^{dijet} (GeV);Gen #eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
                9, 
                bins9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmin9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmax9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi);
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Sumw2();
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted = new THnSparseD("hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted","Title;Gen p_{T}^{dijet} (GeV);Gen #eta^{dijet};Gen #Delta#phi^{dijet} (rad);Gen p_{T}^{Lead} (GeV);Gen #eta^{Lead};Gen #phi^{Lead} (rad);Gen p_{T}^{Sublead} (GeV);Gen #eta^{Sublead};Gen #phi^{Sublead} (rad)",
                9, 
                bins9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmin9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi,
                xmax9D_gen_GenDijetPtEtaDPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi);
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Sumw2();


        hGenPtLeadPtSublead = new TH2D("hGenPtLeadPtSublead","Leading gen jet pT vs subleading gen jet pT;Gen p_{T}^{Leading} (GeV);Gen p_{T}^{Subleading} (GeV)",
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
        hGenPtLeadPtSubleadMcReweight = new TH2D("hGenPtLeadPtSubleadMcReweight","Leading gen jet pT vs subleading gen jet pT (MC reweighted to data);Gen p_{T}^{Leading} (GeV);Gen p_{T}^{Subleading} (GeV)",
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

        for (int i=0; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hGenDijetEta1D[i] = new TH1D(Form("hGenDijetEta1D_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEta1D[i]->Sumw2();
            //hGenDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DWeighted[i] = new TH1D(Form("hGenDijetEta1DWeighted_%d",i), Form("Gen #eta^{dijet} in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEta1DWeighted[i]->Sumw2();
            //hGenDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2D_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead};#eta^{SubLead}",i, ptAveLow, ptAveHi),
                                                            fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DWeighted[i]->Sumw2();
            hGenDijetEtaForward1D[i] = new TH1D(Form("hGenDijetEtaForward1D_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaForward1D[i]->Sumw2();
            //hGenDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
            hGenDijetEtaForward1DWeighted[i] = new TH1D(Form("hGenDijetEtaForward1DWeighted_%d",i), Form("Gen #eta^{dijet} forward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaForward1DWeighted[i]->Sumw2();
            //hGenDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
            hGenDijetEtaBackward1D[i] = new TH1D(Form("hGenDijetEtaBackward1D_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaBackward1D[i]->Sumw2();
            //hGenDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
            hGenDijetEtaBackward1DWeighted[i] = new TH1D(Form("hGenDijetEtaBackward1DWeighted_%d",i), Form("Gen #eta^{dijet} backward in the lab frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaBackward1DWeighted[i]->Sumw2();
            //hGenDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);

            hGenDijetEta1DCM[i] = new TH1D(Form("hGenDijetEta1DCM_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                           2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEta1DCM[i]->Sumw2();
            //hGenDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEta1DCMWeighted[i] = new TH1D(Form("hGenDijetEta1DCMWeighted_%d",i), Form("Gen #eta^{dijet} in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEta1DCMWeighted[i]->Sumw2();
            //hGenDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hGenDijetEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DCM_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hGenDijetEtaLeadVsEtaSubLead2DCMWeighted_%d",i), Form("Gen #eta^{dijet} leading vs subleading #eta in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{Lead}_{CM};#eta^{SubLead}_{CM}",i, ptAveLow, ptAveHi),
                                                            fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
            hGenDijetEtaCMForward1D[i] = new TH1D(Form("hGenDijetEtaCMForward1D_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                  fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaCMForward1D[i]->Sumw2();
            //hGenDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
            hGenDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hGenDijetEtaCMForward1DWeighted_%d",i), Form("Gen #eta^{dijet} forward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaCMForward1DWeighted[i]->Sumw2();
            //hGenDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1D[i] = new TH1D(Form("hGenDijetEtaCMBackward1D_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaCMBackward1D[i]->Sumw2();
            //hGenDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
            hGenDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hGenDijetEtaCMBackward1DWeighted_%d",i), Form("Gen #eta^{dijet} backward in the CM frame in %d for %3.0f<p_{T} (GeV)<%3.0f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                           fEtaBins, 0., fEtaRange[1]);
            hGenDijetEtaCMBackward1DWeighted[i]->Sumw2();
            //hGenDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBVals);
        }

        for (int i=0; i<fPtAveOldBins.size()-1; i++) {

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

        hGenDijetPtEtaDphi = new TH3D("hGenDijetPtEtaDphi","Gen dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                      fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                      fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                      fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphi->Sumw2();
        hGenDijetPtEtaDphiWeighted = new TH3D("hGenDijetPtEtaDphiWeighted","Gen dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                              fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                              fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                              fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphiWeighted->Sumw2();

        hGenDijetEtaCM = new TH1D("hGenDijetEtaCM", "Gen dijet #eta in CM;#eta^{dijet}_{CM}",
                                  fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hGenDijetEtaCM->Sumw2();
        hGenDijetPtEtaDphiCM = new TH3D("hGenDijetPtEtaDphiCM","Gen dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hGenDijetPtEtaDphiCM->Sumw2();
        hGenDijetPtEtaDphiCMWeighted = new TH3D("hGenDijetPtEtaDphiCMWeighted","Gen dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
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
        hRecoInclusiveUnmatchedJetPtVsEta = new TH2D("hRecoInclusiveUnmatchedJetPtVsEta", "Inclusive reco jet unmatched gen pT vs eta;#eta;p_{T} (GeV)",
                                                     fEtaBins, fEtaRange[0], fEtaRange[1],
                                                     fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveUnmatchedJetPtVsEta->Sumw2();

        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed inclusive jets;Reco p_{T, corr}^{Inclusive} (GeV);Reco p_{T, raw}^{Inclusive} (GeV);Ref p_{T}^{Inclusive} (GeV);Reco #eta^{Inclusive};Ref #eta^{Inclusive}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted = new THnSparseD("hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted","Reconstructed inclusive jets weighted;Reco p_{T, corr}^{Inclusive} (GeV);Reco p_{T, raw}^{Inclusive} (GeV);Ref p_{T}^{Inclusive} (GeV);Reco #eta^{Inclusive};Ref #eta^{Inclusive}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Sumw2();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed leading jets;Reco p_{T, corr}^{Leading} (GeV);Reco p_{T, raw}^{Leading} (GeV);Ref p_{T}^{Leading} (GeV);Reco #eta^{Leading};Ref #eta^{Leading}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted = new THnSparseD("hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted","Reconstructed leading jets weighted;Reco p_{T, corr}^{Leading} (GeV);Reco p_{T, raw}^{Leading} (GeV);Ref p_{T}^{Leading} (GeV);Reco #eta^{Leading};Ref #eta^{Leading}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Sumw2();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen = new THnSparseD("hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen","Reconstructed subleading jets;Reco p_{T, corr}^{Subleading} (GeV);Reco p_{T, raw}^{Subleading} (GeV);Ref p_{T}^{Subleading} (GeV);Reco #eta^{Subleading};Ref #eta^{Subleading}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Sumw2();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted = new THnSparseD("hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted","Reconstructed subleading jets weighted;Reco p_{T, corr}^{Subleading} (GeV);Reco p_{T, raw}^{Subleading} (GeV);Ref p_{T}^{Subleading} (GeV);Reco #eta^{Subleading};Ref #eta^{Subleading}",
                5,
                bins5D_jet_PtPtPtEtaEta,
                xmin5D_jet_PtPtPtEtaEta,
                xmax5D_jet_PtPtPtEtaEta);
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Sumw2();

        hJESInclusiveJetPtEtaPhi = new THnSparseD("hJESInclusiveJetPtEtaPhi","JES of inclusive jets;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV);#eta^{reco};#phi^{reco} (rad)",
                4,
                bins4D_jet_JESPtEtaPhi,
                xmin4D_jet_JESPtEtaPhi,
                xmax4D_jet_JESPtEtaPhi);
        hJESInclusiveJetPtEtaPhi->Sumw2();
        hJESInclusiveJetPtEtaPhiWeighted = new THnSparseD("hJESInclusiveJetPtEtaPhiWeighted","JES of inclusive jet weighted;p_{T}^{reco}/p_{T}^{gen};Ref p_{T} (GeV);#eta^{reco};#phi^{reco} (rad)",
                4,
                bins4D_jet_JESPtEtaPhi,
                xmin4D_jet_JESPtEtaPhi,
                xmax4D_jet_JESPtEtaPhi);
        hJESInclusiveJetPtEtaPhiWeighted->Sumw2();


        hRecoInclusiveJetJECFactorVsPtEta = new TH3D("hRecoInclusiveJetJECFactorVsPtEta","JEC factor vs p_{T} and #eta;p_{T}^{corr}/p_{T}^{raw};p_{T}^{gen} (GeV);#eta^{gen};JEC factor",
                                           20, 0., 2., fPtBins, fPtRange[0], fPtRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetJECFactorVsPtEta->Sumw2(); 
        hRecoInclusiveJetJEC2FactorVsPtEta = new TH3D("hRecoInclusiveJetJEC2FactorVsPtEta","JEC factor (my) vs p_{T} and #eta;p_{T}^{my, corr}/p_{T}^{raw};p_{T}^{gen} (GeV);#eta^{gen};JEC2 factor",
                                           20, 0., 2., fPtBins, fPtRange[0], fPtRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetJEC2FactorVsPtEta->Sumw2();
        hRecoInclusiveJetPtRawOverPtRefVsPtEta = new TH3D("hRecoInclusiveJetPtRawOverPtRefVsPtEta","p_{T}^{raw}/p_{T}^{ref} vs p_{T} and #eta;p_{T}^{raw}/p_{T}^{ref};p_{T}^{gen} (GeV);#eta^{gen}",
                                           20, 0., 2., fPtBins, fPtRange[0], fPtRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRecoInclusiveJetPtRawOverPtRefVsPtEta->Sumw2();
        hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning = new TH3D("hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning","p_{T}^{raw}/p_{T}^{ref} vs p_{T} and #eta (std binning);p_{T}^{raw}/p_{T}^{ref};p_{T}^{gen} (GeV);#eta^{gen}",
                                                                    20, 0., 2., 1300, 10., 6510., 104, -5.2, 5.2);
        hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->Sumw2();
        hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->GetZaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);


        hInclusiveJetJESVsPtGen = new TH2D("hInclusiveJetJESVsPtGen","JES vs p_{T}^{gen} for |#eta|<1.4;p_{T}^{gen} (GeV);p_{T}^{reco}/p_{T}^{gen}",
                                           fPtBins, fPtRange[0], fPtRange[1], fJESBins, fJESRange[0], fJESRange[1]);
        hInclusiveJetJESVsPtGen->Draw();
        hInclusiveJetJESGenPtGenEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESGenPtGenEtaPtHatWeighted","JES vs p_{T}^{gen} vs #eta^{gen} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{gen} (GeV);#eta^{gen};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hInclusiveJetJESGenPtGenEtaPtHatWeighted->Sumw2();
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESRecoPtRecoEtaPtHatWeighted","JES vs p_{T}^{reco} vs #eta^{reco} vs #hat{p}_{T} weighted;p_{T}^{reco}/p_{T}^{gen};p_{T}^{reco} (GeV);#eta^{reco};#hat{p}_{T} (GeV)",
                4,
                bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Sumw2();

        hRecoInclusiveMatchedJetPt = new TH1D("hRecoInclusiveMatchedJetPt","Inclusive reco jet that has matching;p_{T} (GeV)",
                                                    fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveMatchedJetPt->Sumw2();
        hRecoInclusiveMatchedJetPtVsEta = new TH2D("hRecoInclusiveMatchedJetPtVsEta","Inclusive reco jet that has matching pT vs eta;#eta;p_{T} (GeV)",
                                                    fEtaBins, fEtaRange[0], fEtaRange[1],
                                                    fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveMatchedJetPtVsEta->Sumw2();



        hRecoLeadJetMatchedPtVsEta = new TH2D("hRecoLeadJetMatchedPtVsEta","Leading jet matched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                              fEtaBins, fEtaRange[0], fEtaRange[1], 
                                              fPtBins, fPtRange[0], fPtRange[1]);
        hRecoLeadJetMatchedPtVsEta->Sumw2();
        hRecoLeadJetUnmatchedPtVsEta = new TH2D("hRecoLeadJetUnmatchedPtVsEta","Leading jet unmatched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                                fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                fPtBins, fPtRange[0], fPtRange[1]);
        hRecoLeadJetUnmatchedPtVsEta->Sumw2();
        hRecoSubLeadJetMatchedPtVsEta = new TH2D("hRecoSubLeadJetMatchedPtVsEta","Subleading jet matched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                                 fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                 fPtBins, fPtRange[0], fPtRange[1]);
        hRecoSubLeadJetMatchedPtVsEta->Sumw2();
        hRecoSubLeadJetUnmatchedPtVsEta = new TH2D("hRecoSubLeadJetUnmatchedPtVsEta","Subleading jet unmatched p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                                   fEtaBins, fEtaRange[0], fEtaRange[1], 
                                                   fPtBins, fPtRange[0], fPtRange[1]);
        hRecoSubLeadJetUnmatchedPtVsEta->Sumw2();


        // Reco dijet with MC
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta = new THnSparseD("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta",
                "Reco to ref correspondence;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV);Ref #eta^{Subleading}",
                12,
                bins12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmin12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmax12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta);
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Sumw2();
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted = new THnSparseD("hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted",
                "Reco to ref correspondence weighted;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV);Ref #eta^{Subleading}",
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

        hRefSelInclusiveJetPt = new TH1D("hRefSelInclusiveJetPt","Ref-selected jet p_{T};Ref p_{T} (GeV);Entries",
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelInclusiveJetPt->Sumw2();
        hRefSelInclusiveJetPtEta = new TH2D("hRefSelInclusiveJetPtEta","Ref-selected jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelInclusiveJetPtEta->Sumw2();
        hRefSelInclusiveJetPtEtaPtHat = new TH3D("hRefSelInclusiveJetPtEtaPtHat","Ref-selected jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefSelInclusiveJetPtEtaPtHat->Sumw2();

        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted = new THnSparseD("hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted",
                "Reco to ref correspondence (via ref selection) weighted;Reco p_{T}^{dijet} (GeV);Reco #eta^{dijet};Reco p_{T}^{Leading} (GeV);Reco #eta^{Leading};Reco p_{T}^{Subleading} (GeV);Reco #eta^{Subleading};Ref p_{T}^{dijet};Ref #eta^{dijet};Ref p_{T}^{Leading} (GeV);Ref #eta^{Leading};Ref p_{T}^{Subleading} (GeV);Ref #eta^{Subleading}",
                12,
                bins12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmin12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta,
                xmax12D_dijet_PtEtaPtEtaPtEtaPtEtaPtEtaPtEta);
        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Sumw2();
        hRefSelDijetEta = new TH1D("hRefSelDijetEta","Ref selected dijets;#eta^{dijet}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEta->Sumw2();
        hRefSelDijetPtEtaDphi = new TH3D("hRefSelDijetPtEtaDphi","RefSel dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphi->Sumw2();
        hRefSelDijetPtEtaDphiWeighted = new TH3D("hRefSelDijetPtEtaDphiWeighted","RefSel dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphiWeighted->Sumw2();
        hRefSelDijetEtaCM = new TH1D("hRefSelDijetEtaCM","Ref selected dijets in CM;#eta^{dijet}_{CM}",
                                    fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefSelDijetEtaCM->Sumw2();
        hRefSelDijetPtEtaDphiCM = new TH3D("hRefSelDijetPtEtaDphiCM","RefSel dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                         fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                         fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                         fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphiCM->Sumw2();
        hRefSelDijetPtEtaDphiCMWeighted = new TH3D("hRefSelDijetPtEtaDphiCMWeighted","RefSel dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                                fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefSelDijetPtEtaDphiCMWeighted->Sumw2();

        // New pT and eta binning
        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hRefSelDijetEta1D[i] = new TH1D(Form("hRefSelDijetEta1D_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                            2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelDijetEta1D[i]->Sumw2();
            //hRefSelDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DWeighted[i] = new TH1D(Form("hRefSelDijetEta1DWeighted_%d",i),Form("Ref selected #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelDijetEta1DWeighted[i]->Sumw2();
            //hRefSelDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1D[i] = new TH1D(Form("hRefSelRecoDijetEta1D_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelRecoDijetEta1D[i]->Sumw2();
            //hRefSelRecoDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelRecoDijetEta1DWeighted[i]->Sumw2();
            //hRefSelRecoDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2D_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DWeighted_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DWeighted[i]->Sumw2();
            hRefSelDijetEtaForward1D[i] = new TH1D(Form("hRefSelDijetEtaForward1D_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                   fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaForward1D[i]->Sumw2();
            //hRefSelDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
            hRefSelDijetEtaForward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaForward1DWeighted_%d",i),Form("Ref selected #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                           fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaForward1DWeighted[i]->Sumw2();
            //hRefSelDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
            hRefSelDijetEtaBackward1D[i] = new TH1D(Form("hRefSelDijetEtaBackward1D_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaBackward1D[i]->Sumw2();
            //hRefSelDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
            hRefSelDijetEtaBackward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaBackward1DWeighted_%d",i),Form("Ref selected #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                            fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaBackward1DWeighted[i]->Sumw2();
            //hRefSelDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);


            hRefSelDijetEta1DCM[i] = new TH1D(Form("hRefSelDijetEta1DCM_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                              2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelDijetEta1DCM[i]->Sumw2();
            //hRefSelDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelDijetEta1DCMWeighted[i] = new TH1D(Form("hRefSelDijetEta1DCMWeighted_%d",i),Form("Ref selected #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                      2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelDijetEta1DCMWeighted[i]->Sumw2();
            //hRefSelDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DCM[i] = new TH1D(Form("hRefSelRecoDijetEta1DCM_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                  2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelRecoDijetEta1DCM[i]->Sumw2();
            //hRefSelRecoDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelRecoDijetEta1DCMWeighted[i] = new TH1D(Form("hRefSelRecoDijetEta1DCMWeighted_%d",i),Form("Ref selected reco #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelRecoDijetEta1DCMWeighted[i]->Sumw2();
            //hRefSelDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefSelEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DCM_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DCMWeighted_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
            hRefSelDijetEtaCMForward1D[i] = new TH1D(Form("hRefSelDijetEtaCMForward1D_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                     fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaCMForward1D[i]->Sumw2();
            //hRefSelDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
            hRefSelDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMForward1DWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaCMForward1DWeighted[i]->Sumw2();
            //hRefSelDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
            hRefSelDijetEtaCMBackward1D[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1D_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaCMBackward1D[i]->Sumw2();
            //hRefSelDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
            hRefSelDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hRefSelDijetEtaCMBackward1DWeighted_%d",i),Form("Ref selected #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                             fEtaBins, 0., fEtaRange[1]);
            hRefSelDijetEtaCMBackward1DWeighted[i]->Sumw2();
            //hRefSelDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBins, dijetEtaFBins);
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
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
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
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                       fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
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
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
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
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                                   fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted_%d",i),Form("Ref selected #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
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


        //
        // Ref jet histograms
        //

        hRefInclusiveJetPt = new TH1D("hRefInclusiveJetPt","Ref jet p_{T};Ref p_{T} (GeV);Entries",
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPt->Sumw2();
        hRefInclusiveJetPtEta = new TH2D("hRefInclusiveJetPtEta","Ref jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPtEta->Sumw2();
        hRefInclusiveJetPtEtaPtHat = new TH3D("hRefInclusiveJetPtEtaPtHat","Ref jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefInclusiveJetPtEtaPtHat->Sumw2();

        // Ref dijets
        hRefPtLeadPtSublead = new TH2D("hRefPtLeadPtSublead","Ref leading vs subleading p_{T};Ref p_{T}^{Leading} (GeV);Ref p_{T}^{Subleading} (GeV)",
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
        hRefPtLeadPtSubleadMcReweight = new TH2D("hRefPtLeadPtSubleadMcReweight","Ref leading vs subleading p_{T} (MC reweighted to data);Ref p_{T}^{Leading} (GeV);Ref p_{T}^{Subleading} (GeV)",
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

        // New pT and eta binning
        for (unsigned int i{0}; i<fPtAveBins.size()-1; i++) {
            double ptAveLow = fPtAveBins.at(i);
            double ptAveHi = fPtAveBins.at(i+1);
            hRefDijetEta1D[i] = new TH1D(Form("hRefDijetEta1D_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                         2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefDijetEta1D[i]->Sumw2();
            //hRefDijetEta1D[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DWeighted[i] = new TH1D(Form("hRefDijetEta1DWeighted_%d",i),Form("Ref #eta^{dijet} in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefDijetEta1DWeighted[i]->Sumw2();
            //hRefDijetEta1DWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefEtaLeadVsEtaSubLead2D[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2D_%d",i),Form("Ref #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2D[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DWeighted_%d",i),Form("Ref #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2D[i] = new TH2D(Form("hRecoVsRefDijetEta2D_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2D[i]->Sumw2();
            hRecoVsRefDijetEta2DWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DWeighted_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2D[i] = new TH2D(Form("hRecoVsRefLeadJetEta2D_%d",i),Form("Reco vs Ref #eta^{Leading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Leading};Ref #eta^{Leading}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2D[i]->Sumw2();
            hRecoVsRefLeadJetEta2DWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DWeighted_%d",i),Form("Reco vs Ref #eta^{Leading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Leading};Ref #eta^{Leading}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2D[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2D_%d",i),Form("Reco vs Ref #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Subleading};Ref #eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2D[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DWeighted_%d",i),Form("Reco vs Ref #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Subleading};Ref #eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DWeighted[i]->Sumw2();
            hRefDijetEtaForward1D[i] = new TH1D(Form("hRefDijetEtaForward1D_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaForward1D[i]->Sumw2();
            //hRefDijetEtaForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaForward1DWeighted[i] = new TH1D(Form("hRefDijetEtaForward1DWeighted_%d",i),Form("Ref #eta^{dijet} forward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaForward1DWeighted[i]->Sumw2();
            //hRefDijetEtaForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1D[i] = new TH1D(Form("hRefDijetEtaBackward1D_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                 fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaBackward1D[i]->Sumw2();
            //hRefDijetEtaBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaBackward1DWeighted[i] = new TH1D(Form("hRefDijetEtaBackward1DWeighted_%d",i),Form("Ref #eta^{dijet} backward in the lab frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet};dN/d#eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaBackward1DWeighted[i]->Sumw2();
            //hRefDijetEtaBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);

            // Ref dijets in CM frame
            hRefDijetEta1DCM[i] = new TH1D(Form("hRefDijetEta1DCM_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                           2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefDijetEta1DCM[i]->Sumw2();
            //hRefDijetEta1DCM[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefDijetEta1DCMWeighted[i] = new TH1D(Form("hRefDijetEta1DCMWeighted_%d",i),Form("Ref #eta^{dijet} in the CM frame for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   2 * fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefDijetEta1DCMWeighted[i]->Sumw2();
            //hRefDijetEta1DCMWeighted[i]->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
            hRefEtaLeadVsEtaSubLead2DCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DCM_%d",i),Form("Ref #eta^{Leading}_{CM} vs #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DCMWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DCMWeighted_%d",i),Form("Ref #eta^{Leading}_{CM} vs #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DCMWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DCM[i]->Sumw2();
            hRecoVsRefDijetEta2DCMWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DCMWeighted_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                    dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DCMWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DCM_%d",i),Form("Reco vs Ref #eta^{Leading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Leading}_{CM};Ref #eta^{Leading}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DCMWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DCMWeighted_%d",i),Form("Reco vs Ref #eta^{Leading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Leading}_{CM};Ref #eta^{Leading}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DCMWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DCM_%d",i),Form("Reco vs Ref #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Subleading}_{CM};Ref #eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DCMWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DCMWeighted_%d",i),Form("Reco vs Ref #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Subleading}_{CM};Ref #eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                    fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DCMWeighted[i]->Sumw2();
            hRefDijetEtaCMForward1D[i] = new TH1D(Form("hRefDijetEtaCMForward1D_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                  fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaCMForward1D[i]->Sumw2();
            //hRefDijetEtaCMForward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMForward1DWeighted[i] = new TH1D(Form("hRefDijetEtaCMForward1DWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} forward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaCMForward1DWeighted[i]->Sumw2();
            //hRefDijetEtaCMForward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1D[i] = new TH1D(Form("hRefDijetEtaCMBackward1D_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                   fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaCMBackward1D[i]->Sumw2();
            //hRefDijetEtaCMBackward1D[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
            hRefDijetEtaCMBackward1DWeighted[i] = new TH1D(Form("hRefDijetEtaCMBackward1DWeighted_%d",i),Form("Ref #eta^{dijet}_{CM} backward for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{dijet}_{CM};dN/d#eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                          fEtaBins, 0., fEtaRange[1]);
            hRefDijetEtaCMBackward1DWeighted[i]->Sumw2();
            //hRefDijetEtaCMBackward1DWeighted[i]->GetXaxis()->Set(dijetEtaFBBins, dijetEtaFBVals);
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
            hRefEtaLeadVsEtaSubLead2DOldPt[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPt_%d",i),Form("Ref #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPt[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtWeighted_%d",i),Form("Ref #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPt[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtWeighted_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{Leading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Leading};Ref #eta^{Leading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPt[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtWeighted_%d",i),Form("Reco vs Ref #eta^{Leading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Leading};Ref #eta^{Leading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPt[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPt_%d",i),Form("Reco vs Ref #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Subleading};Ref #eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPt[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtWeighted_%d",i),Form("Reco vs Ref #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Subleading};Ref #eta^{Subleading}",i, ptAveLow, ptAveHi),
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
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtCM_%d",i),Form("Ref #eta^{Leading}_{CM} vs #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted_%d",i),Form("Ref #eta^{Leading}_{CM} vs #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtCMWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtCMWeighted_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaBins, dijetEtaVals, dijetEtaBins, dijetEtaVals);
            hRecoVsRefDijetEta2DOldPtCMWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{Leading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Leading}_{CM};Ref #eta^{Leading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtCMWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtCMWeighted_%d",i),Form("Reco vs Ref #eta^{Leading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Leading}_{CM};Ref #eta^{Leading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtCMWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtCM_%d",i),Form("Reco vs Ref #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Subleading}_{CM};Ref #eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtCMWeighted_%d",i),Form("Reco vs Ref #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Subleading}_{CM};Ref #eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
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
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinning_%d",i),Form("Ref #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinning[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted_%d",i),Form("Ref #eta^{Leading} vs #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading};#eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningWeighted_%d",i),Form("Reco vs Ref #eta^{dijet} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet};Ref #eta^{dijet}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{Leading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Leading};Ref #eta^{Leading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);                                  
            hRecoVsRefLeadJetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningWeighted_%d",i),Form("Reco vs Ref #eta^{Leading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Leading};Ref #eta^{Leading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinning_%d",i),Form("Reco vs Ref #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Subleading};Ref #eta^{Subleading}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinning[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted_%d",i),Form("Reco vs Ref #eta^{Subleading} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Subleading};Ref #eta^{Subleading}",i, ptAveLow, ptAveHi),
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
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningCM_%d",i),Form("Ref #eta^{Leading}_{CM} vs #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[i]->Sumw2();
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted_%d",i),Form("Ref #eta^{Leading}_{CM} vs #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;#eta^{Leading}_{CM};#eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoVsRefDijetEta2DOldPtBinningCMWeighted_%d",i),Form("Reco vs Ref #eta^{dijet}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}",i, ptAveLow, ptAveHi),
                                                        dijetEtaOldBins, dijetEtaOldVals, dijetEtaOldBins, dijetEtaOldVals);
            hRecoVsRefDijetEta2DOldPtBinningCMWeighted[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{Leading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Leading}_{CM};Ref #eta^{Leading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted_%d",i),Form("Reco vs Ref #eta^{Leading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Leading}_{CM};Ref #eta^{Leading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningCM_%d",i),Form("Reco vs Ref #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f;Reco #eta^{Subleading}_{CM};Ref #eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
                                                        fEtaBins, fEtaRange[0], fEtaRange[1], fEtaBins, fEtaRange[0], fEtaRange[1]);
            hRecoVsRefSubLeadJetEta2DOldPtBinningCM[i]->Sumw2();
            hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[i] = new TH2D(Form("hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted_%d",i),Form("Reco vs Ref #eta^{Subleading}_{CM} for %d in range %3.f<p_{T}^{ave} (GeV)<%3.f weighted;Reco #eta^{Subleading}_{CM};Ref #eta^{Subleading}_{CM}",i, ptAveLow, ptAveHi),
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
        hRefDijetPtEtaDphi = new TH3D("hRefDijetPtEtaDphi","Ref dijet info;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphi->Sumw2();
        hRefDijetPtEtaDphiWeighted = new TH3D("hRefDijetPtEtaDphiWeighted","Ref dijet info weighted;p_{T}^{ave} (GeV);#eta^{dijet};#Delta#phi (rad)",
                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                            fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphiWeighted->Sumw2();

        hRefDijetEtaCM = new TH1D("hRefDijetEtaCM","Ref dijet #eta in CM;Ref #eta^{dijet}_{CM};Entries",
                                fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1]);
        hRefDijetEtaCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV)",
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                            fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Sumw2();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted = new TH3D("hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted", "Reco dijet #eta vs ref dijet #eta vs reco dijet p_{T} weighted in CM;Reco #eta^{dijet}_{CM};Ref #eta^{dijet}_{CM}; Reco p_{T}^{dijet} (GeV)",
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                                                   fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1] );
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Sumw2();
        hRefDijetPtEtaDphiCM = new TH3D("hRefDijetPtEtaDphiCM","Ref dijet info in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
                                        fDijetPtBins, fDijetPtRange[0], fDijetPtRange[1],
                                        fDijetEtaBins, fDijetEtaRange[0], fDijetEtaRange[1],
                                        fDijetDphiBins, fDijetDphiRange[0], fDijetDphiRange[1] );
        hRefDijetPtEtaDphiCM->Sumw2();
        hRefDijetPtEtaDphiCMWeighted = new TH3D("hRefDijetPtEtaDphiCMWeighted","Ref dijet info weighted in CM;p_{T}^{ave} (GeV);#eta^{dijet}_{CM};#Delta#phi (rad)",
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


        //
        // Modify axis binning
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

        //
        // Modify bins of gen dijets
        //
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
        hGenDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
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
    // Modify binning of reco histograms
    //
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->GetAxis(1)->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetEta->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphi->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetEtaCM->GetXaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiCM->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);
    hRecoDijetPtEtaDphiCMWeighted->GetYaxis()->Set(dijetEtaBins, dijetEtaVals);

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
    hRecoInclusiveAllJetPtVsEta->Write();
    hRecoInclusiveJetPtRawVsEta->Write();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Write();
    hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Write();
    hRecoInclusiveJetPt->Write();
    hRecoInclusiveJetPtEta->Write();
    hRecoInclusiveJetPtEtaPtHat->Write();
    hRecoPtLeadPtSublead->Write();
    hRecoEtaLeadEtaSublead->Write();
    hRecoEtaCMLeadEtaCMSublead->Write();
    hRecoPtLeadPtSubleadMcReweight->Write();
    hRecoEtaLeadEtaSubleadMcReweight->Write();
    hRecoDijetPtEta->Write();
    hRecoDijetPtEtaForward->Write();
    hRecoDijetPtEtaBackward->Write();
    hRecoDijetPtEtaCMForward->Write();
    hRecoDijetPtEtaCMBackward->Write();
    hRecoDijetPtEtaForwardWeighted->Write();
    hRecoDijetPtEtaBackwardWeighted->Write();
    hRecoDijetPtEtaCMForwardWeighted->Write();
    hRecoDijetPtEtaCMBackwardWeighted->Write();

    hRecoDijetPtEtaDphi->Write();
    hRecoDijetPtEtaDphiWeighted->Write();
    hRecoDijetPtEtaDphiCM->Write();
    hRecoDijetPtEtaDphiCMWeighted->Write();

    hRecoLeadJetAllPtVsEta->Write();
    hRecoSubLeadJetAllPtVsEta->Write();

    hRecoGoodInclusiveJetEtaLabFrame->Write();
    hRecoGoodInclusiveJetEtaCMFrame->Write();

    hRecoDijetEta->Write();
    hRecoDijetEtaCM->Write();

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
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Write();
        hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Write();
        hGenInclusiveJetPt->Write();
        hGenInclusiveJetPtEta->Write();
        hGenInclusiveJetPtEtaPtHat->Write();
        hGenPtLeadPtSublead->Write();
        hGenEtaLeadEtaSublead->Write();
        hGenEtaCMLeadEtaCMSublead->Write();
        hGenPtLeadPtSubleadMcReweight->Write();
        hGenEtaLeadEtaSubleadMcReweight->Write();
        hGenDijetEta->Write();
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

        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Write();
        hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Write();
        hJESInclusiveJetPtEtaPhi->Write();
        hJESInclusiveJetPtEtaPhiWeighted->Write();

        hRecoInclusiveJetJECFactorVsPtEta->Write();
        hRecoInclusiveJetJEC2FactorVsPtEta->Write();
        hRecoInclusiveJetPtRawOverPtRefVsPtEta->Write();
        hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->Write();

        hInclusiveJetJESVsPtGen->Write();
        hInclusiveJetJESGenPtGenEtaPtHatWeighted->Write();
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Write();

        hRecoInclusiveMatchedJetPt->Write();
        hRecoInclusiveMatchedJetPtVsEta->Write();
        hRecoInclusiveUnmatchedJetPtVsEta->Write();
        hRecoLeadJetMatchedPtVsEta->Write();
        hRecoLeadJetUnmatchedPtVsEta->Write();
        hRecoSubLeadJetMatchedPtVsEta->Write(); 
        hRecoSubLeadJetUnmatchedPtVsEta->Write();

        hRefInclusiveJetPt->Write();
        hRefInclusiveJetPtEta->Write();
        hRefInclusiveJetPtEtaPtHat->Write();

        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Write();
        hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Write();

        hRecoDijetPtEtaRefDijetPtEta->Write();
        hRecoDijetPtEtaRefDijetPtEtaWeighted->Write();

        hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Write();

        hRefDijetEta->Write();
        hRefDijetEtaVsRecoDijetEta->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Write();
        hRefDijetPtEtaDphi->Write();
        hRefDijetPtEtaDphiWeighted->Write();

        hRefDijetPtEtaForward->Write();
        hRefDijetPtEtaBackward->Write();
        hRefDijetPtEtaCMForward->Write();
        hRefDijetPtEtaCMBackward->Write();
        hRefDijetPtEtaForwardWeighted->Write();
        hRefDijetPtEtaBackwardWeighted->Write();
        hRefDijetPtEtaCMForwardWeighted->Write();
        hRefDijetPtEtaCMBackwardWeighted->Write();

        hRefDijetEtaCM->Write();
        hRefDijetPtEtaDphiCM->Write();
        hRefDijetPtEtaDphiCMWeighted->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Write();
        hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Write();

        hRefPtLeadPtSublead->Write();
        hRefEtaLeadEtaSublead->Write();
        hRefEtaCMLeadEtaCMSublead->Write();
        hRefPtLeadPtSubleadMcReweight->Write();
        hRefEtaLeadEtaSubleadMcReweight->Write();

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
        // Ref-seletected histograms
        //

        hRefSelInclusiveJetPt->Sumw2();
        hRefSelInclusiveJetPtEta->Sumw2();
        hRefSelInclusiveJetPtEtaPtHat->Sumw2();

        hRefSelDijetEta->Write();
        hRefSelDijetPtEtaDphi->Write();
        hRefSelDijetPtEtaDphiWeighted->Write();
        hRefSelDijetEtaCM->Write();
        hRefSelDijetPtEtaDphiCM->Write();
        hRefSelDijetPtEtaDphiCMWeighted->Write();

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