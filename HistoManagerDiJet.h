/**
 * @file HistoManagerDiJet.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for dijet studies
 * @version 1.2
 * @date 2025-01-09
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef HistoManagerDiJet_h
#define HistoManagerDiJet_h

// Jet analysis headers
#include "BaseHistoManager.h"

// ROOT headers
#include "TObject.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TMath.h"

// C++ headers
#include <vector>

//________________
class HistoManagerDiJet : public BaseHistoManager {
  public:
    /// @brief Constructor
    HistoManagerDiJet();
    /// @brief Destructor
    virtual ~HistoManagerDiJet();

    /// @brief Initialize and create histograms
    void init();
    /// @brief Use MC histograms
    void setIsMc(const bool& isMc = true) { fIsMc = isMc; }
    /// @brief Write all objects to the output file
    void writeOutput();

    // /// @brief Set number of centrality bins
    // void setCentBins(const int& n = 6) { fCentBins = n; }
    // /// @brief Set centrality range (in percentage)
    // void setCentRange(const double& lo = 0, const double& hi = 60) { fCentRange[0]=lo; fCentRange[1]=hi; }

    /// @brief Set number of jet pT bins
    void setJetPtBins(const int& n = 200) { fPtBins = n; }
    /// @brief Set jet pT range
    void setJetPtRange(const double& lo=0., const double& hi=1000.) { fPtRange[0]=lo; fPtRange[1]=hi; }
    /// @brief Set number of eta bins
    void setJetEtaBins(const int& n=50) { fEtaBins = n; }
    /// @brief Set eta range
    void setJetEtaRange(const double& lo=-2.5, const double& hi=2.5) { fEtaRange[0]=lo; fEtaRange[1]=hi; }
    /// @brief Set number of phi bins
    void setJetPhiBins(const int& n=64) { fPhiBins = n; }
    /// @brief Set phi range
    void setJetPhiRange(const double& lo=-TMath::Pi(), const double& hi=TMath::Pi()) { fPhiRange[0]=lo; fPhiRange[1]=hi; }
    /// @brief Set number of bins for ptHat
    void setPtHatBins(const int& n=10) { fPtHatBins = n; }
    /// @brief Set range of ptHat
    void setPtHatRange(const double& lo=15., const double& hi=215.) { fPtHatRange[0]=lo; fPtHatRange[1]=hi; }

    //
    // Event histograms
    //
    TH1D *hVz;
    TH1D *hVzWeighted;
    TH1D *hPtHat;
    TH1D *hPtHatWeighted;
    TH1D *hHiBin;
    TH1D *hHiBinWeighted;

    TH1D *hVzGenDijetLab;
    TH1D *hVzGenDijetLabWeighted;
    TH1D *hHiBinGenDijetLab;
    TH1D *hHiBinGenDijetLabWeighted;

    TH1D *hVzGenDijetCM;
    TH1D *hVzGenDijetCMWeighted;
    TH1D *hHiBinGenDijetCM;
    TH1D *hHiBinGenDijetCMWeighted;

    TH1D *hVzRecoDijetLab;
    TH1D *hVzRecoDijetLabWeighted;
    TH1D *hHiBinRecoDijetLab;
    TH1D *hHiBinRecoDijetLabWeighted;

    TH1D *hVzRecoDijetCM;
    TH1D *hVzRecoDijetCMWeighted;
    TH1D *hHiBinRecoDijetCM;
    TH1D *hHiBinRecoDijetCMWeighted;

    TH1D *hVzRefSelDijetLab;
    TH1D *hVzRefSelDijetLabWeighted;
    TH1D *hHiBinRefSelDijetLab;
    TH1D *hHiBinRefSelDijetLabWeighted;

    TH1D *hVzRefSelDijetCM;
    TH1D *hVzRefSelDijetCMWeighted;
    TH1D *hHiBinRefSelDijetCM;
    TH1D *hHiBinRefSelDijetCMWeighted;


    //
    // Gen jet histograms
    //

    TH1D *hGenJetCollectionSize;
    TH2D *hGenVsRecoJetCollectionSize;

    TH2D *hGenLeadingJetPtOverPtHatVsLeadingJetPt;
    TH2D *hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted;
    TH2D *hGenDijetPtOverPtHatVsDijetPt;
    TH2D *hGenDijetPtOverPtHatVsDijetPtWeighted;
    TH2D *hGenDijetPtAveOverPtHatVsDijetPtAve;
    TH2D *hGenDijetPtAveOverPtHatVsDijetPtAveWeighted;

    // Dijet pt, dijet eta, dijet dphi, lead pt, lead eta, lead phi, 
    // sublead pt, sublead eta, sublead phi [9]
    THnSparseD *hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    // Dijet pt, dijet eta, dijet dphi, lead pt, lead eta, lead phi, 
    // sublead pt, sublead eta, sublead phi weighted [9]
    THnSparseD *hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    TH1D *hGenInclusiveJetPt;
    TH2D *hGenInclusiveJetPtEta;
    TH3D *hGenInclusiveJetPtEtaPtHat;
    TH3D *hGenLeadJetPtEtaPtHat;
    TH3D *hGenSubLeadJetPtEtaPtHat;
    TH2D *hGenPtLeadPtSublead;
    TH2D *hGenEtaLeadEtaSublead;
    TH2D *hGenEtaCMLeadEtaCMSublead;
    TH2D *hGenPtLeadPtSubleadMcReweight;
    TH2D *hGenEtaLeadEtaSubleadMcReweight;
    TH1D *hGenDijetEta;
    TH3D *hGenDijetPtEtaDphi;
    TH3D *hGenDijetPtEtaDphiWeighted;
    TH1D *hGenDijetEtaCM;
    TH3D *hGenDijetPtEtaDphiCM;
    TH3D *hGenDijetPtEtaDphiCMWeighted;
    TH2D *hGenDijetPtEtaForward;
    TH2D *hGenDijetPtEtaBackward;
    TH2D *hGenDijetPtEtaCMForward;
    TH2D *hGenDijetPtEtaCMBackward;
    TH2D *hGenDijetPtEtaForwardWeighted;
    TH2D *hGenDijetPtEtaBackwardWeighted;
    TH2D *hGenDijetPtEtaCMForwardWeighted;
    TH2D *hGenDijetPtEtaCMBackwardWeighted;
    TH1D *hGenGoodInclusiveJetEtaLabFrame;
    TH1D *hGenGoodInclusiveJetEtaCMFrame;

    TH1D *hGenInclusiveDijetDetaCM;
    TH1D *hGenInclusiveDijetDetaCMWeighted;
    TH2D *hGenInclusiveDijetDetaCMPt;
    TH2D *hGenInclusiveDijetDetaCMPtWeighted;
    TH3D *hGenInclusiveDijetEtaDetaCMPt;
    TH3D *hGenInclusiveDijetEtaDetaCMPtWeighted;
    TH1D *hGenInclusiveDijetXPb;
    TH1D *hGenInclusiveDijetXPbWeighted;
    TH1D *hGenInclusiveDijetXp;
    TH1D *hGenInclusiveDijetXpWeighted;
    TH1D *hGenInclusiveDijetXPbOverXp;
    TH1D *hGenInclusiveDijetXPbOverXpWeighted;
    TH2D *hGenInclusiveDijetXPbOverXpEta;
    TH2D *hGenInclusiveDijetXPbOverXpEtaWeighted;

    TH1D *hGenSelectedDijetDetaCM;
    TH1D *hGenSelectedDijetDetaCMWeighted;
    TH2D *hGenSelectedDijetDetaCMPt;
    TH2D *hGenSelectedDijetDetaCMPtWeighted;
    TH3D *hGenSelectedDijetEtaDetaCMPt;
    TH3D *hGenSelectedDijetEtaDetaCMPtWeighted;
    TH1D *hGenSelectedDijetXPb;
    TH1D *hGenSelectedDijetXPbWeighted;
    TH1D *hGenSelectedDijetXp;
    TH1D *hGenSelectedDijetXpWeighted;
    TH1D *hGenSelectedDijetXPbOverXp;
    TH1D *hGenSelectedDijetXPbOverXpWeighted;
    TH2D *hGenSelectedDijetXPbOverXpEta;
    TH2D *hGenSelectedDijetXPbOverXpEtaWeighted;

    // New ptAve and eta binning
    TH1D *hGenDijetEta1D[16];
    TH1D *hGenDijetEta1DWeighted[16];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2D[16];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DWeighted[16];
    TH1D *hGenDijetEtaForward1D[16];
    TH1D *hGenDijetEtaForward1DWeighted[16];
    TH1D *hGenDijetEtaBackward1D[16];
    TH1D *hGenDijetEtaBackward1DWeighted[16];

    TH1D *hGenDijetEta1DCM[16];
    TH1D *hGenDijetEta1DCMWeighted[16];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DCM[16];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[16];
    TH1D *hGenDijetEtaCMForward1D[16];
    TH1D *hGenDijetEtaCMForward1DWeighted[16];
    TH1D *hGenDijetEtaCMBackward1D[16];
    TH1D *hGenDijetEtaCMBackward1DWeighted[16];

    // Old ptAve and new eta binning
    TH1D *hGenDijetEta1DOldPt[6];
    TH1D *hGenDijetEta1DOldPtWeighted[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPt[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[6];
    TH1D *hGenDijetEtaForward1DOldPt[6];
    TH1D *hGenDijetEtaForward1DOldPtWeighted[6];
    TH1D *hGenDijetEtaBackward1DOldPt[6];
    TH1D *hGenDijetEtaBackward1DOldPtWeighted[6];

    TH1D *hGenDijetEta1DOldPtCM[6];
    TH1D *hGenDijetEta1DOldPtCMWeighted[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[6];
    TH1D *hGenDijetEtaCMForward1DOldPt[6];
    TH1D *hGenDijetEtaCMForward1DOldPtWeighted[6];
    TH1D *hGenDijetEtaCMBackward1DOldPt[6];
    TH1D *hGenDijetEtaCMBackward1DOldPtWeighted[6];

    // Old ptAve and old eta binning
    TH1D *hGenDijetEta1DOldPtBinning[6];
    TH1D *hGenDijetEta1DOldPtBinningWeighted[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[6];
    TH1D *hGenDijetEtaForward1DOldPtBinning[6];
    TH1D *hGenDijetEtaForward1DOldPtBinningWeighted[6];
    TH1D *hGenDijetEtaBackward1DOldPtBinning[6];
    TH1D *hGenDijetEtaBackward1DOldPtBinningWeighted[6];

    TH1D *hGenDijetEta1DOldPtBinningCM[6];
    TH1D *hGenDijetEta1DOldPtBinningCMWeighted[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[6];
    TH2D *hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[6];
    TH1D *hGenDijetEtaCMForward1DOldPtBinning[6];
    TH1D *hGenDijetEtaCMForward1DOldPtBinningWeighted[6];
    TH1D *hGenDijetEtaCMBackward1DOldPtBinning[6];
    TH1D *hGenDijetEtaCMBackward1DOldPtBinningWeighted[6];


    //
    // Reco jet histograms
    //

    TH1D *hRecoJetCollectionSize;

    // Inclusive jets
    TH1D *hRecoInclusiveJetNHF[4];
    TH1D *hRecoInclusiveJetNEmF[4];
    TH1D *hRecoInclusiveJetNumOfConst[4];
    TH1D *hRecoInclusiveJetMUF[4];
    TH1D *hRecoInclusiveJetCHF[4];
    TH1D *hRecoInclusiveJetChargedMult[4];
    TH1D *hRecoInclusiveJetCEmF[4];
    TH1D *hRecoInclusiveJetNumOfNeutPart[4];

    // Reconstructed jet acceptance
    TH2D *hRecoInclusiveAllJetPtRawEta;

    // Reco dijet pt, reco dijet eta, reco dijet dphi,
    // Reco lead pt, reco lead eta, reco lead phi,
    // Reco sublead pt, reco sublead eta, reco sublead phi [9]
    THnSparseD *hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    // Reco dijet pt, reco dijet eta, reco dijet dphi,
    // Reco lead pt, reco lead eta, reco lead phi,
    // Reco sublead pt, reco sublead eta, reco sublead phi weighted [9]
    THnSparseD *hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;


    TH1D *hRecoInclusiveAllJetPt;
    TH2D *hRecoInclusiveAllJetPtEta;
    TH3D *hRecoInclusiveAllJetPtEtaPtHat;
    TH3D *hRecoInclusiveMatchedJetPtEtaPtHat;
    TH3D *hRecoInclusiveUnmatchedJetPtEtaPtHat;

    TH2D *hRecoPtLeadPtSublead;
    TH2D *hRecoEtaLeadEtaSublead;
    TH2D *hRecoEtaCMLeadEtaCMSublead;
    TH2D *hRecoPtLeadPtSubleadMcReweight;
    TH2D *hRecoEtaLeadEtaSubleadMcReweight;
    TH2D *hRecoDijetPtEta;
    TH2D *hRecoDijetPtEtaForward;
    TH2D *hRecoDijetPtEtaBackward;
    TH2D *hRecoDijetPtEtaCMForward;
    TH2D *hRecoDijetPtEtaCMBackward;
    TH2D *hRecoDijetPtEtaForwardWeighted;
    TH2D *hRecoDijetPtEtaBackwardWeighted;
    TH2D *hRecoDijetPtEtaCMForwardWeighted;
    TH2D *hRecoDijetPtEtaCMBackwardWeighted;

    TH3D *hRecoDijetPtEtaDphi;
    TH3D *hRecoDijetPtEtaDphiWeighted;
    TH3D *hRecoDijetPtEtaDphiCM;
    TH3D *hRecoDijetPtEtaDphiCMWeighted;

    // Leading and subleading jet acceptance
    TH2D *hRecoLeadAllJetPtEta;
    TH3D *hRecoLeadAllJetPtEtaPtHat; 
    TH2D *hRecoSubLeadAllJetPtEta;
    TH3D *hRecoSubLeadAllJetPtEtaPtHat;

    TH1D *hRecoGoodInclusiveJetEtaLabFrame;
    TH1D *hRecoGoodInclusiveJetEtaCMFrame;

    // ptAve-integrated dijet pseudorapidity distributions
    TH1D *hRecoDijetEta;
    TH1D *hRecoDijetEtaCM;

    // New ptAve and eta binning
    TH1D *hRecoDijetEta1D[16];
    TH1D *hRecoDijetEta1DWeighted[16];
    TH2D *hRecoDijetEtaLeadVsEtaSubLead2D[16];
    TH2D *hRecoDijetEtaLeadVsEtaSubLead2DWeighted[16];
    TH1D *hRecoDijetEtaForward1D[16];
    TH1D *hRecoDijetEtaForward1DWeighted[16];
    TH1D *hRecoDijetEtaBackward1D[16];
    TH1D *hRecoDijetEtaBackward1DWeighted[16];

    TH1D *hRecoDijetEta1DCM[16];
    TH1D *hRecoDijetEta1DCMWeighted[16];
    TH2D *hRecoEtaLeadVsEtaSubLead2DCM[16];
    TH2D *hRecoEtaLeadVsEtaSubLead2DCMWeighted[16];
    TH1D *hRecoDijetEtaCMForward1D[16];
    TH1D *hRecoDijetEtaCMForward1DWeighted[16];
    TH1D *hRecoDijetEtaCMBackward1D[16];
    TH1D *hRecoDijetEtaCMBackward1DWeighted[16];

    // Old ptAve and new eta binning
    TH1D *hRecoDijetEta1DOldPt[6];
    TH1D *hRecoDijetEta1DOldPtWeighted[6];
    TH2D *hRecoDijetEtaLeadVsEtaSubLead2DOldPt[6];
    TH2D *hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[6];
    TH1D *hRecoDijetEtaForward1DOldPt[6];
    TH1D *hRecoDijetEtaForward1DOldPtWeighted[6];
    TH1D *hRecoDijetEtaBackward1DOldPt[6];
    TH1D *hRecoDijetEtaBackward1DOldPtWeighted[6];

    TH1D *hRecoDijetEta1DOldPtCM[6];
    TH1D *hRecoDijetEta1DOldPtCMWeighted[6];
    TH2D *hRecoEtaLeadVsEtaSubLead2DOldPtCM[6];
    TH2D *hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[6];
    TH1D *hRecoDijetEtaCMForward1DOldPt[6];
    TH1D *hRecoDijetEtaCMForward1DOldPtWeighted[6];
    TH1D *hRecoDijetEtaCMBackward1DOldPt[6];
    TH1D *hRecoDijetEtaCMBackward1DOldPtWeighted[6];

    // Old ptAve and old eta binning
    TH1D *hRecoDijetEta1DOldPtBinning[6];
    TH1D *hRecoDijetEta1DOldPtBinningWeighted[6];
    TH2D *hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[6];
    TH2D *hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[6];
    TH1D *hRecoDijetEtaForward1DOldPtBinning[6];
    TH1D *hRecoDijetEtaForward1DOldPtBinningWeighted[6];
    TH1D *hRecoDijetEtaBackward1DOldPtBinning[6];
    TH1D *hRecoDijetEtaBackward1DOldPtBinningWeighted[6];

    TH1D *hRecoDijetEta1DOldPtBinningCM[6];
    TH1D *hRecoDijetEta1DOldPtBinningCMWeighted[6];
    TH2D *hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[6];
    TH2D *hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[6];
    TH1D *hRecoDijetEtaCMForward1DOldPtBinning[6];
    TH1D *hRecoDijetEtaCMForward1DOldPtBinningWeighted[6];
    TH1D *hRecoDijetEtaCMBackward1DOldPtBinning[6];
    TH1D *hRecoDijetEtaCMBackward1DOldPtBinningWeighted[6];
  
    //
    // Ref jet histograms
    //

    // Inclusive jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen [6]
    THnSparseD *hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
    // Inclusive jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen, weighted [6]
    THnSparseD *hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
    // Leading jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen [6]
    THnSparseD *hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
    // Leading jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen, weighted [6]
    THnSparseD *hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
    // Subleading jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen [6]
    THnSparseD *hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
    // Subleading jet pt corr, pt raw, 
    // pt ref, eta corr, eta gen, weighted [6]
    THnSparseD *hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
    // pt corr / pt gen, gen pt, gen eta, reco phi [4]
    THnSparseD *hJESInclusiveJetPtEtaPhi;
    // pt corr / pt gen, gen pt, gen eta, reco phi [4]
    THnSparseD *hJESInclusiveJetPtEtaPhiWeighted;

    TH2D *hRecoLeadingJetPtOverPtHatVsLeadingJetPt;
    TH2D *hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted;
    TH2D *hRecoDijetPtOverPtHatVsDijetPt;
    TH2D *hRecoDijetPtOverPtHatVsDijetPtWeighted;
    TH2D *hRecoDijetPtAveOverPtHatVsDijetPtAve;
    TH2D *hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted;

    // pt corr / pt raw, gen pt, gen eta
    TH3D *hRecoInclusiveJetJECFactorVsPtEta;
    // pt corr (my) / pt raw, gen pt, gen eta
    TH3D *hRecoInclusiveJetJEC2FactorVsPtEta;
    // pt raw / ref pt, gen pt, gen eta
    TH3D *hRecoInclusiveJetPtRawOverPtRefVsPtEta;
    // pt raw / ref pt, gen pt, gen eta
    TH3D *hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning;
    // pt raw / ref pt, raw pt, reco eta
    TH3D *hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning;

    // pt corr/pt gen vs pt gen at midrapidity -1.4 < eta < 1.4 
    TH2D *hInclusiveJetJESVsPtGen;
    // pt corr/pt, gen pt, gen eta, ptHat
    THnSparseD *hInclusiveJetJESGenPtGenEtaPtHatWeighted;
    // pt corr/pt, reco pt, reco eta, ptHat
    THnSparseD *hInclusiveJetJESRecoPtRecoEtaPtHatWeighted;

    // pt corr/pt gen vs pt gen, gen eta, ptHat
    THnSparseD *hLeadingJetJESGenPtEtaPtHatWeighted;
    // pt corr/pt gen vs pt gen, gen eta, ptHat
    THnSparseD *hSubleadingJetJESGenPtEtaPtHatWeighted;

    // Matched and unmatched jet acceptance
    TH1D *hRecoInclusiveMatchedJetPt;
    TH2D *hRecoInclusiveMatchedJetPtEta;
    TH2D *hRecoInclusiveUnmatchedJetPtEta;

    TH2D *hRecoLeadMatchedJetPtEta;
    TH3D *hRecoLeadMatchedJetPtEtaPtHat;
    TH2D *hRecoLeadUnmatchedJetPtEta;
    TH3D *hRecoLeadUnmatchedJetPtEtaPtHat;

    TH2D *hRecoSubLeadMatchedJetPtEta;
    TH3D *hRecoSubLeadMatchedJetPtEtaPtHat;
    TH2D *hRecoSubLeadUnmatchedJetPtEta;
    TH3D *hRecoSubLeadUnmatchedJetPtEtaPtHat;

    TH1D *hRefInclusiveJetPt;
    TH2D *hRefInclusiveJetPtEta;
    TH3D *hRefInclusiveJetPtEtaPtHat;
    TH3D *hRefLeadJetPtEtaPtHat;
    TH3D *hRefSubLeadJetPtEtaPtHat;

    // Reco dijet pt, dijet eta, 
    // Reco lead pt, lead eta,
    // Reco sublead pt, sublead eta,
    // Ref dijet pt, dijet eta,
    // Ref lead pt, lead eta,
    // Ref sublead pt, sublead eta [12]
    THnSparseD *hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta;
    // Reco dijet pt, dijet eta, 
    // Reco lead pt, lead eta,
    // Reco sublead pt, sublead eta,
    // Ref dijet pt, dijet eta,
    // Ref lead pt, lead eta,
    // Ref sublead pt, sublead eta weighted [12]
    THnSparseD *hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;

    // Reco dijet pt, eta
    // Ref dijet pt, eta
    THnSparseD *hRecoDijetPtEtaRefDijetPtEta;
    THnSparseD *hRecoDijetPtEtaRefDijetPtEtaWeighted;

    // With cut on ref for the selection of dijets !!
    // Reco dijet pt, dijet eta, 
    // Reco lead pt, lead eta,
    // Reco sublead pt, sublead eta,
    // Ref dijet pt, dijet eta,
    // Ref lead pt, lead eta,
    // Ref sublead pt, sublead eta weighted [12]
    THnSparseD *hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;

    TH1D *hRefDijetEta;
    TH2D *hRefDijetEtaVsRecoDijetEta;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted;
    TH3D *hRefDijetPtEtaDphi;
    TH3D *hRefDijetPtEtaDphiWeighted;

    TH2D *hRefDijetPtEtaForward;
    TH2D *hRefDijetPtEtaBackward;
    TH2D *hRefDijetPtEtaCMForward;
    TH2D *hRefDijetPtEtaCMBackward;
    TH2D *hRefDijetPtEtaForwardWeighted;
    TH2D *hRefDijetPtEtaBackwardWeighted;
    TH2D *hRefDijetPtEtaCMForwardWeighted;
    TH2D *hRefDijetPtEtaCMBackwardWeighted;

    TH1D *hRefDijetEtaCM;
    TH3D *hRefDijetPtEtaDphiCM;
    TH3D *hRefDijetPtEtaDphiCMWeighted;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted;

    TH2D *hRefPtLeadPtSublead;
    TH2D *hRefEtaLeadEtaSublead;
    TH2D *hRefEtaCMLeadEtaCMSublead;
    TH2D *hRefPtLeadPtSubleadMcReweight;
    TH2D *hRefEtaLeadEtaSubleadMcReweight;

    // New ptAve and eta binning
    TH1D *hRefDijetEta1D[16];
    TH1D *hRefDijetEta1DWeighted[16];
    TH2D *hRefEtaLeadVsEtaSubLead2D[16];
    TH2D *hRefEtaLeadVsEtaSubLead2DWeighted[16];
    TH2D *hRecoVsRefDijetEta2D[16];
    TH2D *hRecoVsRefDijetEta2DWeighted[16];
    TH2D *hRecoVsRefLeadJetEta2D[16];
    TH2D *hRecoVsRefLeadJetEta2DWeighted[16];
    TH2D *hRecoVsRefSubLeadJetEta2D[16];
    TH2D *hRecoVsRefSubLeadJetEta2DWeighted[16];
    TH1D *hRefDijetEtaForward1D[16];
    TH1D *hRefDijetEtaForward1DWeighted[16];
    TH1D *hRefDijetEtaBackward1D[16];
    TH1D *hRefDijetEtaBackward1DWeighted[16];

    TH1D *hRefDijetEta1DCM[16];
    TH1D *hRefDijetEta1DCMWeighted[16];
    TH2D *hRefEtaLeadVsEtaSubLead2DCM[16];
    TH2D *hRefEtaLeadVsEtaSubLead2DCMWeighted[16];
    TH2D *hRecoVsRefDijetEta2DCM[16];
    TH2D *hRecoVsRefDijetEta2DCMWeighted[16];
    TH2D *hRecoVsRefLeadJetEta2DCM[16];
    TH2D *hRecoVsRefLeadJetEta2DCMWeighted[16];
    TH2D *hRecoVsRefSubLeadJetEta2DCM[16];
    TH2D *hRecoVsRefSubLeadJetEta2DCMWeighted[16];
    TH1D *hRefDijetEtaCMForward1D[16];
    TH1D *hRefDijetEtaCMForward1DWeighted[16];
    TH1D *hRefDijetEtaCMBackward1D[16];
    TH1D *hRefDijetEtaCMBackward1DWeighted[16];

    // Old ptAve and new eta binning
    TH1D *hRefDijetEta1DOldPt[6];
    TH1D *hRefDijetEta1DOldPtWeighted[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPt[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtWeighted[6];
    TH2D *hRecoVsRefDijetEta2DOldPt[6];
    TH2D *hRecoVsRefDijetEta2DOldPtWeighted[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPt[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtWeighted[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPt[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtWeighted[6];
    TH1D *hRefDijetEtaForward1DOldPt[6];
    TH1D *hRefDijetEtaForward1DOldPtWeighted[6];
    TH1D *hRefDijetEtaBackward1DOldPt[6];
    TH1D *hRefDijetEtaBackward1DOldPtWeighted[6];

    TH1D *hRefDijetEta1DOldPtCM[6];
    TH1D *hRefDijetEta1DOldPtCMWeighted[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtCM[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[6];
    TH2D *hRecoVsRefDijetEta2DOldPtCM[6];
    TH2D *hRecoVsRefDijetEta2DOldPtCMWeighted[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtCM[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtCMWeighted[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtCM[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[6];
    TH1D *hRefDijetEtaCMForward1DOldPt[6];
    TH1D *hRefDijetEtaCMForward1DOldPtWeighted[6];
    TH1D *hRefDijetEtaCMBackward1DOldPt[6];
    TH1D *hRefDijetEtaCMBackward1DOldPtWeighted[6];

    // Old ptAve and old eta binning
    TH1D *hRefDijetEta1DOldPtBinning[6];
    TH1D *hRefDijetEta1DOldPtBinningWeighted[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtBinning[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtBinningWeighted[6];
    TH2D *hRecoVsRefDijetEta2DOldPtBinning[6];
    TH2D *hRecoVsRefDijetEta2DOldPtBinningWeighted[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtBinning[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtBinningWeighted[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtBinning[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtBinningWeighted[6];
    TH1D *hRefDijetEtaForward1DOldPtBinning[6];
    TH1D *hRefDijetEtaForward1DOldPtBinningWeighted[6];
    TH1D *hRefDijetEtaBackward1DOldPtBinning[6];
    TH1D *hRefDijetEtaBackward1DOldPtBinningWeighted[6];

    TH1D *hRefDijetEta1DOldPtBinningCM[6];
    TH1D *hRefDijetEta1DOldPtBinningCMWeighted[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[6];
    TH2D *hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[6];
    TH2D *hRecoVsRefDijetEta2DOldPtBinningCM[6];
    TH2D *hRecoVsRefDijetEta2DOldPtBinningCMWeighted[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtBinningCM[6];
    TH2D *hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtBinningCM[6];
    TH2D *hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[6];
    TH1D *hRefDijetEtaCMForward1DOldPtBinning[6];
    TH1D *hRefDijetEtaCMForward1DOldPtBinningWeighted[6];
    TH1D *hRefDijetEtaCMBackward1DOldPtBinning[6];
    TH1D *hRefDijetEtaCMBackward1DOldPtBinningWeighted[6];

    //
    // Ref-selected jet histograms
    //

    TH1D *hRefSelInclusiveJetPt;
    TH2D *hRefSelInclusiveJetPtEta;
    TH3D *hRefSelInclusiveJetPtEtaPtHat;
    TH3D *hRefSelLeadJetPtEtaPtHat;
    TH3D *hRefSelSubLeadJetPtEtaPtHat;

    TH1D *hRefSelDijetEta;
    TH3D *hRefSelDijetPtEtaDphi;
    TH3D *hRefSelDijetPtEtaDphiWeighted;
    TH1D *hRefSelDijetEtaCM;
    TH3D *hRefSelDijetPtEtaDphiCM;
    TH3D *hRefSelDijetPtEtaDphiCMWeighted;

    // New ptAve and eta binning
    TH1D *hRefSelDijetEta1D[16];
    TH1D *hRefSelDijetEta1DWeighted[16];
    TH1D *hRefSelRecoDijetEta1D[16];
    TH1D *hRefSelRecoDijetEta1DWeighted[16];
    TH2D *hRefSelEtaLeadVsEtaSubLead2D[16];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DWeighted[16];
    TH1D *hRefSelDijetEtaForward1D[16];
    TH1D *hRefSelDijetEtaForward1DWeighted[16];
    TH1D *hRefSelDijetEtaBackward1D[16];
    TH1D *hRefSelDijetEtaBackward1DWeighted[16];

    TH1D *hRefSelDijetEta1DCM[16];
    TH1D *hRefSelDijetEta1DCMWeighted[16];
    TH1D *hRefSelRecoDijetEta1DCM[16];
    TH1D *hRefSelRecoDijetEta1DCMWeighted[16];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DCM[16];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DCMWeighted[16];
    TH1D *hRefSelDijetEtaCMForward1D[16];
    TH1D *hRefSelDijetEtaCMForward1DWeighted[16];
    TH1D *hRefSelDijetEtaCMBackward1D[16];
    TH1D *hRefSelDijetEtaCMBackward1DWeighted[16];

    // Old ptAve and new eta binning
    TH1D *hRefSelDijetEta1DOldPt[6];
    TH1D *hRefSelDijetEta1DOldPtWeighted[6];
    TH1D *hRefSelRecoDijetEta1DOldPt[6];
    TH1D *hRefSelRecoDijetEta1DOldPtWeighted[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPt[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[6];
    TH1D *hRefSelDijetEtaForward1DOldPt[6];
    TH1D *hRefSelDijetEtaForward1DOldPtWeighted[6];
    TH1D *hRefSelDijetEtaBackward1DOldPt[6];
    TH1D *hRefSelDijetEtaBackward1DOldPtWeighted[6];

    TH1D *hRefSelDijetEta1DOldPtCM[6];
    TH1D *hRefSelDijetEta1DOldPtCMWeighted[6];
    TH1D *hRefSelRecoDijetEta1DOldPtCM[6];
    TH1D *hRefSelRecoDijetEta1DOldPtCMWeighted[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtCM[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[6];
    TH1D *hRefSelDijetEtaCMForward1DOldPt[6];
    TH1D *hRefSelDijetEtaCMForward1DOldPtWeighted[6];
    TH1D *hRefSelDijetEtaCMBackward1DOldPt[6];
    TH1D *hRefSelDijetEtaCMBackward1DOldPtWeighted[6];

    // Old ptAve and old eta binning
    TH1D *hRefSelDijetEta1DOldPtBinning[6];
    TH1D *hRefSelDijetEta1DOldPtBinningWeighted[6];
    TH1D *hRefSelRecoDijetEta1DOldPtBinning[6];
    TH1D *hRefSelRecoDijetEta1DOldPtBinningWeighted[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[6];
    TH1D *hRefSelDijetEtaForward1DOldPtBinning[6];
    TH1D *hRefSelDijetEtaForward1DOldPtBinningWeighted[6];
    TH1D *hRefSelDijetEtaBackward1DOldPtBinning[6];
    TH1D *hRefSelDijetEtaBackward1DOldPtBinningWeighted[6];

    TH1D *hRefSelDijetEta1DOldPtBinningCM[6];
    TH1D *hRefSelDijetEta1DOldPtBinningCMWeighted[6];
    TH1D *hRefSelRecoDijetEta1DOldPtBinningCM[6];
    TH1D *hRefSelRecoDijetEta1DOldPtBinningCMWeighted[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[6];
    TH2D *hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[6];
    TH1D *hRefSelDijetEtaCMForward1DOldPtBinning[6];
    TH1D *hRefSelDijetEtaCMForward1DOldPtBinningWeighted[6];
    TH1D *hRefSelDijetEtaCMBackward1DOldPtBinning[6];
    TH1D *hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[6];

  private:

    /// @brief Is Monte Carlo
    bool fIsMc;
    // /// @brief  Number of centrality bins
    // int    fCentBins;
    // /// @brief Centrality range
    // double fCentRange[2];
    /// @brief Number of pT bins
    int    fPtBins;
    /// @brief Transverse momentum range
    double fPtRange[2];
    /// @brief Number of pseudorapidity bins
    int    fEtaBins;
    /// @brief Pseudorapidity range
    double fEtaRange[2];
    /// @brief Number of azimuthal angle bins
    int    fPhiBins;
    /// @brief Azimuthal angle range
    double fPhiRange[2];
    /// @brief Number of dijet pT bins
    int    fDijetPtBins;
    /// @brief Dijet pT range
    double fDijetPtRange[2];
    /// @brief Number of dijet eta bins
    int    fDijetEtaBins;
    /// @brief Dijet eta range
    double fDijetEtaRange[2];
    /// @brief Number of dijet delta phi bins
    int    fDijetDphiBins;
    /// @brief Dijet delta phi range
    double fDijetDphiRange[2];
    /// @brief Number of ptHat bins
    int    fPtHatBins;
    /// @brief PtHat range
    double fPtHatRange[2];
    /// @brief Jet type: PF or Calo
    TString  fJetType;
    int    fFracBins;
    double fFracRange[2];
    int   fMultBins;
    double fMultRange[2];

    /// @brief Values for new dijet ptAve binning
    std::vector<double> fPtAveBins;
    /// @brief Values for old dijet ptAve binning
    std::vector<double> fPtAveOldBins;

    ClassDef(HistoManagerDiJet, 0)
};

#endif // #define HistoManagerDiJet_h
