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
    void setJetPhiBins(const int& n=16) { fPhiBins = n; }
    /// @brief Set phi range
    void setJetPhiRange(const double& lo=-TMath::Pi(), const double& hi=TMath::Pi()) { fPhiRange[0]=lo; fPhiRange[1]=hi; }
    /// @brief Number of dijet pT bins
    void setDijetPtBins(const int& n = 196) { fDijetPtBins = n; }
    /// @brief Dijet pT range
    void setDijetPtRange(const double& lo = 20., const double& hi = 1000.) { fDijetPtRange[0] = lo; fDijetPtRange[1] = hi; }
    /// @brief Number of dijet eta bins
    void setDijetEtaBins(const int& n = 32) { fDijetEtaBins = n; }
    /// @brief Dijet eta range
    void setDijetEtaRange(const double& lo = -3.2, const double& hi = 3.2) { fDijetEtaRange[0] = lo; fDijetEtaRange[1] = hi; }
    /// @brief Number of dijet eta forward or backward bins
    void setDijetEtaFBBins(const int& n = 16) { fDijetEtaFBBins = n; }
    /// @brief Dijet eta forward or backward range
    void setDijetEtaFBRange(const double& lo = 0., const double& hi = 3.2) { fDijetEtaFBRange[0] = lo; fDijetEtaFBRange[1] = hi; }
    // /// @brief Number of dijet delta phi bins
    // void setDijetDphiBins(const int& n = 16) { fDijetDphiBins = n; }
    // /// @brief Dijet delta phi range
    // void setDijetDphiRange(const double& lo = -TMath::Pi(), const double& hi = TMath::Pi()) { fDijetDphiRange[0] = lo; fDijetDphiRange[1] = hi; }
    /// @brief Set number of bins for ptHat
    void setPtHatBins(const int& n=10) { fPtHatBins = n; }
    /// @brief Set range of ptHat
    void setPtHatRange(const double& lo=15., const double& hi=215.) { fPtHatRange[0]=lo; fPtHatRange[1]=hi; }
    /// @brief Set variable binning for histograms
    void setUseVariableBinning(const bool& useVariableBinning = true) { fUseVariableBinning = useVariableBinning; }

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

    TH2D *hGenLeadJetPtOverPtHatVsLeadJetPt;
    TH2D *hGenLeadJetPtOverPtHatVsLeadJetPtWeighted;
    TH2D *hGenDijetPtOverPtHatVsDijetPt;
    TH2D *hGenDijetPtOverPtHatVsDijetPtWeighted;
    TH2D *hGenDijetPtAveOverPtHatVsDijetPtAve;
    TH2D *hGenDijetPtAveOverPtHatVsDijetPtAveWeighted;

    // Gen dijet info [10 dimensions]
    // 0 - dijet pt ave, 1 - dijet eta lab, 2 - dijet eta cm, 3 - dijet delta phi,
    // 4 - lead pt, 5 - lead eta lab, 6 - lead eta cm,
    // 7 - sublead pt, 8 - sublead eta lab, 9 - sublead eta cm
    // THnSparseD *hGenDijetInfo;

    TH1D *hGenInclusiveJetPt;
    TH1D *hGenInclusiveJetEta;
    TH1D *hGenInclusiveJetEtaUnweighted;
    TH2D *hGenInclusiveJetPtEta;
    TH2D *hGenInclusiveJetPtEtaCM;
    TH3D *hGenInclusiveJetPtEtaPtHat;
    TH2D *hGenLeadJetPtEta;
    TH2D *hGenLeadJetPtEtaCM;
    TH3D *hGenLeadJetPtEtaPtHat;
    TH2D *hGenSubLeadJetPtEta;
    TH2D *hGenSubLeadJetPtEtaCM;
    TH3D *hGenSubLeadJetPtEtaPtHat;
    TH2D *hGenPtLeadPtSublead;
    TH2D *hGenEtaLeadEtaSublead;
    TH2D *hGenEtaCMLeadEtaCMSublead;
    TH2D *hGenPtLeadPtSubleadMcReweight;
    TH2D *hGenEtaLeadEtaSubleadMcReweight;

    TH1D *hGenDijetEta;
    TH2D *hGenDijetPtEta;
    TH2D *hGenDijetPtEtaWeighted;
    TH2D *hGenDijetPtEtaCMInLab;
    TH1D *hGenDijetEtaCM;
    TH2D *hGenDijetPtEtaCM;
    TH2D *hGenDijetPtEtaCMWeighted;
    TH2D *hGenDijetPtEtaLabInCM;
    TH2D *hGenDijetPtEtaForward;
    TH2D *hGenDijetPtEtaBackward;
    TH2D *hGenDijetPtEtaForwardCMInLab;
    TH2D *hGenDijetPtEtaBackwardCMInLab;
    TH2D *hGenDijetPtEtaCMForward;
    TH2D *hGenDijetPtEtaCMBackward;
    TH2D *hGenDijetPtEtaForwardLabInCM;
    TH2D *hGenDijetPtEtaBackwardLabInCM;
    TH2D *hGenDijetPtEtaForwardWeighted;
    TH2D *hGenDijetPtEtaBackwardWeighted;
    TH2D *hGenDijetPtEtaCMForwardWeighted;
    TH2D *hGenDijetPtEtaCMBackwardWeighted;

    TH3D *hGenDijetPtAveLeadPtSubLeadPt;
    TH3D *hGenDijetPtAveLeadPtSubLeadPtCM;
    TH3D *hGenDijetPtAveLeadEtaSubLeadEta;
    TH3D *hGenDijetPtAveLeadEtaSubLeadEtaCM;
    TH3D *hGenDijetEtaLeadEtaSubLeadEta;
    TH3D *hGenDijetEtaLeadEtaSubLeadEtaCM;

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

    // Reco dijet info [10 dimensions]
    // 0 - dijet pt ave, 1 - dijet eta lab, 2 - dijet eta cm, 3 - dijet delta phi
    // 4 - lead pt, 5 - lead eta, 6 - lead eta cm,
    // 7 - sublead pt, 8 - sublead eta, 9 - sublead eta cm
    // THnSparseD *hRecoDijetInfo;

    TH1D *hRecoInclusiveAllJetPt;
    TH1D *hRecoInclusiveAllJetEta;
    TH1D *hRecoInclusiveAllJetEtaUnweighted;
    TH2D *hRecoInclusiveAllJetPtEta;
    TH2D *hRecoInclusiveAllJetPtEtaCM;
    TH3D *hRecoInclusiveAllJetPtEtaPtHat;
    TH3D *hRecoInclusiveAllJetPtRawEtaPtHat;
    TH2D *hRecoInclusiveMatchedAllJetPtEta;
    TH3D *hRecoInclusiveMatchedJetPtEtaPtHat;
    TH2D *hRecoInclusiveUnmatchedAllJetPtEta;
    TH3D *hRecoInclusiveUnmatchedJetPtEtaPtHat;
    
    TH2D *hRecoPtLeadPtSublead;
    TH2D *hRecoEtaLeadEtaSublead;
    TH2D *hRecoEtaCMLeadEtaCMSublead;
    TH2D *hRecoPtLeadPtSubleadMcReweight;
    TH2D *hRecoEtaLeadEtaSubleadMcReweight;

    TH2D *hRecoDijetPtEtaForward;
    TH2D *hRecoDijetPtEtaBackward;
    TH2D *hRecoDijetPtEtaForwardCMInLab;
    TH2D *hRecoDijetPtEtaBackwardCMInLab;
    TH2D *hRecoDijetPtEtaCMForward;
    TH2D *hRecoDijetPtEtaCMBackward;
    TH2D *hRecoDijetPtEtaForwardLabInCM;
    TH2D *hRecoDijetPtEtaBackwardLabInCM;
    TH2D *hRecoDijetPtEtaForwardWeighted;
    TH2D *hRecoDijetPtEtaBackwardWeighted;
    TH2D *hRecoDijetPtEtaCMForwardWeighted;
    TH2D *hRecoDijetPtEtaCMBackwardWeighted;

    TH2D *hRecoDijetPtEta;
    TH2D *hRecoDijetPtEtaWeighted;
    TH2D *hRecoDijetPtEtaCMInLab;
    TH2D *hRecoDijetPtEtaCM;
    TH2D *hRecoDijetPtEtaCMWeighted;
    TH2D *hRecoDijetPtEtaLabInCM;
    TH2D *hRecoDijetPtEtaMatched;
    TH2D *hRecoDijetPtEtaCMMatched;

    TH3D *hRecoDijetPtAveLeadPtSubLeadPt;
    TH3D *hRecoDijetPtAveLeadPtSubLeadPtCM;
    TH3D *hRecoDijetPtAveLeadEtaSubLeadEta;
    TH3D *hRecoDijetPtAveLeadEtaSubLeadEtaCM;
    TH3D *hRecoDijetEtaLeadEtaSubLeadEta;
    TH3D *hRecoDijetEtaLeadEtaSubLeadEtaCM;

    // Lead and SubLead jet acceptance
    TH2D *hRecoLeadAllJetPtEta;
    TH2D *hRecoLeadAllJetPtEtaCM;
    TH3D *hRecoLeadAllJetPtEtaPtHat; 
    TH2D *hRecoSubLeadAllJetPtEta;
    TH2D *hRecoSubLeadAllJetPtEtaCM;
    TH3D *hRecoSubLeadAllJetPtEtaPtHat;

    TH1D *hRecoGoodInclusiveJetEtaLabFrame;
    TH1D *hRecoGoodInclusiveJetEtaCMFrame;

    // ptAve-integrated dijet pseudorapidity distributions
    TH1D *hRecoDijetEta;
    TH1D *hRecoDijetEtaCM;

 
    //
    // Ref jet histograms
    //

    // Inclusive jet correlation between reco and ref [5 dimensions]
    // 0 - corrected reco pt, 1 - raw reco pt, 2 - ref pt, 
    // 3 - reco eta, 4 - ref eta
    THnSparseD *hRecoInclusiveJetReco2Ref;
    // Leading jet correlation between reco and ref [5 dimensions]
    // 0 - corrected reco pt corr, 1 - raw reco pt, 2 - ref pt, 
    // 3 - reco eta, 4 - ref eta
    THnSparseD *hRecoLeadJetReco2Ref;
    // SubLeading jet correlation between reco and ref [5 dimensions]
    // 0 - corrected reco pt, 1 - raw reco pt, 2 - ref pt, 
    // 3 - reco eta, 4 - ref eta
    THnSparseD *hRecoSubLeadJetReco2Ref;
    // Jet energy scale (JES) vs ref parameters [4 dimensions]
    // 0 - corrected reco pt/gen pt, 1 - gen pt, 2 - gen eta, 3 - gen phi
    THnSparseD *hJESInclusiveJetPtEtaPhi;

    TH2D *hRecoLeadJetPtOverPtHatVsLeadJetPt;
    TH2D *hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted;
    TH2D *hRecoDijetPtOverPtHatVsDijetPt;
    TH2D *hRecoDijetPtOverPtHatVsDijetPtWeighted;
    TH2D *hRecoDijetPtAveOverPtHatVsDijetPtAve;
    TH2D *hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted;

    // // pt corr / pt raw, gen pt, gen eta
    // TH3D *hRecoInclusiveJetJECFactorVsPtEta;
    // // pt corr (my) / pt raw, gen pt, gen eta
    // TH3D *hRecoInclusiveJetJEC2FactorVsPtEta;
    // // pt raw / ref pt, gen pt, gen eta
    // TH3D *hRecoInclusiveJetPtRawOverPtRefVsPtEta;
    // // pt raw / ref pt, gen pt, gen eta
    // TH3D *hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning;
    // // pt raw / ref pt, raw pt, reco eta
    // TH3D *hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning;

    // pt corr/pt gen vs pt gen at midrapidity -1.4 < eta < 1.4 
    TH2D *hInclusiveJetJESVsPtGen;
    // pt corr/pt, gen pt, gen eta, ptHat
    THnSparseD *hInclusiveJetJESGenPtGenEtaPtHatWeighted;
    // pt corr/pt, gen pt, gen eta, ptHat
    THnSparseD *hInclusiveJetJESGenPtGenEtaCMPtHatWeighted;
    // pt corr/pt, reco pt, reco eta, ptHat
    THnSparseD *hInclusiveJetJESRecoPtRecoEtaPtHatWeighted;

    // pt corr/pt gen vs pt gen, gen eta, ptHat
    THnSparseD *hLeadJetJESGenPtEtaPtHatWeighted;
    THnSparseD *hLeadJetJESGenPtEtaCMPtHatWeighted;
    // pt corr/pt gen vs pt gen, gen eta, ptHat
    THnSparseD *hSubLeadJetJESGenPtEtaPtHatWeighted;
    THnSparseD *hSubLeadJetJESGenPtEtaCMPtHatWeighted;

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
    TH1D *hRefInclusiveJetEta;
    TH1D *hRefInclusiveJetEtaUnweighted;
    TH2D *hRefInclusiveJetPtEta;
    TH2D *hRefInclusiveJetPtEtaCM;
    TH3D *hRefInclusiveJetPtEtaPtHat;
    TH2D *hRefLeadJetPtEta;
    TH2D *hRefLeadJetPtEtaCM;
    TH3D *hRefLeadJetPtEtaPtHat;
    TH2D *hRefLeadUnswappedJetPtEta;
    TH3D *hRefLeadUnswappedJetPtEtaPtHat;
    TH2D *hRefSubLeadJetPtEta;
    TH2D *hRefSubLeadJetPtEtaCM;
    TH3D *hRefSubLeadJetPtEtaPtHat;
    TH2D *hRefSubLeadUnswappedJetPtEta;
    TH3D *hRefSubLeadUnswappedJetPtEtaPtHat;

    // Reco 2 Ref parameters [12 dimensions]
    // 0 - reco dijet ptAve, 1 - dijet eta,
    // 2 - reco lead pt, 3 - lead eta,
    // 4 - reco sublead pt, 5 - sublead eta,
    // 6 - ref dijet ptAve, 7 - dijet eta,
    // 8 - ref lead pt, 9 - lead eta,
    // 10 - ref sublead pt, 11 - sublead eta
    // THnSparseD *hReco2RefFull;

    // Reco dijet pt, eta
    // Ref dijet pt, eta
    THnSparseD *hRecoDijetPtEtaRefDijetPtEta;
    THnSparseD *hRecoDijetPtEtaRefDijetPtEtaWeighted;

    // Reco 2 Ref parameters (RefSelected) [12 dimensions]
    // 0 - reco dijet ptAve, 1 - dijet eta,
    // 2 - reco lead pt, 3 - lead eta,
    // 4 - reco sublead pt, 5 - sublead eta,
    // 6 - ref dijet ptAve, 7 - dijet eta,
    // 8 - ref lead pt, 9 - lead eta,
    // 10 - ref sublead pt, 11 - sublead eta
    // THnSparseD *hRefSel2RecoFull;

    TH1D *hRefDijetEta;
    TH2D *hRefDijetEtaVsRecoDijetEta;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted;
    TH2D *hRefDijetPtEta;
    TH2D *hRefDijetPtEtaWeighted;
    TH2D *hRefDijetPtEtaCMInLab;

    TH2D *hRefDijetPtEtaForward;
    TH2D *hRefDijetPtEtaBackward;
    TH2D *hRefDijetPtEtaForwardCMInLab;
    TH2D *hRefDijetPtEtaBackwardCMInLab;
    TH2D *hRefDijetPtEtaCMForward;
    TH2D *hRefDijetPtEtaCMBackward;
    TH2D *hRefDijetPtEtaForwardLabInCM;
    TH2D *hRefDijetPtEtaBackwardLabInCM;
    TH2D *hRefDijetPtEtaForwardWeighted;
    TH2D *hRefDijetPtEtaBackwardWeighted;
    TH2D *hRefDijetPtEtaCMForwardWeighted;
    TH2D *hRefDijetPtEtaCMBackwardWeighted;

    TH1D *hRefDijetEtaCM;
    TH2D *hRefDijetPtEtaCM;
    TH2D *hRefDijetPtEtaCMWeighted;
    TH2D *hRefDijetPtEtaLabInCM;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted;

    TH2D *hRefPtLeadPtSublead;
    TH2D *hRefEtaLeadEtaSublead;
    TH2D *hRefEtaCMLeadEtaCMSublead;
    TH2D *hRefPtLeadPtSubleadMcReweight;
    TH2D *hRefEtaLeadEtaSubleadMcReweight;

    //
    // Ref-selected jet histograms
    //

    TH1D *hRefSelInclusiveJetPt;
    TH1D *hRefSelInclusiveJetEta;
    TH1D *hRefSelInclusiveJetEtaUnweighted;
    TH2D *hRefSelInclusiveJetPtEta;
    TH3D *hRefSelInclusiveJetPtEtaPtHat;
    TH2D *hRefSelLeadJetPtEta;
    TH3D *hRefSelLeadJetPtEtaPtHat;
    TH2D *hRefSelSubLeadJetPtEta;
    TH3D *hRefSelSubLeadJetPtEtaPtHat;

    TH1D *hRefSelDijetEta;
    TH2D *hRefSelDijetPtEta;
    TH2D *hRefSelDijetPtEtaWeighted;
    TH2D *hRefSelDijetPtEtaCMInLab;
    TH1D *hRefSelDijetEtaCM;
    TH2D *hRefSelDijetPtEtaCM;
    TH2D *hRefSelDijetPtEtaCMWeighted;
    TH2D *hRefSelDijetPtEtaLabInCM;

  private:

    /// @brief Is Monte Carlo
    bool fIsMc;
    /// @brief Use variable binning (default is true)
    bool fUseVariableBinning;
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
    /// @brief Number of dijet eta bins in forward or backward directions
    int    fDijetEtaFBBins;
    /// @brief Dijet eta range in forward or backward directions (typically, 0., fDijetEtaRange[1])
    double fDijetEtaFBRange[2];
    // /// @brief Number of dijet delta phi bins
    // int    fDijetDphiBins;
    // /// @brief Dijet delta phi range
    // double fDijetDphiRange[2];
    /// @brief Number of ptHat bins
    int    fPtHatBins;
    /// @brief PtHat range
    double fPtHatRange[2];
    /// @brief Number of bins for fraction calculations
    int    fFracBins;
    /// @brief Fraction range for fraction calculations
    double fFracRange[2];
    /// @brief Number of multiplicity bins
    int   fMultBins;
    /// @brief Multiplicity range
    double fMultRange[2];

    // ClassDef(HistoManagerDiJet, 0)
};

#endif // #define HistoManagerDiJet_h
