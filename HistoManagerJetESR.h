/**
 * @file HistoManagerJetESR.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for JES and JER studies
 * @version 1.01
 * @date 2025-01-06
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef HistoManagerJetESR_h
#define HistoManagerJetESR_h

// Jet analysis headers
#include "BaseHistoManager.h"

// ROOT headers
#include "TObject.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TMath.h"

//________________
class HistoManagerJetESR : public BaseHistoManager {
  public:
    /// @brief Constructor
    HistoManagerJetESR();
    /// @brief Destructor
    virtual ~HistoManagerJetESR();

    /// @brief Initialize and create histograms
    void init();
    /// @brief Write all objects to the output file
    void writeOutput();

    /// @brief Set verbose mode
    void setVerbose(const bool& v = false) { fVerbose = v; }
    /// @brief Set if is Monte Carlo
    void setIsMc(const bool& isMc = true)  { fIsMc = isMc; }

    /// @brief Set number of centrality bins
    void setCentBins(const int& n = 6) { fCentBins = n; }
    /// @brief Set centrality range (in percentage)
    void setCentRange(const double& lo = 0, const double& hi = 60) { fCentRange[0]=lo; fCentRange[1]=hi; }

    /// @brief Set number of jet pT bins
    void setJetPtBins(const int& n = 200) { fPtBins = n; }
    /// @brief Set jet pT range
    void setJetPtRange(const double& lo=15., const double& hi=1015.) { fPtRange[0]=lo; fPtRange[1]=hi; }
    /// @brief Set number of eta bins
    void setJetEtaBins(const int& n=50) { fEtaBins = n; }
    /// @brief Set eta range
    void setJetEtaRange(const double& lo=-2.5, const double& hi=2.5) { fEtaRange[0]=lo; fEtaRange[1]=hi; }
    /// @brief Set number of phi bins
    void setJetPhiBins(const int& n=64) { fPhiBins = n; }
    /// @brief Set phi range
    void setJetPhiRange(const double& lo=-TMath::Pi(), const double& hi=TMath::Pi()) { fPhiRange[0]=lo; fPhiRange[1]=hi; }
    /// @brief Set number of JES bins
    void setJESBins(const int& n=500) { fJESBins = n; }
    /// @brief Set JES range
    void setJESRange(const double& lo=0., const double& hi=5.) { fJESRange[0]=lo; fJESRange[1]=hi; }
    /// @brief Set number of bins for ptHat
    void setPtHatBins(const int& n=10) { fPtHatBins = n; }
    /// @brief Set range of ptHat
    void setPtHatRange(const double& lo=15., const double& hi=215.) { fPtHatRange[0]=lo; fPtHatRange[1]=hi; }

    //
    // Event histograms
    //
    TH1D *hVz;
    TH1D *hVzCentWeighted;
    TH1D *hVzPtHatWeighted;
    TH1D *hVzWeighted;
    TH1D *hPtHat;
    TH1D *hPtHatWeighted;
    TH1D *hHiBin;
    TH1D *hHiBinWeighted;
    TH1D *hPtHatWeight;
    TH1D *hPtHatWeightWeighted;
    TH1D *hCentrality;
    TH1D *hCentralityWeighted;
    TH3D *hVzPtHatCent;
    TH3D *hVzPtHatCentWeighted;

    //
    // Gen jet histograms
    //

    TH1D *hGenJetCollectionSize;
    TH2D *hGenVsRecoJetCollectionSize;

    // Histograms to check event overweighting
    TH2D *hGenLeadingJetPtOverPtHatVsLeadingJetPt;
    TH2D *hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted;
    TH2D *hGenDijetPtOverPtHatVsDijetPt;
    TH2D *hGenDijetPtOverPtHatVsDijetPtWeighted;
    TH2D *hGenDijetPtAveOverPtHatVsDijetPtAve;
    TH2D *hGenDijetPtAveOverPtHatVsDijetPtAveWeighted;

    // Single-jet distributions (inclusive, leading, subleading)
    TH1D *hGenInclusiveJetPt;
    TH1D *hGenInclusiveJetEta;
    TH1D *hGenInclusiveJetEtaUnweighted;
    TH2D *hGenInclusiveJetPtEta;
    TH3D *hGenInclusiveJetPtEtaPtHat;
    TH2D *hGenLeadJetPtEta;
    TH3D *hGenLeadJetPtEtaPtHat;
    TH2D *hGenSubLeadJetPtEta;
    TH3D *hGenSubLeadJetPtEtaPtHat;
    TH1D *hGenGoodInclusiveJetEtaLabFrame;
    TH1D *hGenGoodInclusiveJetEtaCMFrame;

    //
    // Reco jet histograms
    //

    TH1D *hRecoJetCollectionSize;

    TH2D *hRecoLeadingJetPtOverPtHatVsLeadingJetPt;
    TH2D *hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted;
    TH2D *hRecoDijetPtOverPtHatVsDijetPt;
    TH2D *hRecoDijetPtOverPtHatVsDijetPtWeighted;
    TH2D *hRecoDijetPtAveOverPtHatVsDijetPtAve;
    TH2D *hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted;

    // JetId histograms
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

    // Single-jet distributions (inclusive, leading, subleading)
    TH1D *hRecoInclusiveAllJetPt;
    TH1D *hRecoInclusiveAllJetEta;
    TH1D *hRecoInclusiveAllJetEtaUnweighted;
    TH2D *hRecoInclusiveAllJetPtEta;
    TH3D *hRecoInclusiveAllJetPtEtaPtHat;
    TH2D *hRecoLeadAllJetPtEta;
    TH3D *hRecoLeadAllJetPtEtaPtHat; 
    TH2D *hRecoSubLeadAllJetPtEta;
    TH3D *hRecoSubLeadAllJetPtEtaPtHat;
    TH1D *hRecoGoodInclusiveJetEtaLabFrame;
    TH1D *hRecoGoodInclusiveJetEtaCMFrame;

    //
    // Ref jet histograms
    //

    // Single-jet distributions (inclusive, leading, subleading)
    TH1D *hRefInclusiveJetPt;
    TH1D *hRefInclusiveJetEta;
    TH1D *hRefInclusiveJetEtaUnweighted;
    TH2D *hRefInclusiveJetPtEta;
    TH3D *hRefInclusiveJetPtEtaPtHat;
    TH2D *hRefLeadJetPtEta;
    TH3D *hRefLeadJetPtEtaPtHat;
    TH2D *hRefLeadUnswappedJetPtEta;
    TH3D *hRefLeadUnswappedJetPtEtaPtHat;
    TH2D *hRefSubLeadJetPtEta;
    TH3D *hRefSubLeadJetPtEtaPtHat;
    TH2D *hRefSubLeadUnswappedJetPtEta;
    TH3D *hRefSubLeadUnswappedJetPtEtaPtHat;

    // Matched and unmatched jet acceptance

    TH1D *hRecoInclusiveMatchedJetPt;
    TH2D *hRecoInclusiveMatchedJetPtEta;
    TH3D *hRecoInclusiveMatchedJetPtEtaPtHat;

    TH1D *hRecoInclusiveUnmatchedJetPt;
    TH2D *hRecoInclusiveUnmatchedJetPtEta;
    TH3D *hRecoInclusiveUnmatchedJetPtEtaPtHat;

    TH2D *hRecoLeadMatchedJetPtEta;
    TH3D *hRecoLeadMatchedJetPtEtaPtHat;
    TH2D *hRecoLeadUnmatchedJetPtEta;
    TH3D *hRecoLeadUnmatchedJetPtEtaPtHat;

    TH2D *hRecoSubLeadMatchedJetPtEta;
    TH3D *hRecoSubLeadMatchedJetPtEtaPtHat;
    TH2D *hRecoSubLeadUnmatchedJetPtEta;
    TH3D *hRecoSubLeadUnmatchedJetPtEtaPtHat;

    // JES and JER for inclusive, leading, and subleading jets

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

    // Reco 2 Ref correlations for inclusive jets
    // 0 - reco ptCorr, 1- reco eta, 2 - reco phi,
    // 3 - ref pt, 4 - ref eta, 5 - ref phi, 6 - ptHat
    THnSparseD *hInclusiveReco2RefJetPtEtaPhiPtHat;
    THnSparseD *hLeadReco2RefJetPtEtaPhiPtHat;
    THnSparseD *hSubLeadReco2RefJetPtEtaPhiPtHat;

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

  private:
    
    /// @brief Number of ptHat bins (for plotting)
    int    fPtHatBins;
    /// @brief Range of ptHat (for plotting)
    double fPtHatRange[2];
    /// @brief Number of centrality bins (for plotting)
    int    fCentBins;
    /// @brief Range of centrality (for plotting)
    double fCentRange[2];
    /// @brief Number of pT bins (for plotting)
    int    fPtBins;
    /// @brief Range of pT (for plotting)
    double fPtRange[2];
    /// @brief Number of eta bins (for plotting)
    int    fEtaBins;
    /// @brief Range of eta (for plotting)
    double fEtaRange[2];
    /// @brief Number of azimuthal angle bins (for plotting)
    int    fPhiBins;
    /// @brief Range of azimuthal angle (for plotting)
    double fPhiRange[2];
    /// @brief Number of bins for JES (for plotting)
    int    fJESBins;
    /// @brief Range of JES (for plotting)
    double fJESRange[2];
    /// @brief Verbose mode (default: false)
    bool   fVerbose;
    /// @brief Is this a Monte Carlo sample
    bool   fIsMc;

    ClassDef(HistoManagerJetESR, 0)
};

#endif // #define HistoManagerJetESR_h
