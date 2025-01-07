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
#include "TProfile.h"
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
    void init(const bool& isMc = false);
    /// @brief Write all objects to the output file
    void writeOutput();

    /// @brief Set verbose mode
    void setVerbose(const bool& v = false) { fVerbose = v; }

    /// @brief Set number of centrality bins
    void setCentBins(const int& n = 6) { fCentBins = n; }
    /// @brief Set centrality range (in percentage)
    void setCentRange(const double& lo = 0, const double& hi = 60) { fCentRange[0]=lo; fCentRange[1]=hi; }

    /// @brief Set number of jet pT bins
    void setJetPtBins(const int& n = 200) { fJetPtBins = n; }
    /// @brief Set jet pT range
    void setJetPtRange(const double& lo=15., const double& hi=1015.) { fJetPtRange[0]=lo; fJetPtRange[1]=hi; }
    /// @brief Set number of eta bins
    void setJetEtaBins(const int& n=50) { fJetEtaBins = n; }
    /// @brief Set eta range
    void setJetEtaRange(const double& lo=-2.5, const double& hi=2.5) { fJetEtaRange[0]=lo; fJetEtaRange[1]=hi; }
    /// @brief Set number of phi bins
    void setJetPhiBins(const int& n=64) { fJetPhiBins = n; }
    /// @brief Set phi range
    void setJetPhiRange(const double& lo=-TMath::Pi(), const double& hi=TMath::Pi()) { fJetPhiRange[0]=lo; fJetPhiRange[1]=hi; }
    /// @brief Set number of JES bins
    void setJESBins(const int& n=500) { fJESBins = n; }
    /// @brief Set JES range
    void setJESRange(const double& lo=0., const double& hi=5.) { fJESRange[0]=lo; fJESRange[1]=hi; }
    /// @brief Set number of flavorForB bins
    void setFlavorForBBins(const int& n=14) { fJetFlavorForBBins = n; }
    /// @brief Set range for flavorForB
    void setFlavorForBRange(const double& lo=-6.5, const double& hi=6.5) { fJetFlavorForBRange[0]=lo; fJetFlavorForBRange[1]=hi; }
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

    TH1D *hHiBin;
    TH1D *hHiBinPtHatWeighted;
    TH1D *hHiBinWeighted;

    TH1D *hPtHat;
    TH1D *hPtHatPtHatWeighted;
    TH1D *hPtHatCentWeighted;
    TH1D *hPtHatWeighted;

    TH1D *hPtHatWeight;
    TH1D *hPtHatWeightWeighted;

    TH1D *hCentrality;
    TH1D *hCentralityPtHatWeighted;
    TH1D *hCentralityWeighted;

    TH3D *hVzPtHatCent;
    TH3D *hVzPtHatCentPtHatWeighted;
    TH3D *hVzPtHatCentWeighted;

    //
    // Gen jet histograms
    //

    // Inclusive jets
    TH1D *hGenInclusiveJetPt;
    TH1D *hGenInclusiveJetPtWeighted;
    TH2D *hGenInclusiveJetEtaPt;
    TH2D *hGenInclusiveJetEtaPtWeighted;
    THnSparseD *hGenInclusiveJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted;

    // Leading and subleading jets
    TH1D *hGenInclusiveLeadJetPt;
    TH1D *hGenInclusiveLeadJetPtWeighted;
    TH2D *hGenInclusiveLeadJetEtaPt;
    TH2D *hGenInclusiveLeadJetEtaPtWeighted;
    THnSparseD *hGenInclusiveLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hGenInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted;

    TH1D *hGenInclusiveSubLeadJetPt;
    TH1D *hGenInclusiveSubLeadJetPtWeighted;
    TH2D *hGenInclusiveSubLeadJetEtaPt;
    TH2D *hGenInclusiveSubLeadJetEtaPtWeighted;
    THnSparseD *hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    // Inclusive dijets
    TH1D *hGenInclusiveDijetDphi;
    TH1D *hGenInclusiveDijetDphiWeighted;
    TH2D *hGenInclusiveDijetEtaPt;
    TH2D *hGenInclusiveDijetEtaPtWeighted;
    TH3D *hGenInclusiveDijetEtaPtDphi;
    TH3D *hGenInclusiveDijetEtaPtDphiWeighted;
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

    // Selected leading and subleading jets
    TH1D *hGenSelectedLeadJetPt;
    TH1D *hGenSelectedLeadJetPtWeighted;
    TH2D *hGenSelectedLeadJetEtaPt;
    TH2D *hGenSelectedLeadJetEtaPtWeighted;
    THnSparseD *hGenSelectedLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted;

    TH1D *hGenSelectedSubLeadJetPt;
    TH1D *hGenSelectedSubLeadJetPtWeighted;
    TH2D *hGenSelectedSubLeadJetEtaPt;
    TH2D *hGenSelectedSubLeadJetEtaPtWeighted;
    THnSparseD *hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    // Selected dijets
    TH1D *hGenSelectedDijetPt;
    TH1D *hGenSelectedDijetPtWeighted;
    TH1D *hGenSelectedDijetEta;
    TH1D *hGenSelectedDijetEtaWeighted;
    TH2D *hGenSelectedDijetEtaPt;
    TH2D *hGenSelectedDijetEtaPtWeighted;
    TH3D *hGenSelectedDijetEtaPtDphi;
    TH3D *hGenSelectedDijetEtaPtDphiWeighted;
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

    // Inclusive jets
    TH1D *hRecoInclusiveJetNHF[4];
    TH1D *hRecoInclusiveJetNEmF[4];
    TH1D *hRecoInclusiveJetNumOfConst[4];
    TH1D *hRecoInclusiveJetMUF[4];
    TH1D *hRecoInclusiveJetCHF[4];
    TH1D *hRecoInclusiveJetChargedMult[4];
    TH1D *hRecoInclusiveJetCEmF[4];
    TH1D *hRecoInclusiveJetNumOfNeutPart[4];
    // 0 - pt, 1 - eta, 2 - phi, 3 - dummyIter, 4 - NHF, 5 - NEF, 6 - num of constituents, 
    // 7 - MUF, 8 - CHF, 9 - charged mult, 10 - CEF, 11 - neutral mult, 12 - ptHat, 13 - centrality
    THnSparseD *hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent;

    TH1D *hRecoInclusiveJetPt;
    TH1D *hRecoInclusiveJetPtWeighted;
    TH2D *hRecoInclusiveJetEtaPt;
    TH2D *hRecoInclusiveJetEtaPtWeighted;
    TH3D *hRecoInclusiveJetPtRawPtCorrEta;
    TH3D *hRecoInclusiveJetPtRawPtCorrEtaWeighted;

    // Inclusive unmatched jets
    TH1D *hRecoInclusiveUnmatchedJetPt;
    TH1D *hRecoInclusiveUnmatchedJetPtWeighted;
    TH2D *hRecoInclusiveUnmatchedJetEtaPt;
    TH2D *hRecoInclusiveUnmatchedJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted;

    // Inclusive matched jets
    TH1D *hRecoInclusiveMatchedJetPt;
    TH1D *hRecoInclusiveMatchedJetPtWeighted;
    TH2D *hRecoInclusiveMatchedJetEtaPt;
    TH2D *hRecoInclusiveMatchedJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted;
    // 0 - JES, 1 - pt, 2 - eta, 3 - phi, 4 - centrality
    THnSparseD *hRecoInclusiveJetJESPtEtaPhiCent;
    THnSparseD *hRecoInclusiveJetJESPtEtaPhiCentWeighted;

    // Correlation between reco and ref for inclusive jets
    TH2D *hReco2RefInclusiveJetPt;
    TH2D *hReco2RefInclusiveJetPtWeighted;
    TH2D *hReco2RefInclusiveJetEta;
    TH2D *hReco2RefInclusiveJetEtaWeighted;
    TH2D *hReco2RefInclusiveJetPhi;
    TH2D *hReco2RefInclusiveJetPhiWeighted;
    THnSparseD *hReco2RefInclusiveJetPtEtaPhiPtHatCentrality;
    THnSparseD *hReco2RefInclusiveJetPtEtaPhiPtHatCentralityWeighted;

    // Inclusive leading jets
    TH1D *hRecoInclusiveLeadJetPt;
    TH1D *hRecoInclusiveLeadJetPtWeighted;
    TH2D *hRecoInclusiveLeadJetEtaPt;
    TH2D *hRecoInclusiveLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted;
    THnSparseD *hRecoInclusiveLeadJetJESPtEtaPhiCent;
    THnSparseD *hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted;

    // Correlation between inclusive reco and ref for leading jets
    TH2D *hReco2RefInclusiveLeadJetPt;
    TH2D *hReco2RefInclusiveLeadJetPtWeighted;
    TH2D *hReco2RefInclusiveLeadJetEta;
    TH2D *hReco2RefInclusiveLeadJetEtaWeighted;
    TH2D *hReco2RefInclusiveLeadJetPhi;
    TH2D *hReco2RefInclusiveLeadJetPhiWeighted;
    // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
    THnSparseD *hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality;
    THnSparseD *hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentralityWeighted;

    // Inclusive subleading jets
    TH1D *hRecoInclusiveSubLeadJetPt;
    TH1D *hRecoInclusiveSubLeadJetPtWeighted;
    TH2D *hRecoInclusiveSubLeadJetEtaPt;
    TH2D *hRecoInclusiveSubLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted;
    THnSparseD *hRecoInclusiveSubLeadJetJESPtEtaPhiCent;
    THnSparseD *hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted;

    // Correlation between inclusive reco and ref for subleading jets
    TH2D *hReco2RefInclusiveSubLeadJetPt;
    TH2D *hReco2RefInclusiveSubLeadJetPtWeighted;
    TH2D *hReco2RefInclusiveSubLeadJetEta;
    TH2D *hReco2RefInclusiveSubLeadJetEtaWeighted;
    TH2D *hReco2RefInclusiveSubLeadJetPhi;
    TH2D *hReco2RefInclusiveSubLeadJetPhiWeighted;
    // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
    THnSparseD *hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality;
    THnSparseD *hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted;

    // Inclusive dijets
    TH1D *hRecoInclusiveDijetDphi;
    TH1D *hRecoInclusiveDijetDphiWeighted;
    TH2D *hRecoInclusiveDijetEtaPt;
    TH2D *hRecoInclusiveDijetEtaPtWeighted;
    TH3D *hRecoInclusiveDijetEtaPtDphi;
    TH3D *hRecoInclusiveDijetEtaPtDphiWeighted;
    TH1D *hRecoInclusiveDijetDetaCM;
    TH1D *hRecoInclusiveDijetDetaCMWeighted;
    TH2D *hRecoInclusiveDijetDetaCMPt;
    TH2D *hRecoInclusiveDijetDetaCMPtWeighted;
    TH3D *hRecoInclusiveDijetEtaDetaCMPt;
    TH3D *hRecoInclusiveDijetEtaDetaCMPtWeighted;
    THnSparseD *hRecoInclusiveDijetJESPtEtaDphiCent;
    THnSparseD *hRecoInclusiveDijetJESPtEtaDphiCentWeighted;

    // Selected leading and subleading jets
    TH1D *hRecoSelectedLeadJetPt;
    TH1D *hRecoSelectedLeadJetPtWeighted;
    TH2D *hRecoSelectedLeadJetEtaPt;
    TH2D *hRecoSelectedLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted;

    TH1D *hRecoSelectedSubLeadJetPt;
    TH1D *hRecoSelectedSubLeadJetPtWeighted;
    TH2D *hRecoSelectedSubLeadJetEtaPt;
    TH2D *hRecoSelectedSubLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    // Selected dijets
    TH1D *hRecoSelectedDijetPt;
    TH1D *hRecoSelectedDijetPtWeighted;
    TH1D *hRecoSelectedDijetEta;
    TH1D *hRecoSelectedDijetEtaWeighted;
    TH2D *hRecoSelectedDijetEtaPt;
    TH2D *hRecoSelectedDijetEtaPtWeighted;
    TH3D *hRecoSelectedDijetEtaPtDphi;
    TH3D *hRecoSelectedDijetEtaPtDphiWeighted;
    TH1D *hRecoSelectedDijetDetaCM;
    TH1D *hRecoSelectedDijetDetaCMWeighted;
    TH2D *hRecoSelectedDijetDetaCMPt;
    TH2D *hRecoSelectedDijetDetaCMPtWeighted;
    TH3D *hRecoSelectedDijetEtaDetaCMPt;
    TH3D *hRecoSelectedDijetEtaDetaCMPtWeighted;

    // Correlation between reco and ref for selected dijets
    // 0 - reco dijetpt, 1 - reco dijeteta, 
    // 2 - reco leading jet pt, 3 - reco leading jet eta,
    // 4 - reco subleading jet pt, 5 - reco subleading jet eta, 
    // 6 - ref dijet pt, 7 - ref dijet eta,
    // 8 - ref leading jet pt, 9 - ref leading jet eta,
    // 10 - ref subleading jet pt, 11 - ref subleading jet eta
    THnSparseD *hReco2RefSelectedDijetPtEtaFull;
    THnSparseD *hReco2RefSelectedDijetPtEtaFullWeighted;


    //
    // Ref jet histograms
    //

    // Inclusive jets
    TH1D *hRefInclusiveJetPt;
    TH1D *hRefInclusiveJetPtWeighted;
    TH2D *hRefInclusiveJetEtaPt;
    TH2D *hRefInclusiveJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRefInclusiveJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted;

    // Inclusive leading jets
    TH1D *hRefInclusiveLeadJetPt;
    TH1D *hRefInclusiveLeadJetPtWeighted;
    TH2D *hRefInclusiveLeadJetEtaPt;
    TH2D *hRefInclusiveLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted;

    // Inclusive subleading jets
    TH1D *hRefInclusiveSubLeadJetPt;
    TH1D *hRefInclusiveSubLeadJetPtWeighted;
    TH2D *hRefInclusiveSubLeadJetEtaPt;
    TH2D *hRefInclusiveSubLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality 
    THnSparseD *hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    // Inclusive dijets
    TH1D *hRefInclusiveDijetDphi;
    TH1D *hRefInclusiveDijetDphiWeighted;
    TH1D *hRefInclusiveDijetPt;
    TH1D *hRefInclusiveDijetPtWeighted;
    TH1D *hRefInclusiveDijetEta;
    TH1D *hRefInclusiveDijetEtaWeighted;
    TH2D *hRefInclusiveDijetEtaPt;
    TH2D *hRefInclusiveDijetEtaPtWeighted;
    TH3D *hRefInclusiveDijetEtaPtDphi;
    TH3D *hRefInclusiveDijetEtaPtDphiWeighted;

    // Selected leading and subleading jets
    TH1D *hRefSelectedLeadJetPt;
    TH1D *hRefSelectedLeadJetPtWeighted;
    TH2D *hRefSelectedLeadJetEtaPt;
    TH2D *hRefSelectedLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRefSelectedLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted;

    TH1D *hRefSelectedSubLeadJetPt;
    TH1D *hRefSelectedSubLeadJetPtWeighted;
    TH2D *hRefSelectedSubLeadJetEtaPt;
    TH2D *hRefSelectedSubLeadJetEtaPtWeighted;
    // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
    THnSparseD *hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent;
    THnSparseD *hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted;

    // Selected dijets
    TH1D *hRefSelectedDijetPt;
    TH1D *hRefSelectedDijetPtWeighted;
    TH1D *hRefSelectedDijetEta;
    TH1D *hRefSelectedDijetEtaWeighted;
    TH2D *hRefSelectedDijetEtaPt;
    TH2D *hRefSelectedDijetEtaPtWeighted;
    TH3D *hRefSelectedDijetEtaPtDphi;
    TH3D *hRefSelectedDijetEtaPtDphiWeighted;

  private:

    bool   fVerbose;
    int    fPtHatBins;
    double fPtHatRange[2];
    int    fCentBins;
    double fCentRange[2];
    int    fJetPtBins;
    double fJetPtRange[2];
    int    fJetEtaBins;
    double fJetEtaRange[2];
    int    fJetPhiBins;
    double fJetPhiRange[2];
    int    fJetFlavorForBBins;
    double fJetFlavorForBRange[2];
    int    fJESBins;
    double fJESRange[2];
    int    fDijetPtBins;
    double fDijetPtRange[2];
    int    fDijetEtaBins;
    double fDijetEtaRange[2];
    int    fDijetDphiBins;
    double fDijetDphiRange[2];

    ClassDef(HistoManagerJetESR, 0)
};

#endif // #define HistoManagerJetESR_h