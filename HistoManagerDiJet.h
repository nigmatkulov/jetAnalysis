/**
 * @file HistoManagerDiJet.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for dijet studies
 * @version 0.1
 * @date 2024-01-10
 * 
 * @copyright Copyright (c) 2023
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

//________________
class HistoManagerDiJet : public BaseHistoManager {
  public:
    /// @brief Constructor
    HistoManagerDiJet();
    /// @brief Destructor
    virtual ~HistoManagerDiJet();

    /// @brief Initialize and create histograms
    void init(const Bool_t& isMc = kFALSE);
    /// @brief Use MC histograms
    void setIsMc(const Bool_t& isMc = kTRUE) { fIsMc = isMc; }
    /// @brief Write all objects to the output file
    void writeOutput();

    // /// @brief Set number of centrality bins
    // void setCentBins(const Int_t& n = 6) { fCentBins = n; }
    // /// @brief Set centrality range (in percentage)
    // void setCentRange(const Double_t& lo = 0, const Double_t& hi = 60) { fCentRange[0]=lo; fCentRange[1]=hi; }

    /// @brief Set number of jet pT bins
    void setJetPtBins(const Int_t& n = 200) { fPtBins = n; }
    /// @brief Set jet pT range
    void setJetPtRange(const Double_t& lo=0., const Double_t& hi=1000.) { fPtRange[0]=lo; fPtRange[1]=hi; }
    /// @brief Set number of eta bins
    void setJetEtaBins(const Int_t& n=50) { fEtaBins = n; }
    /// @brief Set eta range
    void setJetEtaRange(const Double_t& lo=-2.5, const Double_t& hi=2.5) { fEtaRange[0]=lo; fEtaRange[1]=hi; }
    /// @brief Set number of phi bins
    void setJetPhiBins(const Int_t& n=64) { fPhiBins = n; }
    /// @brief Set phi range
    void setJetPhiRange(const Double_t& lo=-TMath::Pi(), const Double_t& hi=TMath::Pi()) { fPhiRange[0]=lo; fPhiRange[1]=hi; }
    /// @brief Set number of bins for ptHat
    void setPtHatBins(const Int_t& n=10) { fPtHatBins = n; }
    /// @brief Set range of ptHat
    void setPtHatRange(const Double_t& lo=15., const Double_t& hi=215.) { fPtHatRange[0]=lo; fPtHatRange[1]=hi; }

    //
    // Event histograms
    //
    TH1D *hVz;
    TH1D *hVzWeighted;
    TH1D *hMult;
    TH1D *hHiBin;
    TH1D *hHiBinWeighted;
    TH1D *hPtHat;
    TH1D *hPtHatWeighted;
    TH1D *hPtHatWeight;
    // TH1D *hCentrality;
    // TH1D *hCentralityWeighted;
    THnSparseD *hVzPtHat;
    THnSparseD *hVzPtHatWeighted;

    TH1D *hNHF[4];
    TH1D *hNEmF[4];
    TH1D *hNumOfConst[4];
    TH1D *hMUF[4];
    TH1D *hCHF[4];
    TH1D *hChargedMult[4];
    TH1D *hCEmF[4];
    TH1D *hNumOfNeutPart[4];

    //
    // Gen jet histograms
    //

    // Dijet pt, dijet eta, dijet dphi, lead pt, lead eta, lead phi, 
    // sublead pt, sublead eta, sublead phi [9]
    THnSparseD *hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    // Dijet pt, dijet eta, dijet dphi, lead pt, lead eta, lead phi, 
    // sublead pt, sublead eta, sublead phi weighted [9]
    THnSparseD *hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    TH1D *hGenInclusiveJetPt;
    TH2D *hGenInclusiveJetPtEta;
    TH2D *hGenPtLeadPtSublead;
    TH2D *hGenEtaLeadEtaSublead;
    TH1D *hGenDijetEta;
    TH3D *hGenDijetPtEtaDphi;

    //
    // Reco jet histograms
    //


    // Single jets

    // Inclusive jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen [5]
    THnSparseD *hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen;
    // Inclusive jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen, weighted [5]
    THnSparseD *hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
    // Leading jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen [5]
    THnSparseD *hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
    // Leading jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen, weighted [5]
    THnSparseD *hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
    // Subleading jet pt corr, pt raw, pt ref, 
    // eta corr, eta gen [5]
    THnSparseD *hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen;
    // Subleading jet pt corr, pt raw, 
    // pt ref, eta corr, eta gen, weighted [5]
    THnSparseD *hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted;
    // pt corr / pt gen, gen pt, reco eta, reco phi [4]
    THnSparseD *hJESInclusiveJetPtEtaPhi;
    // pt corr / pt gen, gen pt, reco eta, reco phi [4]
    THnSparseD *hJESInclusiveJetPtEtaPhiWeighted;
    TH2D *hRecoMatchedPtEta;

    // Dijets

    // Experiment

    // Reco dijet pt, reco dijet eta, reco dijet dphi,
    // Reco lead pt, reco lead eta, reco lead phi,
    // Reco sublead pt, reco sublead eta, reco sublead phi [9]
    THnSparseD *hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi;
    // Reco dijet pt, reco dijet eta, reco dijet dphi,
    // Reco lead pt, reco lead eta, reco lead phi,
    // Reco sublead pt, reco sublead eta, reco sublead phi weighted [9]
    THnSparseD *hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted;
    TH1D *hRecoInclusiveJetPt;
    TH2D *hRecoPtLeadPtSublead;
    TH2D *hRecoEtaLeadEtaSublead;
    TH3D *hRecoDijetPtEtaDphi;

    TH3D *hRecoDijetPtEtaDphiJetId;

    TH2D *hRecoInclusiveAllJetPtVsEta;
    TH2D *hRecoInclusiveMatchedJetPtVsEta;
    TH2D *hRecoInclusiveUnmatchedJetPtVsEta;

    // Jet selection algo check
    TH2D       *hRecoInclusiveJetPtVsEtaKineCut;
    THnSparseD *hRecoInclusiveJetJESPtEtaPhiKineCut;
    TH3D       *hRecoInclusiveJetDEtaPtEtaKineCut;
    TH2D       *hRecoInclusiveMatchedJetPtVsEtaKineCut;
    TH2D       *hRecoInclusiveUnmatchedJetPtVsEtaKineCut;
    TH2D       *hRecoInclusiveJetRefPtVsEtaKineCut;

    TH2D       *hRecoInclusiveJetPtVsEtaTrkMaxCut;
    THnSparseD *hRecoInclusiveJetJESPtEtaPhiTrkMaxCut;
    TH3D       *hRecoInclusiveJetDEtaPtEtaTrkMaxCut;
    TH2D       *hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut;
    TH2D       *hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut;
    TH2D       *hRecoInclusiveJetRefPtVsEtaTrkMaxCut;

    TH2D       *hRecoInclusiveJetPtVsEtaJetIdCut;
    THnSparseD *hRecoInclusiveJetJESPtEtaPhiJetIdCut;
    TH3D       *hRecoInclusiveJetDEtaPtEtaJetIdCut;
    TH2D       *hRecoInclusiveMatchedJetPtVsEtaJetIdCut;
    TH2D       *hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut;
    TH2D       *hRecoInclusiveJetRefPtVsEtaJetIdCut;

    TH1D       *hRecoTrkMaxToJetIdDijetMatching;
    
    // MC

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

    // With cut on ref for the selection of dijets !!
    // Reco dijet pt, dijet eta, 
    // Reco lead pt, lead eta,
    // Reco sublead pt, sublead eta,
    // Ref dijet pt, dijet eta,
    // Ref lead pt, lead eta,
    // Ref sublead pt, sublead eta weighted [12]
    THnSparseD *hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted;


    TH1D *hRecoDijetEta;
    TH1D *hRefDijetEta;
    TH2D *hRefDijetEtaVsRecoDijetEta;
    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt;
    TH3D *hRefDijetPtEtaDphi;
    TH1D *hRefSelDijetEta;
    TH3D *hRefSelDijetPtEtaDphi;

    TH3D *hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId;
    TH3D *hRefDijetPtEtaDphiJetId;

    TH1D *hRefInclusiveJetPt;
    TH2D *hRefInclusiveJetPtEta;
    TH2D *hRefPtLeadPtSublead;
    TH2D *hRefEtaLeadEtaSublead;

  private:

    /// @brief Is Monte Carlo
    Bool_t fIsMc;
    // /// @brief  Number of centrality bins
    // Int_t    fCentBins;
    // /// @brief Centrality range
    // Double_t fCentRange[2];
    /// @brief Number of pT bins
    Int_t    fPtBins;
    /// @brief Transverse momentum range
    Double_t fPtRange[2];
    /// @brief Number of pseudorapidity bins
    Int_t    fEtaBins;
    /// @brief Pseudorapidity range
    Double_t fEtaRange[2];
    /// @brief Number of azimuthal angle bins
    Int_t    fPhiBins;
    /// @brief Azimuthal angle range
    Double_t fPhiRange[2];
    /// @brief Number of dijet pT bins
    Int_t    fDijetPtBins;
    /// @brief Dijet pT range
    Double_t fDijetPtRange[2];
    /// @brief Number of dijet eta bins
    Int_t    fDijetEtaBins;
    /// @brief Dijet eta range
    Double_t fDijetEtaRange[2];
    /// @brief Number of dijet delta phi bins
    Int_t    fDijetDphiBins;
    /// @brief Dijet delta phi range
    Double_t fDijetDphiRange[2];
    /// @brief Number of ptHat bins
    Int_t    fPtHatBins;
    /// @brief PtHat range
    Double_t fPtHatRange[2];
    /// @brief Jet type: PF or Calo
    TString  fJetType;
    Int_t    fFracBins;
    Double_t fFracRange[2];
    Int_t   fMultBins;
    Double_t fMultRange[2];

    ClassDef(HistoManagerDiJet, 0)
};

#endif // #define HistoManagerDiJet_h