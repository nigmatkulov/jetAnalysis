#ifndef BASICHISTOMANAGER_H
#define BASICHISTOMANAGER_H

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
class BasicHistoManager : public BaseHistoManager {
  public:
    /// @brief Constructor
    BasicHistoManager();
    /// @brief Destructor
    virtual ~BasicHistoManager();

    /// @brief Initialize and create histograms
    void init(const Bool_t& isMc = kFALSE);
    /// @brief Use MC histograms
    void setIsMc(const Bool_t& isMc = kTRUE) { fIsMc = isMc; }
    /// @brief Write all objects to the output file
    void writeOutput();

    /// @brief Set number of centrality bins
    void setCentBins(const Int_t& n = 6) { fCentBins = n; }
    /// @brief Set centrality range (in percentage)
    void setCentRange(const Double_t& lo = 0 , const Double_t& hi = 60) { fCentRange[0]=lo; fCentRange[1]=hi; }

    /// @brief Set number of jet pT bins
    void setJetPtBins(const Int_t& n = 200) { fJetPtBins = n; }
    /// @brief Set jet pT range
    void setJetPtRange(const Double_t& lo=0., const Double_t& hi=1000.) { fJetPtRange[0]=lo; fJetPtRange[1]=hi; }
    /// @brief Set number of eta bins
    void setJetEtaBins(const Int_t& n=50) { fJetEtaBins = n; }
    /// @brief Set eta range
    void setJetEtaRange(const Double_t& lo=-2.5, const Double_t& hi=2.5) { fJetEtaRange[0]=lo; fJetEtaRange[1]=hi; }
    /// @brief Set number of phi bins
    void setJetPhiBins(const Int_t& n=64) { fJetPhiBins = n; }
    /// @brief Set phi range
    void setJetPhiRange(const Double_t& lo=-TMath::Pi(), const Double_t& hi=TMath::Pi()) { fJetPhiRange[0]=lo; fJetPhiRange[1]=hi; }
    /// @brief Set number of JES bins
    void setJESBins(const Int_t& n=500) { fJESBins = n; }
    /// @brief Set JES range
    void setJESRange(const Double_t& lo=0., const Double_t& hi=5.) { fJESRange[0]=lo; fJESRange[1]=hi; }
    /// @brief Set number of JER bins
    void setJERBins(const Int_t& n=400) { fJERBins = n; }
    /// @brief Set JER range
    void setJERRange(const Double_t& lo=-2., const Double_t& hi=2.) { fJERRange[0]=lo; fJERRange[1]=hi; }


    //
    // Event histograms
    //
    TH1D *hVz;
    TH1D *hVzWeighted;
    TH1D *hMult;
    TH1D *hHiBin;
    TH1D *hHiBinWieghted;
    TH1D *hPtHat;
    TH1D *hPtHatWeighted;
    TH1D *hPtHatWeight;
    TH1D *hCentrality;

    //
    // Reco jet histograms
    //
    TH1D* hNRecoJets;
    TH1D* hRecoJetPtRaw;
    TH2F* hRecoJetPtCorrVsPtRaw;
    TH1D* hRecoJetPt; // after applying JEC
    TH1D* hRecoJetPtWeighted;
    TH1D* hRecoJetEta;
    TH1D* hRecoJetPhi;
    TH2F* hRecoJetPtVsEta;
    TH2F* hRecoJetPhiVsPt;
    TH2F* hRecoJetEtaVsPhi;

    THnSparseD *hRecoJet;
    THnSparseD *hRecoJetCorr;
    THnSparseD *hRecoJetCorrWeighted;
    
    //
    // Ref jet histograms
    //
    TH1D* hNRefJets;
    TH1D* hRefJetPt;
    TH1D* hRefJetPtWeighted;
    TH1D* hRefJetEta;
    TH1D* hRefJetPhi;
    TH2F* hRefJetPtVsEta;

    THnSparseD *hRefJet;
    THnSparseD *hRefJetWeighted;

    // Jet Energy Scale
    THnSparseD *hJESRaw;
    THnSparseD *hJESRawWeighted;
    THnSparseD *hJESReco;
    THnSparseD *hJESRecoWeighted;

    // Jet Energy Resolution
    THnSparseD *hJERRaw;
    THnSparseD *hJERRawWeighted;
    THnSparseD *hJERReco;
    THnSparseD *hJERRecoWeighted;

  private:

    /// @brief Is Monte Carlo
    Bool_t fIsMc;

    Int_t    fCentBins;
    Double_t fCentRange[2];
    Int_t    fJetPtBins;
    Double_t fJetPtRange[2];
    Int_t    fJetEtaBins;
    Double_t fJetEtaRange[2];
    Int_t    fJetPhiBins;
    Double_t fJetPhiRange[2];
    Int_t    fJESBins;
    Double_t fJESRange[2];
    Int_t    fJERBins;
    Double_t fJERRange[2];

    ClassDef(BasicHistoManager, 0)
};

#endif // #define BASICHISTOMANAGER_H