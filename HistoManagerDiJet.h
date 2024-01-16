/**
 * @file HistoManagerDiJetR.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for JES and JER studies
 * @version 0.1
 * @date 2024-01-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef HistoManagerDiJetR_h
#define HistoManagerDiJetR_h

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
class HistoManagerDiJetR : public BaseHistoManager {
  public:
    /// @brief Constructor
    HistoManagerDiJetR();
    /// @brief Destructor
    virtual ~HistoManagerDiJetR();

    /// @brief Initialize and create histograms
    void init(const Bool_t& isMc = kFALSE);
    /// @brief Use MC histograms
    void setIsMc(const Bool_t& isMc = kTRUE) { fIsMc = isMc; }
    /// @brief Write all objects to the output file
    void writeOutput();

    /// @brief Set number of centrality bins
    void setCentBins(const Int_t& n = 6) { fCentBins = n; }
    /// @brief Set centrality range (in percentage)
    void setCentRange(const Double_t& lo = 0, const Double_t& hi = 60) { fCentRange[0]=lo; fCentRange[1]=hi; }

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
    TH1D *hCentrality;
    TH1D *hCentralityWeighted;
    TH1D *hNBadJets[5];   // pThat >0, >20, >40, >60, >80
    THnSparseD *hVzPtHatCent;
    THnSparseD *hVzPtHatCentWeighted;

    //
    // Gen jet histograms
    //
    THnSparseD *hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent;
    THnSparseD *hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted;
    THnSparseD *hGenDiJetPtEtaDeltaPhiCent;
    THnSparseD *hGenDiJetPtEtaDeltaPhiCentWeighted;

  private:

    /// @brief Is Monte Carlo
    Bool_t fIsMc;

    /// @brief  Number of centrality bins
    Int_t    fCentBins;
    /// @brief Centrality range
    Double_t fCentRange[2];
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
    /// @brief Number of ptHat bins
    Int_t    fPtHatBins;

    Double_t fPtHatRange[2];
    /// @brief Jet type: PF or Calo
    TString  fJetType;

    ClassDef(HistoManagerDiJetR, 0)
};

#endif // #define HistoManagerDiJetR_h