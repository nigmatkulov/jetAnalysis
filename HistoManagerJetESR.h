/**
 * @file HistoManagerJetESR.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for JES and JER studies
 * @version 0.1
 * @date 2023-10-24
 * 
 * @copyright Copyright (c) 2023
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
    /// @brief Set number of JES bins
    void setJESBins(const Int_t& n=500) { fJESBins = n; }
    /// @brief Set JES range
    void setJESRange(const Double_t& lo=0., const Double_t& hi=5.) { fJESRange[0]=lo; fJESRange[1]=hi; }
    /// @brief Set number of flavorForB bins
    void setFlavorForBBins(const Int_t& n=14) { fFlavorForBBins = n; }
    /// @brief Set range for flavorForB
    void setFlavorForBRange(const Double_t& lo=-6.5, const Double_t& hi=6.5) { fFlavorForBRange[0]=lo; fFlavorForBRange[1]=hi; }
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

    //
    // Gen jet histograms
    //
    TH1D *hNGenJets[5];   // for jets with pT cuts: >0, >20, >50, >80, >120 GeV
    THnSparseD *hGenJetPtEtaPhiCent;              // pt, eta, phi, centrality
    THnSparseD *hGenJetPtEtaPhiCentWeighted;      // pt, eta, phi, centrality
    THnSparseD *hGenJetPtFlavPtHatCent;         // pt, flavorForB, ptHat, centrality
    THnSparseD *hGenJetPtFlavPtHatCentWeighted; // pt, flavorForB, ptHat, centrality

    //
    // Reco jet histograms
    //
    TH1D *hNRecoJets[5];  // for jets with pT cuts: >0, >20, >50, >80, >120 GeV
    THnSparseD *hRecoJetRawPtEtaPhiCent;       // ptRaw, eta, phi, centrality
    THnSparseD *hRecoJetPtEtaPhiCent;          // ptCorr, eta, phi, centrality
    THnSparseD *hRecoJetPtEtaPhiCentWeighted;  // ptCorr, eta, phi, centrality
    THnSparseD *hRecoJetPtFlavPtHatCent;          // ptCorr, flavorForB, ptHat, centrality
    THnSparseD *hRecoJetPtFlavPtHatCentWeighted;  // ptCorr, flavorForB, ptHat, centrality
    THnSparseD *hRecoJetPtFlavPtHatCentInclusive; // (matched + unmatched) ptCorr, flavorForB, ptHat, centrality
    THnSparseD *hRecoJetPtFlavPtHatCentInclusiveWeighted; // (matched + unmatched) ptCorr, flavorForB, ptHat, centrality
    THnSparseD *hRecoJetDeltaRPtCent;           // deltaR=sqrt((eta-WTAeta)^2+(phi-WTAphi)^2),ptCorr,centrality
    THnSparseD *hRecoUmnatchedJetPtFlavPtHatCent; // ptCorr of unmatched jets, flavorForB, ptHat, centrality
    THnSparseD *hRecoUmnatchedJetPtFlavPtHatCentWeighted; // ptCorr of unmatched jets, flavorForB, ptHat, centrality
    THnSparseD *hRecoLeadJetPtFlavPtHatCent;    // Leading jet ptCorr, ptHat, centrality
    THnSparseD *hRecoLeadJetPtFlavPtHatCentWeighted;    // Leading jet ptCorr, ptHat, centrality

    THnSparseD *hRecoJetRawPtCorrPtGenPtCent;   // Reconstructed jet raw pT, corrected pT, gen pT, centrality
    
    //
    // Ref jet histograms
    //
    TH1D *hNRefJets[5];   // for jets with pT cuts: >0, >20, >50, >80, >120 GeV
    THnSparseD *hRefJetPtEtaPhiCent;               // pt, eta, phi, centrality
    THnSparseD *hRefJetPtEtaPhiCentWeighted;       // pt, eta, phi, centrality
    THnSparseD *hRefJetPtFlavPtHatCent;            // ptCorr, flavorForB, ptHat, centrality
    THnSparseD *hRefJetPtFlavPtHatCentWeighted;    // ptCorr, flavorForB, ptHat, centrality

    THnSparseD *hRecoJetPtRefJetPtPtHatCent;         // reco pt corr, gen pt, ptHat, centrality
    THnSparseD *hRecoJetPtRefJetPtPtHatCentWeighted; // reco pt corr, gen pt, ptHat, centrality

    //
    // Jet energy scale and resolution
    //
    THnSparseD *hJESPtEtaPhiCent;               // JES via pt corr, gen pt, flavorForB, ptHat, centrality
    THnSparseD *hJESPtEtaPhiCentWeighted;       // JES via pt corr, gen pt, flavorForB, ptHat, centrality
    THnSparseD *hJESRawPtFlavPtHatCent;         // JES via pt raw,  gen pt, flavorForB, ptHat, centrality
    THnSparseD *hJESRawPtFlavPtHatCentWeighted; // JES via pt raw,  gen pt, flavorForB, ptHat, centrality
    THnSparseD *hJESPtFlavPtHatCent;            // JES via pt corr, gen pt, flavorForB, ptHat, centrality
    THnSparseD *hJESPtFlavPtHatCentWeighted;    // JES via pt corr, gen pt, flavorForB, ptHat, centrality

  private:

    /// @brief Is Monte Carlo
    Bool_t fIsMc;

    Int_t    fCentBins;
    Double_t fCentRange[2];
    Int_t    fPtBins;
    Double_t fPtRange[2];
    Int_t    fEtaBins;
    Double_t fEtaRange[2];
    Int_t    fPhiBins;
    Double_t fPhiRange[2];
    Int_t    fJESBins;
    Double_t fJESRange[2];
    Int_t    fPtHatBins;
    Double_t fPtHatRange[2];
    Int_t    fFlavorForBBins;
    Double_t fFlavorForBRange[2];
    /// @brief Jet type: PF or Calo
    TString  fJetType;

    ClassDef(HistoManagerJetESR, 0)
};

#endif // #define HistoManagerJetESR_h