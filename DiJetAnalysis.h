/**
 * @file DiJetAnalysis.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Dijet analysis
 * @version 0.1
 * @date 2024-01-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef DiJetAnalysis_h
#define DiJetAnalysis_h

// Load ROOT libraries
#include "TObject.h"
#include "TString.h"
#include "Rtypes.h"
#include "TChain.h"

// Jet analysis headers
#include "BaseAnalysis.h"
#include "HistoManagerDiJet.h"
#include "Event.h"

//________________
class DiJetAnalysis : public BaseAnalysis {
  public:
    /// @brief Default constructor
    DiJetAnalysis();
    /// @brief Destructor
    virtual ~DiJetAnalysis();

    /// @brief Initialize variables and functions
    void init();
    /// @brief Process event
    void processEvent(const Event* ev);
    /// @brief Finish analysis
    void finish();

    /// @brief Returns reports of all cuts applied and correlation functions being done
    virtual void report();
    /// @brief Return a TList of objects to be written as output
    virtual TList* getOutputList();

    /// @brief Set debug information
    void setDebug(const Bool_t& debug) { fDebug = debug; }
    /// @brief Add histogram manager to the analysis
    void addHistoManager(HistoManagerDiJet *hm) { fHM = hm; }
    /// @brief Add lorentz shift
    void setEtaShift(const Double_t& shift) { fEtaShift = shift; }
    /// @brief Set dataset to be MC
    void setIsMc(const Bool_t& isMc) { fIsMc = isMc; }
    /// @brief Is pPb dataset
    void setIsPPb()                  { fIsPPb = kTRUE; }
    /// @brief Set cut on the ptHat of the event (for MC in pPb only due to the xsection matching)
    void setPtHatRange(const Double_t& lo, const Double_t& hi) { fPtHatRange[0] = lo; fPtHatRange[1] = hi; }
    /// @brief Cut on the lowest momentum of leading jet
    void setLeadJetPtLow(const Double_t& lo) { fLeadJetPtLow = lo; }
    /// @brief Cut on the lowest momentum of subleading jet
    void setSubLeadJetPtLow(const Double_t& lo) { fSubleadJetPtLow = lo; }
    /// @brief Cut on angle between leading and subleading jet
    void setDijetPhiCut(const Double_t& cut) { fDijetPhiCut = cut; }
    /// @brief Set the direction of Pb-going ion
    void setPbGoing()                        { fIsPbGoingDir = {kTRUE}; }
    /// @brief Set the direction of p-going ion
    void setpGoing()                         { fIsPbGoingDir = {kFALSE}; }
    /// @brief Set verbose mode
    void setVerbose()                        { fVerbose = {kTRUE}; }
    /// @brief Set number of events in the embedding sample (for the given ptHat)
    void setNEventsInSample(const Int_t& n)  { fNEventsInSample = n; }
    /// @brief Print DiJetAnalysis setup
    void print();

  private:

    /// Calculate event weight
    Double_t eventWeight(const Bool_t& isMc, const Bool_t& isPPb, const Double_t& ptHat, const Double_t& vz);
    /// Process gen jets
    void processGenJets(const Event* event, Double_t ptHatW);
    /// Process reco jets
    void processRecoJets(const Event* event, Double_t ptHatW);
    /// Process ref jets
    void processRefJets(const Event* event, Double_t ptHatW);
    /// Dijet selection
    Bool_t isGoodDijet(const Double_t& ptLead, const Double_t& ptSublead, const Double_t& dphi);
    /// Calculate delta phi between two jets in the range [-pi, pi]
    Double_t deltaPhi(const Double_t& phi1, const Double_t phi2);

    /// @brief Print debug information
    Bool_t   fDebug;
    /// @brief Centrality weight
    Bool_t   fUseCentralityWeight;
    /// @brief Histogram manager
    HistoManagerDiJet *fHM;
    /// @brief  Pseudorapidity shift for asymmetric collisions (pPb)
    Double_t fEtaShift;
    /// @brief Is MC sample (needed for event weight corrections)
    Bool_t   fIsMc;
    /// @brief  Is pPb dataset
    Bool_t   fIsPPb;
    /// @brief ptHat range for the generated events (must cut events on this one)
    Double_t fPtHatRange[2];

    /// @brief Momentum selection of the leading jet
    Double_t fLeadJetPtLow;
    /// @brief Momentum selection of the subleading jet
    Double_t fSubleadJetPtLow;
    /// @brief Angular selection of dijet
    Double_t fDijetPhiCut;
    /// @brief Lead going direction for pPb collisions
    Bool_t   fIsPbGoingDir;
    /// @brief Verbose mode
    Bool_t   fVerbose;
    /// @brief Number of events in the embedding sample
    Int_t    fNEventsInSample;

  ClassDef(DiJetAnalysis, 0)
};

#endif // #define DiJetAnalysis_h