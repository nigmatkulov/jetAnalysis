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
#include "TF1.h"

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

    /// @brief Add histogram manager to the analysis
    void addHistoManager(HistoManagerDiJet *hm) { fHM = hm; }
    /// @brief Add lorentz shift
    void setEtaShift(const Double_t& shift)  { fEtaShift = shift; }
    /// @brief Set dataset to be MC
    void setIsMc(const Bool_t& isMc)         { fIsMc = isMc; }
    /// @brief Is pPb dataset
    void setIsPPb()                          { fIsPPb = kTRUE; }
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
    /// @brief Set loose jetId cut
    void setLooseJetIdCut()                  { fIsLooseJetIdCut = {kTRUE}; }
    /// @brief Select inclusive jets in the center-of-mass frame
    void selectJetsInCMFrame()               { fSelectJetsInCMFrame = {kTRUE}; }
    /// @brief Set eta range to select jets in the lab frame
    void setJetEtaLabRange(const Double_t& lo, const Double_t& hi) { fJetEtaLab[0]=lo; fJetEtaLab[1]=hi; }
    /// @brief Set eta range to select jets in the center-of-mass frame
    void setJetEtaCMRange(const Double_t& lo, const Double_t& hi) { fJetEtaCM[0]=lo; fJetEtaCM[1]=hi; }

    /// @brief Reweight MC to data (trigger-dependent): 
    /// 0 - do not reweight (default)
    /// 1 - MB
    /// 2 - Jet60
    /// 3 - Jet80
    /// 4 - Jet100
    void setUseMcReweighting(const Int_t& w = 0) { fUseMcReweighting = (Short_t)w; }

    void findMcWeight(const Double_t& ptLead, const Double_t& ptSublead);
    /// @brief Print DiJetAnalysis setup
    void print();

  private:

    /// @brief Calculate event weight
    Double_t eventWeight(const Bool_t& isMc, const Bool_t& isPPb, const Double_t& ptHat, const Double_t& vz);
    /// @brief Process gen jets
    void processGenJets(const Event* event, Double_t ptHatW);
    /// @brief Process reco jets
    void processRecoJets(const Event* event, Double_t ptHatW);
    /// @brief Process ref jets
    void processRefJets(const Event* event, Double_t ptHatW);
    /// @brief Dijet selection
    Bool_t isGoodDijet(const Double_t& ptLead, const Double_t& ptSublead, const Double_t& dphi);
    /// @brief Calculate delta phi between two jets in the range [-pi, pi]
    Double_t deltaPhi(const Double_t& phi1, const Double_t &phi2);
    /// @brief Single gen/ref jet selection criteria
    Bool_t isGoodGenJet(const GenJet* jet);
    /// @brief Single reco jet selection criteria
    Bool_t isGoodRecoJet(const RecoJet* jet);
    /// @brief Check if jet passes jetId requirements
    Bool_t isGoodJetId(const RecoJet* jet);
    /// @brief Check if good track max cut
    Bool_t isGoodTrkMax(const RecoJet* jet);
    /// @brief Boost eta to the center-of-mass frame
    Double_t boostEta2CM(const Double_t &etaLab);
    /// @brief Get proper eta in the lab frame depending on beam direction 
    Double_t etaLab(const Double_t &eta);
    /// @brief Dijet eta calculation
    Double_t dijetEtaInFrame(const Double_t& eta1, const Double_t& eta2, Bool_t isCM = kFALSE);
    
    /// @brief Pass pt of the jet and check if it is leading or subleading jet
    void findLeadSubleadJets(const double &pt, const int &counter, double &ptLead, double &ptSublead, 
                             int &idLead, int &idSubLead);

    /// @brief Vz weight to match MC to data
    TF1 *fVzWeight;
    /// @brief Dijet ptAve weight (to match PYTHIA 2 pp data)
    TF1 *fDijetPtAveWeight;
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
    /// @brief Is loose/tight jetId cut (default: false = tight)
    Bool_t   fIsLooseJetIdCut;
    /// @brief Check if dijet passed trkMax cut is found
    Bool_t   fIsDijetFound;
    /// @brief Check if dijet passed jetId cut is found
    Bool_t   fIsDijetJetIdFound;
    /// @brief Select jets in the center-of-mass frame (default: false)
    Bool_t   fSelectJetsInCMFrame;
    /// @brief Reweight MC to data (trigger-dependent): 
    /// 0 - do not reweight (default)
    /// 1 - MB
    /// 2 - Jet60
    /// 3 - Jet80
    /// 4 - Jet100
    Short_t   fUseMcReweighting;
    Int_t     fJetPtBins;
    Double_t  fJetPtLow;
    Double_t  fJetPtHi;
    Double_t  fJetPtStep;
    Double_t  fJetPtLeadPtSubleadReweightMatrix[75][75];
    Double_t  fMcReweight;
    /// Range of eta selection in the lab frame
    Double_t  fJetEtaLab[2];
    /// Range of eta selection in the center-of-mass frame
    Double_t  fJetEtaCM[2];

    Int_t     fEventCounter;
    Int_t     fCycleCounter;
    UInt_t    fTotalCounter;

  ClassDef(DiJetAnalysis, 0)
};

#endif // #define DiJetAnalysis_h