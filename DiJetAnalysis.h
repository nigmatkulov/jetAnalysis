/**
 * @file DiJetAnalysis.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Dijet analysis
 * @version 1.1
 * @date 2025-01-07
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

// C++ headers
#include <vector>

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
    void setEtaShift(const double& shift)  { fEtaShift = shift; }
    /// @brief Set dataset to be MC
    void setIsMc(const bool& isMc)         { fIsMc = isMc; }
    /// @brief Is pPb dataset
    void setIsPPb()                          { fIsPPb = true; }
    /// @brief Set cut on the ptHat of the event (for MC in pPb only due to the xsection matching)
    void setPtHatRange(const double& lo, const double& hi) { fPtHatRange[0] = lo; fPtHatRange[1] = hi; }
    /// @brief Cut on the lowest momentum of leading jet
    void setLeadJetPtLow(const double& lo) { fLeadJetPtLow = lo; }
    /// @brief Cut on the lowest momentum of subleading jet
    void setSubLeadJetPtLow(const double& lo) { fSubleadJetPtLow = lo; }
    /// @brief Cut on angle between leading and subleading jet
    void setDijetPhiCut(const double& cut) { fDijetPhiCut = cut; }
    /// @brief Set the direction of Pb-going ion
    void setPbGoing()                        { fIsPbGoingDir = {true}; }
    /// @brief Set the direction of p-going ion
    void setpGoing()                         { fIsPbGoingDir = {false}; }
    /// @brief Set verbose mode
    void setVerbose()                        { fVerbose = {true}; }
    /// @brief Set number of events in the embedding sample (for the given ptHat)
    void setNEventsInSample(const int& n)  { fNEventsInSample = n; }
    /// @brief Set loose jetId cut
    void setLooseJetIdCut()                  { fIsLooseJetIdCut = {true}; }
    /// @brief Select inclusive jets in the center-of-mass frame
    void selectJetsInCMFrame()               { fSelectJetsInCMFrame = {true}; }
    /// @brief Set eta range to select jets in the lab frame
    void setJetEtaLabRange(const double& lo, const double& hi) { fJetEtaLab[0]=lo; fJetEtaLab[1]=hi; }
    /// @brief Set eta range to select jets in the center-of-mass frame
    void setJetEtaCMRange(const double& lo, const double& hi) { fJetEtaCM[0]=lo; fJetEtaCM[1]=hi; }

    /// @brief Reweight MC to data (trigger-dependent): 
    /// 0 - do not reweight (default)
    /// 1 - MB
    /// 2 - Jet60
    /// 3 - Jet80
    /// 4 - Jet100
    void setUseMcReweighting(const int& w = 0) { fUseMcReweighting = (Short_t)w; }

    /// @brief Set jetId selection of the jets (default: trkMax)
    void useJetIdSelection()                       { fUseJetIdSelection = {true}; }

    void findMcWeight(const double& ptLead, const double& ptSublead);
    /// @brief Print DiJetAnalysis setup
    void print();

  private:

    /// @brief Calculate event weight
    double eventWeight(const bool& isMc, const bool& isPPb, const double& ptHat, const double& vz);
    /// @brief Process gen jets
    void processGenJets(const Event* event, double ptHatW);
    /// @brief Process reco jets
    void processRecoJets(const Event* event, double ptHatW);
    /// @brief Process ref jets
    void processRefJets(const Event* event, double ptHatW);
    /// @brief Dijet selection
    bool isGoodDijet(const double& ptLead, const double& ptSublead, const double& dphi);

    /// @brief Dijet selection
    bool isGoodDijet(const double& ptLead, const double& etaLead, const double& ptSubLead, 
                     const double& etaSubLead, const double& dphi, const bool& isCM = false);


    /// @brief Calculate delta phi between two jets in the range [-pi, pi]
    double deltaPhi(const double& phi1, const double &phi2);
    /// @brief Single gen/ref jet selection criteria
    bool isGoodGenJet(const GenJet* jet);
    /// @brief Single reco jet selection criteria
    bool isGoodRecoJet(const RecoJet* jet);
    /// @brief Check if jet passes jetId requirements
    bool isGoodJetId(const RecoJet* jet);
    /// @brief Check if good track max cut
    bool isGoodTrkMax(const RecoJet* jet);
    /// @brief Boost eta to the center-of-mass frame
    double boostEta2CM(const double &etaLab);
    /// @brief Get proper eta in the lab frame depending on beam direction 
    double etaLab(const double &eta);
    /// @brief Dijet eta calculation
    double dijetEtaInFrame(const double& eta1, const double& eta2, bool isCM = false);

    /// @brief Pass pt of the jet and check if it is leading or subleading jet
    void findLeadSubleadJets(const double &pt, const int &counter, double &ptLead, double &ptSublead, 
                             int &idLead, int &idSubLead);
    /// @brief Find dijet ptAve bin
    int  findDijetPtAveBin(const double& pt);
    /// @brief Find dijet ptAve bin (old binning)
    int  findDijetPtAveOldBin(const double& pt);

    /// @brief Vz weight to match MC to data
    TF1 *fVzWeight;
    /// @brief Dijet ptAve weight (to match PYTHIA 2 pp data)
    TF1 *fDijetPtAveWeight;
    /// @brief Centrality weight
    bool   fUseCentralityWeight;
    /// @brief Histogram manager
    HistoManagerDiJet *fHM;
    /// @brief  Pseudorapidity shift for asymmetric collisions (pPb)
    double fEtaShift;
    /// @brief Is MC sample (needed for event weight corrections)
    bool   fIsMc;
    /// @brief  Is pPb dataset
    bool   fIsPPb;
    /// @brief ptHat range for the generated events (must cut events on this one)
    double fPtHatRange[2];

    /// @brief Momentum selection of the leading jet
    double fLeadJetPtLow;
    
    /// @brief Momentum selection of the subleading jet
    double fSubleadJetPtLow;
    /// @brief Angular selection of dijet
    double fDijetPhiCut;
    /// @brief Lead going direction for pPb collisions
    bool   fIsPbGoingDir;
    /// @brief Verbose mode
    bool   fVerbose;
    /// @brief Number of events in the embedding sample
    int    fNEventsInSample;
    /// @brief Use jetId selection (default - false, i.e. trkMax)
    bool   fUseJetIdSelection;
    /// @brief Is loose/tight jetId cut (default: false = tight)
    bool   fIsLooseJetIdCut;
    /// @brief Check if dijet passed trkMax cut is found
    bool   fIsDijetFound;
    /// @brief Check if dijet passed jetId cut is found
    bool   fIsDijetJetIdFound;
    /// @brief Select jets in the center-of-mass frame (default: false)
    bool   fSelectJetsInCMFrame;
    /// @brief Reweight MC to data (trigger-dependent): 
    /// 0 - do not reweight (default)
    /// 1 - MB
    /// 2 - Jet60
    /// 3 - Jet80
    /// 4 - Jet100
    Short_t   fUseMcReweighting;
    int     fJetPtBins;
    double  fJetPtLow;
    double  fJetPtHi;
    double  fJetPtStep;
    double  fJetPtLeadPtSubleadReweightMatrix[75][75];
    double  fMcReweight;
    /// Range of eta selection in the lab frame
    double  fJetEtaLab[2];
    /// Range of eta selection in the center-of-mass frame
    double  fJetEtaCM[2];

    int     fEventCounter;
    int     fCycleCounter;
    int     fTotalCounter;

    std::vector<double> fPtAveBins;
    std::vector<double> fPtAveOldBins;

  ClassDef(DiJetAnalysis, 0)
};

#endif // #define DiJetAnalysis_h