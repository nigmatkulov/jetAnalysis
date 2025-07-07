/**
 * @file JetESRAnalysis.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet energy scale and resolution analysis
 * @version 1.2
 * @date 2025-06-02
 * 
 * @copyright Copyright (c) 2025
 * 
 */

#ifndef JetESRAnalysis_h
#define JetESRAnalysis_h

// Load ROOT libraries
#include "TObject.h"
#include "TString.h"
#include "Rtypes.h"
#include "TChain.h"
#include "TF1.h"

// Jet analysis headers
#include "BaseAnalysis.h"
#include "HistoManagerJetESR.h"
#include "Event.h"

// Forward declarations
class JetCut;

//________________
class JetESRAnalysis : public BaseAnalysis {
  public:
    /// @brief Default constructor
    JetESRAnalysis();
    /// @brief Destructor
    virtual ~JetESRAnalysis();

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

    /// @brief Set centrality weight
    void useCentralityWeight()                     { fUseCentralityWeight = {true}; }
    /// @brief Set debug information
    void setVerbose()                              { fVerbose = {true}; }
    /// @brief Add histogram manager to the analysis
    void addHistoManager(HistoManagerJetESR *hm)   { fHM = hm; }
    /// @brief Set collision system
    void setCollisionSystem(const int& syst)       { fCollisionSystem = syst; }
    /// @brief Set collision energy in GeV (default: 8160 GeV)
    void setCollisionEnergyInGeV(const int& en)    { fCollisionEnergy = en; }
    /// @brief Set lead going direction for pPb collisions
    void setPbGoing()                              { fIsPbGoingDir = {true}; }
    /// @brief Set cut on the ptHat of the event (for MC in pPb only due to the xsection matching)
    void setPtHatRange(const double& lo, const double& hi) { fPtHatRange[0] = lo; fPtHatRange[1] = hi; }
    /// @brief Set if is Monte Carlo (default: true)
    void setIsMc(const bool& isMc = true)           { fIsMc = isMc; }
    /// @brief Add lorentz shift
    void setEtaShift(const float& shift)            { fEtaShift = shift; }
    /// @brief Set reco jet cut
    void setRecoJetCut(JetCut *cut)                 { fRecoJetCut = cut; }
    /// @brief Set gen jet cut
    void setGenJetCut(JetCut *cut)                  { fGenJetCut = cut; }

    /// @brief Set number of events in the embedding sample (for cross section recovery in pPb embedding)
    void setNEventsInSample(const int& n)           { fNEventsInSample = n; }

    /// @brief Return collision system name 
    TString collisionSystem() const;

    /// @brief Print JetESRAnalysis setup
    void print();

  private:

    /// @brief Loop over reco, gen and ref-selected reco jets and save jet indices in pT-sorted vectors
    void makePtSortedJetVectors(const Event* event);

    /// Check if event is overweighted in MC
    bool isOverweightedEvent(const Event* event, const double& weight);
    /// @brief  Check if event is overweighted
    bool isOverweighted(const float& ptLead, const float& dijetPtAve, const float& ptHat);

      /// Loop over reco, gen and ref jets and search for leading and subleading jets
    void processInclusiveJets(const Event* event, const double& weight);
    /// @brief Calculate event weight
    double eventWeight(const double& ptHat, const double& vz, const double& centWeight, const double& ptHatW);
    /// @brief Process gen jets
    void processGenJets(const Event* event, const double &weight);
    /// @brief Process reco jets
    void processRecoJets(const Event* event, const double &weight);
    /// @brief Process ref jets
    void processRefJets(const Event* event, const double &weight);
    /// @brief Boost eta to the center-of-mass frame
    float boostEta2CM(const float &etaLab);
    /// @brief Get proper eta in the lab frame depending on beam direction 
    float etaLab(const float &eta);

    /// @brief Initialize vz weight function
    void initVzWeightFunction();

    /// @brief Print debug information
    bool   fVerbose;

    /// @brief Histogram manager
    HistoManagerJetESR *fHM;
    /// @brief Vz weight to match MC to data
    TF1 *fVzWeight;

    /// @brief Centrality weight
    bool   fUseCentralityWeight;
    /// @brief  Pseudorapidity shift for asymmetric collisions (pPb)
    double fEtaShift;
    /// @brief Collision energy in GeV
    int    fCollisionEnergy;
    /// @brief  Type of collision system: 0 - pp, 1 - pPb, 2 - PbPb. Default is PbPb
    int    fCollisionSystem;
    /// @brief Lead going direction for pPb collisions
    bool   fIsPbGoingDir;
    /// @brief ptHat range for the generated events (must cut events on this one) in case of pPb
    double fPtHatRange[2];
    /// @brief Set if is Monte Carlo
    bool   fIsMc;

    /// @brief Number of events in the sample (for pPb cross section MC recovery)
    int    fNEventsInSample;

    // Indices of the leading and subleading jets (at the beginning of the event processing must be set to -1)
    int    fRecoIdLead;
    int    fRecoIdSubLead;
    int    fGenIdLead;
    int    fGenIdSubLead;
    int    fRefSelRecoIdLead;
    int    fRefSelRecoIdSubLead;
    std::vector<int> fRecoPtSortedJetIds;
    std::vector<int> fGenPtSortedJetIds;
    std::vector<int> fRefSelRecoPtSortedJetIds;

    JetCut *fRecoJetCut;
    JetCut *fGenJetCut;

    ClassDef(JetESRAnalysis, 0)
};

#endif // #define JetESRAnalysis_h
