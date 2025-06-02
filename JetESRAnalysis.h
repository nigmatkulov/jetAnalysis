/**
 * @file JetESRAnalysis.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet energy scale and resolution analysis
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
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
    /// @brief Set jetId selection of the jets (default: trkMax)
    void useJetIdSelection()                       { fUseJetIdSelection = {true}; }
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

    /// @brief Set leading jet pt low cut
    void setLeadJetPtLowCut(const double& pt)       { fLeadJetPtLow = pt; }
    /// @brief Set leading jet eta range
    void setLeadJetEtaCut(const double& lo, const double& hi) { fLeadJetEta[0] = lo; fLeadJetEta[1] = hi; }
    /// @brief Set subleading jet pt low cut
    void setSubleadJetPtLowCut(const double& pt)    { fSubleadJetPtLow = pt; }
    /// @brief Set subleading jet eta range
    void setSubleadJetEtaCut(const double& lo, const double& hi) { fSubleadJetEta[0] = lo; fSubleadJetEta[1] = hi; }
    /// @brief Set dijet dPhi cut
    void setDijetDPhiCut(const double& cut)         { fDijetDPhiCut = cut; }
    /// @brief Set number of events in the embedding sample (for cross section recovery in pPb embedding)
    void setNEventsInSample(const int& n)           { fNEventsInSample = n; }

    /// @brief Set loose jetId cut (default: tight = false)
    void setLooseJetIdCut()                         { fIsLooseJetIdCut = {true}; }

    /// @brief Return collision system name 
    TString collisionSystem() const;

    /// @brief Print JetESRAnalysis setup
    void print();

  private:

    /// @brief Check if dijet is good
    bool isGoodDijet(const double& ptLead, const double& etaLead, const double& ptSubLead, 
                     const double& etaSubLead, const double& dphi);
    /// @brief Calculate event weight
    double eventWeight(const double& ptHat, const double& vz, const double& centWeight, const double& ptHatW);
    /// @brief Process gen jets
    void processGenJets(const Event* event, double weight);
    /// @brief Process reco jets
    void processRecoJets(const Event* event, double weight);
    /// @brief Process ref jets
    //void processRefJets(const Event* event, double weight);
    /// @brief Pass pt of the jet and check if it is leading or subleading jet
    void findLeadSubleadJets(const double &pt, const int &counter, double &ptLead, double &ptSublead, 
                             int &idLead, int &idSubLead);
    /// @brief Check if jet passes jetId requirements
    bool isGoodJetId(const RecoJet* jet);
    /// @brief Check if good track max cut
    bool isGoodTrkMax(const RecoJet* jet);

    /// @brief Initialize vz weight function
    void initVzWeightFunction();

    /// @brief Check if gen jet is good
    bool isGoodGenJet(const GenJet* jet);
    /// @brief Check if reco jet is good
    bool isGoodRecoJet(const RecoJet* jet);
    /// @brief Angle between two jets
    double deltaPhi(const double& phi1, const double &phi2);

    /// @brief Pring debug information
    bool   fVerbose;
    /// @brief Histogram manager
    HistoManagerJetESR *fHM;

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

    /// @brief Use jetId selection (default - false, i.e. trkMax)
    bool   fUseJetIdSelection;
    /// @brief Is loose/tight jetId cut (default: false = tight)
    bool   fIsLooseJetIdCut;

    /// @brief Leading jet pt low cut
    double fLeadJetPtLow;
    /// @brief Leading jet eta range
    double fLeadJetEta[2];
    /// @brief Subleading jet pt low cut
    double fSubleadJetPtLow;
    /// @brief Subleading jet eta range
    double fSubleadJetEta[2];
    /// @brief Dijet phi cut
    double fDijetDPhiCut;
    /// @brief Number of events in the sample (for pPb cross section MC recovery)
    int    fNEventsInSample;

    /// @brief Vz weight to match MC to data
    TF1 *fVzWeight;
};

#endif // #define JetESRAnalysis_h
