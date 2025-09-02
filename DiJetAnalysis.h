/**
 * @file DiJetAnalysis.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Dijet analysis
 * @version 1.1
 * @date 2025-01-09
 * 
 * @copyright Copyright (c) 2025
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
#include "DiJet.h"

// C++ headers
#include <vector>
#include <map>

// Forward declarations
class JetCut;
class DiJetCut;

//_________________
struct MixEvent {
    /// @brief Leading and subleading jets in the current event
    DiJet fGenDijetCM; 
};

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
    void setEtaShift(const float& shift)  { fEtaShift = shift; }
    /// @brief Set dataset to be MC
    void setIsMc()         { fIsMc = {true}; }
    /// @brief Set collision system: 0 - pp, 1 - pPb, 2 - PbPb
    void setCollisionSystem(const int& syst) { fCollisionSystem = syst; }
    /// @brief Set collision energy in GeV (default: 8160 GeV)
    void setCollisionEnergyInGeV(const int& en)    { fCollisionEnergy = en; }
    /// @brief Set cut on the ptHat of the event (for MC in pPb only due to the xsection matching)
    void setPtHatRange(const float& lo, const float& hi) { fPtHatRange[0] = lo; fPtHatRange[1] = hi; }
    /// @brief Set the direction of Pb-going ion
    void setPbGoing()                        { fIsPbGoingDir = {true}; }
    /// @brief Set the direction of p-going ion
    void setpGoing()                         { fIsPbGoingDir = {false}; }
    /// @brief Set verbose mode
    void setVerbose()                        { fVerbose = {true}; }
    /// @brief Set number of events in the embedding sample (for the given ptHat)
    void setNEventsInSample(const int& n)  { fNEventsInSample = n; }

    /// @brief Set reco jet cut
    void setRecoJetCut(JetCut *cut) { fRecoJetCut = cut; }
    /// @brief Set gen jet cut
    void setGenJetCut(JetCut *cut)  { fGenJetCut = cut; }
    /// @brief Set dijet cut
    void setDiJetCut(DiJetCut *cut)   { fDiJetCut = cut; }

    /// @brief Reweight MC to data (trigger-dependent): 
    /// 0 - do not reweight (default)
    /// 1 - MB
    /// 2 - Jet60
    /// 3 - Jet80
    /// 4 - Jet100
    void setUseMcReweighting(const int& w = 0) { fUseMcReweighting = (short)w; }

    void findMcWeight(const float& ptLead, const float& ptSublead);
    /// @brief Print DiJetAnalysis setup
    void print();

    /// @brief Return collision system name 
    TString collisionSystem() const;

    /// @brief Set the size of the mixing buffer
    void setMixBufferSize(const int& size) { (size > 0) ? fMixBufferSize = size : fMixBufferSize = 10; }

  private:

    /// @brief Loop over reco, gen and ref-selected reco jets and save jet indices in pT-sorted vectors
    void makePtSortedJetVectors(const Event* event);

    /// Check if event is overweighted in MC
    bool isOverweightedEvent(const Event* event, const double& weight);
    /// @brief  Check if event is overweighted
    bool isOverweighted(const float& ptLead, const float& dijetPtAve, const float& ptHat);
    /// @brief Calculate event weight
    double eventWeight(const float& ptHat, const float& vz, const float& centWeight, const float& ptHatW);

    /// Loop over reco, gen and ref jets and search for leading and subleading jets
    void processInclusiveJets(const Event* event, const double& weight);
    /// @brief Process gen jets
    void processGenJets(const Event* event, const double &weight);
    /// @brief Process reco jets
    void processRecoJets(const Event* event, const double &weight);
    /// @brief Process ref jets
    void processRefJets(const Event* event, const double &weight);

    /// @brief Process dijets
    void processDijets(const Event* event, const double &weight);
    /// @brief Process gen dijets
    void processGenDijets(const Event* event, const double &weight);
    /// @brief Process reco dijets
    void processRecoDijets(const Event* event, const double &weight);
    /// @brief Process ref dijets
    void processRefDijets(const Event* event, const double &weight);

    /// @brief Boost eta to the center-of-mass frame
    float boostEta2CM(const float &etaLab);
    /// @brief Get proper eta in the lab frame depending on beam direction 
    float etaLab(const float &eta);
    /// @brief Dijet eta calculation
    float dijetEtaInFrame(const float& eta1, const float& eta2, bool isCM = false);

    /// @brief Initialize vz weight function
    void initVzWeightFunction();

    /// @brief Add event to mixing buffer (Olga)
    void addEventToMixBufferOlga(const double &ptAve, const DiJet& dijet);
    /// @brief Add event to mixing buffer
    void addEventToMixBuffer(const double &vz, const DiJet& dijet);

    /// @brief Vz weight to match MC to data
    TF1 *fVzWeight;
    /// @brief Dijet ptAve weight (to match PYTHIA 2 pp data)
    TF1 *fDijetPtAveWeight;
    /// @brief Centrality weight
    bool   fUseCentralityWeight;
    /// @brief Histogram manager
    HistoManagerDiJet *fHM;
    /// @brief  Pseudorapidity shift for asymmetric collisions (pPb)
    float fEtaShift;
    /// @brief Is MC sample (needed for event weight corrections)
    bool   fIsMc;
    /// @brief  Type of collision system: 0 - pp, 1 - pPb, 2 - PbPb. Default is PbPb
    int    fCollisionSystem;
    /// @brief Collision energy in GeV
    int    fCollisionEnergy;
    /// @brief ptHat range for the generated events (must cut events on this one)
    float fPtHatRange[2];

    /// @brief Lead going direction for pPb collisions
    bool   fIsPbGoingDir;
    /// @brief Verbose mode
    bool   fVerbose;
    /// @brief Number of events in the embedding sample
    int    fNEventsInSample;

    /// @brief Gen dijet in the lab frame found (default: false)
    bool   fIsGenDijetLabFound;
    /// @brief Gen dijet in the center-of-mass frame found
    bool   fIsGenDijetCMFound;
    /// @brief Reco dijet in the lab frame found
    bool   fIsRecoDijetLabFound;
    /// @brief Reco dijet in the center-of-mass frame found
    bool   fIsRecoDijetCMFound;
    /// @brief Ref-selected dijet in the lab frame found
    bool   fIsRefSelDijetLabFound;
    /// @brief Ref-selected dijet in the center-of-mass frame found
    bool   fIsRefSelDijetCMFound;

    /// @brief Reweight MC to data (trigger-dependent): 
    /// 0 - do not reweight (default)
    /// 1 - MB
    /// 2 - Jet60
    /// 3 - Jet80
    /// 4 - Jet100
    short   fUseMcReweighting;
    int     fJetPtBins;
    float   fJetPtLeadPtSubleadReweightMatrix[75][75];
    double  fMcReweight;

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

    DiJet *fRecoDijet;
    DiJet *fGenDijet;
    DiJet *fRefDijet;

    JetCut *fRecoJetCut;
    JetCut *fGenJetCut;
    DiJetCut *fDiJetCut;

    /// @brief Size of the mixing buffer (default: 10)
    int fMixBufferSize = 10; 
    /// @brief Map of events for mixing (key: ptAve bin, value: MixEvent)
    std::map<int, MixEvent> fMixBufferOlga;
    /// @brief Map of events for mixing (key: vz, value: MixEvent)
    std::map<int, MixEvent> fMixBuffer;

  ClassDef(DiJetAnalysis, 0)
};

#endif // #define DiJetAnalysis_h
