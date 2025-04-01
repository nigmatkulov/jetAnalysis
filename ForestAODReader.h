/**
 * @file ForestAODReader.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief CMS ForestAOD reader
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef ForestAODReader_h
#define ForestAODReader_h

// ROOT headers
#include "Rtypes.h"
#include "TChain.h"
#include "TString.h"
#include "TF1.h"
#include "TRandom3.h"

// JetAnalysis headers
#include "BaseReader.h"
#include "Event.h"
#include "EventCut.h"
#include "JetCut.h"
#include "JetCorrector.h"
#include "JetUncertainty.h"

// C++ headers
#include <list>

//_________________
class ForestAODReader : public BaseReader {

  public:
    /// @brief Default constructor
    ForestAODReader();
    /// @brief Constructor for ForestAODReader
    /// @param inputStream Input file (.root) or list of ROOT files that contain CMS ForestAOD
    ForestAODReader(const char* inputStream, 
                    const bool& useHltBranch = true, const bool& useSkimmingBranch = true, 
                    const bool& useRecoJetBranch = true, 
                    const bool& useTrackBranch = false, const bool& useGenTrackBranch = false, 
                    const bool& isMc = false );
    /// @brief Destructor
    virtual ~ForestAODReader();

    /// @brief  Initialize input
    int init();
    /// @brief Finish (print final information)
    void finish();
    /// Read event and fill objects
    Event* returnEvent();
    /// @brief Report event from reader
    void report();

    /// Turn-on HLT branch to be read
    void useHltBranch()         { fUseHltBranch = {true}; }

    void useSkimmingBranch()    { fUseSkimmingBranch = {true}; }
    /// Turn-on particle flow branch to be read
    void useRecoJetBranch() { fUseRecoJetBranch = {true}; }
    /// @brief Set particle flow jet branch name
    void setRecoJetBranchName(const char *name = "akCs4PFJetAnalyzer") { fRecoJetTreeName = name; }
    /// Turn-on track branch to be read
    void useTrackBranch()       { fUseTrackBranch = {true}; }

    /// @brief Set colliding system
    void setCollidingSystem(const char *sys = "PbPb") { fCollidingSystem = sys; }
    /// @brief Set colliding energy
    void setCollidingEnergy(const int& ene = 5020)    { fCollidingEnergyGeV = {ene}; }
    /// @brief Set year of data taking
    void setYearOfDataTaking(const int& year = 2018)  { fYearOfDataTaking = {year}; }
    /// @brief Path to jetAnalysis directory (or any folder) that contains aux_files/... with JEC corrections
    void setPath2JetAnalysis(const char *name = "../") { fJECPath = name; }
    /// @brief Add JEC file name to the list of JEC files
    void addJECFile(const char *name = "Autumn18_HI_V8_MC_L2Relative_AK4PF.txt") { fJECFiles.push_back(name); }
    /// @brief Use JEU: 0 - not use (default), 1 = +JEU, -1 = -JEU
    void setUseJEU(const int &use = 0)  { fUseJEU = (int)use; }
    /// @brief Add JEU file to the list of JEU files
    void setJEUFileName(const char *name = "Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF.txt") { fJEUInputFileName = name; }
    /// @brief Apply jet pT-smearing
    void setJetPtSmearing(const bool& smear = false) { fDoJetPtSmearing = {smear}; }
    /// @brief Set event cut
    void setEventCut(EventCut *cut) { fEventCut = {cut}; }
    /// @brief Set jet cut
    void setJetCut(JetCut *cut) { fJetCut = {cut}; }
    /// @brief Is the dataset from MC
    void setIsMc() { fIsMc = {true}; }
    /// @brief Correct centrality for MC events
    void setCorrectCentMC() { fCorrectCentMC = {true}; }
    /// @brief Fix jet arrays
    void fixJetArrays() { fFixJetArrays = {true}; }
    /// @brief Use extra JEC correction (in pPb 8160 it is needed for ak4cs)
    void useExtraJECCorrForConstSubtraction() { fUseExtraJECforAk4Cs = {true}; }
    /// @brief Set verbose mode
    void setVerbose()       { fVerbose = {true}; }
    /// @brief Use JES systematics variation: 0 - default, +1 - JES+, -1 - JES-
    void useJERSystematics(const int &syst)  { fUseJERSystematics = syst; }
    /// @brief Set parameters of JER fit with sqrt(a*a + b*b/x)
    void setJERFitParams(const double &a = 0.0415552, const double &b = 0.960013) { fAlphaJER = a; fBetaJER = b; }
    /// @brief Set default parameters of JER for systematics
    void setJERSystParams();
    /// @brief Set use manually calculate JEC correction function
    void setUseManualJEC() { fUseManualJEC = true; }
    /// @brief Set Pb-direction (default is true)
    void setPbGoingDir(const bool &pb = true) { fIsPbGoing = pb; }

    /// @brief Return amount of events to read
    Long64_t nEventsTotal() const { return fEvents2Read; }

  private:

    /// @brief Setup input stream (either single file or a list of files)
    void setInputStream(const char *inputStream) { fInFileName = {inputStream}; }
    /// @brief Instantiate extra JEC scale correction factor (for pPb 2016 for ak4cs)
    void createExtraJECScaleCorrFunction();

    /// Setup input all input
    void setupInput(TString input, TChain *hltChain, TChain *eveChain, TChain *jetChain, 
                    TChain *trkChain, bool useMC, TChain *genTrkChain);
    /// Setup chains to be filled
    int setupChains();
    /// Setup branches
    void setupBranches();

    /// @brief Clear variables for reading
    void clearVariables();
    /// @brief Fix jet arrays
    void fixIndices();

    /// @brief Call read event
    void readEvent();

    /// @brief Setup JEC
    void setupJEC();
    /// @brief Setup JEU
    void setupJEU();
    /// @brief Calculate and return smearing factor
    double  extraJERCorr(const double &pt, const double& eta);
    /// @brief Find resolution factor from JERSyst values for the given eta
    double  retrieveResolutionFactor(const double& eta);

    /// @brief Calculate centrality weight
    double evalCentralityWeight(const double& hiBin);

    /// @brief Find bin index for the given value in the array
    int findBinIndex(const double &val, double *array, int nBins);
    /// @brief Manual JEC correction
    double jecManualCorrection(const double &pt, const double &eta);

    /// @brief Event with jets and other variables
    Event *fEvent;

    /// @brief Input filename (name.root) or file with list of ROOT files 
    const char *fInFileName;
    /// @brief Number of events to process from input file(s)
    Long64_t fEvents2Read;
    /// @brief How many events were processed
    Long64_t fEventsProcessed;

    /// @brief Is file with MC information
    bool fIsMc;
    /// @brief Use to apply hiBin shift for centrality correction
    bool fCorrectCentMC;
    /// @brief Use JEC computed manually (default false)
    bool fUseManualJEC;
    /// @brief Pb-direction (default true)
    bool fIsPbGoing;

    /// @brief Switch HLT branch ON
    bool fUseHltBranch;
    /// @brief Switch skimming branch ON
    bool fUseSkimmingBranch;

    /// @brief Switch particle flow jet branch ON
    bool fUseRecoJetBranch;
    /// @brief Switch track branch ON
    bool fUseTrackBranch;
    /// @brief Switch MC track branch ON
    bool fUseGenTrackBranch;

    /// @brief Chain conaining HLT information (used to friend other trees)
    TChain *fHltTree;
    /// @brief Chain containing skimming information
    TChain *fSkimTree;
    /// @brief Chain containing event information
    TChain *fEventTree;
    /// @brief Chain containing particle flow jets
    TChain *fRecoJetTree;
    /// @brief Chain containing tracks
    TChain *fTrkTree;
    /// @brief Chain containing Monte Carlo tracks
    TChain *fGenTrkTree;

    /// @brief Name of the reconstructed jet tree (e.g. akCs4PFJetAnalyzer for PbPb or ak4PFJetAnalyzer for pp)
    TString fRecoJetTreeName;

    //
    // Variables to store information from TTree
    //

    //
    // Event quantities
    //

    /// @brief Run ID
    unsigned int fRunId;
    /// @brief Event ID
    ULong64_t fEventId;
    /// @brief Luminosity
    unsigned int    fLumi;
    /// @brief Vertex position z
    float fVertexZ; 
    /// @brief Collision centrality (0-200)
    int   fHiBin;
    /// @brief Event weight (prescale from generator level)
    float fPtHatWeight; 
    /// @brief ptHat (initial parton pT) - from generator level
    float fPtHat;

    //
    // Trigger and skimming information
    //

    int fHLT_HIAK4CaloJet60_v1;
    int fHLT_HIAK4CaloJet80_v1;
    int fHLT_PAAK4CaloJet60_Eta5p1_v3;
    int fHLT_PAAK4CaloJet80_Eta5p1_v3;
    int fHLT_PAAK4CaloJet100_Eta5p1_v3;
    int fHLT_PAAK4PFJet60_Eta5p1_v4;
    int fHLT_PAAK4PFJet80_Eta5p1_v3;
    int fHLT_PAAK4PFJet100_Eta5p1_v3;
    int fHLT_PAAK4PFJet120_Eta5p1_v2;

    int fHLT_HIAK4PFJet15_v1;
    int fHLT_HIAK4PFJet15_v1_Prescl;
    int fHLT_HIAK4PFJet30_v1;
    int fHLT_HIAK4PFJet30_v1_Prescl;
    int fHLT_HIAK4PFJet40_v1;
    int fHLT_HIAK4PFJet40_v1_Prescl;
    int fHLT_HIAK4PFJet60_v1;
    int fHLT_HIAK4PFJet60_v1_Prescl;
    int fHLT_HIAK4PFJet80_v1;
    int fHLT_HIAK4PFJet80_v1_Prescl;
    int fHLT_HIAK4PFJet120_v1;
    int fHLT_HIAK4PFJet120_v1_Prescl;

    int fHLT_HIAK8PFJet15_v1;
    int fHLT_HIAK8PFJet15_v1_Prescl;
    int fHLT_HIAK8PFJet25_v1;
    int fHLT_HIAK8PFJet25_v1_Prescl;
    int fHLT_HIAK8PFJet40_v1;
    int fHLT_HIAK8PFJet40_v1_Prescl;
    int fHLT_HIAK8PFJet60_v1;
    int fHLT_HIAK8PFJet60_v1_Prescl;
    int fHLT_HIAK8PFJet80_v1;
    int fHLT_HIAK8PFJet80_v1_Prescl;
    int fHLT_HIAK8PFJet140_v1;
    int fHLT_HIAK8PFJet140_v1_Prescl;

    int fHLT_HIPFJet25_v1;
    int fHLT_HIPFJet25_v1_Prescl;
    int fHLT_HIPFJet140_v1;
    int fHLT_HIPFJet140_v1_Prescl;

    int fHLT_HIPuAK4CaloJet80Eta5p1_v1;
    int fHLT_HIPuAK4CaloJet100Eta5p1_v1;


    // Skimanalysis part
    int fHBHENoiseFilterResultRun2Loose;
    int fHBHENoiseFilterResultRun2Tight;
    int fHBHEIsoNoiseFilterResult;
    int fCollisionEventSelectionAODv2;
    int fPhfCoincFilter2Th4;
    int fPPAprimaryVertexFilter;
    int fPBeamScrapingFilter;
    int fPprimaryVertexFilter;
    int fPVertexFilterCutG;
    int fPVertexFilterCutGloose;
    int fPVertexFilterCutGtight;
    int fPVertexFilterCutE;
    int fPVertexFilterCutEandG;
    int fPClusterCompatibilityFilter;

    int fPhfCoincFilter;
    int fPVertexFilterCutdz1p0;
    int fPVertexFilterCutGplus;
    int fPVertexFilterCutVtx1;

    //
    // Jet information
    //

    /// @brief Number of reconstructed jets
    int   fNRecoJets;
    /// @brief Reconstructed jet transverse momentum (without JEC)
    float fRecoJetPt[10000];
    /// @brief Pseudorapidity of reconstructed jet
    float fRecoJetEta[10000];
    /// @brief Azimuthal angle of reconstructed jet
    float fRecoJetPhi[10000];
    /// @brief WTA eta of reconstructed jet
    float fRecoJetWTAEta[10000];
    /// @brief WTA phi of reconstructed jet
    float fRecoJetWTAPhi[10000];
    /// @brief Track with maximum pT in reconstructed jet
    float fRecoJetTrackMax[10000];

    float fRecoJtPfNHF[10000];
    float fRecoJtPfNEF[10000];
    float fRecoJtPfCHF[10000];
    float fRecoJtPfMUF[10000];
    float fRecoJtPfCEF[10000];
    int fRecoJtPfCHM[10000];
    int fRecoJtPfCEM[10000];
    int fRecoJtPfNHM[10000];
    int fRecoJtPfNEM[10000];
    int fRecoJtPfMUM[10000];

    /// @brief Transverse momentum of generated jet that was matched with reconstructed jet
    float fRefJetPt[10000];
    /// @brief /// @brief Pseudorapidity of generated jet that was matched with reconstructed jet
    float fRefJetEta[10000];
    /// @brief Azimuthal angle of generated jet that was matched with reconstructed jet
    float fRefJetPhi[10000];
    /// @brief WTA eta of generated jet that was matched with reconstructed jet
    float fRefJetWTAEta[10000];
    /// @brief WTA phi of generated jet that was matched with reconstructed jet
    float fRefJetWTAPhi[10000];
    /// @brief Parton flavor of generated jet that was matched with reconstructed jet
    int   fRefJetPartonFlavor[10000];
    /// @brief Parton flavor for B of generated jet that was matched with reconstructed jet
    int   fRefJetPartonFlavorForB[10000];

    /// @brief Number of generated jets
    int   fNGenJets;
    /// @brief Generated jet transverse momentum
    float fGenJetPt[10000];
    /// @brief Pseudorapidity of generated jet
    float fGenJetEta[10000];
    /// @brief Azimuthal angle of generated jet
    float fGenJetPhi[10000];
    /// @brief WTA eta of generated jet
    float fGenJetWTAEta[10000];
    /// @brief WTA phi of generated jet
    float fGenJetWTAPhi[10000];

    //
    // Reconstructed tracks
    //

    /// @brief Number of tracks
    int   fNTracks;
    /// @brief Track transverse momentum
    float fTrackPt[20000];
    /// @brief Track pseudorapidity
    float fTrackEta[20000];
    /// @brief Track azimuthal angle
    float fTrackPhi[20000];
    /// @brief Track pT error (uncertainty)
    float fTrackPtErr[20000];
    /// @brief Track distance of closest approach in transverse plane (XY)
    float fTrackDcaXY[20000];
    /// @brief Track distance of closest approach in beam direction (z)
    float fTrackDcaZ[20000];
    /// @brief Track distance of closest approach error in transverse plane (XY)
    float fTrackDcaXYErr[20000];
    /// @brief Track distance of closest approach error in beam direction (z)
    float fTrackDcaZErr[20000];
    /// @brief Track fitting (reconstruction) chi2
    float fTrackChi2[20000];
    /// @brief Track number of degrees of freedom in the fitting 
    unsigned char fTrackNDOF[20000];    
    /// @brief Particle flow energy deposited in ECAL from the given track
    float fTrackPartFlowEcal[20000];
    /// @brief Particle flow energy deposited in HCAL from the given track
    float fTrackPartFlowHcal[20000];
    /// @brief Track MVA for each step
    float fTrackMVA[20000];
    /// @brief Track algorithm/step
    unsigned char fTrackAlgo[20000];
    /// @brief Track charge
    int   fTrackCharge[20000];
    /// @brief Number of hits in the tracker
    unsigned char fTrackNHits[20000];
    /// @brief Number of layers with measurement in the tracker
    unsigned char fTrackNLayers[20000];
    /// @brief Tracker steps MVA selection
    bool  fTrackHighPurity[20000];

    //
    // Monte Carlo tracks
    //

    /// @brief Generated particle transverse momentum
    std::vector<float> fGenTrackPt;
    /// @brief Generated particle pseudorapidity
    std::vector<float> fGenTrackEta;
    /// @brief Generated particle azimuthal angle
    std::vector<float> fGenTrackPhi;
    /// @brief Generated particle charge
    std::vector<int>   fGenTrackCharge;
    /// @brief Generated particle PID
    std::vector<int>   fGenTrackPid;
    /// @brief Generated particle sube (?)
    std::vector<int>   fGenTrackSube;

    /// @brief Jet Energy Corrector instance
    JetCorrector *fJEC;
    /// @brief List of files with JEC
    std::vector< std::string > fJECFiles;
    /// @brief Path to jetAnalysis directory
    TString fJECPath;
    /// @brief Jet Energy Uncertainty instance
    JetUncertainty *fJEU;
    /// @brief Input file name with JEU correction
    TString fJEUInputFileName;

    /// @brief Colliding system: pp, pPb or PbPb
    TString fCollidingSystem;
    /// @brief Colliding energy
    int fCollidingEnergyGeV;
    /// @brief Year of data taking
    int fYearOfDataTaking;
    /// @brief Apply jet pT-smearing 
    bool fDoJetPtSmearing;

    /// @brief Fix indices
    bool fFixJetArrays;

    /// @brief Event cut
    EventCut *fEventCut;
    /// @brief Jet cut
    JetCut *fJetCut;

    /// @brief Vector that contains indices of generated jets that matched to the reconsructed 
    /// particle flow jet (should be of the reco/red size)
    std::vector<int> fRecoJet2GenJetId;
    /// @brief Vector that contains indices of the reconstructed particle flow jets that 
    /// macthed to generated jet
    std::vector<int> fGenJet2RecoJet;

    /// @brief Use extra correction for JEC (for 8160 pPb year 2016 ak4pf is default. extra needed for ak4cs)
    bool  fUseExtraJECforAk4Cs;
    /// @brief JEC extra correction
    TF1    *fJECScaleCorr;
    /// @brief Use JEU correction: 0 - not use (default), 1 - +JEU, -1 - -JEU
    int   fUseJEU;
    /// @brief  For MC use extra correction to address imperfection of MC w.r.t. data.
    /// 0 - default, 1 - JES+, -1 - JES-
    int   fUseJERSystematics;
    std::vector<double> fJerEtaLow;
    std::vector<double> fJerEtaHi;
    std::vector<double> fJerDef;
    std::vector<double> fJerLow;
    std::vector<double> fJerHi;
    /// @brief First parameter of JER fit with sqrt(a*a + b*b/x)
    double fAlphaJER;
    /// @brief Second parameter of JER fit with sqrt(a*a + b*b/x)
    double fBetaJER;
    /// @brief Gaussian distribution that is used to smear JER
    TF1 *fJERSmearFunc;
    /// @brief Random number generator
    TRandom3 *fRndm;

    /// @brief  Verbose mode
    bool  fVerbose;

    ClassDef(ForestAODReader, 1)
};

#endif // #define ForestAODReader_h