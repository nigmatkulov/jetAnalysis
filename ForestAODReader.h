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
    ForestAODReader(const Char_t* inputStream, 
                    const Bool_t& useHltBranch = kTRUE, const Bool_t& useSkimmingBranch = kTRUE, 
                    const Bool_t& usePFJetBranch = kTRUE, const Bool_t& useCaloJetBranch = kFALSE, 
                    const Bool_t& useTrackBranch = kFALSE, const Bool_t& useGenTrackBranch = kFALSE, 
                    const Bool_t& isMc = kFALSE );
    /// @brief Destructor
    virtual ~ForestAODReader();

    /// @brief  Initialize input
    Int_t init();
    /// @brief Finish (print final information)
    void finish();
    /// Read event and fill objects
    Event* returnEvent();
    /// @brief Report event from reader
    void report();

    /// Turn-on HLT branch to be read
    void useHltBranch()         { fUseHltBranch = {kTRUE}; }

    void useSkimmingBranch()    { fUseSkimmingBranch = {kTRUE}; }
    /// Turn-on particle flow branch to be read
    void usePartFlowJetBranch() { fUsePartFlowJetBranch = {kTRUE}; }
    /// @brief Set particle flow jet branch name
    void setPartFlowJetBranchName(const Char_t *name = "akCs4PFJetAnalyzer") { fPFTreeName = name; }
    /// Turn-on calorimeter jet branch to be read
    void useCaloJetBranch()     { fUseCaloJetBranch = {kTRUE}; }
    /// @brief Set calorimeter jet branch name
    void setCaloJetBranchName(const Char_t *name = "akPu4CaloJetAnalyzer") { fCaloTreeName = name; }
    /// Turn-on calorimeter jet branch to be read
    void useTrackBranch()       { fUseTrackBranch = {kTRUE}; }

    /// @brief Set colliding system
    void setCollidingSystem(const Char_t *sys = "PbPb") { fCollidingSystem = sys; }
    /// @brief Set colliding energy
    void setCollidingEnergy(const Int_t& ene = 5020)    { fCollidingEnergyGeV = {ene}; }
    /// @brief Set year of data taking
    void setYearOfDataTaking(const Int_t& year = 2018)  { fYearOfDataTaking = {year}; }
    /// @brief Path to jetAnalysis directory (or any folder) that contains aux_files/... with JEC corrections
    void setPath2JetAnalysis(const Char_t *name = "/Users/gnigmat/work/cms/soft/jetAnalysis") { fJECPath = name; }
    /// @brief Set input JEC file name
    void setJECFileName(const Char_t *name = "Autumn18_HI_V8_MC_L2Relative_AK4PF.txt") { fJECInputFileName = name; }
    /// @brief Set residual JEC file name (for experimental data)
    void setJECFileDataName(const Char_t *name = "Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt") { fJECInputFileDataName = name; }
    /// @brief Apply jet pT-smearing
    void setJetPtSmearing(const Bool_t& smear = kFALSE) { fDoJetPtSmearing = {smear}; }
    /// @brief Set event cut
    void setEventCut(EventCut *cut) { fEventCut = {cut}; }
    /// @brief Set jet cut
    void setJetCut(JetCut *cut) { fJetCut = {cut}; }
    /// @brief Is the dataset from MC
    void setIsMc() { fIsMc = {kTRUE}; }
    /// @brief Correct centrality for MC events
    void setCorrectCentMC() { fCorrectCentMC = {kTRUE}; }
    /// @brief Fix jet arrays
    void fixJetArrays() { fFixJetArrays = {kTRUE}; }
    /// @brief Use extra JEC correction (in pPb 8160 it is needed for ak4cs)
    void useExtraJECCorr() { fUseExtraJEC = {kTRUE}; }

    /// @brief Return amount of events to read
    Long64_t nEventsTotal() const { return fEvents2Read; }

  private:

    /// @brief Setup input stream (either single file or a list of files)
    void setInputStream(const Char_t *inputStream) { fInFileName = {inputStream}; }
    /// @brief Instantiate extra JEC scale correction factor (for pPb 2016 for ak4cs)
    void createExtraJECScaleCorrFunction();

    /// Setup input all input
    void setupInput(TString input, TChain *hltChain, TChain *eveChain, TChain *partFlowChain, 
                    TChain *trkChain, Bool_t useMC, TChain *genTrkChain);
    /// Setup chains to be filled
    Int_t setupChains();
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

    /// @brief Calculate event weight
    /// @param isMC Is MC sample
    /// @param use_centrality If use centrality instead of multiplicity
    /// @param system Colliding system (pp, pPb or PbPb)
    /// @param year Year of the data taking
    /// @param energy Colliding energy (in GeV)
    /// @param vz Vertex Z position
    /// @param mult Charge particle multiplicity
    /// @param weighttree PtHat of the event (for MC only)
    /// @param leadjetpt Leading jet pT
    Float_t eventWeight(const Bool_t& isMC, const Bool_t& use_centrality, const std::string& system, 
                        const Int_t& year, const Int_t& energy, const Float_t& vz, 
                        const Int_t mult, const Float_t& weighttree, const Float_t& leadjetpt) const;
    /// @brief For compatibility between MC reco and data
    /// @param isMC 
    /// @param system 
    /// @param year 
    /// @param energy 
    /// @param jetpt 
    Float_t jetPtWeight(const Bool_t& isMC, const std::string& system, const Int_t& year, 
                        const Int_t& energy, float jetpt) const;
    /// @brief For compatibility between MC RECO and Data
    /// @param isMC Is MC sample
    /// @param system Colliding system (pp, pPb or PbPb)
    /// @param year Year of the data taking
    /// @param energy Colliding energy (in GeV)
    /// @param leadjetpt Leading jet pT
    Float_t leadJetPtWeight(const Bool_t& isMC, const std::string& system, const Int_t& year, 
                            const Int_t& energy, const Float_t& leadjetpt) const;
    /// @brief For compatibility between MC RECO and Data
    /// @param isMC Is MC sample
    /// @param system Colliding system (pp, pPb or PbPb)
    /// @param year Year of the data taking
    /// @param energy Colliding energy (in GeV)
    /// @param subleadjetpt Leading jet pT
    Float_t subleadJetPtWeight(const Bool_t& isMC, const std::string& system, const Int_t& year, 
                               const Int_t& energy, const Float_t& subleadjetpt);
    /// @brief Jet smearing resolution effect
    /// @param isMC True for MC and false for Data
    /// @param system Colliding system (pp, pPb or PbPb)
    /// @param year Year of the data taking
    /// @param energy Colliding energy (in GeV)
    /// @param jetpt Jet pT weight
    /// @param dosmearing Apply smearing
    /// @param resolutionfactor Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
    Float_t jetPtSmeringWeight(const Bool_t& isMC, const std::string& system, const Int_t& year, 
                               const Int_t& energy, const Float_t& jetpt, const Bool_t& dosmearing, 
                               const Float_t resolutionfactor) const;
    /// @brief Track mixing effect (Seagull)
    /// @param isMC True for MC and false for Data
    /// @param system Colliding system (pp, pPb or PbPb)
    /// @param year Year of the data taking
    /// @param energy Colliding energy (in GeV)
    /// @param trketa Track eta
    /// @param reco Is reco track
    Float_t trkEtaMixWeight(const Bool_t& isMC, const std::string& system, const Int_t& year, 
                            const Int_t& energy, const Float_t& trketa, const Bool_t& reco) const;

    /// @brief Calculate centrality weight
    Double_t evalCentralityWeight(const Double_t& hiBin);

    /// @brief Event with jets and other variables
    Event *fEvent;

    /// @brief Input filename (name.root) or file with list of ROOT files 
    const Char_t *fInFileName;
    /// @brief Number of events to process from input file(s)
    Long64_t fEvents2Read;
    /// @brief How many events were processed
    Long64_t fEventsProcessed;

    /// @brief Is file with MC information
    Bool_t fIsMc;
    /// @brief Use to apply hiBin shift for centrality correction
    Bool_t fCorrectCentMC;

    /// @brief Switch HLT branch ON
    Bool_t fUseHltBranch;
    /// @brief Switch skimming branch ON
    Bool_t fUseSkimmingBranch;

    /// @brief Switch particle flow jet branch ON
    Bool_t fUsePartFlowJetBranch;
    /// @brief Switch calorimeter-reconstructed jet branch ON
    Bool_t fUseCaloJetBranch;
    /// @brief Switch track branch ON
    Bool_t fUseTrackBranch;
    /// @brief Switch MC track branch ON
    Bool_t fUseGenTrackBranch;

    /// @brief Chain conaining HLT information (used to friend other trees)
    TChain *fHltTree;
    /// @brief Chain containing skimming information
    TChain *fSkimTree;
    /// @brief Chain containing event information
    TChain *fEventTree;
    /// @brief Chain containing jets from calorimeters
    TChain *fCaloJetTree;
    /// @brief Chain containing particle flow jets
    TChain *fPartFlowJetTree;
    /// @brief Chain containing tracks
    TChain *fTrkTree;
    /// @brief Chain containing Monte Carlo tracks
    TChain *fGenTrkTree;

    /// @brief Name of the particle flow tree (e.g. akCs4PFJetAnalyzer for PbPb or ak4PFJetAnalyzer for pp)
    TString fPFTreeName;
    /// @brief Name of the particle flow tree (e.g. akPu4CaloJetAnalyzer for PbPb)
    TString fCaloTreeName;

    //
    // Variables to store information from TTree
    //

    //
    // Event quantities
    //

    /// @brief Run ID
    UInt_t    fRunId;
    /// @brief Event ID
    ULong64_t fEventId;
    /// @brief Luminosity
    UInt_t    fLumi;
    /// @brief Vertex position z
    Float_t fVertexZ; 
    /// @brief Collision centrality (0-200)
    Int_t   fHiBin;
    /// @brief Event weight (prescale from generator level)
    Float_t fPtHatWeight; 
    /// @brief ptHat (initial parton pT) - from generator level
    Float_t fPtHat;

    //
    // Trigger and skimming information
    //

    Int_t fHLT_PAAK4CaloJet60_Eta5p1_v3;
    Int_t fHLT_PAAK4CaloJet80_Eta5p1_v3;
    Int_t fHLT_PAAK4CaloJet100_Eta5p1_v3;
    Int_t fHLT_PAAK4PFJet60_Eta5p1_v4;
    Int_t fHLT_PAAK4PFJet80_Eta5p1_v3;
    Int_t fHLT_PAAK4PFJet100_Eta5p1_v3;
    Int_t fHLT_PAAK4PFJet120_Eta5p1_v2;

    Int_t fHLT_HIAK4PFJet15_v1;
    Int_t fHLT_HIAK4PFJet15_v1_Prescl;
    Int_t fHLT_HIAK4PFJet30_v1;
    Int_t fHLT_HIAK4PFJet30_v1_Prescl;
    Int_t fHLT_HIAK4PFJet40_v1;
    Int_t fHLT_HIAK4PFJet40_v1_Prescl;
    Int_t fHLT_HIAK4PFJet60_v1;
    Int_t fHLT_HIAK4PFJet60_v1_Prescl;
    Int_t fHLT_HIAK4PFJet80_v1;
    Int_t fHLT_HIAK4PFJet80_v1_Prescl;
    Int_t fHLT_HIAK4PFJet120_v1;
    Int_t fHLT_HIAK4PFJet120_v1_Prescl;

    Int_t fHLT_HIAK8PFJet15_v1;
    Int_t fHLT_HIAK8PFJet15_v1_Prescl;
    Int_t fHLT_HIAK8PFJet25_v1;
    Int_t fHLT_HIAK8PFJet25_v1_Prescl;
    Int_t fHLT_HIAK8PFJet40_v1;
    Int_t fHLT_HIAK8PFJet40_v1_Prescl;
    Int_t fHLT_HIAK8PFJet60_v1;
    Int_t fHLT_HIAK8PFJet60_v1_Prescl;
    Int_t fHLT_HIAK8PFJet80_v1;
    Int_t fHLT_HIAK8PFJet80_v1_Prescl;
    Int_t fHLT_HIAK8PFJet140_v1;
    Int_t fHLT_HIAK8PFJet140_v1_Prescl;

    Int_t fHLT_HIPFJet25_v1;
    Int_t fHLT_HIPFJet25_v1_Prescl;
    Int_t fHLT_HIPFJet140_v1;
    Int_t fHLT_HIPFJet140_v1_Prescl;

    Int_t fHLT_HIPuAK4CaloJet80Eta5p1_v1;
    Int_t fHLT_HIPuAK4CaloJet100Eta5p1_v1;


    // Skimanalysis part
    Int_t fHBHENoiseFilterResultRun2Loose;
    Int_t fHBHENoiseFilterResultRun2Tight;
    Int_t fHBHEIsoNoiseFilterResult;
    Int_t fCollisionEventSelectionAODv2;
    Int_t fPhfCoincFilter2Th4;
    Int_t fPPAprimaryVertexFilter;
    Int_t fPBeamScrapingFilter;
    Int_t fPprimaryVertexFilter;
    Int_t fPVertexFilterCutG;
    Int_t fPVertexFilterCutGloose;
    Int_t fPVertexFilterCutGtight;
    Int_t fPVertexFilterCutE;
    Int_t fPVertexFilterCutEandG;
    Int_t fPClusterCompatibilityFilter;

    Int_t fPhfCoincFilter;
    Int_t fPVertexFilterCutdz1p0;
    Int_t fPVertexFilterCutGplus;
    Int_t fPVertexFilterCutVtx1;

    //
    // Jet information
    //

    /// @brief Number of reconstructed jets
    Int_t   fNPFRecoJets;
    /// @brief Reconstructed jet transverse momentum (without JEC)
    Float_t fPFRecoJetPt[100];
    /// @brief Pseudorapidity of reconstructed jet
    Float_t fPFRecoJetEta[100];
    /// @brief Azimuthal angle of reconstructed jet
    Float_t fPFRecoJetPhi[100];
    /// @brief WTA eta of reconstructed jet
    Float_t fPFRecoJetWTAEta[100];
    /// @brief WTA phi of reconstructed jet
    Float_t fPFRecoJetWTAPhi[100];
    /// @brief Track with maximum pT in reconstructed jet
    Float_t fPFRecoJetTrackMax[100];

    Float_t fPFRecoJtPfNHF[100];
    Float_t fPFRecoJtPfNEF[100];
    Float_t fPFRecoJtPfCHF[100];
    Float_t fPFRecoJtPfMUF[100];
    Float_t fPFRecoJtPfCEF[100];
    Int_t fPFRecoJtPfCHM[100];
    Int_t fPFRecoJtPfCEM[100];
    Int_t fPFRecoJtPfNHM[100];
    Int_t fPFRecoJtPfNEM[100];
    Int_t fPFRecoJtPfMUM[100];

    /// @brief Transverse momentum of generated jet that was matched with reconstructed jet
    Float_t fPFRefJetPt[100];
    /// @brief /// @brief Pseudorapidity of generated jet that was matched with reconstructed jet
    Float_t fPFRefJetEta[100];
    /// @brief Azimuthal angle of generated jet that was matched with reconstructed jet
    Float_t fPFRefJetPhi[100];
    /// @brief WTA eta of generated jet that was matched with reconstructed jet
    Float_t fPFRefJetWTAEta[100];
    /// @brief WTA phi of generated jet that was matched with reconstructed jet
    Float_t fPFRefJetWTAPhi[100];
    /// @brief Parton flavor of generated jet that was matched with reconstructed jet
    Int_t   fPFRefJetPartonFlavor[100];
    /// @brief Parton flavor for B of generated jet that was matched with reconstructed jet
    Int_t   fPFRefJetPartonFlavorForB[100];

    /// @brief Number of generated jets
    Int_t   fNPFGenJets;
    /// @brief Generated jet transverse momentum
    Float_t fPFGenJetPt[100];
    /// @brief Pseudorapidity of generated jet
    Float_t fPFGenJetEta[100];
    /// @brief Azimuthal angle of generated jet
    Float_t fPFGenJetPhi[100];
    /// @brief WTA eta of generated jet
    Float_t fPFGenJetWTAEta[100];
    /// @brief WTA phi of generated jet
    Float_t fPFGenJetWTAPhi[100];

    /// @brief Number of reconstructed jets
    Int_t   fNCaloRecoJets;
    /// @brief Reconstructed jet transverse momentum (without JEC)
    Float_t fCaloRecoJetPt[100];
    /// @brief Pseudorapidity of reconstructed jet
    Float_t fCaloRecoJetEta[100];
    /// @brief Azimuthal angle of reconstructed jet
    Float_t fCaloRecoJetPhi[100];
    /// @brief WTA eta of reconstructed jet
    Float_t fCaloRecoJetWTAEta[100];
    /// @brief WTA eta of reconstructed jet
    Float_t fCaloRecoJetWTAPhi[100];
    /// @brief Track with maximum pT in reconstructed jet
    Float_t fCaloRecoJetTrackMax[100];

    /// @brief Transverse momentum of generated jet that was matched with reconstructed jet
    Float_t fCaloRefJetPt[100];
    /// @brief /// @brief Pseudorapidity of generated jet that was matched with reconstructed jet
    Float_t fCaloRefJetEta[100];
    /// @brief Azimuthal angle of generated jet that was matched with reconstructed jet
    Float_t fCaloRefJetPhi[100];
    /// @brief WTA eta of generated jet that was matched with reconstructed jet
    Float_t fCaloRefJetWTAEta[100];
    /// @brief WTA eta of generated jet that was matched with reconstructed jet
    Float_t fCaloRefJetWTAPhi[100];
    /// @brief Parton flavor of generated jet that was matched with reconstructed jet
    Int_t   fCaloRefJetPartonFlavor[100];
    /// @brief Parton flavor for B of generated jet that was matched with reconstructed jet
    Int_t   fCaloRefJetPartonFlavorForB[100];

    /// @brief Number of generated jets
    Int_t   fNCaloGenJets;
    /// @brief Generated jet transverse momentum
    Float_t fCaloGenJetPt[100];
    /// @brief Pseudorapidity of generated jet
    Float_t fCaloGenJetEta[100];
    /// @brief Azimuthal angle of generated jet
    Float_t fCaloGenJetPhi[100];
    /// @brief WTA eta of generated jet
    Float_t fCaloGenJetWTAEta[100];
    /// @brief WTA phi of generated jet
    Float_t fCaloGenJetWTAPhi[100];

    //
    // Reconstructed tracks
    //

    /// @brief Number of tracks
    Int_t   fNTracks;
    /// @brief Track transverse momentum
    Float_t fTrackPt[9999];
    /// @brief Track pseudorapidity
    Float_t fTrackEta[9999];
    /// @brief Track azimuthal angle
    Float_t fTrackPhi[9999];
    /// @brief Track pT error (uncertainty)
    Float_t fTrackPtErr[9999];
    /// @brief Track distance of closest approach in transverse plane (XY)
    Float_t fTrackDcaXY[9999];
    /// @brief Track distance of closest approach in beam direction (z)
    Float_t fTrackDcaZ[9999];
    /// @brief Track distance of closest approach error in transverse plane (XY)
    Float_t fTrackDcaXYErr[9999];
    /// @brief Track distance of closest approach error in beam direction (z)
    Float_t fTrackDcaZErr[9999];
    /// @brief Track fitting (reconstruction) chi2
    Float_t fTrackChi2[9999];
    /// @brief Track number of degrees of freedom in the fitting 
    UChar_t fTrackNDOF[9999];    
    /// @brief Particle flow energy deposited in ECAL from the given track
    Float_t fTrackPartFlowEcal[9999];
    /// @brief Particle flow energy deposited in HCAL from the given track
    Float_t fTrackPartFlowHcal[9999];
    /// @brief Track MVA for each step
    Float_t fTrackMVA[9999];
    /// @brief Track algorithm/step
    UChar_t fTrackAlgo[9999];
    /// @brief Track charge
    Int_t   fTrackCharge[9999];
    /// @brief Number of hits in the tracker
    UChar_t fTrackNHits[9999];
    /// @brief Number of layers with measurement in the tracker
    UChar_t fTrackNLayers[9999];
    /// @brief Tracker steps MVA selection
    Bool_t  fTrackHighPurity[9999];

    //
    // Monte Carlo tracks
    //

    /// @brief Generated particle transverse momentum
    std::vector<Float_t> fGenTrackPt;
    /// @brief Generated particle pseudorapidity
    std::vector<Float_t> fGenTrackEta;
    /// @brief Generated particle azimuthal angle
    std::vector<Float_t> fGenTrackPhi;
    /// @brief Generated particle charge
    std::vector<Int_t>   fGenTrackCharge;
    /// @brief Generated particle PID
    std::vector<Int_t>   fGenTrackPid;
    /// @brief Generated particle sube (?)
    std::vector<Int_t>   fGenTrackSube;

    /// @brief Jet Energy Corrector instance
    JetCorrector *fJEC;
    /// @brief List of files with JEC
    std::vector< std::string > fJECFiles;
    /// @brief Path to jetAnalysis directory
    TString fJECPath;
    /// @brief Input file with JEC corrections
    TString fJECInputFileName;
    /// @brief Input file for residual corrections (first apply fJECInputFileName and 
    ///        the residual (for experimental data only))
    TString fJECInputFileDataName;
    /// @brief Jet Energy Uncertainty instance
    JetUncertainty *fJEU;
    /// @brief List of files with JEU
    std::vector< std::string > fJEUFiles;

    /// @brief Colliding system: pp, pPb or PbPb
    TString fCollidingSystem;
    /// @brief Colliding energy
    Int_t fCollidingEnergyGeV;
    /// @brief Year of data taking
    Int_t fYearOfDataTaking;
    /// @brief Apply jet pT-smearing 
    Bool_t fDoJetPtSmearing;

    /// @brief Fix indices
    Bool_t fFixJetArrays;

    /// @brief Event cut
    EventCut *fEventCut;
    /// @brief Jet cut
    JetCut *fJetCut;

    /// @brief Vector that contains indices of generated jets that matched to the reconsructed 
    /// particle flow jet (should be of the reco/red size)
    std::vector<Int_t> fRecoPFJet2GenJetId;
    /// @brief Vector that contains indices of the reconstructed particle flow jets that 
    /// macthed to generated jet
    std::vector<Int_t> fGenJet2RecoPFJet;

    /// @brief Vector that contains indices of generated jets that matched to the reconsructed 
    /// calorimeter jet (should be of the reco/red size)
    std::vector<Int_t> fRecoCaloJet2GenJetId;
    /// @brief Vector that contains indices of the reconstructed calorimeter jets that 
    /// macthed to generated jet
    std::vector<Int_t> fGenJet2RecoCaloJet;
    /// @brief Use extra correction for JEC (for 8160 pPb year 2016 ak4pf is default. extra needed for ak4cs)
    Bool_t  fUseExtraJEC;
    /// @brief JEC extra correction
    TF1    *fJECScaleCorr;

    ClassDef(ForestAODReader, 1)
};

#endif // #define ForestAODReader_h