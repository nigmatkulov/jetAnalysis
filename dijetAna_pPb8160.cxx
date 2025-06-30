// C++ headers
#include <iostream>

// Jet analysis headers
#include "Manager.h"
#include "ForestAODReader.h"
#include "DiJetAnalysis.h"
#include "HistoManagerDiJet.h"
#include "EventCut.h"
#include "JetCut.h"
#include "DiJetCut.h"

// ROOT headers
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "dijetAna_pPb8160.h"

//________________
void usage() {
    std::cout << "./programName inputFileList oFileName isMc isPbGoingDir ptHatLow ptHatHi jeuSyst jerSyst triggerId" << std::endl;
    std::cout << "isMc: 1 (embedding), 0 (data)" << std::endl;
    std::cout << "isPbGoingDir: 1 (Pb-going), 0 (p-going)" << std::endl;
    std::cout << "ptHatLow: Low ptHat cut (for embedding)" << std::endl;
    std::cout << "ptHatHi: High ptHat cut (for embedding)" << std::endl;
    std::cout << "jeuSyst: 0 (default), 1 (JEU+), -1 (JEU-)" << std::endl;
    std::cout << "jerSyst: 0 (default), 1 (JER+), -1 (JER-), other - only JEC is applied" << std::endl;
    std::cout << "triggerId: 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100" << std::endl;
}

//________________
EventCut *createEventCut(const bool &isMc, const int &triggerId, const float *ptHatCut) {
    //
    // Create event cut object
    //
    EventCut *eventCut = new EventCut{};

    //
    // Set vertex position limits
    //
    eventCut->setVz(-15.0f, 15.0f);

    //
    // Skim
    //
    eventCut->usePBeamScrapingFilter();
    eventCut->usePPAprimaryVertexFilter();
    eventCut->useHBHENoiseFilterResultRun2Loose();
    eventCut->usePhfCoincFilter();

    eventCut->usePVertexFilterCutdz1p0();    // Default cut
    // eventCut->usePVertexFilterCutGplus(); // Pile-up systematics
    // Pile-up systematics
    //eventCut->usePVertexFilterCutVtx1();

    //
    // Trigger
    //
    if ( !isMc ) {
        if ( triggerId == 1 ) {            
            eventCut->useHLT_PAAK4PFJet60_Eta5p1_v4();
        }
        else if ( triggerId == 2 ) {
            eventCut->useHLT_PAAK4PFJet80_Eta5p1_v3();
        }
        else if ( triggerId == 3 ) {
            eventCut->useHLT_PAAK4PFJet100_Eta5p1_v3();
        }
        else {
            // No trigger = MB trigger
        }
    }

    // Set ptHat cut for embedding
    if ( isMc ) {
        eventCut->setPtHat(ptHatCut[0], ptHatCut[1]);
    }

    //
    // Select specific run IDs
    //
    // eventCut->addRunIdToSelect( 285480 ); // PU 0.04
    // eventCut->addRunIdToSelect( 285505 ); // PU 0.25
    // eventCut->addRunIdToSelect( 285517 ); // PU 0.1
    // eventCut->addRunIdToSelect( 285832 ); // PU 0.004
    // eventCut->addRunIdToSelect( 285993 ); // PU 0.2
    
    //
    // Set verbose mode
    //
    //eventCut->setVerbose();

    return eventCut;
}

//________________
JetCut *createRecoJetCut(int collEnergyGeV) {
    JetCut *jetCut = new JetCut{};
    jetCut->setPt(1.0f, static_cast<float>(collEnergyGeV));
    jetCut->setEtaLab(-5.2f, 5.2f);
    jetCut->setEtaCM(-5.2f, 5.2f);
    jetCut->setSelectionMethod(1); // 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId
    // jetCut->setVerbose();
    return jetCut;
}

//________________
JetCut *createGenJetCut(int collEnergyGeV) {
    JetCut *jetCut = new JetCut{};
    jetCut->setPt(1.0f, static_cast<float>(collEnergyGeV));
    jetCut->setEtaLab(-5.2f, 5.2f);
    jetCut->setEtaCM(-5.2f, 5.2f);
    // jetCut->setVerbose();
    return jetCut;
}

//________________
DiJetCut *createDiJetCut() {
    DiJetCut *dijetCut = new DiJetCut{};
    
    dijetCut->setLeadJetPtMinimum(50.0f);
    dijetCut->setLeadJetEtaLab(-2.5f, 2.5f);
    dijetCut->setLeadJetEtaCM(-2.0f, 2.0f);

    dijetCut->setSubLeadJetPtMinimum(40.0f);
    dijetCut->setSubLeadJetEtaLab(-2.5f, 2.5f);
    dijetCut->setSubLeadJetEtaCM(-2.0f, 2.0f);

    dijetCut->setDijetDPhi(TMath::TwoPi() / 3); // 120 degrees

    // dijetCut->setVerbose();

    return dijetCut;
}

//________________
ForestAODReader *createForestAODReader(const TString &inFileName, const bool &isMc, const bool &isCentWeightCalc, const bool &isPbGoingDir, 
                                       const TString &recoJetBranchName, const TString &collisionSystemName, const int &collEnergyGeV, 
                                       const int &collYear, const TString &path2JEC, const TString &JECFileName, const TString &JECFileDataName, 
                                       const TString &JEUFileName, const int &useJEUSyst, const int &useJERSyst, 
                                       EventCut *eventCut = nullptr, JetCut *jetCut = nullptr) {

    // Create ForestAODReader object
    ForestAODReader *forestReader = new ForestAODReader{inFileName};

    // Define collision system parameters
    forestReader->setCollidingSystem( collisionSystemName.Data() );
    forestReader->setCollidingEnergy( collEnergyGeV ) ;
    forestReader->setYearOfDataTaking( collYear );
    forestReader->setPbGoingDir( isPbGoingDir );

    // Set isMc flag anc centrality correction flag
    if (isMc) {
        // If is MC
        forestReader->setIsMc();
        if ( isCentWeightCalc ) {
            // Apply hiBin shift and centrality weight calculation
            forestReader->setCorrectCentMC();
        }
    }

    // Set branches to read
    forestReader->useHltBranch();
    forestReader->useSkimmingBranch();
    forestReader->useRecoJetBranch();
    forestReader->setRecoJetBranchName( recoJetBranchName.Data() );
    if ( recoJetBranchName.CompareTo("akcs4pfjetanalyzer", TString::kIgnoreCase) == 0 ) {
        std::cout << "Extra correction will be used for JEC" << std::endl;
        forestReader->useExtraJECCorrForConstSubtraction();
    }

    // Perform jet manual jet matching
    forestReader->fixJetArrays();

    // Set path to jet analysis directory (then will automatically add path to aux_files)
    forestReader->setPath2JetAnalysis( path2JEC.Data() );
    forestReader->addJECFile( JECFileName.Data() );
    if ( !isMc ) {
        forestReader->setUseJEU( useJEUSyst );
        forestReader->addJECFile( JECFileDataName.Data() );
        forestReader->setJEUFileName( JEUFileName );
    }
    if ( isMc ) {
        forestReader->useJERSystematics( useJERSyst ); // 0-default, 1-JER+, -1-JER-, other - not use
        if ( isPbGoingDir ) {
            forestReader->setJERFitParams(0.0018, 0.9352); // in |eta|<1.6
        }
        else {
            forestReader->setJERFitParams(0.0018, 0.9352); // in |eta|<1.6
        }
        forestReader->setJERSystParams();
    }
    // If want to use manual JEC
    // forestReader->setUseManualJEC();


    // Set event cut
    if ( eventCut ) forestReader->setEventCut(eventCut);
    if ( jetCut ) forestReader->setJetCut(jetCut);

    // Set verbose mode
    // forestReader->setVerbose();

    return forestReader;
}

//________________
DiJetAnalysis *createDiJetAnalysis(const int &collisionSystem, const int &collEnergyGeV, const bool &isMc, 
                                   const bool &isPbGoingDir, const float *ptHatCut, 
                                   JetCut *recoJetCut, JetCut *genJetCut, DiJetCut *dijetCut, 
                                   const float &etaShift) {

    // Create DiJetAnalysis object
    DiJetAnalysis *analysis = new DiJetAnalysis{};

    analysis->setCollisionSystem( collisionSystem );
    analysis->setCollisionEnergyInGeV( collEnergyGeV );
    analysis->setIsMc( isMc );
    if ( isMc ) {
        analysis->setPtHatRange(ptHatCut[0], ptHatCut[1]);
    }
    if ( isPbGoingDir ) {
        analysis->setPbGoing();
    }
    analysis->setRecoJetCut( recoJetCut );
    analysis->setGenJetCut( genJetCut );
    analysis->setDiJetCut( dijetCut );
    analysis->setEtaShift( etaShift );

    if ( isMc ) {
        analysis->setUseMcReweighting(0); // 0 - no reweighting, 1 - reweight to MB, 2 - reweight to Jet60, 3 - reweight to Jet80, 4 - reweight to Jet100
    }

    //analysis->setVerbose();

    return analysis;
}

//________________
/// @brief The prorgram that launches the physics analysis
/// @param argc Number of arguments
/// @param argv Argument list
/// @return 0 in case of OKAY
int main(int argc, char const *argv[]) {

    // Set default values for arguments
    bool isMc{true};
    bool isCentWeightCalc{false};
    bool isPbGoingDir{};
    TString inFileName{};
    int   collEnergyGeV{8160};
    int   collisionSystem{1}; // 0 - pp, 1 -pPb, 2 - PbPb 
    TString collisionSystemName{"pPb"};
    int   collYear{2016};
    //TString recoJetBranchName{"akCs4PFJetAnalyzer"};
    TString recoJetBranchName{"ak4PFJetAnalyzer"};
    TString oFileName{};
    TString JECFileName;
    TString JECFileDataName;
    TString JEUFileName;
    TString path2JEC = "..";
    // For debugging purposes (using VS Code)
    // TString path2JEC = "/Users/gnigmat/work/cms/soft/jetAnalysis";
    float ptHatCut[2] {15., 30.};
    int   useJEUSyst{0};     // 0-default, 1-JEU+, -1-JEU-
    int   useJERSyst{-99};   // 0-default, 1-JER+, -1-JER-, other - only JEC is applied
    float etaShift = 0.465;
    int   triggerId{0};     // 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100

    // Sequence of command line arguments:
    //
    // inputFileList (or forest.root) - input file list with forest
    // outputFileName.root            - output file name
    // isMc                           - 1 (embedding), 0 (data)
    // isPbGoingDir                   - 1 (Pb-going), 0 (p-going)
    // ptHatLow                       - Low ptHat cut (for embedding)
    // ptHatHi                        - High ptHat cut (for embedding)
    // useJEUSyst                     - 0 (default), 1 (JEU+), -1 (JEU-)
    // useJERSyst                     - 0 (default), 1 (JER+), -1 (JER-)
    // triggerId                      - 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100

    // Read input argument list 
    if (argc <= 1) {
        std::cout << "Too few arguments passed. Terminating" << std::endl;
        usage();
		return -1;
    }
    else {
        std::cout << "Number of arguments passed: " << argc << std::endl;
        inFileName   = argv[1];
        oFileName    = argv[2];
        isMc         = atoi( argv[3] );
        isPbGoingDir = atoi( argv[4] );
        ptHatCut[0]  = atoi( argv[5] );
        ptHatCut[1]  = atoi( argv[6] );
        useJEUSyst   = atoi( argv[7] );
        useJERSyst   = atoi( argv[8] );
        if (argc <= 9 ) {
            triggerId = 0;
        }
        else {
            triggerId = atoi( argv[9] );
        }
    }

    std::cout << "Arguments passed:\n"
              << "Input file name                        : " << inFileName << std::endl
              << "Output file name                       : " << oFileName << std::endl
              << "Is MC                                  : " << isMc << std::endl
              << "Collision system (0-pp, 1-pPb, 2-PbPb) : " << collisionSystem << std::endl
              << "Collision system name                  : " << collisionSystemName.Data() << std::endl
              << "Is Pb-going direction                  : " << isPbGoingDir << std::endl
              << "ptHat range                            : " << ptHatCut[0] << "-" << ptHatCut[1] << std::endl
              << "Use centrality weight                  : " << isCentWeightCalc << std::endl
              << "Use JEU systematics                    : " << useJEUSyst << std::endl
              << "Use JER systematics                    : " << useJERSyst << std::endl
              << "Trigger ID                             : " << triggerId << std::endl
              << std::endl;

    if (isMc) {
        if (isPbGoingDir) {
            // PYTHIA+EPOS
            JECFileName = "Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt";
            // PYTHIA
            // JECFileName = "Autumn16_HI_pPb_Pbgoing_Unembedded_MC_L2Relative_AK4PF.txt";
        }
        else {
            // PYTHIA+EPOS
            JECFileName = "Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt";
            // PYTHIA
            // JECFileName = "Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF";
        }
    }
    else {
        if (isPbGoingDir) { // Remember to flip to p-going for data
            JECFileName = "Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt";
            
        }
        else {
            JECFileName = "Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt";
        }
        JECFileDataName = "Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt";
        JEUFileName = "Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF.txt";
    } // else

    //
    // Initialize package manager
    //
    Manager *manager = new Manager{};

    //
    // Initialize event cut
    //
    EventCut *eventCut = createEventCut(isMc, triggerId, ptHatCut);

    //
    // Initialize reco jet cut 
    //
    JetCut *recoJetCut = createRecoJetCut(collEnergyGeV);

    //
    // Initialize gen jet cut
    //
    JetCut *genJetCut = createGenJetCut(collEnergyGeV);

    //
    // Initialize dijet cut
    //
    DiJetCut *dijetCut = createDiJetCut();


    //
    // Initialize event reader
    //
    ForestAODReader *reader = createForestAODReader(inFileName, isMc, isCentWeightCalc, isPbGoingDir, 
                                                    recoJetBranchName, collisionSystemName, collEnergyGeV, 
                                                    collYear, path2JEC, JECFileName, JECFileDataName, 
                                                    JEUFileName, useJEUSyst, useJERSyst, eventCut);

    // Pass reader to the manager
    manager->setEventReader(reader);

    //
    // Initialize analysis
    //
    DiJetAnalysis *analysis = createDiJetAnalysis(collisionSystem, collEnergyGeV, isMc, isPbGoingDir, 
                                                  ptHatCut, recoJetCut, genJetCut, dijetCut, etaShift);
    
    //
    // Initialize histogram manager
    //
    HistoManagerDiJet *hm = new HistoManagerDiJet{};
    hm->setIsMc( isMc );
    hm->init();

    // Add histogram manager to analysis
    analysis->addHistoManager( hm );

    // Add analysis to manager
    manager->addAnalysis( analysis );

    // Run chain of analyses
    manager->init();

    // Important for pPb8160 embedding reweightening
    if ( isMc ) {
        analysis->setNEventsInSample( reader->nEventsTotal() );
    }

    manager->performAnalysis();
    manager->finish();

    // Create output file and store results of calculations
    TFile* oFile = new TFile(oFileName, "recreate");
    hm->writeOutput();
    oFile->Close();

    return 0;
}
