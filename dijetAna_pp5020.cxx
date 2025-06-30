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

//________________
void usage() {
    std::cout << "./programName inputFileList oFileName isMc isPbGoingDir ptHatLow ptHatHi jeuSyst jerSyst triggerId recoJetSelMethod" << std::endl;
    std::cout << "isMc: 1 (embedding), 0 (data)" << std::endl;
    std::cout << "isPbGoingDir: 1 (Pb-going), 0 (p-going)" << std::endl;
    std::cout << "ptHatLow: Low ptHat cut (for embedding)" << std::endl;
    std::cout << "ptHatHi: High ptHat cut (for embedding)" << std::endl;
    std::cout << "jeuSyst: 0 (default), 1 (JEU+), -1 (JEU-)" << std::endl;
    std::cout << "jerSyst: 0 (default), 1 (JER+), -1 (JER-), other - only JEC is applied" << std::endl;
    std::cout << "triggerId: 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100" << std::endl;
    std::cout << "recoJetSelMethod: 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId" << std::endl;
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
    //eventCut->usePhfCoincFilter();

    // Default cut
    // eventCut->usePVertexFilterCutdz1p0();
    // Pile-up systematics
    // eventCut->usePVertexFilterCutGplus();
    // Pile-up systematics
    //eventCut->usePVertexFilterCutVtx1();
    
    // Triggers from Vipul
    if ( !isMc ) {
        eventCut->useHLT_HIAK4PFJet60_v1();
        // eventCut->useHLT_HIAK4PFJet80_v1();

        //eventCut->useHLT_HIAK4CaloJet60_v1();
        // eventCut->useHLT_HIAK4CaloJet80_v1();
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
JetCut *createRecoJetCut(int collEnergyGeV, int recoJetSelMethod) {
    JetCut *jetCut = new JetCut{};
    jetCut->setPt(1.0f, static_cast<float>(collEnergyGeV));
    jetCut->setEtaLab(-5.2f, 5.2f);
    jetCut->setEtaCM(-5.2f, 5.2f);
    jetCut->setSelectionMethod(recoJetSelMethod); // 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId
    jetCut->setLooseJetIdCut(true); // true = loose, false = tight
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
ForestAODReader *createForestAODReader(const TString &inFileName, const bool &isMc, 
    const bool &isCentWeightCalc, const bool &isPbGoingDir, const TString &recoJetBranchName, 
    const TString &collisionSystemName, const int &collisionSystem, const int &collEnergyGeV, 
    const int &collYear, const float& etaShift, const TString &path2JEC, const TString &JECFileName, 
    const TString &JECFileDataName, const TString &JEUFileName, const int &useJEUSyst, 
    const int &useJERSyst, EventCut *eventCut = nullptr, JetCut *jetCut = nullptr) {

    // Create ForestAODReader object
    ForestAODReader *forestReader = new ForestAODReader{inFileName};

    // Define collision system parameters
    forestReader->setCollisionSystemName( collisionSystemName.Data() );
    forestReader->setCollisionSystem( collisionSystem );
    forestReader->setCollisionEnergyInGeV( collEnergyGeV );
    forestReader->setYearOfDataTaking( collYear );
    forestReader->setPbGoingDir( isPbGoingDir );
    // Set eta shift for the analysis
    forestReader->setEtaShift( etaShift );

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
        forestReader->setJERFitParams(0.0415552, 0.960013); // in |eta|<1.6
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
    if ( isMc ) {
        analysis->setIsMc();
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

    std::cout << "Starting dijetAna_pp5020 program" << std::endl;

    // Set default values for arguments
    bool isMc{true};
    bool isCentWeightCalc{false};
    bool isPbGoingDir{};
    TString inFileName{};
    int   collEnergyGeV{5020}; // Keep this value intentionally to use the corrections
    int   collisionSystem{0};  // 0 - pp, 1 -pPb, 2 - PbPb 
    TString collisionSystemName;
    if (collisionSystem == 0) {
        collisionSystemName = "pp";
    } 
    else if (collisionSystem == 1) {
        collisionSystemName = "pPb";
    } 
    else if (collisionSystem == 2) {
        collisionSystemName = "PbPb";
    }
    else {
        collisionSystemName = "pPb";
    }
    int   collYear{2016};
    // TString recoJetBranchName{"akCs4PFJetAnalyzer"};
    TString recoJetBranchName{"ak4PFJetAnalyzer"};
    TString oFileName{};
    TString JECFileName;
    TString JECFileDataName;
    TString JEUFileName;
    TString path2JEC = "..";
    float ptHatCut[2] {-100000000, 100000000};
    int   useJEUSyst{0};       // 0-default, 1-JEU+, -1-JEU-
    int   useJERSyst{-99};     // 0-default, 1-JER+, -1-JER-, other - not use extra JER smearing
    float etaShift = 0.465;
    int   triggerId{0};        // 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100
    int   recoJetSelMethod{1}; // 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId

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
    // recoJetSelMethod               - 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetId

    // Read input argument list 
    if (argc <= 1) {
        std::cout << "Too few arguments passed. Terminating" << std::endl;
        usage();
		return -1;
    }
    else {
        inFileName   = argv[1];
        oFileName    = argv[2];
        isMc         = atoi(argv[3]);
        isPbGoingDir = atoi(argv[4]);
        ptHatCut[0]  = atoi(argv[5]);
        ptHatCut[1]  = atoi(argv[6]);
        useJEUSyst   = atoi( argv[7] );
        useJERSyst   = atoi( argv[8] );
        if (argc <= 9 ) {
            triggerId = 0;
        }
        else {
            triggerId = atoi( argv[9] );
        }
        if (argc <= 10 ) {
            recoJetSelMethod = 1; // Default is trkMaxPt/RawPt
        }
        else {
            recoJetSelMethod = atoi( argv[10] );
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
              << "Reco jet selection method              : " << recoJetSelMethod << std::endl
              << std::endl;

    if (isMc) {
        JECFileName = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt";
    }
    else {
        JECFileName = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt";
        JECFileDataName = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt";
        JEUFileName = "Spring18_ppRef5TeV_V6_DATA_Uncertainty_AK4PF.txt";
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
    JetCut *recoJetCut = createRecoJetCut(collEnergyGeV, recoJetSelMethod);

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
                                                    recoJetBranchName, collisionSystemName, collisionSystem, collEnergyGeV, 
                                                    collYear, etaShift, path2JEC, JECFileName, JECFileDataName, 
                                                    JEUFileName, useJEUSyst, useJERSyst, eventCut, nullptr);

    // Pass reader to the manager
    manager->setEventReader(reader);

    //
    // Initialize analysis
    //
    DiJetAnalysis *analysis = createDiJetAnalysis(collisionSystem, collEnergyGeV, isMc, isPbGoingDir, 
                                                  ptHatCut, recoJetCut, genJetCut, dijetCut, etaShift);

    
    // Initialize histogram manager
    HistoManagerDiJet *hm = new HistoManagerDiJet{};
    hm->setIsMc( isMc );
    hm->init(); // true stands up for use MC; need to FIX

    // Add histogram manager to analysis
    analysis->addHistoManager(hm);

    // Add analysis to manager
    manager->addAnalysis(analysis);

    // Run chain of analyses
    manager->init();

    // Important for embedding reweightening
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
