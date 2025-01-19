// C++ headers
#include <iostream>

// Jet analysis headers
#include "Manager.h"
#include "ForestAODReader.h"
#include "DiJetAnalysis.h"
#include "HistoManagerDiJet.h"
#include "EventCut.h"
#include "JetCut.h"

// ROOT headers
#include "TFile.h"
#include "TMath.h"
#include "TString.h"

//________________
void usage() {
    std::cout << "./programName inputFileList oFileName isMc isPbGoingDir ptHatLow ptHatHi jeuSyst jerSyst" << std::endl;
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
    int   collEnergyGeV{5020}; // Keep this value intentionally to use the corrections
    int   collisionSystem{0};  // 0 - pp, 1 -pPb, 2 - PbPb 
    TString collisionSystemName{"pp"};
    int   collYear{2016};
    // TString recoJetBranchName{"akCs4PFJetAnalyzer"};
    TString recoJetBranchName{"ak4PFJetAnalyzer"};
    TString oFileName{};
    TString JECFileName;
    TString JECFileDataName;
    TString JEUFileName;
    TString path2JEC = "..";
    double ptHatCut[2] {-100000000, 100000000};
    int   useJEUSyst{0};  // 0-default, 1-JEU+, -1-JEU-
    int   useJERSyst{-99};  // 0-default, 1-JER+, -1-JER-, other - not use extra JER smearing
    double etaShift = 0.465;

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
              << std::endl;

    if (isMc) {
        JECFileName = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt";
    }
    else {
        JECFileName = "Spring18_ppRef5TeV_V6_MC_L2Relative_AK4PF.txt";
        JECFileDataName = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt";
        JEUFileName = "Spring18_ppRef5TeV_V6_DATA_Uncertainty_AK4PF.txt";
    } // else

    // Initialize package manager
    Manager *manager = new Manager{};

    // Initialize event cut
    EventCut *eventCut = new EventCut{};
    eventCut->setVz(-15., 15.);
    // Skim
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
    //eventCut->setVerbose();

    // Initialize jet cut 
    JetCut *jetCut = new JetCut{};
    //jetCut->setMustHaveGenMathing();
    jetCut->setPt(20., 1500.);
    jetCut->setEta(-5.1, 5.1);
    //jetCut->setVerbose();

    // Initialize event reader
    ForestAODReader *reader = new ForestAODReader(inFileName);
    if (isMc) {
        // If is MC
        reader->setIsMc();
        if ( isCentWeightCalc ) {
            // Apply hiBin shift and centrality weight calculation
            reader->setCorrectCentMC();
        }
    }
    reader->useHltBranch();
    reader->useSkimmingBranch();
    reader->useRecoJetBranch();
    reader->setRecoJetBranchName( recoJetBranchName.Data() );

    if ( recoJetBranchName.CompareTo("akcs4pfjetanalyzer", TString::kIgnoreCase) == 0 ) {
        std::cout << "Extra correction will be used for JEC" << std::endl;
        reader->useExtraJECCorrForConstSubtraction();
    }

    reader->setCollidingSystem( collisionSystemName.Data() );
    reader->setCollidingEnergy( collEnergyGeV ) ;
    reader->setYearOfDataTaking( collYear );
    reader->setEventCut(eventCut);
    reader->fixJetArrays();

    // Set path to jet analysis (then will automatically add path to aux_files)
    reader->setPath2JetAnalysis( path2JEC.Data() );
    reader->addJECFile( JECFileName.Data() ); // Default JEC file for both data and MC
    if ( !isMc ) {
        reader->addJECFile( JECFileDataName.Data() );
        reader->setUseJEU( useJEUSyst );
        reader->setJEUFileName( JEUFileName );
    }
    if ( isMc ) {
        reader->useJERSystematics( useJERSyst ); // 0-default, 1-JER+, -1-JER-, other - not use (to check JES and JER)
        reader->setJERFitParams(0.0415552, 0.960013);
        reader->setJERSystParams();
    }

    //reader->setVerbose();

    // Pass reader to the manager
    manager->setEventReader(reader);

    // Initialize analysis
    DiJetAnalysis *analysis = new DiJetAnalysis{};
    analysis->setCollisionSystem( collisionSystem );
    analysis->setCollisionEnergyInGeV( collEnergyGeV );
    analysis->setIsMc(isMc);
    if (isMc) {
        analysis->setPtHatRange(ptHatCut[0], ptHatCut[1]);
    }
    analysis->setEtaShift( etaShift );
    analysis->setLeadJetPtLow( 30. );
    analysis->setSubLeadJetPtLow( 20. );
    analysis->setJetEtaLabRange( -3., 3. );
    analysis->setJetEtaCMRange( -2.5, 2.5 );
    analysis->setDijetPhiCut( 2. * TMath::Pi() / 3 );
    if ( isMc ) {
        analysis->setUseMcReweighting(0); // 0 - no reweighting, 1 - reweight to MB, 2 - reweight to Jet60, 3 - reweight to Jet80, 4 - reweight to Jet100
    }
    //analysis->setVerbose();
    
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
