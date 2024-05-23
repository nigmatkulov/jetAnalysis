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

//________________
void usage() {
    std::cout << "./programName inputFileList oFileName isMc isPbGoingDir ptHatLow ptHatHi" << std::endl;
}

//________________
/// @brief The prorgram that launches the physics analysis
/// @param argc Number of arguments
/// @param argv Argument list
/// @return 0 in case of OKAY
int main(int argc, char const *argv[]) {

    // Set default values for arguments
    Bool_t isMc{kTRUE};
    Bool_t isCentWeightCalc{kFALSE};
    Bool_t isPbGoingDir{};
    TString inFileName{};
    Int_t   collEnergyGeV{8160};
    TString collSystem{"pPb"};
    Int_t   collYear{2016};
    TString pfBranchName{"akCs4PFJetAnalyzer"};
    //TString pfBranchName{"ak4PFJetAnalyzer"};
    TString oFileName{};
    TString JECFileName;
    TString JECFileDataName;
    TString path2JEC = "../";
    Double_t ptHatCut[2] {15., 30.};

    // Sequence of command line arguments:
    //
    // inputFileList (or forest.root) - input file list with forest
    // outputFileName.root            - output file name
    // isMc                           - 1 (embedding), 0 (data)
    // isPbGoingDir                   - 1 (Pb-going), 0 (p-going)
    // ptHatLow                       - Low ptHat cut (for embedding)
    // ptHatHi                        - High ptHat cut (for embedding)

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
    }

    std::cout << "Arguments passed:\n"
              << "Input file name       : " << inFileName << std::endl
              << "Output file name      : " << oFileName << std::endl
              << "Is MC                 : " << isMc << std::endl
              << "Is Pb-going direction : " << isPbGoingDir << std::endl
              << "ptHat range           : " << ptHatCut[0] << "-" << ptHatCut[1] << std::endl
              << "Use centrality weight : " << isCentWeightCalc << std::endl
              << std::endl;

    if (isMc) {
        if (isPbGoingDir) {
            JECFileName = "Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt";
        }
        else {
            JECFileName = "Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt";
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
    eventCut->usePhfCoincFilter();
    eventCut->usePVertexFilterCutdz1p0();
    // Trigger
    //eventCut->useHLT_PAAK4PFJet80_Eta5p1_v3();

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
    reader->usePartFlowJetBranch();
    reader->setPartFlowJetBranchName( pfBranchName.Data() );

    if ( pfBranchName.CompareTo("akcs4pfjetanalyzer", TString::kIgnoreCase) == 0 ) {
        std::cout << "Extra correction will be used for JEC" << std::endl;
        reader->useExtraJECCorr();
    }

    //reader->useCaloJetBranch();
    reader->setCollidingSystem( collSystem.Data() );
    reader->setCollidingEnergy( collEnergyGeV ) ;
    reader->setYearOfDataTaking( collYear );
    reader->setEventCut(eventCut);
    //reader->setJetCut(jetCut);
    reader->fixJetArrays();

    // Set path to jet analysis (then will automatically add path to aux_files)
    reader->setPath2JetAnalysis( path2JEC.Data() );
    reader->setJECFileName( JECFileName.Data() );
    if ( !isMc ) {
        reader->setJECFileDataName( JECFileDataName.Data() );
    }

    // Pass reader to the manager
    manager->setEventReader(reader);

    // Initialize analysis
    DiJetAnalysis *analysis = new DiJetAnalysis{};
    analysis->setIsMc(isMc);
    if (isMc) {
        analysis->setPtHatRange(ptHatCut[0], ptHatCut[1]);
    }
    analysis->setIsPPb();
    if ( isPbGoingDir ) {
        analysis->setPbGoing();
    }
    analysis->setEtaShift( 0.4654094531 );
    analysis->setLeadJetPtLow( 50. );
    analysis->setSubLeadJetPtLow( 30. );
    analysis->setDijetPhiCut( TMath::TwoPi() / 3 );
    //analysis->useEtaShiftAndSignSwap();
    //analysis->setVerbose();
    
    // Initialize histogram manager
    HistoManagerDiJet *hm = new HistoManagerDiJet{};
    hm->setIsMc(kTRUE);
    hm->init(kTRUE); // kTRUE stands up for use MC; need to FIX

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
