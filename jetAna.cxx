// C++ headers
#include <iostream>

// Jet analysis headers
#include "Manager.h"
#include "ForestAODReader.h"
#include "JetESRAnalysis.h"
#include "HistoManagerJetESR.h"
#include "EventCut.h"
#include "JetCut.h"

// ROOT headers
#include "TFile.h"

//________________
void usage() {
    std::cout << "./programName inputFileList oFileName" << std::endl;
}

//________________
/// @brief The prorgram that launches the physics analysis
/// @param argc Number of arguments
/// @param argv Argument list
/// @return 0 in case of OKAY
int main(int argc, char const *argv[]) {
    /* code */

    // Set default values for arguments
    
    //"../../../data/HiForestAOD_PbPbMC2018skim_10.root"

    // Sequence of command line arguments:
    //
    // inputFileList (or forest.root) - input file list with forest
    // outputFileName.root            - output file name
    //

    Bool_t isPbPb{kTRUE};
    Bool_t isMc{kTRUE};
    Bool_t isCentWeightCalc{kTRUE};
    TString inFileName{};
    Int_t   collEnergyGeV{5020};
    TString collSystem{};
    Int_t   collisionsSystem{2}; // 0 - pp, 1 - pPb, 2 - PbPb
    Int_t   collYear{2018};
    TString recoJetBranchName{};
    TString oFileName{};
    TString JECFileName;
    //TString path2JEC = "/Users/gnigmat/work/cms/soft/jetAnalysis";
    TString path2JEC = "../";
    if ( isPbPb ) {
        inFileName = "../../../data/HiForestAOD_PbPb_sim.list";
        //inFileName = "../../../data/HiForestAOD_PbPb_exp.list";
        collSystem = "PbPb";
        recoJetBranchName = "akCs4PFJetAnalyzer";
        oFileName = "oTestReadForest_PbPb.root";
        JECFileName = "Autumn18_HI_V8_MC_L2Relative_AK4PF.txt";
    }
    else { // pp
        //inFileName = "../../../data/pp/HiForestAOD_1113.root";
        inFileName = "../../../data/HiForestAOD_pp.list";
        collSystem = "pp";
        recoJetBranchName = "ak4PFJetAnalyzer";
        oFileName = "oTestReadForest_pp.root";
        JECFileName = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt";
    }

    if (argc <= 1) {
        std::cout << "Too few arguments passed. Terminating" << std::endl;
        usage();
		return -1;
    }
    else {
        // Read input argument list 
        inFileName = argv[1];
        oFileName = argv[2];
    }

    std::cout << "Arguments passed:\n"
              << "Input file name       : " << inFileName << std::endl
              << "Output file name      : " << oFileName << std::endl
              << "Is MC                 : " << isMc << std::endl
              << "Use centrality weight : " << isCentWeightCalc << std::endl
              << std::endl;


    Long64_t nEventsToRead = 500;

    Manager *manager = new Manager{};
    EventCut *eventCut = new EventCut{};
    eventCut->setVz(-15., 15.);
    if ( isPbPb ) {
        eventCut->usePPrimaryVertexFilter();
        eventCut->useHBHENoiseFilterResultRun2Loose();
        eventCut->useCollisionEventSelectionAODv2();
        eventCut->usePhfCoincFilter2Th4();
        eventCut->usePClusterCompatibilityFilter();
        // Trigger
        eventCut->useHLT_HIPuAK4CaloJet80Eta5p1_v1();
    }
    else { // pp case
        eventCut->useHBHENoiseFilterResultRun2Loose();
        eventCut->usePPAprimaryVertexFilter();
        eventCut->usePBeamScrapingFilter();
    }

    if ( isMc ) {
        //eventCut->setPtHat(50, 1e6);
        eventCut->setPtHat(15, 1e6);
    }
    //eventCut->setVerbose();

    JetCut *jetCut = new JetCut{};
    //jetCut->setMustHaveGenMathing();
    //jetCut->setPt(50., 1500.);
    jetCut->setPt(20., 1500.);
    //jetCut->setEta(-1.6, 1.6);
    //jetCut->setVerbose();

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
    reader->setCollidingSystem( collSystem.Data() );
    reader->setCollidingEnergy( collEnergyGeV ) ;
    reader->setYearOfDataTaking( collYear );
    reader->setEventCut(eventCut);
    reader->setJetCut(jetCut);
    reader->fixJetArrays();

    // Set path to jet analysis (then will automatically add path to aux_files)
    reader->setPath2JetAnalysis( path2JEC.Data() );
    reader->addJECFile( JECFileName.Data() );
    
    manager->setEventReader(reader);

    JetESRAnalysis *analysis = new JetESRAnalysis{};
    analysis->setCollisionSystem(collisionsSystem);
    // analysis->useCentralityWeight();  // For PbPb
    // analysis->setPtHatRange(15., 1e6);
    analysis->setLeadJetPtLowCut(50.);
    analysis->setSubleadJetPtLowCut(40.);
    analysis->setDijetDPhiCut( 2. * TMath::Pi() / 3. );
    analysis->setLeadJetEtaCut(-3., 3.);
    analysis->setSubleadJetEtaCut(-3., 3.);


    HistoManagerJetESR *hm = new HistoManagerJetESR{};
    // hm->setIsMc(kTRUE);
    hm->init(true);
    analysis->addHistoManager(hm);
    manager->addAnalysis(analysis);

    manager->init();
    manager->performAnalysis();
    manager->finish();

    TFile* oFile = new TFile(oFileName, "recreate");
    hm->writeOutput();
    oFile->Close();
    
    return 0;
}
