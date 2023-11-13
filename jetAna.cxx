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
/// @brief The prorgram that launches the physics analysis
/// @param argc Number of arguments
/// @param argv Argument list
/// @return 0 in case of OKAY
int main(int argc, char const *argv[]) {
    /* code */

    // Set default values for arguments
    
    //"../../../data/HiForestAOD_PbPbMC2018skim_10.root"

    Bool_t isPbPb = kTRUE;
    TString inFileName{};
    Int_t   collEnergyGeV{};
    TString collSystem{};
    Int_t   collYear{};
    TString pfBranchName{};
    TString oFileName{};
    TString JECFileName;
    if ( isPbPb ) {
        inFileName = "../../../data/HiForestAOD_PbPb.list";
        collEnergyGeV = {5020};
        collSystem = "PbPb";
        collYear = 2018;
        pfBranchName = "akCs4PFJetAnalyzer";
        oFileName = "oTestReadForest_PbPb.root";
    }
    else { // pp
        //inFileName = "../../../data/pp/HiForestAOD_1113.root";
        inFileName = "../../../data/HiForestAOD_pp.list";
        collEnergyGeV = {5020};
        collSystem = "pp";
        collYear = 2018;
        pfBranchName = "ak4PFJetAnalyzer";
        oFileName = "oTestReadForest_pp.root";
        JECFileName = "Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK4PF.txt";
    }

    Long64_t nEventsToRead = 500;

    // Read input argument list 
    if (argc > 1) inFileName = argv[1];
    if (argc > 2) oFileName = argv[2];

    Manager *manager = new Manager{};
    EventCut *eventCut = new EventCut{};
    eventCut->setVz(-15., 15.);
    if ( isPbPb ) {
        eventCut->usePPrimaryVertexFilter();
        eventCut->useHBHENoiseFilterResultRun2Loose();
        eventCut->useCollisionEventSelectionAODv2();
        eventCut->usePhfCoincFilter2Th4();
        eventCut->usePClusterCompatibilityFilter();
    }
    else { // pp case
        eventCut->useHBHENoiseFilterResultRun2Loose();
        eventCut->usePPAprimaryVertexFilter();
        eventCut->usePBeamScrapingFilter();
    }
    eventCut->setPtHat(15, 1e6);

    JetCut *jetCut = new JetCut{};
    //jetCut->setMustHaveGenMathing();
    jetCut->setPt(20., 1500.);
    jetCut->setEta(-1.6, 1.6);

    //jetCut->setVerbose();
    ForestAODReader *reader = new ForestAODReader(inFileName);
    reader->setIsMc(kTRUE);
    reader->useHltBranch();
    reader->useSkimmingBranch();
    reader->usePartFlowJetBranch();
    reader->setPartFlowJetBranchName( pfBranchName.Data() );
    //reader->useCaloJetBranch();
    reader->setCollidingSystem( collSystem.Data() );
    reader->setCollidingEnergy( collEnergyGeV) ;
    reader->setYearOfDataTaking( collYear );
    reader->setEventCut(eventCut);
    reader->setJetCut(jetCut);
    reader->fixJetArrays();
    if ( !isPbPb ) {
        reader->setJECFileName(JECFileName.Data());
    }
    manager->setEventReader(reader);

    JetESRAnalysis *analysis = new JetESRAnalysis{};
    HistoManagerJetESR *hm = new HistoManagerJetESR{};
    hm->setIsMc(kTRUE);
    hm->init(kTRUE); // kTRUE stands up for use MC; need to FIX
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
