/**
 * @file Manager.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Program manager - main entry point to the framework
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "Manager.h"
#include "Event.h"

ClassImp(Manager)

//________________
Manager::Manager() : 
    fAnalysisCollection{nullptr}, fEventReader{nullptr},
    fEventsInChain{0} {
    fAnalysisCollection = new AnalysisCollection;
}

//________________
Manager::~Manager() {
    AnalysisIterator iter;
    for (iter = fAnalysisCollection->begin(); 
         iter != fAnalysisCollection->end(); 
         iter++ ) {
        delete *iter;
        *iter = nullptr;
    }
    if (fEventReader) delete fEventReader;
}

//________________
void Manager::init() {
    if (fEventReader) {
        fEventReader->init();
        fEventReader->report();
    }
    
    fEventsInChain = fEventReader->nEventsTotal();

    AnalysisIterator anaIter;
    for ( anaIter = fAnalysisCollection->begin();
          anaIter != fAnalysisCollection->end();
          anaIter++ ) {
        (*anaIter)->init();
    }
}

//________________
void Manager::finish() {
    if (fEventReader) {
        fEventReader->finish();
    }

    AnalysisIterator anaIter;
    for ( anaIter = fAnalysisCollection->begin();
          anaIter != fAnalysisCollection->end();
          anaIter++ ) {
        (*anaIter)->finish();
    }
}

//________________
void Manager::performAnalysis() {

    //std::cout << "Manager::performAnalysis - Number of events in chain: " << fEventsInChain << std::endl;

    int nEventsPerCycle{50000};
    double progress{0.0};

    // Loop over all events available
    for (Long64_t iEvent=0; iEvent<fEventsInChain; iEvent++) {

        // Print progress
        if ( iEvent % nEventsPerCycle == 0 ) {
            progress = 100.0 * iEvent / fEventsInChain;
            std::cout << "Processed " << iEvent << " events. Progress: " << progress << "%" << std::endl;
        }

        //std::cout << "=================================" << std::endl;
        //std::cout << "Manager::performAnalysis - Processing event: " << iEvent << std::endl;
        Event *currentEvent = fEventReader->returnEvent();

        if ( !currentEvent ) {
            if ( fEventReader->status() != 0) {
                std::cout << "Reader returned status: " 
                        << fEventReader->status() 
                        << ". Terminating\n";
            }
        } // if ( !currentEvent)
        else {
            // Perform data processing by all analyses
            AnalysisIterator anaIter;
            for ( anaIter = fAnalysisCollection->begin();
                anaIter != fAnalysisCollection->end();
                anaIter++ ) {

                (*anaIter)->processEvent( currentEvent );
            }
        }

        if ( currentEvent ) {
            delete currentEvent;
            currentEvent = nullptr;
        }
    } // for (Long64_t iEvent=0; iEvent<fEventsInChain; iEvent++)
}

//________________
void Manager::report() {
    std::cout << "Reporting from analysis manager" << std::endl;
}