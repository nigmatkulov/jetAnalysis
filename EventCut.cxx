// Jet analysis headers
#include "EventCut.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

// C++ headers
#include <iostream>

ClassImp(EventCut)

//________________
EventCut::EventCut() : fVx{-1e9, 1e9}, fVy{-1e9, 1e9}, fVz{-1e9, 1e9},
    fShiftVx{0}, fShiftVy{0}, fVR{1e9}, 
    fHiBin{-1000, 1000}, fCentVal{-1000., 1000.},
    fPtHat{-1e9, 1e9}, fWeight{-1e9, 1e9}, fVerbose{kFALSE},
    fEventsPassed{0}, fEventsFailed{0} {
    fLumi[0] = 0;
    fLumi[1] = std::numeric_limits<unsigned int>::max();
}

//________________
EventCut::~EventCut() {
    /* Empty */
}

//________________
void EventCut::report() {
    TString report = "\nReporting from EventCut";
    report += TString::Format( "Vx              :\t %f - %f\n", fVx[0], fVx[1] );
    report += TString::Format( "Vy              :\t %f - %f\n", fVy[0], fVy[1] );
    report += TString::Format( "Vz              :\t %f - %f\n", fVz[0], fVz[1] );
    report += TString::Format( "HiBin           :\t %d - %d\n", fHiBin[0], fHiBin[1] );
    report += TString::Format( "Centrality      :\t %f - %f\n", fCentVal[0], fCentVal[1] );
    report += TString::Format( "pThat           :\t %f - %f\n", fPtHat[0], fPtHat[1] );
    report += TString::Format( "eventWeight     :\t %f - %f\n", fWeight[0], fWeight[1] );
    report += TString::Format( "Events passed   :\t %lld\n", fEventsPassed );
    report += TString::Format( "Events failed   :\t %lld\n", fEventsFailed );
    std::cout << report.Data() << std::endl;
}

//________________
Bool_t EventCut::pass(const Event* ev) {
    
    if (fVerbose) {
        std::cout << "\n----- Event cut -----\n";
    }

    const Bool_t goodVx = kTRUE;
    const Bool_t goodVy = kTRUE;

    const Bool_t goodVz = ( fVz[0] <= ev->vz() ) &&
                          ( ev->vz() < fVz[1] );
    if (fVerbose) {
        std::cout << Form("vz        : %5.2f <= %5.2f < %5.2f \t %s \n",
                          fVz[0], ev->vz(), fVz[1], ( goodVz ) ? "true" : "false" );
    }

    const Bool_t goodHiBin = ( fHiBin[0] <= ev->hiBin() ) &&
                             ( ev->hiBin() < fHiBin[1] );
    if (fVerbose) {
        std::cout << Form("hiBin        : %d <= %d < %d \t %s \n",
                          fHiBin[0], ev->hiBin(), fHiBin[1], ( goodHiBin ) ? "true" : "false" );
    }
    const Bool_t goodCent = ( fCentVal[0] <=  ev->centrality() ) &&
                            ( ev->centrality() < fCentVal[1] );

    if (fVerbose) {
        std::cout << Form("centrality   : %5.2f <= %5.2f < %5.2f \t %s \n",
                          fCentVal[0], ev->centrality(), fCentVal[1], ( goodCent ) ? "true" : "false" );
    }

    const Bool_t goodPtHat = ( fPtHat[0] <= ev->ptHat() ) &&
                             ( ev->ptHat() < fPtHat[1] );

    if (fVerbose) {
        std::cout << Form("ptHat        : %9.2f <= %9.2f < %9.2f \t %s \n",
                          fPtHat[0], ev->ptHat(), fPtHat[1], ( goodPtHat ) ? "true" : "false" );
    }

    const Bool_t goodEventWeight = ( fWeight[0] <= ev->weight() ) &&
                                   ( ev->weight() < fWeight[1] );
                
    if (fVerbose) {
        std::cout << Form("eventWeight  : %7.2f <= %7.2f < %7.2f \t %s \n",
                          fWeight[0], ev->weight(), fWeight[1], ( goodEventWeight ) ? "true" : "false" );
    }    

    Bool_t passEvent = goodVx && goodVy && goodVz && goodHiBin &&
                       goodCent && goodPtHat && goodEventWeight;
    ( passEvent ) ? fEventsPassed++ : fEventsFailed++;
    
    return passEvent;
}