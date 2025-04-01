/**
 * @file EventCut.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Example of the event cut
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "EventCut.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

// C++ headers
#include <iostream>

//________________
EventCut::EventCut() : fVx{-1e9, 1e9}, fVy{-1e9, 1e9}, fVz{-1e9, 1e9},
    fShiftVx{0}, fShiftVy{0}, fVR{1e9}, 
    fHiBin{-1000, 1000}, fCentVal{-1000., 1000.},
    fPtHat{-1e9, 1e9}, fPtHatWeight{-1e9, 1e9}, fVerbose{false},
    fPPrimaryVertexFilter{false},
    fHBHENoiseFilterResultRun2Loose{false},
    fCollisionEventSelectionAODc2{false},
    fPhfCoincFilter2Th4{false},
    fPPAprimaryVertexFilter{false},
    fPBeamScrapingFilter{false},
    fPClusterCompatibilityFilter{false},
    fPhfCoincFilter{false},
    fPVertexFilterCutdz1p0{false},
    fPVertexFilterCutGplus{false},
    fPVertexFilterCutVtx1{false},
    fHLT_HIAK4CaloJet60_v1{false},
    fHLT_HIAK4CaloJet80_v1{false},
    fHLT_HIPuAK4CaloJet80Eta5p1_v1{false},
    fHLT_PAAK4PFJet60_Eta5p1_v4{false},
    fHLT_PAAK4PFJet80_Eta5p1_v3{false},
    fHLT_PAAK4PFJet100_Eta5p1_v3{false},
    fHLT_PAAK4PFJet120_Eta5p1_v2{false},
    fHLT_HIAK4PFJet60_v1{false},
    fHLT_HIAK4PFJet80_v1{false},
    fRunIdsToSelect{},
    fRunIdsToExclude{},
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
    report += TString::Format( "ptHatWeight     :\t %f - %f\n", fPtHatWeight[0], fPtHatWeight[1] );
    report += TString::Format( "Events passed   :\t %lld\n", fEventsPassed );
    report += TString::Format( "Events failed   :\t %lld\n", fEventsFailed );
    std::cout << report.Data() << std::endl;
}

//________________
bool EventCut::pass(const Event* ev) {
    
    if (fVerbose) {
        std::cout << "\n----- Event cut -----\n";
    }

    const bool goodVx = true;
    const bool goodVy = true;

    const bool goodVz = ( fVz[0] <= ev->vz() ) &&
                          ( ev->vz() < fVz[1] );
    if (fVerbose) {
        std::cout << Form("vz        : %5.2f <= %5.2f < %5.2f \t %s \n",
                          fVz[0], ev->vz(), fVz[1], ( goodVz ) ? "true" : "false" );
    }

    const bool goodHiBin = ( fHiBin[0] <= ev->hiBin() ) &&
                             ( ev->hiBin() < fHiBin[1] );
    if (fVerbose) {
        std::cout << Form("hiBin        : %d <= %d < %d \t %s \n",
                          fHiBin[0], ev->hiBin(), fHiBin[1], ( goodHiBin ) ? "true" : "false" );
    }
    const bool goodCent = ( fCentVal[0] <=  ev->centrality() ) &&
                            ( ev->centrality() < fCentVal[1] );

    if (fVerbose) {
        std::cout << Form("centrality   : %5.2f <= %5.2f < %5.2f \t %s \n",
                          fCentVal[0], ev->centrality(), fCentVal[1], ( goodCent ) ? "true" : "false" );
    }

    const bool goodPtHat = ( fPtHat[0] <= ev->ptHat() ) &&
                             ( ev->ptHat() < fPtHat[1] );

    if (fVerbose) {
        std::cout << Form("ptHat        : %9.2f <= %9.2f < %9.2f \t %s \n",
                          fPtHat[0], ev->ptHat(), fPtHat[1], ( goodPtHat ) ? "true" : "false" );
    }

    const bool goodPtHatWeight = ( fPtHatWeight[0] <= ev->ptHatWeight() ) &&
                                   ( ev->ptHatWeight() < fPtHatWeight[1] );
                
    if (fVerbose) {
        std::cout << Form("eventWeight  : %7.2f <= %7.2f < %7.2f \t %s \n",
                          fPtHatWeight[0], ev->ptHatWeight(), fPtHatWeight[1], ( goodPtHatWeight ) ? "true" : "false" );
    }

    bool goodFilters{true};
    if ( fPPrimaryVertexFilter ) {
        if ( ev->trigAndSkim()->pprimaryVertexFilter() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pprimaryVertexFilter\n");
            }
        }

    }
    if ( fHBHENoiseFilterResultRun2Loose ) {
        if ( ev->trigAndSkim()->HBHENoiseFilterResultRun2Loose() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad HBHENoiseFilterResultRun2Loos\n");
            }
        }
    }
    if ( fCollisionEventSelectionAODc2 ) {
        if ( ev->trigAndSkim()->collisionEventSelectionAODv2() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad collisionEventSelectionAODv2\n");
            }
        }
    }
    if ( fPhfCoincFilter2Th4 ) {
        if ( ev->trigAndSkim()->phfCoincFilter2Th4() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad phfCoincFilter2Th4\n");
            }
        }
    }
    if ( fPPAprimaryVertexFilter ) {
        if ( ev->trigAndSkim()->pPAprimaryVertexFilter() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pPAprimaryVertexFilter\n");
            }
        }
    }
    if ( fPBeamScrapingFilter ) {
        if ( ev->trigAndSkim()->pBeamScrapingFilter() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pBeamScrapingFilter\n");
            }
        }
    }
    if ( fPClusterCompatibilityFilter ) {
        if ( ev->trigAndSkim()->pClusterCompatibilityFilter() == 0 ) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pClusterCompatibilityFilter\n");
            }
        }
    }
    if ( fPhfCoincFilter ) {
        if ( ev->trigAndSkim()->phfCoincFilter() == 0) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad phfCoincFilter\n");
            }
        }
    }
    if ( fPVertexFilterCutdz1p0 ) {
        if ( ev->trigAndSkim()->pVertexFilterCutdz1p0() == 0) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pVertexFilterCutdz1p0\n");
            }
        }
    }
    if ( fPVertexFilterCutGplus ) {
        if ( ev->trigAndSkim()->pVertexFilterCutGplus() == 0) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pVertexFilterCutGplus\n");
            }
        }
    }
    if ( fPVertexFilterCutVtx1 ) {
        if ( ev->trigAndSkim()->pVertexFilterCutVtx1() == 0) {
            goodFilters = {false};
            if ( fVerbose ) {
                std::cout << Form("Bad pVertexFilterCutVtx1\n");
            }
        }
    }
    if (fVerbose) {
        std::cout << Form("Event filters passed: %s\n", (goodFilters) ? "true" : "false");
    }

    bool goodTrigger{true};
    if ( fHLT_HIAK4CaloJet60_v1 ) {
        if ( ev->trigAndSkim()->HLT_HIAK4CaloJet60_v1() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_HIAK4CaloJet60_v1\n");
            }
        }
    }
    if ( fHLT_HIAK4CaloJet80_v1 ) {
        if ( ev->trigAndSkim()->HLT_HIAK4CaloJet80_v1() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_HIAK4CaloJet80_v1\n");
            }
        }
    }
    if ( fHLT_HIPuAK4CaloJet80Eta5p1_v1 ) {
        if ( ev->trigAndSkim()->HLT_HIPuAK4CaloJet80Eta5p1_v1() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_HIPuAK4CaloJet80Eta5p1_v1\n");
            }
        }
    }
    if ( fHLT_PAAK4PFJet60_Eta5p1_v4 ) {
        if ( ev->trigAndSkim()->HLT_PAAK4PFJet60_Eta5p1_v4() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_PAAK4PFJet60_Eta5p1_v4\n");
            }
        }
    }
    if ( fHLT_PAAK4PFJet80_Eta5p1_v3 ) {
        if ( ev->trigAndSkim()->HLT_PAAK4PFJet80_Eta5p1_v3() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_PAAK4PFJet80_Eta5p1_v3\n");
            }
        }
    }
    if ( fHLT_PAAK4PFJet100_Eta5p1_v3 ) {
        if ( ev->trigAndSkim()->HLT_PAAK4PFJet100_Eta5p1_v3() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_PAAK4PFJet100_Eta5p1_v3\n");
            }
        }
    }
    if ( fHLT_PAAK4PFJet120_Eta5p1_v2 ) {
        if ( ev->trigAndSkim()->HLT_PAAK4PFJet120_Eta5p1_v2() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_PAAK4PFJet120_Eta5p1_v2\n");
            }
        }
    }
    if ( fHLT_HIAK4PFJet60_v1) {
        if ( ev->trigAndSkim()->HLT_HIAK4PFJet60_v1() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_HIAK4PFJet60_v1\n");
            }
        }
    }
    if ( fHLT_HIAK4PFJet80_v1) {
        if ( ev->trigAndSkim()->HLT_HIAK4PFJet80_v1() == 0 ) {
            goodTrigger = { false };
            if ( fVerbose ) {
                std::cout << Form("Bad trigger: HLT_HIAK4PFJet80_v1\n");
            }
        }
    }
    if ( fVerbose ) {
        std::cout << Form("Event triggers passed: %s\n", (goodTrigger) ? "true" : "false");
    }

    bool isRunIdToSelect = ( fRunIdsToSelect.size() == 0 ) ? true : false;
    for ( auto runId : fRunIdsToSelect ) {
        if ( ev->runId() == runId ) {
            isRunIdToSelect = true;
            break;
        }
    }
    if ( fVerbose ) {
        std::cout << Form("RunId to select: %s\n", (isRunIdToSelect) ? "true" : "false");
    }

    bool isBadRunId = ( fRunIdsToExclude.size() == 0 ) ? false : true;
    for ( auto runId : fRunIdsToExclude ) {
        if ( ev->runId() == runId ) {
            isBadRunId = true;
            break;
        }
    }
    if ( fVerbose ) {
        std::cout << Form("RunId to exclude: %s\n", (isBadRunId) ? "true" : "false");
    }

    bool passEvent = goodVx && goodVy && goodVz && goodHiBin && goodFilters &&
                     goodCent && goodPtHat && goodPtHatWeight && goodTrigger && 
                     isRunIdToSelect && !isBadRunId;
    ( passEvent ) ? fEventsPassed++ : fEventsFailed++;
    
    return passEvent;
}