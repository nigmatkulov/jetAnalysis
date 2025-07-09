/**
 * @file JetESRAnalysis.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet energy scale and resolution analysis
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// ROOT headers
#include "TF1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"

// C++ headers
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <cmath>

// Jet analysis headers
#include "JetESRAnalysis.h"
#include "JetCut.h"

//________________
JetESRAnalysis::JetESRAnalysis() : BaseAnalysis(), 
    fVerbose{false}, fHM{nullptr}, fVzWeight{nullptr}, 
    fUseCentralityWeight{false}, fEtaShift{0},
    fCollisionEnergy{8160}, fCollisionSystem{1}, 
    fIsPbGoingDir{false}, fPtHatRange{0., 10000000.},
    fIsMc{true}, fNEventsInSample{1000000}, 
    fRecoIdLead{-1}, fRecoIdSubLead{-1},
    fGenIdLead{-1}, fGenIdSubLead{-1},
    fRefSelRecoIdLead{-1}, fRefSelRecoIdSubLead{-1},
    fRecoPtSortedJetIds{},
    fGenPtSortedJetIds{},
    fRefSelRecoPtSortedJetIds{},
    fRecoJetCut{nullptr}, fGenJetCut{nullptr} {
    /* Empty */
}

//________________
JetESRAnalysis::~JetESRAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
    if (fVzWeight) { delete fVzWeight; fVzWeight = nullptr; }
    if (fRecoJetCut) { delete fRecoJetCut; fRecoJetCut = nullptr; }
    if (fGenJetCut) { delete fGenJetCut; fGenJetCut = nullptr; }
}

//________________
void JetESRAnalysis::init() {
    // Initialize analysis
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::init - begin" << std::endl;
    }

    // Print initial setup of the analysis
    print();

    if ( fIsMc ) {
        // Initialize vz weight function
        initVzWeightFunction();
    }

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::init - end" << std::endl;
    }
}

//________________
void JetESRAnalysis::initVzWeightFunction() {

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::initVzWeightFunction - begin" << std::endl;
    }

    // Check if Vz weight function exists
    if ( !fVzWeight ) {
        if ( fCollisionSystem == 0 ) { // Assume pp 5020
            fVzWeight = new TF1("fVzWeight", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x", -20., 20.);
            fVzWeight->SetParameters(0.973941, 0.00310622, 0.000711664, -1.83098e-06, 6.9346e-07, 0., 0.);
        }
        else if ( fCollisionSystem == 1 ) { // Assume pPb 8160
            fVzWeight = new TF1("fVzWeight", "pol8", -15.1, 15.1);
            fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);
        }
        else if ( fCollisionSystem == 2 ) { // Assume PbPb 5020
            fVzWeight = new TF1("fVzWeight", "pol0", -15.1, 15.1);
            fVzWeight->SetParameter(0, 1.);
        }
        else { // Unknown collision system
            fVzWeight = new TF1("fVzWeight", "pol0", -200.1, 200.1);
            fVzWeight->SetParameter(0, 1.); 
        }
    }

    if ( fVerbose ) {
        std::cout << "Vz weight function: ";
        fVzWeight->Print();

        std::cout << "JetESRAnalysis::initVzWeightFunction - end" << std::endl;
    }
}

//________________
TString JetESRAnalysis::collisionSystem() const {
    TString collSys = "PbPb";
    if ( fCollisionSystem == 0 ) {
        collSys = "pp";
    }
    else if ( fCollisionSystem == 1 ) {
        collSys = "pPb";
    }
    else if ( fCollisionSystem == 2 ) {
        collSys = "PbPb";
    }
    else {
        collSys = "Unknown";
    }

    return collSys;
}

//________________
void JetESRAnalysis::print() {

    if ( fCollisionSystem)
    std::cout << "----------------------------------------\n";
    std::cout << "JetESRAnalysis initial parameters:\n";
    std::cout << "Use centrality weight       : " << fUseCentralityWeight << std::endl
              << "Histogram manager           : " << fHM << std::endl
              << "Is MC                       : " << fIsMc << std::endl
              << "Collisions system           : " << collisionSystem().Data() << std::endl
              << "Is Pb-going direction       : " << fIsPbGoingDir << std::endl
              << "eta shift                   : " << fEtaShift << std::endl
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl;
              if ( fRecoJetCut ) {
                  std::cout << "Reco jet cut parameters     : " << std::endl;
                  fRecoJetCut->report();
              }
              if ( fGenJetCut ) {
                  std::cout << "Gen jet cut parameters      : " << std::endl;
                  fGenJetCut->report();
              }
    std::cout << "----------------------------------------\n";
}

//________________
double JetESRAnalysis::eventWeight(const double& ptHat, const double& vz, 
                                   const double& centWeight, const double& ptHatW) {

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::eventWeight - begin" << std::endl;
    }

    // Calculate event weight
    double weight{1.};
    double genWeight{1.};
    double vzWeight{1.};

    if ( fIsMc ) {
        // In case of pPb (assumed to be pPb8160)
        if ( fCollisionSystem == 0 ) {     // Assuming pp 5020
            weight = ptHatW;
            if ( fVzWeight ) {
                vzWeight = fVzWeight->Eval( vz );
            }
            weight *= vzWeight;
        }
        else if ( fCollisionSystem == 1 ) { // Assuming pPb 8160

            // Magic numbers are (cross section x Nevents generated)
            if (ptHat > 15.0 && ptHat <= 30.)       { genWeight = 1.0404701e-06 * 961104 ; }
            else if (ptHat > 30. && ptHat <= 50.)   { genWeight = 7.7966624e-08 * 952110 ; }
            else if (ptHat > 50. && ptHat <= 80.)   { genWeight = 1.0016052e-08 * 952554 ; }
            else if (ptHat > 80. && ptHat <= 120.)  { genWeight = 1.3018269e-09 * 996844 ; }
            else if (ptHat > 120.&& ptHat <= 170.)  { genWeight = 2.2648493e-10 * 964681 ; }
            else if (ptHat > 170. && ptHat <= 220.) { genWeight = 4.0879112e-11 * 999260 ; }
            else if (ptHat > 220. && ptHat <= 280.) { genWeight = 1.1898939e-11 * 964336 ; }
            else if (ptHat > 280. && ptHat <= 370.) { genWeight = 3.3364433e-12 * 995036 ; }
            else if (ptHat > 370. && ptHat <= 460.) { genWeight = 7.6612402e-13 * 958160 ; }
            else if (ptHat > 460. && ptHat <= 540.) { genWeight = 2.1341026e-13 * 981427 ; }
            else if (ptHat > 540.)                  { genWeight = 7.9191586e-14 * 1000000; }
            genWeight /= fNEventsInSample;
            
            // Vz weighting
            if ( fVzWeight ) {
                vzWeight = fVzWeight->Eval( vz );
            }
            vzWeight = 1. / vzWeight;

            weight = genWeight * vzWeight;
        } 
        else if ( fCollisionSystem == 2 ) { // Assuming PbPb 5020
            weight = ptHatW;
            if ( fUseCentralityWeight ) {
                weight *= centWeight;
            }
        }
        else {
            weight = 1.;
        }
    } // if ( fIsMc )

    if ( fVerbose) {
        std::cout << Form("Input parameters:\nptHat: %5.2f vz: %5.2f centWeight: %5.2f ptHatW: %5.2f\n", 
                          ptHat, vz, centWeight, ptHatW);
        std::cout << Form("Calculated parameters:\nweight: %5.2f genWeight: %5.2f vzWeight: %5.2f\n", weight, genWeight, vzWeight);
        if ( fCollisionSystem == 1 ) {
            std::cout << Form("fNEventsInSample: %d\n", fNEventsInSample);
        }
        std::cout << "JetESRAnalysis::eventWeight - end" << std::endl;
    }

    return weight;
}

//________________
bool JetESRAnalysis::isOverweightedEvent(const Event* event, const double& weight) {
    
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::isOverweightedEvent -- begin" << std::endl;
    }

    // Clear the state
    bool isRecoOverweightedEvent = {false};
    bool isGenOverweightedEvent = {false};

    // Event ptHat value
    float ptHat = event->ptHat();

    //
    // Process reco jets in order to remove x-jets
    //
    if ( fRecoIdLead >= 0 && fRecoIdSubLead >= 0 ) {
        RecoJet *recoLeadJet = event->recoJetCollection()->at( fRecoIdLead );
        RecoJet *recoSubLeadJet = event->recoJetCollection()->at( fRecoIdSubLead );
        
        float ptRecoLead = recoLeadJet->pt();
        float ptRecoSubLead = recoSubLeadJet->pt();
        float dijetRecoPt = ptRecoLead + ptRecoSubLead;
        float dijetRecoPtAve = 0.5 * dijetRecoPt;

        fHM->hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Fill( ptRecoLead, ptRecoLead/ptHat, 1. );
        fHM->hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Fill( ptRecoLead, ptRecoLead/ptHat, weight );
        
        fHM->hRecoDijetPtOverPtHatVsDijetPt->Fill(dijetRecoPt, dijetRecoPt/ptHat, 1.);
        fHM->hRecoDijetPtOverPtHatVsDijetPtWeighted->Fill(dijetRecoPt, dijetRecoPt/ptHat, weight);
        fHM->hRecoDijetPtAveOverPtHatVsDijetPtAve->Fill(dijetRecoPtAve, dijetRecoPtAve/ptHat, 1.);
        fHM->hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Fill(dijetRecoPtAve, dijetRecoPtAve/ptHat, weight);

        if ( isOverweighted( ptRecoLead, dijetRecoPtAve, ptHat ) ) {
            if ( fVerbose ) {
                std::cout << Form("Overweighted event [reco]: ptLead/ptHat = %3.2f ptAve/ptHat = %3.2f", ptRecoLead/ptHat, dijetRecoPtAve/ptHat) << std::endl;
            }
            isRecoOverweightedEvent = {true};
        } // if ( isOverweightedEvent( ptRecoLead, ptHat ) )
    } // if ( fRecoIdLead >= 0 && fRecoIdSubLead >= 0 )
    else {
        // Skip events with less than 2 jets
        isRecoOverweightedEvent = {true};
    }

    //
    // Process gen jets in order to remove overweighted gen events
    //

    if ( fGenIdLead >= 0 && fGenIdSubLead >= 0 ) {

        GenJet* leadJet = event->genJetCollection()->at( fGenIdLead );
        GenJet* subLeadJet = event->genJetCollection()->at( fGenIdSubLead );
        float ptGenLead = leadJet->pt();
        float ptGenSubLead = subLeadJet->pt();
        float dijetGenPt = ptGenLead + ptGenSubLead;
        float dijetGenPtAve = 0.5 * dijetGenPt;

        fHM->hGenLeadingJetPtOverPtHatVsLeadingJetPt->Fill( ptGenLead, ptGenLead/ptHat, 1. );
        fHM->hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Fill( ptGenLead, ptGenLead/ptHat, weight );

        fHM->hGenDijetPtOverPtHatVsDijetPt->Fill(dijetGenPt, dijetGenPt/ptHat, 1.);
        fHM->hGenDijetPtOverPtHatVsDijetPtWeighted->Fill(dijetGenPt, dijetGenPt/ptHat, weight);
        fHM->hGenDijetPtAveOverPtHatVsDijetPtAve->Fill(dijetGenPtAve, dijetGenPtAve/ptHat, 1.);
        fHM->hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Fill(dijetGenPtAve, dijetGenPtAve/ptHat, weight);

        if ( isOverweighted( ptGenLead, dijetGenPtAve, ptHat ) ) {
            if ( fVerbose ) {
                std::cout << Form("Overweighted event [Gen]: ptLead/ptHat = %3.2f ptAve/ptHat = %3.2f", ptGenLead/ptHat, dijetGenPtAve/ptHat) << std::endl;
            }
            isGenOverweightedEvent = {true};
        } // if ( isOverweightedEvent( ptGenLead, ptHat ) )
    } // if ( fGenIdLead >= 0 && fGenIdSubLead >= 0 )
    else {
        // Skip events with less than 2 jets
        isGenOverweightedEvent = {true};
    }


    bool overweightedEvent = isRecoOverweightedEvent || isGenOverweightedEvent;

    if ( fVerbose ) {
        std::cout << Form("Event overweighted: %s Reco overweighted: %s Gen overweighted: %s", 
                          ((overweightedEvent) ? "[true]" : "[false]"),
                          ((isRecoOverweightedEvent) ? "[true]" : "[false]"), 
                          ((isGenOverweightedEvent) ? "[true]" : "[false]")) << std::endl;
        std::cout << "JetESRAnalysis::isOverweightedEvent -- end" << std::endl;
    }
    return overweightedEvent;
}

//________________
bool JetESRAnalysis::isOverweighted(const float& ptLead, const float& dijetPtAve, const float& ptHat) {
    return (  ( ( ptLead / ptHat ) > 2.5) || ( ( dijetPtAve / ptHat ) > 1.7) );
}

//________________
float JetESRAnalysis::boostEta2CM(const float &eta) {
    // if ( fVerbose ) {
    //     std::cout << "JetESRAnalysis::boostEta2CM -- begin" << std::endl;
    // }
    float etaCM = eta;

    // Apply lab frame boost to CM
    if ( fCollisionSystem == 0 ) { // pp
        // For pp do nothing. Already in the CM frame
    }

    else if ( fCollisionSystem == 1 ) { // pPb
        if ( fIsMc ) { // For embedding: Pb goes to negative, p goes to positive
            if ( fIsPbGoingDir ) {
                etaCM += fEtaShift;
                etaCM = -etaCM;
            }
            else {
                etaCM -= fEtaShift;
            }
        }
        else { // For data: p goes to negative, Pb goes to positive
            if ( fIsPbGoingDir ) {
                etaCM -= fEtaShift;
            }
            else {
                etaCM += fEtaShift;
                etaCM = -etaCM;
            }
        }
    } 
    else if ( fCollisionSystem == 2 ) { // PbPb
        // For PbPb do nothing. Already in the CM frame
    }
    else {
        // Unknown collision system
        // Do nothing
    }

    // if ( fVerbose ) {
    //     std::cout << Form("eta: %5.2f  ->  etaCM: %5.2f", eta, etaCM) << std::endl;
    //     std::cout << "JetESRAnalysis::boostEta2CM -- end" << std::endl;
    // }
    return etaCM;
}

//________________
float JetESRAnalysis::etaLab(const float &eta) {
    // if ( fVerbose ) {
    //     std::cout << "JetESRAnalysis::etaLab -- begin" << std::endl;
    // }

    float etaL = eta;
    // Check collision system
    if ( fCollisionSystem == 0 ) { 
        // For pp apply eta shift (to move from CM to lab frame, to match pPb)
        etaL += fEtaShift;
    }
    else if ( fCollisionSystem == 1 ) { 
        // For pPb we already in the lab frame. Just need to properly address
        // beam direction
        if ( fIsMc ) { // For embedding: Pb goes to negative, p goes to positive
            if (fIsPbGoingDir) {
                etaL = -etaL;
            }
        }
        else { // For data: p goes to negative, Pb goes to positive
            if (fIsPbGoingDir) {
            }
            else {
                etaL = -etaL;
            }
        }
    }
    else if ( fCollisionSystem == 2 ) { 
        // For PbPb apply eta shift (to move from CM to lab frame, to match pPb)
        etaL += fEtaShift;
    }
    else {
        // Unknown collision system
        // Do nothing
    }

    // if ( fVerbose ) {
    //     std::cout << Form("eta: %5.2f  ->  etaLab: %5.2f", eta, etaL) << std::endl;
    //     std::cout << "DiJetAnalysis::etaLab -- end" << std::endl;
    // }
    return etaL;
}

//________________
void JetESRAnalysis::makePtSortedJetVectors(const Event* event) {
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::makePtSortedJetVectors -- begin" << std::endl;
    }

    // Clear vectors
    fRecoPtSortedJetIds.clear();
    if ( fIsMc ) {
        fGenPtSortedJetIds.clear();
        fRefSelRecoPtSortedJetIds.clear();
    }

    // Reco jet iterators
    RecoJetIterator recoJetIter;
    // Gen jet iterators
    GenJetIterator genJetIter;

    //
    // Reco jets
    //

    // Jet counter
    int recoJetCounter{0};
    // std::cout << "Reco jet collection size: " << event->recoJetCollection()->size() << std::endl;
    // Loop over reconstructed jets and store indices of good jets
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        recoJetCounter++;

        // Check if reconstructed jet passes selection criteria (*jet, isCM, isMC, requireMatching)
        if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, true, false) ) continue; 

        fRecoPtSortedJetIds.push_back( recoJetCounter-1 );
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    // Sort indices based on the jet corrected pT (from high to low)
        std::sort( fRecoPtSortedJetIds.begin(), fRecoPtSortedJetIds.end(), [&](int i, int j) { 
            return event->recoJetCollection()->at(i)->ptJECCorr() > event->recoJetCollection()->at(j)->ptJECCorr(); 
    } );

    if ( fVerbose ) {
        // Print sorted indices and corresponding jet pT    
        for (const auto& id : fRecoPtSortedJetIds) {
            std::cout << Form("Sorted reco jet index: %d | pT: %5.1f eta: %3.2f\n", id, event->recoJetCollection()->at(id)->ptJECCorr(), etaLab(event->recoJetCollection()->at(id)->eta()));
        }
    }

    if ( fRecoPtSortedJetIds.size() >= 2) {
        // Indices of leading and subleading reco jets
        fRecoIdLead = fRecoPtSortedJetIds.at(0);
        fRecoIdSubLead = fRecoPtSortedJetIds.at(1);
    }
    else if ( fRecoPtSortedJetIds.size() == 1 ) {
        fRecoIdLead = fRecoPtSortedJetIds.at(0);
        fRecoIdSubLead = -1;
    }

    //
    // Pure Monte Carlo part
    //

    if ( fIsMc ) {

        //
        // Gen jets
        //

        int genJetCounter{0};
        // Loop over generated jets and store indices of good jets
        for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {
            genJetCounter++;

            // Check gen jet passes the selection criteria (*genJet, isCM)
            if ( fGenJetCut && !fGenJetCut->pass(*genJetIter, false) ) continue;

            fGenPtSortedJetIds.push_back( genJetCounter-1 );
        } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )

        // Sort indices based on the jet pT (from high to low)
        std::sort( fGenPtSortedJetIds.begin(), fGenPtSortedJetIds.end(), [&](int i, int j) { 
            return event->genJetCollection()->at(i)->pt() > event->genJetCollection()->at(j)->pt(); 
        } );

        if ( fVerbose ) {
            // Print sorted indices and corresponding jet pT
            for (const auto& id : fGenPtSortedJetIds) {
                std::cout << Form("Sorted gen jet index: %d | pT: %5.1f eta: %3.2f\n", id, event->genJetCollection()->at(id)->pt(), etaLab(event->genJetCollection()->at(id)->eta()));
            }
        }

        if ( fGenPtSortedJetIds.size() >= 2) {
            // Indices of leading and subleading gen jets
            fGenIdLead = fGenPtSortedJetIds.at(0);
            fGenIdSubLead = fGenPtSortedJetIds.at(1);
        }
        else if ( fGenPtSortedJetIds.size() == 1 ) {
            fGenIdLead = fGenPtSortedJetIds.at(0);
            fGenIdSubLead = -1;
        }

        //
        // Ref-selected reco jets
        //

        int refSelRecoJetCounter{0};
        // Loop over reconstructed jets and select those only that have matching gen jets
        for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {
            refSelRecoJetCounter++;

            // Check selection criteria (*recoJet, isCM, isMC, requireMatching)
            if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, true, true) ) continue;

            fRefSelRecoPtSortedJetIds.push_back( refSelRecoJetCounter-1 );
        } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

        // Sort indices based on the gen jet pT from the reco matched jets (from high to low)
        std::sort( fRefSelRecoPtSortedJetIds.begin(), fRefSelRecoPtSortedJetIds.end(), [&](int i, int j) { 
            return event->genJetCollection()->at( event->recoJetCollection()->at(i)->genJetId() )->pt() > event->genJetCollection()->at( event->recoJetCollection()->at(j)->genJetId() )->pt(); 
        } );

        if ( fVerbose ) {
            // Print sorted indices and corresponding gen jet pT
            for (const auto& id : fRefSelRecoPtSortedJetIds) {
                std::cout << Form("Sorted ref-selected reco jet index: %d | gen pT: %5.1f | reco pT: %5.1f\n", id, event->genJetCollection()->at( event->recoJetCollection()->at(id)->genJetId() )->pt(), event->recoJetCollection()->at(id)->ptJECCorr());
            }
        }

        if ( fRefSelRecoPtSortedJetIds.size() >= 2) {
            // Indices of leading and subleading ref-selected reco jets
            fRefSelRecoIdLead = fRefSelRecoPtSortedJetIds.at(0);
            fRefSelRecoIdSubLead = fRefSelRecoPtSortedJetIds.at(1);
        }
        else if ( fRefSelRecoPtSortedJetIds.size() == 1 ) {
            fRefSelRecoIdLead = fRefSelRecoPtSortedJetIds.at(0);
            fRefSelRecoIdSubLead = -1;
        }

    } // if ( fIsMc )

    if ( fVerbose ) {
        std::cout << Form("Reco leadId: %d, subleadId: %d\nGen leadId: %d, subleadId: %d\n", fRecoIdLead, fRecoIdSubLead, fRefSelRecoIdLead, fRefSelRecoIdSubLead) << std::endl;
        std::cout << "JetESRAnalysis::makePtSortedJetVectors -- end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processInclusiveJets(const Event* event, const double& weight) {
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processInclusiveJets -- begin" << std::endl;
    }

    fHM->hRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), 1. );
    processRecoJets( event, weight );

    if ( fIsMc ) {
        fHM->hGenJetCollectionSize->Fill( event->genJetCollection()->size(), 1. );
        fHM->hGenVsRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), event->genJetCollection()->size(), 1. );
        processGenJets( event, weight );
        processRefJets( event, weight );
    } // if ( fIsMc )

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processInclusiveJets -- end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processGenJets(const Event* event, const double &weight) {

    if ( weight <= 0. ) {
        std::cout << "JetESRAnalysis::processGenJets -- weight is zero or negative. Skip processing." << std::endl;
        return;
    }
    // ptHat value
    float ptHat = event->ptHat();

    // Jet iterators
    GenJetIterator genJetIter;

    // Jet counter
    int genJetCounter{0};
    if ( event->genJetCollection()->size() > 0 ) {

        // Loop over generated jets and search for leading and subleading jets
        for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

            float pt = (*genJetIter)->pt();
            float eta = etaLab( (*genJetIter)->eta() );
            float phi = (*genJetIter)->phi();
    
            if ( fVerbose ) {
                std::cout << Form("Gen jet #%d: pt: %5.2f eta: %5.2f phi: %5.2f\n", genJetCounter, pt, eta, phi);
            }

            genJetCounter++;

            // Selection criteria (*genJet, isCM)
            if ( fGenJetCut && !fGenJetCut->pass(*genJetIter, false) ) continue;

            // Inclusive gen jet
            fHM->hGenInclusiveJetPt->Fill(pt, weight);
            fHM->hGenInclusiveJetEta->Fill(eta, weight);
            fHM->hGenInclusiveJetEtaUnweighted->Fill(eta, 1.);
            fHM->hGenInclusiveJetPtEta->Fill(eta, pt, weight);
            fHM->hGenInclusiveJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

            // Leading gen jet
            if ( (fGenIdLead >= 0) && ((genJetCounter - 1) ==  fGenIdLead) ) {
                fHM->hGenLeadJetPtEta->Fill(eta, pt, weight);
                fHM->hGenLeadJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
            }

            // Subleading gen jet
            if ( (fGenIdSubLead >= 0) && ((genJetCounter - 1) ==  fGenIdSubLead) ) {
                fHM->hGenSubLeadJetPtEta->Fill(eta, pt, weight);
                fHM->hGenSubLeadJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
            }

            if ( pt > 30. ) {
                fHM->hGenGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), weight);
                fHM->hGenGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), weight );
            }

        } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )
    } // if ( event->genJetCollection()->size() > 0 )

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processGenJets -- end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processRecoJets(const Event* event, const double &weight) {

    if ( weight <= 0. ) {
        std::cout << "JetESRAnalysis::processRecoJets -- weight is zero or negative. Skip processing." << std::endl;
        return;
    }

    // ptHat value
    float ptHat = event->ptHat();

    // Jet iterators
    RecoJetIterator recoJetIter;

    // Jet counter
    int recoJetCounter{0};

    // Loop over reconstructed jets
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        float pt = (*recoJetIter)->ptJECCorr();
        float eta = etaLab( (*recoJetIter)->eta() );
        float phi = (*recoJetIter)->phi();
        float ptRaw = (*recoJetIter)->rawPt();

        if ( fVerbose ) {
            if ( !fIsMc ) {
                std::cout << Form("Reco jet #%d pt: %5.1f eta: %3.2f phi: %3.2f ptRaw: %5.1f", recoJetCounter, pt, eta, phi, ptRaw) << std::endl;
            }
            else {
                std::cout << Form("Reco jet #%d pt: %5.1f eta: %3.2f phi: %3.2f ptRaw: %5.1f genJetId: %d", recoJetCounter, pt, eta, phi, ptRaw, (*recoJetIter)->genJetId() ) << std::endl;
            }
        }

        recoJetCounter++;

        // JetId parameters
        int chargedMult = (*recoJetIter)->jtPfCHM() + (*recoJetIter)->jtPfCEM() + (*recoJetIter)->jtPfMUM();
        int neutralMult = (*recoJetIter)->jtPfNHM() + (*recoJetIter)->jtPfNEM();
        int numberOfConstituents = chargedMult + neutralMult;
        
        int dummyIter{0};
        if ( fabs( eta ) <= 2.4 ) { dummyIter = {0}; }
        else if ( fabs( eta ) <= 2.7 ) { dummyIter = {1}; }
        else if ( fabs( eta ) <= 3.0 ) { dummyIter = {2}; }
        else { dummyIter = {3}; }

        // JetId histograms
        fHM->hRecoInclusiveJetNHF[dummyIter]->Fill( (*recoJetIter)->jtPfNHF(), weight );
        fHM->hRecoInclusiveJetNEmF[dummyIter]->Fill( (*recoJetIter)->jtPfNEF(), weight );
        fHM->hRecoInclusiveJetNumOfConst[dummyIter]->Fill( numberOfConstituents, weight );
        fHM->hRecoInclusiveJetMUF[dummyIter]->Fill( (*recoJetIter)->jtPfMUF(), weight );
        fHM->hRecoInclusiveJetCHF[dummyIter]->Fill( (*recoJetIter)->jtPfCHF(), weight );
        fHM->hRecoInclusiveJetChargedMult[dummyIter]->Fill( chargedMult, weight );
        fHM->hRecoInclusiveJetCEmF[dummyIter]->Fill( (*recoJetIter)->jtPfCEF(), weight );
        fHM->hRecoInclusiveJetNumOfNeutPart[dummyIter]->Fill( neutralMult, weight );

        // Check if jet passes the selection criteria (*recoJet, isCM, isMC, requireMatching)
        if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, false, false) ) continue;

        //Inclusive (matched+unmatched) jets
        fHM->hRecoInclusiveAllJetPt->Fill(pt, weight);
        fHM->hRecoInclusiveAllJetEta->Fill(eta, weight);
        fHM->hRecoInclusiveAllJetEtaUnweighted->Fill(eta, 1.);
        fHM->hRecoInclusiveAllJetPtEta->Fill(eta, pt, weight);
        fHM->hRecoInclusiveAllJetPtRawEta->Fill(eta, ptRaw, weight);
        fHM->hRecoInclusiveAllJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

        // Inclusive (matched+unmatched) leading jet
        if ( (fRecoIdLead >= 0) && ((recoJetCounter - 1) == fRecoIdLead) ) {
            fHM->hRecoLeadAllJetPtEta->Fill( eta, pt, weight );
            fHM->hRecoLeadAllJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
        }

        // Inclusive (matched+unmatched) subleading jet
        if ( (fRecoIdSubLead >= 0) && ((recoJetCounter - 1) == fRecoIdSubLead) ) {
            fHM->hRecoSubLeadAllJetPtEta->Fill( eta, pt, weight );
            fHM->hRecoSubLeadAllJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
        }

        if ( pt > 30. ) {
            fHM->hRecoGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), weight );
            fHM->hRecoGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), weight );
        }

        // For the MC, check and study matching to gen reconstructed jets
        if ( fIsMc ) {

            // If reco jet has matching to gen jet
            if ( (*recoJetIter)->hasMatching() ) {

                GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
                if ( !matchedJet ) {
                    std::cerr << Form("Cannot retrieve gen jet with id: %d", (*recoJetIter)->genJetId() ) << std::endl;
                    continue;
                }
                
                float genPt = matchedJet->pt();
                float genEta = etaLab( matchedJet->eta() );
                float genPhi = matchedJet->phi();

                if ( fabs( genPt ) < 1e-6 ) {
                    std::cerr << "JetESRAnalysis::processRecoJets -- ref jet pt is zero. Skip processing." << std::endl;
                    continue;
                }

                if ( fabs( genPt ) > std::numeric_limits<float>::max() ||
                     fabs( genEta ) > std::numeric_limits<float>::max() ||
                     fabs( genPhi )  > std::numeric_limits<float>::max() ) {
                    std::cerr << "JetESRAnalysis::processRecoJets -- ref jet properties are crazy. Skip processing." << std::endl;
                    continue;
                }

                float JES = pt/genPt;

                double res1[4] = { JES, genPt, genEta, ptHat };
                double res2[4] = { JES, pt, eta, ptHat };

                // Reco 2 Ref correlations for inclusive jets
                // 0 - reco ptCorr, 1- reco eta, 2 - reco phi,
                // 3 - ref pt, 4 - ref eta, 5 - ref phi, 6 - ptHat
                double reco2refJetPtEtaPhiPtHat[7] = { pt, eta, phi, genPt, genEta, genPhi, ptHat };

                // Inclusive ref jets
                fHM->hRefInclusiveJetPt->Fill( genPt, weight );
                fHM->hRefInclusiveJetEta->Fill( genEta, weight );
                fHM->hRefInclusiveJetEtaUnweighted->Fill( genEta, 1. );
                fHM->hRefInclusiveJetPtEta->Fill( genEta, genPt, weight );
                fHM->hRefInclusiveJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );

                // Inclusive matched reco jets
                fHM->hRecoInclusiveMatchedJetPt->Fill( pt, weight );
                fHM->hRecoInclusiveMatchedJetPtEta->Fill( eta, pt, weight );
                fHM->hRecoInclusiveMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );

                fHM->hInclusiveReco2RefJetPtEtaPhiPtHat->Fill( reco2refJetPtEtaPhiPtHat, weight );

                // Fill JES vs pt for |eta| < 1.4 (midrapidity)
                if ( fabs( genEta ) < 1.4 ) {
                    fHM->hInclusiveJetJESVsPtGen->Fill( genPt, JES, weight );
                }
                fHM->hInclusiveJetJESGenPtGenEtaPtHatWeighted->Fill( res1, weight );
                fHM->hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Fill( res2, weight );

                // Leading jet
                if ( (fRecoIdLead >= 0) && ((recoJetCounter-1) == fRecoIdLead) ) {

                    fHM->hLeadingJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                    // Leading matched reco jet
                    fHM->hRecoLeadMatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoLeadMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                    // Leading ref jet
                    fHM->hRefLeadJetPtEta->Fill( genEta, genPt, weight );
                    fHM->hRefLeadJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );
                    fHM->hLeadReco2RefJetPtEtaPhiPtHat->Fill( reco2refJetPtEtaPhiPtHat, weight );
                    if ( (*recoJetIter)->genJetId() == fGenIdLead ) {
                        fHM->hRefLeadUnswappedJetPtEta->Fill( genEta, genPt, weight );
                        fHM->hRefLeadUnswappedJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );
                    }

                }

                // Subleading jet
                if ( (fRecoIdSubLead >= 0) && ((recoJetCounter-1) == fRecoIdSubLead) ) {
                    fHM->hSubleadingJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                    // Subleading matched reco jet
                    fHM->hRecoSubLeadMatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoSubLeadMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                    // Subleading ref jet
                    fHM->hRefSubLeadJetPtEta->Fill( genEta, genPt, weight );
                    fHM->hRefSubLeadJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );

                    fHM->hSubLeadReco2RefJetPtEtaPhiPtHat->Fill( reco2refJetPtEtaPhiPtHat, weight );

                    if ( (*recoJetIter)->genJetId() == fGenIdSubLead ) {
                        fHM->hRefSubLeadUnswappedJetPtEta->Fill( genEta, genPt, weight );
                        fHM->hRefSubLeadUnswappedJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );
                    }
                }
            } // if ( (*recoJetIter)->hasMatching() )
            else {

                // Fill unmatched reco jets
                fHM->hRecoInclusiveUnmatchedJetPt->Fill(pt, weight);
                fHM->hRecoInclusiveUnmatchedJetPtEta->Fill(eta, pt, weight);
                fHM->hRecoInclusiveUnmatchedJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

                // Leading unmatched reco jet
                if ( (fRecoIdLead >= 0) && ((recoJetCounter-1) == fRecoIdLead) ) {
                    fHM->hRecoLeadUnmatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoLeadUnmatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                }

                // Subleading unmatched reco jet
                if ( (fRecoIdSubLead >= 0) && ((recoJetCounter-1) == fRecoIdSubLead) ) {
                    fHM->hRecoSubLeadUnmatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoSubLeadUnmatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );

                }
            } // else
        } // if ( fIsMc )
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processRecoJets -- end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processRefJets(const Event* event, const double &weight) {

    if ( weight <= 0. ) {
        std::cout << "JetESRAnalysis::processRefJets -- weight is zero or negative. Skip processing." << std::endl;
        return;
    }

    // ptHat value
    float ptHat = event->ptHat();

    // Jet iterators
    RecoJetIterator recoJetIter;

    // Jet counter
    int refSelJetCounter{0};
    if ( event->recoJetCollection()->size() > 0 ) {

        // Loop over reconstructed jets and search for leading and subleading gen-matched jets
        for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

            refSelJetCounter++;

            // if ( fabs( (*recoJetIter)->phi() ) > std::numeric_limits<float>::max() ) continue;

            // Check selection criteria (*recoJet, isCM, isMC, requireMatching)
            if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, true, true) ) continue; 

            // Retrieve matched gen jet
            GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
            float genPt = matchedJet->pt();
            float genEta = etaLab( matchedJet->eta() );
            float genPhi = matchedJet->phi();

            if ( fVerbose ) {
                std::cout << Form("Ref jet #%d pt: %5.2f eta: %5.2f phi: %5.2f --> Reco jet #%d pt: %5.2f eta: %5.2f phi: %5.2f\n", 
                                  (*recoJetIter)->genJetId(), genPt, genEta, genPhi, refSelJetCounter-1, (*recoJetIter)->ptJECCorr(), (*recoJetIter)->eta(), (*recoJetIter)->phi());
            }

            // Inclusive ref-selected ref jet
            fHM->hRefSelInclusiveJetPt->Fill( genPt, weight );
            fHM->hRefSelInclusiveJetEta->Fill( genEta, weight );
            fHM->hRefSelInclusiveJetEtaUnweighted->Fill( genEta, 1. );
            fHM->hRefSelInclusiveJetPtEta->Fill(genEta, genPt, weight);
            fHM->hRefSelInclusiveJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight);

            // Leading ref-selected ref jet
            if ( (fRefSelRecoIdLead >= 0) && ((refSelJetCounter - 1) == fRefSelRecoIdLead) ) {
                fHM->hRefSelLeadJetPtEta->Fill(genEta, genPt, weight);
                fHM->hRefSelLeadJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight);
            }

            // Subleading ref-selected ref jet
            if ( (fRefSelRecoIdSubLead >= 0) && ((refSelJetCounter - 1) == fRefSelRecoIdSubLead) ) {
                fHM->hRefSelSubLeadJetPtEta->Fill(genEta, genPt, weight);
                fHM->hRefSelSubLeadJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight);
            }
        } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )
    } // if ( event->recoJetCollection()->size() > 0 )
    
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processRefJets -- end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processEvent(const Event* event) {

    // Perform the analysis
    if ( fVerbose ) {
        std::cout << "\n\n++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "JetESRAnalysis::processEvent -- begin" << std::endl;
    }

    if ( !fHM ) {
        std::cout << "[Warning] No histogram manager connected to the JetESRAnalysis\n";
        return;
    }

    // Must be flushed for each event !!!!
    fRecoIdLead = {-1}; 
    fRecoIdSubLead = {-1}; 
    fGenIdLead = {-1}; 
    fGenIdSubLead = {-1}; 
    fRefSelRecoIdLead = {-1}; 
    fRefSelRecoIdSubLead = {-1};
    
    fRecoPtSortedJetIds.clear();
    fGenPtSortedJetIds.clear();
    fRefSelRecoPtSortedJetIds.clear();

    //
    // Event quantities
    //

    // ptHat
    float ptHat = event->ptHat();
    // Vertex z position
    float vz = event->vz();
    // ptHat weight 
    float ptHatW = event->ptHatWeight();
    // Centrality
    float centrality = event->centrality();
    // Centrality weight
    float centW = event->centralityWeight();
    if ( fCollisionSystem != 2 ) { // Apply centrality weight only for PbPb
        centW = 1.;
    }
    // Final weight
    double weight{1.};

    // Check correctness of MC sample for pPb 8160
    if ( fCollisionSystem == 1 ) { 
        // Skip events with ptHat that is outside the ranged embedded
        if ( ptHat <= fPtHatRange[0] || ptHat > fPtHatRange[1] ) {
            if ( fVerbose ) {
                std::cout << Form("[WARNING] Bad ptHat value: %4.1f < %4.1f <= %4.1f\n", fPtHatRange[0], ptHat, fPtHatRange[1]);
            }
            return;
        }

        // For MC we need to flip the direction of Pb-going in order to properly reweight distributions
        if ( fIsPbGoingDir ) {
            vz = -vz;
        }
    } // if ( fCollisionSystem == 1 )

    // Calculate event weight
    weight = eventWeight(ptHat, vz, centW, ptHatW);

    //
    // Create pT-sorted jet indices
    //
    makePtSortedJetVectors( event );

    // Check for event overweight in both reco and gen jets 
    // (x-jets and purelly overweighted gen jets) in Monte Carlo
    if ( fIsMc ) {
        bool overweight = isOverweightedEvent( event, weight );
        if ( overweight ) {
            if ( fVerbose ) {
                std::cout << "Overweighted event. Skip it." << std::endl;
            }
            return;
        }
    } // if ( fIsMc )

    //
    // Loop over jet collections (reco, gen, ref)
    //
    processInclusiveJets(event, weight);

    // Fill event histograms
    fHM->hVz->Fill( vz,  1. );
    fHM->hVzPtHatWeighted->Fill( vz, ptHatW);
    fHM->hVzWeighted->Fill( vz, weight );

    // fHM->hVzCentWeighted->Fill( vz, centW );

    fHM->hPtHat->Fill( ptHat, 1. );
    fHM->hPtHatWeighted->Fill( ptHat, weight );

    fHM->hHiBin->Fill( event->hiBin(), 1. );
    fHM->hHiBinWeighted->Fill( event->hiBin(), weight );

    fHM->hCentrality->Fill( centrality, 1. );
    // fHM->hCentralityPtHatWeighted->Fill( centrality, ptHatW );
    fHM->hCentralityWeighted->Fill( centrality, weight );

    fHM->hPtHatWeight->Fill( ptHatW, 1. );
    fHM->hPtHatWeightWeighted->Fill( ptHatW, weight );

    // Fill histograms with vz, ptHat and centrality
    fHM->hVzPtHatCent->Fill( vz, ptHat, centrality, 1. );
    // fHM->hVzPtHatCentPtHatWeighted->Fill( vz, ptHat, centrality, ptHatW );
    fHM->hVzPtHatCentWeighted->Fill( vz, ptHat, centrality, weight );



    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processEvent - end" << std::endl;
    }
}

//________________
void JetESRAnalysis::finish() {
    std::cout << "JetESRAnalysis::finish" << std::endl;
}

//________________
void JetESRAnalysis::report() {
    // Force to report everyone
}

//________________
TList* JetESRAnalysis::getOutputList() {
    TList *outputList = new TList();

    // Add list of settings for cuts
    return outputList;
}
