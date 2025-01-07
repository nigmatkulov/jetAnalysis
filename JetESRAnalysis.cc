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

// Jet analysis headers
#include "JetESRAnalysis.h"

//________________
JetESRAnalysis::JetESRAnalysis() : BaseAnalysis(), 
    fVerbose{false}, fHM{nullptr}, fUseCentralityWeight{false}, fEtaShift{0},
    fCollisionEnergy{8160}, fCollisionSystem{2}, fIsPbGoingDir{false}, fPtHatRange{15., 30.},
    fUseJetIdSelection{false}, fIsLooseJetIdCut{false}, 
    fLeadJetPtLow{50.}, fLeadJetEta{-3., 3.}, fSubleadJetPtLow{40.}, fSubleadJetEta{-3., 3.},
    fNEventsInSample{1000000}, fVzWeight{nullptr} {
    /* Empty */
}

//________________
JetESRAnalysis::~JetESRAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void JetESRAnalysis::init() {
    // Initialize analysis

    // Initialize vz weight function
    initVzWeightFunction();

    // Print initial setup of the analysis
    print();
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
              << "Collisions system           : " << collisionSystem().Data() << std::endl
              << "Is Pb-going direction       : " << fIsPbGoingDir << std::endl
              << "eta shift                   : " << fEtaShift << std::endl
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl
              << "Use jetId selection         : " << fUseJetIdSelection << std::endl
              << "Use loose jetId cut         : " << fIsLooseJetIdCut << std::endl
              << "Leading jet pT              : " << fLeadJetPtLow << std::endl
              << "Leading jet eta range       : " << fLeadJetEta[0] << "-" << fLeadJetEta[1] << std::endl
              << "SubLeading jet pT           : " << fSubleadJetPtLow << std::endl
              << "SubLeading jet eta range    : " << fSubleadJetEta[0] << "-" << fSubleadJetEta[1] << std::endl
              << "Dijet phi cut               : " << fDijetDPhiCut << std::endl;
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
void JetESRAnalysis::findLeadSubleadJets(const double &pt, const int &counter,
                                         double &ptLead, double &ptSublead,
                                         int &idLead, int &idSubLead) {
    // Find leading and subleading jets
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::findLeadSubleadJets - begin" << std::endl;
    }

    if ( pt > ptLead ) {
        ptSublead = ptLead;
        idSubLead = idLead;
        ptLead = pt;
        idLead = counter;
    }
    else if ( pt > ptSublead ) {
        ptSublead = pt;
        idSubLead = counter;
    }

    if ( fVerbose ) {
        std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f Lead id: %d SubLead id: %d\n", 
                          ptLead, ptSublead, idLead, idSubLead);
        std::cout << "JetESRAnalysis::findLeadSubleadJets - end" << std::endl;
    }
}

//________________
bool JetESRAnalysis::isGoodGenJet(const GenJet* jet) {
    bool goodJet{false};
    double etaCut[2] { fSubleadJetEta[0], fSubleadJetEta[1] }; 
    double eta = jet->eta();

    // if ( fSelectJetsInCMFrame ) {

    //     eta = boostEta2CM( eta );

    //     etaCut[0] = fJetEtaCM[0];
    //     etaCut[1] = fJetEtaCM[1];
    // }
    // else {
    //     eta = etaLab( eta );
    // }

    if ( jet->pt() > fSubleadJetPtLow && etaCut[0] < eta && eta < etaCut[1] ) {
        goodJet = {true};
    }
    
    if ( fVerbose ) {
        std::cout << Form("Gen jet cut %s\n", goodJet ? "\t[passed]" : "\t[failed]" );
    }

    return goodJet;
}

//________________
bool JetESRAnalysis::isGoodRecoJet(const RecoJet* jet) {
    bool goodJet{false};
    bool goodKine{false};
    bool hasMatching{false};

    double etaCut[2] { fSubleadJetEta[0], fSubleadJetEta[1] };
    double eta = jet->eta();

    // if ( fSelectJetsInCMFrame ) {

    //     eta = boostEta2CM( eta );

    //     etaCut[0] = fJetEtaCM[0]; 
    //     etaCut[1] = fJetEtaCM[1];
    // }
    // else {
    //     eta = etaLab( eta );
    // }

    if ( jet->ptJECCorr() > 15. && etaCut[0] < eta && eta < etaCut[1] ) {
        goodKine = {true};
    }

    if ( jet->hasMatching() ) {
        hasMatching = {true};
    }

    goodJet = goodKine && hasMatching;

    if ( fVerbose ) {
        std::cout << Form("Reco jet cut %s", goodJet ? "\t[passed]" : "\t[failed]"); 
        if ( goodJet ) {
            std::cout << std::endl;
        }
        else {
            std::cout << Form("\t goodKine: %d hasMatching: %d\n", goodKine, hasMatching);
        }
    } // if ( fVerbose )

    return goodJet;
}

//________________
double JetESRAnalysis::deltaPhi(const double& phi1, const double &phi2) {
    double dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ) dphi -= TMath::TwoPi();
    if ( dphi < -TMath::Pi() ) dphi += TMath::TwoPi();
    return dphi;
}

//________________
bool JetESRAnalysis::isGoodDijet(const double& ptLead, const double& etaLead, 
                                 const double& ptSubLead, const double& etaSubLead,
                                 const double& dphi) {
    bool isGood = ( ptLead > fLeadJetPtLow &&
                    fLeadJetEta[0] < etaLead && etaLead < fLeadJetEta[1] &&
                    ptSubLead > fSubleadJetPtLow && 
                    fSubleadJetEta[0] < etaSubLead && etaSubLead < fSubleadJetEta[1] &&
                    dphi > fDijetDPhiCut );
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::isGoodDijet " << ( (isGood) ? "[true]" : "[false]");
        std::cout << Form("\tpTlead: %5.2f etaLead: %5.2f pTsub: %5.2f etaSub: %5.2f dphi: %4.2f\n", 
                          ptLead, etaLead, ptSubLead, etaSubLead, dphi);
    }
    return isGood;
}

//________________
bool JetESRAnalysis::isGoodTrkMax(const RecoJet* jet) {
    bool goodTrackMax = {true};
    double rawPt = jet->rawPt();
    double trackMaxPt = jet->trackMaxPt();
    if ( TMath::Abs( jet->eta() ) < 2.4 && 
         ( trackMaxPt/rawPt < 0.01 ||
           trackMaxPt/rawPt > 0.98) ) {
        goodTrackMax = {kFALSE};
    }

    if ( fVerbose ) {
        std::cout << "TrackMaxPt/rawPt: " << trackMaxPt/rawPt << ( (goodTrackMax) ? " [passed]" : " [failed]" ) 
                  << ( (trackMaxPt/rawPt < 0.01) ? " too low value " : "" ) << ( (trackMaxPt/rawPt > 0.98) ? " too large value " : "" )
                  << std::endl;
    }

    return goodTrackMax;
}

//_________________
bool JetESRAnalysis::isGoodJetId(const RecoJet* jet) {

    bool passJetId = {false};

    int chm = jet->jtPfCHM();
    int cem = jet->jtPfCEM();
    int mum = jet->jtPfMUM();
    int nhm = jet->jtPfNHM();
    int nem = jet->jtPfNEM();

    float chf = jet->jtPfCHF();
    float cef = jet->jtPfCEF();
    float nhf = jet->jtPfNHF();
    float nef = jet->jtPfNEF();
    float muf = jet->jtPfMUF();

    int chargedMult = chm + cem + mum;
    int neutralMult = nhm + nem;
    int numberOfConstituents = chargedMult + neutralMult;

    float eta = jet->eta();

    float chargedEmFracCut{1.}, neutFracCut{1.};
    if ( !fIsLooseJetIdCut ) {
        chargedEmFracCut = {0.9};
        neutFracCut = {0.9};
    }
    else {
        chargedEmFracCut = {0.99};
        neutFracCut = {0.99};
    }

    bool passNHF{false};
    bool passNEF{false};
    bool passNumOfConstituents{true};
    bool passMuonFrac{true};
    bool passChargedFrac{true};
    bool passChargedMult{true};
    bool passChargedEmFrac{true};
    bool passNeutralMult{true};

    // Check cuts depending on jet pseudorapdity
    if ( TMath::Abs( eta ) <= 2.7 ) {
        
        passNHF = ( nhf < neutFracCut ) ? true : false;
        passNEF = ( nhf < neutFracCut ) ? true : false;
        passNumOfConstituents = ( numberOfConstituents > 1 ) ? true : false;

        if ( !fIsLooseJetIdCut ) { 
            passMuonFrac = ( muf < 0.8 ) ? true : false; 
        } // if ( !fIsLooseJetIdCut )

        if( TMath::Abs( eta ) <= 2.4 ) {
            passChargedFrac = ( chf > 0 ) ? true : false;
            passChargedMult = ( chargedMult > 0 ) ? true : false;
            passChargedEmFrac = ( cef < chargedEmFracCut ) ? true : false;
        } // if( TMath::Abs( eta ) <= 2.4 )

    } // if ( TMath::Abs( eta ) <= 2.7 )
    else if ( TMath::Abs( eta ) <= 3.0) {

        passNEF = ( nef > 0.01 ) ? true : false;
        passNHF = ( nhf < 0.98 ) ? true : false;
        passNeutralMult = ( neutralMult > 2 ) ? true : false;

    } // else if ( TMath::Abs( eta ) <= 3.0)
    else  {
        passNEF = ( nef < 0.9 ) ? true : false;
        passNeutralMult = (neutralMult > 10 ) ? true : false; // CAUTION: JET MET says it should be >10
    } // else 

    passJetId = passNHF && passNEF && passNumOfConstituents && passMuonFrac && 
                passChargedFrac && passChargedMult && passChargedEmFrac && passNeutralMult;

    if ( fVerbose ) {
        std::cout << "JetId selection results: " << ( (passJetId) ? "[passed]" : "[failed]" ) << " Reasons ";
        std::cout << Form("passNHF: %d \tpassNEF: %d \tpassNumConst: %d \tpassMuonFrac: %d \tpassChFrac: %d \tpassChMult: %d \tpassChEmFrac: %d \tpassNeutMult: %d\n", 
                            passNHF, passNEF, passNumOfConstituents, passMuonFrac, passChargedFrac, 
                            passChargedMult , passChargedEmFrac , passNeutralMult);
    }
        
    return passJetId;
}


//________________
void JetESRAnalysis::processGenJets(const Event* event, double weight) {
    // Process and analyze gen jets
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processGenJets - begin" << std::endl;
    }

    // Variables for dijets
    double ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
           phiLead{0.}, phiSubLead{0.};
    int idLead{-1}, idSubLead{-1};

    GenJetIterator genJetIter;
    int counter{0};

    int centrality = event->centrality(); // -1 if no centrality available
    double ptHat = event->ptHat();

    // Loop over generated jets
    for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

        if ( fVerbose ) {
            std::cout << "Gen jet #" << counter << " ";
            (*genJetIter)->print();
        }

        double pt = (*genJetIter)->pt();
        double eta = (*genJetIter)->eta();
        double phi = (*genJetIter)->phi();
        int flavB{-6};
        switch ( (*genJetIter)->flavorForB() )
        {
        case -99:
            flavB = {-6};
            break;
        case -5:
            flavB = {-5};
            break;
        case -4:
            flavB = {-4};
            break;
        case -3:
            flavB = {-3};
            break;
        case -2:
            flavB = {-2};
            break;
        case -1:
            flavB = {-1};
            break;
        case 0:
            flavB = {0};
            break;
        case 1:
            flavB = {1};
            break;
        case 2:
            flavB = {2};
            break;
        case 3:
            flavB = {3};
            break;
        case 4:
            flavB = {4};
            break;
        case 5:
            flavB = {5};
            break;
        case 21:
            flavB = {6};
            break;
        default:
            flavB = {-6};
            break;
        }

        // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
        double jetPtEtaPhiFlavPtHatCent[6] = { pt, eta, phi, (double)flavB, ptHat, (double)centrality };

        // Fill inclusive jet histograms
        fHM->hGenInclusiveJetPt->Fill( pt, 1. );
        fHM->hGenInclusiveJetPtWeighted->Fill( pt, weight );
        fHM->hGenInclusiveJetEtaPt->Fill( eta, pt, 1. );
        fHM->hGenInclusiveJetEtaPtWeighted->Fill( eta, pt, weight );
        fHM->hGenInclusiveJetPtEtaPhiFlavPtHatCent->Fill( jetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Fill( jetPtEtaPhiFlavPtHatCent, weight );

        // Search for leading and subleading jets
        findLeadSubleadJets( pt, counter, ptLead, ptSubLead, idLead, idSubLead );

        counter++;
    } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )

    // Check if leading and subleading jets were found
    if ( idLead >= 0 && idSubLead >= 0 ) {
        if ( fVerbose ) {
            std::cout << "Gen dijet found: [true]" << std::endl;
        }

        // Retrieve leading jet
        GenJet *leadJet = event->genJetCollection()->at( idLead );
        
        ptLead = leadJet->pt();
        etaLead = leadJet->eta();
        phiLead = leadJet->phi();
        int flavBLead{-6};
        switch ( leadJet->flavorForB() )
        {
        case -99:
            flavBLead = {-6};
            break;
        case -5:
            flavBLead = {-5};
            break;
        case -4:
            flavBLead = {-4};
            break;
        case -3:
            flavBLead = {-3};
            break;
        case -2:
            flavBLead = {-2};
            break;
        case -1:
            flavBLead = {-1};
            break;
        case 0:
            flavBLead = {0};
            break;
        case 1:
            flavBLead = {1};
            break;
        case 2:
            flavBLead = {2};
            break;
        case 3:
            flavBLead = {3};
            break;
        case 4:
            flavBLead = {4};
            break;
        case 5:
            flavBLead = {5};
            break;
        case 21:
            flavBLead = {6};
            break;
        default:
            flavBLead = {-6};
            break;
        }
        
        // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
        double leadJetPtEtaPhiFlavPtHatCent[6] = { ptLead, etaLead, phiLead, (double)flavBLead, ptHat, (double)centrality };
        fHM->hGenInclusiveLeadJetPt->Fill( ptLead, 1. );
        fHM->hGenInclusiveLeadJetPtWeighted->Fill( ptLead, weight );
        fHM->hGenInclusiveLeadJetEtaPt->Fill( etaLead, ptLead, 1. );
        fHM->hGenInclusiveLeadJetEtaPtWeighted->Fill( etaLead, ptLead, weight );
        fHM->hGenInclusiveLeadJetPtEtaPhiFlavPtHatCent->Fill( leadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( leadJetPtEtaPhiFlavPtHatCent, weight );

        // Retrieve subleading jet
        GenJet *subLeadJet = event->genJetCollection()->at( idSubLead );

        ptSubLead = subLeadJet->pt();
        etaSubLead = subLeadJet->eta();
        phiSubLead = subLeadJet->phi();
        int flavBSubLead{-6};
        switch ( subLeadJet->flavorForB() )
        {
        case -99:
            flavBSubLead = {-6};
            break;
        case -5:
            flavBSubLead = {-5};
            break;
        case -4:
            flavBSubLead = {-4};
            break;
        case -3:
            flavBSubLead = {-3};
            break;
        case -2:
            flavBSubLead = {-2};
            break;
        case -1:
            flavBSubLead = {-1};
            break;
        case 0:
            flavBSubLead = {0};
            break;
        case 1:
            flavBSubLead = {1};
            break;
        case 2:
            flavBSubLead = {2};
            break;
        case 3:
            flavBSubLead = {3};
            break;
        case 4:
            flavBSubLead = {4};
            break;
        case 5:
            flavBSubLead = {5};
            break;
        case 21:
            flavBSubLead = {6};
            break;
        default:
            flavBSubLead = {-6};
            break;
        }

        // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
        double subLeadJetPtEtaPhiFlavPtHatCent[6] = { ptSubLead, etaSubLead, phiSubLead, (double)flavBSubLead, ptHat, (double)centrality };
        fHM->hGenInclusiveSubLeadJetPt->Fill( ptSubLead, 1. );
        fHM->hGenInclusiveSubLeadJetPtWeighted->Fill( ptSubLead, weight );
        fHM->hGenInclusiveSubLeadJetEtaPt->Fill( etaSubLead, ptSubLead, 1. );
        fHM->hGenInclusiveSubLeadJetEtaPtWeighted->Fill( etaSubLead, ptSubLead, weight );
        fHM->hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Fill( subLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( subLeadJetPtEtaPhiFlavPtHatCent, weight );

        // Calculate dijet quantities
        double dijetDphi = deltaPhi( phiLead, phiSubLead );
        double dijetEta = 0.5 * ( etaLead + etaSubLead );
        double dijetPt = 0.5 * ( ptLead + ptSubLead );
        double dijetDetaCM = 0.5 * ( etaLead - etaSubLead );
        double x_Pb = 2. * dijetPt / fCollisionEnergy * TMath::Exp( -1. * dijetDetaCM ) * TMath::CosH( dijetDetaCM );
        double x_p = 2. * dijetPt / fCollisionEnergy * TMath::Exp( dijetDetaCM ) * TMath::CosH( dijetDetaCM );
        double xPbOverXp = x_Pb / x_p;

        fHM->hGenInclusiveDijetDphi->Fill( dijetDphi, 1. );
        fHM->hGenInclusiveDijetDphiWeighted->Fill( dijetDphi, weight );
        fHM->hGenInclusiveDijetEtaPt->Fill( dijetEta, dijetPt, 1. );
        fHM->hGenInclusiveDijetEtaPtWeighted->Fill( dijetEta, dijetPt, weight );
        fHM->hGenInclusiveDijetEtaPtDphi->Fill( dijetEta, dijetPt, dijetDphi, 1. );
        fHM->hGenInclusiveDijetEtaPtDphiWeighted->Fill( dijetEta, dijetPt, dijetDphi, weight );
        fHM->hGenInclusiveDijetDetaCM->Fill( dijetDetaCM, 1. );
        fHM->hGenInclusiveDijetDetaCMWeighted->Fill( dijetDetaCM, weight );
        fHM->hGenInclusiveDijetDetaCMPt->Fill( dijetDetaCM, dijetPt, 1. );
        fHM->hGenInclusiveDijetDetaCMPtWeighted->Fill( dijetDetaCM, dijetPt, weight );
        fHM->hGenInclusiveDijetEtaDetaCMPt->Fill( dijetEta, dijetDetaCM, dijetPt, 1. );
        fHM->hGenInclusiveDijetEtaDetaCMPtWeighted->Fill( dijetEta, dijetDetaCM, dijetPt, weight );
        fHM->hGenInclusiveDijetXPb->Fill( x_Pb, 1. );
        fHM->hGenInclusiveDijetXPbWeighted->Fill( x_Pb, weight );
        fHM->hGenInclusiveDijetXp->Fill( x_p, 1. );
        fHM->hGenInclusiveDijetXpWeighted->Fill( x_p, weight );
        fHM->hGenInclusiveDijetXPbOverXp->Fill( xPbOverXp, 1. );
        fHM->hGenInclusiveDijetXPbOverXpWeighted->Fill( xPbOverXp, weight );
        fHM->hGenInclusiveDijetXPbOverXpEta->Fill( xPbOverXp, dijetEta, 1. );
        fHM->hGenInclusiveDijetXPbOverXpEtaWeighted->Fill( xPbOverXp, dijetEta, weight );

        if ( !isGoodDijet( ptLead, etaLead, ptSubLead, etaSubLead, dijetDphi ) ) {
            return;
        }

        // Fill selected leading jet histograms
        fHM->hGenSelectedLeadJetPt->Fill( ptLead, 1. );
        fHM->hGenSelectedLeadJetPtWeighted->Fill( ptLead, weight );
        fHM->hGenSelectedLeadJetEtaPt->Fill( etaLead, ptLead, 1. );
        fHM->hGenSelectedLeadJetEtaPtWeighted->Fill( etaLead, ptLead, weight );
        fHM->hGenSelectedLeadJetPtEtaPhiFlavPtHatCent->Fill( leadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( leadJetPtEtaPhiFlavPtHatCent, weight );

        // Fill selected subleading jet histograms
        fHM->hGenSelectedSubLeadJetPt->Fill( ptSubLead, 1. );
        fHM->hGenSelectedSubLeadJetPtWeighted->Fill( ptSubLead, weight );
        fHM->hGenSelectedSubLeadJetEtaPt->Fill( etaSubLead, ptSubLead, 1. );
        fHM->hGenSelectedSubLeadJetEtaPtWeighted->Fill( etaSubLead, ptSubLead, weight );
        fHM->hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Fill( subLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( subLeadJetPtEtaPhiFlavPtHatCent, weight );

        // Fill selected dijet histograms
        fHM->hGenSelectedDijetPt->Fill( dijetPt, 1. );
        fHM->hGenSelectedDijetPtWeighted->Fill( dijetPt, weight );
        fHM->hGenSelectedDijetEta->Fill( dijetEta, 1. );
        fHM->hGenSelectedDijetEtaWeighted->Fill( dijetEta, weight );
        fHM->hGenSelectedDijetEtaPt->Fill( dijetEta, dijetPt, 1. );
        fHM->hGenSelectedDijetEtaPtWeighted->Fill( dijetEta, dijetPt, weight );
        fHM->hGenSelectedDijetEtaPtDphi->Fill( dijetEta, dijetPt, dijetDphi, 1. );
        fHM->hGenSelectedDijetEtaPtDphiWeighted->Fill( dijetEta, dijetPt, dijetDphi, weight );
        fHM->hGenSelectedDijetDetaCM->Fill( dijetDetaCM, 1. );
        fHM->hGenSelectedDijetDetaCMWeighted->Fill( dijetDetaCM, weight );
        fHM->hGenSelectedDijetDetaCMPt->Fill( dijetDetaCM, dijetPt, 1. );
        fHM->hGenSelectedDijetDetaCMPtWeighted->Fill( dijetDetaCM, dijetPt, weight );
        fHM->hGenSelectedDijetEtaDetaCMPt->Fill( dijetEta, dijetDetaCM, dijetPt, 1. );
        fHM->hGenSelectedDijetEtaDetaCMPtWeighted->Fill( dijetEta, dijetDetaCM, dijetPt, weight );
        fHM->hGenSelectedDijetXPb->Fill( x_Pb, 1. );
        fHM->hGenSelectedDijetXPbWeighted->Fill( x_Pb, weight );
        fHM->hGenSelectedDijetXp->Fill( x_p, 1. );
        fHM->hGenSelectedDijetXpWeighted->Fill( x_p, weight );
        fHM->hGenSelectedDijetXPbOverXp->Fill( xPbOverXp, 1. );
        fHM->hGenSelectedDijetXPbOverXpWeighted->Fill( xPbOverXp, weight );
        fHM->hGenSelectedDijetXPbOverXpEta->Fill( xPbOverXp, dijetEta, 1. );
        fHM->hGenSelectedDijetXPbOverXpEtaWeighted->Fill( xPbOverXp, dijetEta, weight );

    } // if ( idLead >= 0 && idSubLead >= 0 )
    else {
        if ( fVerbose ) {
            std::cout << "Gen dijet found: [false]" << std::endl;
        }
    } // else

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processGenJets - end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processRecoJets(const Event* event, double weight) {
    // Process and analyze reco jets
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processRecoJets - begin" << std::endl;
    }

    // Define variables
    double ptRecoLead{-1.}, ptRecoSubLead{-1.},
           ptRawRecoLead{-1.}, ptRawRecoSubLead{-1.},
           etaRecoLead{0.}, etaRecoSubLead{0.},
           phiRecoLead{0.},  phiRecoSubLead{0.}, 
           ptRefLead{-1.}, ptRefSubLead{-1.},
           etaRefLead{0.}, etaRefSubLead{0.},
           phiRefLead{0.}, phiRefSubLead{0.};
    int idRecoLead{-1}, idRecoSubLead{-1};

    int centrality = event->centrality(); // -1 if no centrality available
    double ptHat = event->ptHat();

    // Loop over reconstructed jets
    RecoJetIterator recoJetIter;
    int counter{0};
    double pt{-999.}, eta{-999.}, phi{-999.}, ptRaw{-999.};
    int chargedMult{-1}, neutralMult{-1}, numberOfConstituents{-1};

    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        pt = (*recoJetIter)->ptJECCorr();
        eta = (*recoJetIter)->eta();
        phi = (*recoJetIter)->phi();
        ptRaw = (*recoJetIter)->pt();

        if ( fVerbose ) {
            std::cout << "Reco jet #" << counter << " ";
            (*recoJetIter)->print();
        }

        // JetId parameters and histograms
        chargedMult = (*recoJetIter)->jtPfCHM() + (*recoJetIter)->jtPfCEM() + (*recoJetIter)->jtPfMUM();
        neutralMult = (*recoJetIter)->jtPfNHM() + (*recoJetIter)->jtPfNEM();
        numberOfConstituents = chargedMult + neutralMult;

        // 4 cases for jet eta: 0 - |eta| < 2.4, 1 - 2.4 < |eta| < 2.7, 2 - 2.7 < |eta| < 3.0, 3 - 3.0 < |eta|
        int dummyIter{0};
        if ( TMath::Abs( eta ) <= 2.4 ) { dummyIter = {0}; }
        else if ( TMath::Abs( eta ) <= 2.7 ) { dummyIter = {1}; }
        else if ( TMath::Abs( eta ) <= 3.0 ) { dummyIter = {2}; }
        else { dummyIter = {3}; }

        // Fill jet with 0 - pt, 1 - eta, 2 - phi, 3 - dummyIter, 4 - NHF, 5 - NEF, 6 - num of constituents, 
        // 7 - MUF, 8 - CHF, 9 - charged mult, 10 - CEF, 11 - neutral mult, 12 - ptHat, 13 - centrality
        double jetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent[14] = 
        { pt, eta, phi, (double)dummyIter, 
          (*recoJetIter)->jtPfNHF(), 
          (*recoJetIter)->jtPfNEF(), 
          (double)numberOfConstituents, 
          (*recoJetIter)->jtPfMUF(), 
          (*recoJetIter)->jtPfCHF(), 
          (double)chargedMult, 
          (*recoJetIter)->jtPfCEF(), 
          (double)neutralMult, 
          ptHat, (double)centrality 
        };
        
        fHM->hRecoInclusiveJetNHF[dummyIter]->Fill( (*recoJetIter)->jtPfNHF(), weight );
        fHM->hRecoInclusiveJetNEmF[dummyIter]->Fill( (*recoJetIter)->jtPfNEF(), weight );
        fHM->hRecoInclusiveJetNumOfConst[dummyIter]->Fill( numberOfConstituents, weight );
        fHM->hRecoInclusiveJetMUF[dummyIter]->Fill( (*recoJetIter)->jtPfMUF(), weight );
        fHM->hRecoInclusiveJetCHF[dummyIter]->Fill( (*recoJetIter)->jtPfCHF(), weight );
        fHM->hRecoInclusiveJetChargedMult[dummyIter]->Fill( chargedMult, weight );
        fHM->hRecoInclusiveJetCEmF[dummyIter]->Fill( (*recoJetIter)->jtPfCEF(), weight );
        fHM->hRecoInclusiveJetNumOfNeutPart[dummyIter]->Fill( neutralMult, weight );
        fHM->hRecoInclusiveJetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent->Fill( jetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent, weight );

        // Check jet selection, e.g. jetId or trackMaxPt/rawPt
        bool goodTrackMax = isGoodTrkMax( (*recoJetIter) );
        bool goodJetId = isGoodJetId( (*recoJetIter) );
        if ( fUseJetIdSelection && !goodJetId ) { if ( fVerbose ) { std::cout << "JetId selection failed. Skip jet" << std::endl; } continue; }
        if ( !fUseJetIdSelection && !goodTrackMax ) { if ( fVerbose ) { std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; } continue; }

        // Fill inclusive jet histograms
        fHM->hRecoInclusiveJetPt->Fill( pt, 1. );
        fHM->hRecoInclusiveJetPtWeighted->Fill( pt, weight );
        fHM->hRecoInclusiveJetEtaPt->Fill( eta, pt, 1. );
        fHM->hRecoInclusiveJetEtaPtWeighted->Fill( eta, pt, weight );
        fHM->hRecoInclusiveJetPtRawPtCorrEta->Fill( ptRaw, pt, eta, 1. );
        fHM->hRecoInclusiveJetPtRawPtCorrEtaWeighted->Fill( ptRaw, pt, eta, weight );

        // Look at the unmatched jets
        if ( !(*recoJetIter)->hasMatching() ) {
            // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
            double jetPtEtaPhiFlavPtHatCent[6] = { pt, eta, phi, (double)-6, ptHat, (double)centrality };

            fHM->hRecoInclusiveUnmatchedJetPt->Fill( pt, 1. );
            fHM->hRecoInclusiveUnmatchedJetPtWeighted->Fill( pt, weight );
            fHM->hRecoInclusiveUnmatchedJetEtaPt->Fill( eta, pt, 1. );
            fHM->hRecoInclusiveUnmatchedJetEtaPtWeighted->Fill( eta, pt, weight );
            fHM->hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCent->Fill( jetPtEtaPhiFlavPtHatCent, 1. );
            fHM->hRecoInclusiveUnmatchedJetPtEtaPhiFlavPtHatCentWeighted->Fill( jetPtEtaPhiFlavPtHatCent, weight );
        }
        else { // Matched jets
            // Retrieve matched gen jet a.k.a. ref jet
            GenJet *refJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
            double refPt = refJet->pt();
            double refEta = refJet->eta();
            double refPhi = refJet->phi();
            int flavB{-6};
            switch ( refJet->flavorForB() )
            {
            case -99:
                flavB = {-6};
                break;
            case -5:
                flavB = {-5};
                break;
            case -4:
                flavB = {-4};
                break;
            case -3:
                flavB = {-3};
                break;
            case -2:
                flavB = {-2};
                break;
            case -1:
                flavB = {-1};
                break;
            case 0:
                flavB = {0};
                break;
            case 1:
                flavB = {1};
                break;
            case 2:
                flavB = {2};
                break;
            case 3:
                flavB = {3};
                break;
            case 4:
                flavB = {4};
                break;
            case 5:
                flavB = {5};
                break;
            case 21:
                flavB = {6};
                break;
            default:
                flavB = {-6};
                break;
            }

            // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
            double recoJetPtEtaPhiFlavPtHatCent[6] = { pt, eta, phi, (double)flavB, ptHat, (double)centrality };
            double refJetPtEtaPhiFlavPtHatCent[6] = { refPt, refEta, refPhi, (double)flavB, ptHat, (double)centrality };

            fHM->hRecoInclusiveMatchedJetPt->Fill( pt, 1. );
            fHM->hRecoInclusiveMatchedJetPtWeighted->Fill( pt, weight );
            fHM->hRecoInclusiveMatchedJetEtaPt->Fill( eta, pt, 1. );
            fHM->hRecoInclusiveMatchedJetEtaPtWeighted->Fill( eta, pt, weight );
            fHM->hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCent->Fill( recoJetPtEtaPhiFlavPtHatCent, 1. );
            fHM->hRecoInclusiveMatchedJetPtEtaPhiFlavPtHatCentWeighted->Fill( recoJetPtEtaPhiFlavPtHatCent, weight );

            // Jet energy scale
            double JES = pt / refPt;
            double jetJESPtEtaPhiCent[5] = { JES, refPt, refEta, refPhi, (double)centrality };
            fHM->hRecoInclusiveJetJESPtEtaPhiCent->Fill( jetJESPtEtaPhiCent, 1. );
            fHM->hRecoInclusiveJetJESPtEtaPhiCentWeighted->Fill( jetJESPtEtaPhiCent, weight );

            fHM->hRefInclusiveJetPt->Fill( refPt, 1. );
            fHM->hRefInclusiveJetPtWeighted->Fill( refPt, weight );
            fHM->hRefInclusiveJetEtaPt->Fill( refEta, refPt, 1. );
            fHM->hRefInclusiveJetEtaPtWeighted->Fill( refEta, refPt, weight );
            fHM->hRefInclusiveJetPtEtaPhiFlavPtHatCent->Fill( refJetPtEtaPhiFlavPtHatCent, 1. );
            fHM->hRefInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Fill( refJetPtEtaPhiFlavPtHatCent, weight );

            // Correlations between reco and ref jets

            // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
            double reco2refPtEtaPhiPtHatCentrality[8] = { pt, refPt, eta, refEta, phi, refPhi, ptHat, (double)centrality };
            fHM->hReco2RefInclusiveJetPt->Fill( refPt, pt, 1. );
            fHM->hReco2RefInclusiveJetPtWeighted->Fill( refPt, pt, weight );
            fHM->hReco2RefInclusiveJetEta->Fill( refEta, eta, 1. );
            fHM->hReco2RefInclusiveJetEtaWeighted->Fill( refEta, eta, weight );
            fHM->hReco2RefInclusiveJetPhi->Fill( refPhi, phi, 1. );
            fHM->hReco2RefInclusiveJetPhiWeighted->Fill( refPhi, phi, weight );
            fHM->hReco2RefInclusiveJetPtEtaPhiPtHatCentrality->Fill( reco2refPtEtaPhiPtHatCentrality, 1. );
            fHM->hReco2RefInclusiveJetPtEtaPhiPtHatCentralityWeighted->Fill( reco2refPtEtaPhiPtHatCentrality, weight );
        } // else 

        // Search for leading and subleading jets
        findLeadSubleadJets( pt, counter, ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead );

        counter++;
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    // Check if leading and subleading jets were found
    if ( idRecoLead >= 0 && idRecoSubLead >= 0 ) {
        if ( fVerbose ) {
            std::cout << "Reco dijet found: [true]" << std::endl;
        }

        // Retrieve leading jet
        RecoJet *leadJet = event->recoJetCollection()->at( idRecoLead );
        GenJet  *refLeadJet = { nullptr };
        
        ptRecoLead = leadJet->ptJECCorr();
        etaRecoLead = leadJet->eta();
        phiRecoLead = leadJet->phi();
        
        // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
        double leadJetPtEtaPhiFlavPtHatCent[6] = { ptRecoLead, etaRecoLead, phiRecoLead, (double)-6, ptHat, (double)centrality };
        fHM->hRecoInclusiveLeadJetPt->Fill( ptRecoLead, 1. );
        fHM->hRecoInclusiveLeadJetPtWeighted->Fill( ptRecoLead, weight );
        fHM->hRecoInclusiveLeadJetEtaPt->Fill( etaRecoLead, ptRecoLead, 1. );
        fHM->hRecoInclusiveLeadJetEtaPtWeighted->Fill( etaRecoLead, ptRecoLead, weight );
        fHM->hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCent->Fill( leadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hRecoInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( leadJetPtEtaPhiFlavPtHatCent, weight );

        // Plot leading jet JES and correlations between reco and ref jet
        if ( leadJet->hasMatching() ) {
            refLeadJet = event->genJetCollection()->at( leadJet->genJetId() );
            ptRefLead = refLeadJet->pt();
            etaRefLead = refLeadJet->eta();
            phiRefLead = refLeadJet->phi();
            int flavBLead{-6};
            switch ( refLeadJet->flavorForB() )
            {
            case -99:
                flavBLead = {-6};
                break;
            case -5:
                flavBLead = {-5};
                break;
            case -4:
                flavBLead = {-4};
                break;
            case -3:
                flavBLead = {-3};
                break;
            case -2:
                flavBLead = {-2};
                break;
            case -1:
                flavBLead = {-1};
                break;
            case 0:
                flavBLead = {0};
                break;
            case 1:
                flavBLead = {1};
                break;
            case 2:
                flavBLead = {2};
                break;
            case 3:
                flavBLead = {3};
                break;
            case 4:
                flavBLead = {4};
                break;
            case 5:
                flavBLead = {5};
                break;
            case 21:
                flavBLead = {6};
                break;
            default:
                flavBLead = {-6};
                break;
            }

            // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
            double refLeadJetPtEtaPhiFlavPtHatCent[6] = { ptRefLead, etaRefLead, phiRefLead, (double)flavBLead, ptHat, (double)centrality };

            fHM->hRefInclusiveLeadJetPt->Fill( ptRefLead, 1. );
            fHM->hRefInclusiveLeadJetPtWeighted->Fill( ptRefLead, weight );
            fHM->hRefInclusiveLeadJetEtaPt->Fill( etaRefLead, ptRefLead, 1. );
            fHM->hRefInclusiveLeadJetEtaPtWeighted->Fill( etaRefLead, ptRefLead, weight );
            fHM->hRefInclusiveLeadJetPtEtaPhiFlavPtHatCent->Fill( refLeadJetPtEtaPhiFlavPtHatCent, 1. );
            fHM->hRefInclusiveLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( refLeadJetPtEtaPhiFlavPtHatCent, weight );
            
            double JES = ptRecoLead / ptRefLead;
            double jetJESPtEtaPhiCent[5] = { JES, ptRefLead, etaRefLead, phiRefLead, (double)centrality };
            fHM->hRecoInclusiveLeadJetJESPtEtaPhiCent->Fill( jetJESPtEtaPhiCent, 1. );
            fHM->hRecoInclusiveLeadJetJESPtEtaPhiCentWeighted->Fill( jetJESPtEtaPhiCent, weight );

            // Correlations between reco and ref leading jets
            // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
            double reco2refPtEtaPhiPtHatCentrality[8] = { ptRecoLead, ptRefLead, etaRecoLead, etaRefLead, phiRecoLead, phiRefLead, ptHat, (double)centrality };
            fHM->hReco2RefInclusiveLeadJetPt->Fill( ptRefLead, ptRecoLead, 1. );
            fHM->hReco2RefInclusiveLeadJetPtWeighted->Fill( ptRefLead, ptRecoLead, weight );
            fHM->hReco2RefInclusiveLeadJetEta->Fill( etaRefLead, etaRecoLead, 1. );
            fHM->hReco2RefInclusiveLeadJetEtaWeighted->Fill( etaRefLead, etaRecoLead, weight );
            fHM->hReco2RefInclusiveLeadJetPhi->Fill( phiRefLead, phiRecoLead, 1. );
            fHM->hReco2RefInclusiveLeadJetPhiWeighted->Fill( phiRefLead, phiRecoLead, weight );
            fHM->hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentrality->Fill( reco2refPtEtaPhiPtHatCentrality, 1. );
            fHM->hReco2RefInclusiveLeadJetPtEtaPhiPtHatCentralityWeighted->Fill( reco2refPtEtaPhiPtHatCentrality, weight );
        } // if ( leadJet->hasMatching() )

        // Retrieve subleading jet
        RecoJet *subLeadJet = event->recoJetCollection()->at( idRecoSubLead );
        GenJet  *refSubLeadJet = { nullptr };

        ptRecoSubLead = subLeadJet->ptJECCorr();
        etaRecoSubLead = subLeadJet->eta();
        phiRecoSubLead = subLeadJet->phi();
        
        // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
        double subLeadJetPtEtaPhiFlavPtHatCent[6] = { ptRecoSubLead, etaRecoSubLead, phiRecoSubLead, (double)-6, ptHat, (double)centrality };
        fHM->hRecoInclusiveSubLeadJetPt->Fill( ptRecoSubLead, 1. );
        fHM->hRecoInclusiveSubLeadJetPtWeighted->Fill( ptRecoSubLead, weight );
        fHM->hRecoInclusiveSubLeadJetEtaPt->Fill( etaRecoSubLead, ptRecoSubLead, 1. );
        fHM->hRecoInclusiveSubLeadJetEtaPtWeighted->Fill( etaRecoSubLead, ptRecoSubLead, weight );
        fHM->hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Fill( subLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hRecoInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( subLeadJetPtEtaPhiFlavPtHatCent, weight );

        // Plot subleading jet JES and correlations between reco and ref jet
        if ( subLeadJet->hasMatching() ) {
            refSubLeadJet = event->genJetCollection()->at( subLeadJet->genJetId() );
            ptRefSubLead = refSubLeadJet->pt();
            etaRefSubLead = refSubLeadJet->eta();
            phiRefSubLead = refSubLeadJet->phi();
            int flavBSubLead{-6};
            switch ( refSubLeadJet->flavorForB() )
            {
            case -99:
                flavBSubLead = {-6};
                break;
            case -5:
                flavBSubLead = {-5};
                break;
            case -4:
                flavBSubLead = {-4};
                break;
            case -3:
                flavBSubLead = {-3};
                break;
            case -2:
                flavBSubLead = {-2};
                break;
            case -1:
                flavBSubLead = {-1};
                break;
            case 0:
                flavBSubLead = {0};
                break;
            case 1:
                flavBSubLead = {1};
                break;
            case 2:
                flavBSubLead = {2};
                break;
            case 3:
                flavBSubLead = {3};
                break;
            case 4:
                flavBSubLead = {4};
                break;
            case 5:
                flavBSubLead = {5};
                break;
            case 21:
                flavBSubLead = {6};
                break;
            default:
                flavBSubLead = {-6};
                break;
            }

            // 0 - pt, 1 - eta, 2 - phi, 3 - flavorForB, 4 - ptHat, 5 - centrality
            double refSubLeadJetPtEtaPhiFlavPtHatCent[6] = { ptRefSubLead, etaRefSubLead, phiRefSubLead, (double)flavBSubLead, ptHat, (double)centrality };

            fHM->hRefInclusiveSubLeadJetPt->Fill( ptRefSubLead, 1. );
            fHM->hRefInclusiveSubLeadJetPtWeighted->Fill( ptRefSubLead, weight );
            fHM->hRefInclusiveSubLeadJetEtaPt->Fill( etaRefSubLead, ptRefSubLead, 1. );
            fHM->hRefInclusiveSubLeadJetEtaPtWeighted->Fill( etaRefSubLead, ptRefSubLead, weight );
            fHM->hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCent->Fill( refSubLeadJetPtEtaPhiFlavPtHatCent, 1. );
            fHM->hRefInclusiveSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( refSubLeadJetPtEtaPhiFlavPtHatCent, weight );

            double JES = ptRecoSubLead / ptRefSubLead;
            double jetJESPtEtaPhiCent[5] = { JES, ptRefSubLead, etaRefSubLead, phiRefSubLead, (double)centrality };
            fHM->hRecoInclusiveSubLeadJetJESPtEtaPhiCent->Fill( jetJESPtEtaPhiCent, 1. );
            fHM->hRecoInclusiveSubLeadJetJESPtEtaPhiCentWeighted->Fill( jetJESPtEtaPhiCent, weight );
            
            // Correlations between reco and ref subleading jets
            // 0 - pt, 1 - refPt, 2 - eta, 3 - refEta, 4 - phi, 5 - refPhi, 6 - ptHat, 7 - centrality
            double reco2refPtEtaPhiPtHatCentrality[8] = { ptRecoSubLead, ptRefSubLead, etaRecoSubLead, etaRefSubLead, phiRecoSubLead, phiRefSubLead, ptHat, (double)centrality };
            fHM->hReco2RefInclusiveSubLeadJetPt->Fill( ptRefSubLead, ptRecoSubLead, 1. );
            fHM->hReco2RefInclusiveSubLeadJetPtWeighted->Fill( ptRefSubLead, ptRecoSubLead, weight );
            fHM->hReco2RefInclusiveSubLeadJetEta->Fill( etaRefSubLead, etaRecoSubLead, 1. );
            fHM->hReco2RefInclusiveSubLeadJetEtaWeighted->Fill( etaRefSubLead, etaRecoSubLead, weight );
            fHM->hReco2RefInclusiveSubLeadJetPhi->Fill( phiRefSubLead, phiRecoSubLead, 1. );
            fHM->hReco2RefInclusiveSubLeadJetPhiWeighted->Fill( phiRefSubLead, phiRecoSubLead, weight );
            fHM->hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentrality->Fill( reco2refPtEtaPhiPtHatCentrality, 1. );
            fHM->hReco2RefInclusiveSubLeadJetPtEtaPhiPtHatCentralityWeighted->Fill( reco2refPtEtaPhiPtHatCentrality, weight );
        } // if ( subLeadJet->hasMatching() )

        // Calculate dijet quantities
        double dijetRecoDphi = deltaPhi( phiRecoLead, phiRecoSubLead );
        double dijetRecoEta = 0.5 * ( etaRecoLead + etaRecoSubLead );
        double dijetRecoPt = 0.5 * ( ptRecoLead + ptRecoSubLead );
        double dijetRecoDetaCM = 0.5 * ( etaRecoLead - etaRecoSubLead );

        fHM->hRecoInclusiveDijetDphi->Fill( dijetRecoDphi, 1. );
        fHM->hRecoInclusiveDijetDphiWeighted->Fill( dijetRecoDphi, weight );
        fHM->hRecoInclusiveDijetEtaPt->Fill( dijetRecoEta, dijetRecoPt, 1. );
        fHM->hRecoInclusiveDijetEtaPtWeighted->Fill( dijetRecoEta, dijetRecoPt, weight );
        fHM->hRecoInclusiveDijetEtaPtDphi->Fill( dijetRecoEta, dijetRecoPt, dijetRecoDphi, 1. );
        fHM->hRecoInclusiveDijetEtaPtDphiWeighted->Fill( dijetRecoEta, dijetRecoPt, dijetRecoDphi, weight );
        fHM->hRecoInclusiveDijetDetaCM->Fill( dijetRecoDetaCM, 1. );
        fHM->hRecoInclusiveDijetDetaCMWeighted->Fill( dijetRecoDetaCM, weight );
        fHM->hRecoInclusiveDijetDetaCMPt->Fill( dijetRecoDetaCM, dijetRecoPt, 1. );
        fHM->hRecoInclusiveDijetDetaCMPtWeighted->Fill( dijetRecoDetaCM, dijetRecoPt, weight );
        fHM->hRecoInclusiveDijetEtaDetaCMPt->Fill( dijetRecoEta, dijetRecoDetaCM, dijetRecoPt, 1. );
        fHM->hRecoInclusiveDijetEtaDetaCMPtWeighted->Fill( dijetRecoEta, dijetRecoDetaCM, dijetRecoPt, weight );

        // Dijet JES
        if ( leadJet->hasMatching() && subLeadJet->hasMatching() ) {
            refLeadJet = event->genJetCollection()->at( leadJet->genJetId() );
            refSubLeadJet = event->genJetCollection()->at( subLeadJet->genJetId() );

            double ptRefLead = refLeadJet->pt();
            double etaRefLead = refLeadJet->eta();
            double phiRefLead = refLeadJet->phi();

            double ptRefSubLead = refSubLeadJet->pt();
            double etaRefSubLead = refSubLeadJet->eta();
            double phiRefSubLead = refSubLeadJet->phi();

            double dijetRefPt = 0.5 * ( ptRefLead + ptRefSubLead );
            double dijetRefEta = 0.5 * ( etaRefLead + etaRefSubLead );
            double dijetRefDphi = deltaPhi( phiRefLead, phiRefSubLead );

            double JES = dijetRecoPt / dijetRefPt;
            double dijetJESPtEtaDphiCent[4] = { JES, dijetRefPt, dijetRefEta, dijetRefDphi };
            fHM->hRecoInclusiveDijetJESPtEtaDphiCent->Fill( dijetJESPtEtaDphiCent, 1. );
            fHM->hRecoInclusiveDijetJESPtEtaDphiCentWeighted->Fill( dijetJESPtEtaDphiCent, weight );

            fHM->hRefInclusiveDijetPt->Fill( dijetRefPt, 1. );
            fHM->hRefInclusiveDijetPtWeighted->Fill( dijetRefPt, weight );
            fHM->hRefInclusiveDijetEta->Fill( dijetRefEta, 1. );
            fHM->hRefInclusiveDijetEtaWeighted->Fill( dijetRefEta, weight );
            fHM->hRefInclusiveDijetEtaPt->Fill( dijetRefEta, dijetRefPt, 1. );
            fHM->hRefInclusiveDijetEtaPtWeighted->Fill( dijetRefEta, dijetRefPt, weight );
            fHM->hRefInclusiveDijetDphi->Fill( dijetRefDphi, 1. );
            fHM->hRefInclusiveDijetDphiWeighted->Fill( dijetRefDphi, weight );
            fHM->hRefInclusiveDijetEtaPtDphi->Fill( dijetRefEta, dijetRefPt, dijetRefDphi, 1. );
            fHM->hRefInclusiveDijetEtaPtDphiWeighted->Fill( dijetRefEta, dijetRefPt, dijetRefDphi, weight );
        } // if ( leadJet->hasMatching() && subLeadJet->hasMatching() )


        // Check if dijet is good (on MC level, both leading and subleading 
        // jets must have matching ref partners)

        if ( !leadJet->hasMatching() || !subLeadJet->hasMatching() ) {
            if ( fVerbose ) {
                std::cout << "Reco dijet do not have matching ref jets" << std::endl;
            }
            return;
        }

        // Check dijet kinematics
        if ( !isGoodDijet( ptRecoLead, etaRecoLead, ptRecoSubLead, etaRecoSubLead, dijetRecoDphi ) ) {
            return;
        }

        refLeadJet = event->genJetCollection()->at( leadJet->genJetId() );
        refSubLeadJet = event->genJetCollection()->at( subLeadJet->genJetId() );

        double ptRefLead = refLeadJet->pt();
        double etaRefLead = refLeadJet->eta();
        double phiRefLead = refLeadJet->phi();
        double refLeadJetPtEtaPhiFlavPtHatCent[6] = { ptRefLead, etaRefLead, phiRefLead, (double)-6, ptHat, (double)centrality };

        double ptRefSubLead = refSubLeadJet->pt();
        double etaRefSubLead = refSubLeadJet->eta();
        double phiRefSubLead = refSubLeadJet->phi();
        double refSubLeadJetPtEtaPhiFlavPtHatCent[6] = { ptRefSubLead, etaRefSubLead, phiRefSubLead, (double)-6, ptHat, (double)centrality };

        double dijetRefPt = 0.5 * ( ptRefLead + ptRefSubLead );
        double dijetRefEta = 0.5 * ( etaRefLead + etaRefSubLead );
        double dijetRefDphi = deltaPhi( phiRefLead, phiRefSubLead );


        // Fill selected leading jet histograms
        fHM->hRecoSelectedLeadJetPt->Fill( ptRecoLead, 1. );
        fHM->hRecoSelectedLeadJetPtWeighted->Fill( ptRecoLead, weight );
        fHM->hRecoSelectedLeadJetEtaPt->Fill( etaRecoLead, ptRecoLead, 1. );
        fHM->hRecoSelectedLeadJetEtaPtWeighted->Fill( etaRecoLead, ptRecoLead, weight );
        fHM->hRecoSelectedLeadJetPtEtaPhiFlavPtHatCent->Fill( leadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hRecoSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( leadJetPtEtaPhiFlavPtHatCent, weight );

        fHM->hRefSelectedLeadJetPt->Fill( ptRefLead, 1. );
        fHM->hRefSelectedLeadJetPtWeighted->Fill( ptRefLead, weight );
        fHM->hRefSelectedLeadJetEtaPt->Fill( etaRefLead, ptRefLead, 1. );
        fHM->hRefSelectedLeadJetEtaPtWeighted->Fill( etaRefLead, ptRefLead, weight );
        fHM->hRefSelectedLeadJetPtEtaPhiFlavPtHatCent->Fill( refLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hRefSelectedLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( refLeadJetPtEtaPhiFlavPtHatCent, weight );

        // Fill selected subleading jet histograms
        fHM->hRecoSelectedSubLeadJetPt->Fill( ptRecoSubLead, 1. );
        fHM->hRecoSelectedSubLeadJetPtWeighted->Fill( ptRecoSubLead, weight );
        fHM->hRecoSelectedSubLeadJetEtaPt->Fill( etaRecoSubLead, ptRecoSubLead, 1. );
        fHM->hRecoSelectedSubLeadJetEtaPtWeighted->Fill( etaRecoSubLead, ptRecoSubLead, weight );
        fHM->hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Fill( subLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hRecoSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( subLeadJetPtEtaPhiFlavPtHatCent, weight );

        fHM->hRefSelectedSubLeadJetPt->Fill( ptRefSubLead, 1. );
        fHM->hRefSelectedSubLeadJetPtWeighted->Fill( ptRefSubLead, weight );
        fHM->hRefSelectedSubLeadJetEtaPt->Fill( etaRefSubLead, ptRefSubLead, 1. );
        fHM->hRefSelectedSubLeadJetEtaPtWeighted->Fill( etaRefSubLead, ptRefSubLead, weight );
        fHM->hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCent->Fill( refSubLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hRefSelectedSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( refSubLeadJetPtEtaPhiFlavPtHatCent, weight );

        // Fill selected dijet histograms
        fHM->hRecoSelectedDijetPt->Fill( dijetRecoPt, 1. );
        fHM->hRecoSelectedDijetPtWeighted->Fill( dijetRecoPt, weight );
        fHM->hRecoSelectedDijetEta->Fill( dijetRecoEta, 1. );
        fHM->hRecoSelectedDijetEtaWeighted->Fill( dijetRecoEta, weight );
        fHM->hRecoSelectedDijetEtaPt->Fill( dijetRecoEta, dijetRecoPt, 1. );
        fHM->hRecoSelectedDijetEtaPtWeighted->Fill( dijetRecoEta, dijetRecoPt, weight );
        fHM->hRecoSelectedDijetEtaPtDphi->Fill( dijetRecoEta, dijetRecoPt, dijetRecoDphi, 1. );
        fHM->hRecoSelectedDijetEtaPtDphiWeighted->Fill( dijetRecoEta, dijetRecoPt, dijetRecoDphi, weight );
        fHM->hRecoSelectedDijetDetaCM->Fill( dijetRecoDetaCM, 1. );
        fHM->hRecoSelectedDijetDetaCMWeighted->Fill( dijetRecoDetaCM, weight );
        fHM->hRecoSelectedDijetDetaCMPt->Fill( dijetRecoDetaCM, dijetRecoPt, 1. );
        fHM->hRecoSelectedDijetDetaCMPtWeighted->Fill( dijetRecoDetaCM, dijetRecoPt, weight );
        fHM->hRecoSelectedDijetEtaDetaCMPt->Fill( dijetRecoEta, dijetRecoDetaCM, dijetRecoPt, 1. );
        fHM->hRecoSelectedDijetEtaDetaCMPtWeighted->Fill( dijetRecoEta, dijetRecoDetaCM, dijetRecoPt, weight );

        fHM->hRefSelectedDijetPt->Fill( dijetRefPt, 1. );
        fHM->hRefSelectedDijetPtWeighted->Fill( dijetRefPt, weight );
        fHM->hRefSelectedDijetEta->Fill( dijetRefEta, 1. );
        fHM->hRefSelectedDijetEtaWeighted->Fill( dijetRefEta, weight );
        fHM->hRefSelectedDijetEtaPt->Fill( dijetRefEta, dijetRefPt, 1. );
        fHM->hRefSelectedDijetEtaPtWeighted->Fill( dijetRefEta, dijetRefPt, weight );
        fHM->hRefSelectedDijetEtaPtDphi->Fill( dijetRefEta, dijetRefPt, dijetRefDphi, 1. );
        fHM->hRefSelectedDijetEtaPtDphiWeighted->Fill( dijetRefEta, dijetRefPt, dijetRefDphi, weight );

        // Fill full correlation between pt and eta of reco and ref dijets 
        // and correspoinding leading and subleading jets
        
        // 0 - reco dijetpt, 1 - reco dijeteta, 
        // 2 - reco leading jet pt, 3 - reco leading jet eta,
        // 4 - reco subleading jet pt, 5 - reco subleading jet eta, 
        // 6 - ref dijet pt, 7 - ref dijet eta,
        // 8 - ref leading jet pt, 9 - ref leading jet eta,
        // 10 - ref subleading jet pt, 11 - ref subleading jet eta
        double reco2refDijetPtEtaPtEta[12] = { dijetRecoPt, dijetRecoEta, 
                                               ptRecoLead, etaRecoLead, 
                                               ptRecoSubLead, etaRecoSubLead, 
                                               dijetRefPt, dijetRefEta, 
                                               ptRefLead, etaRefLead, 
                                               ptRefSubLead, etaRefSubLead };

        fHM->hReco2RefSelectedDijetPtEtaFull->Fill( reco2refDijetPtEtaPtEta, 1. );
        fHM->hReco2RefSelectedDijetPtEtaFullWeighted->Fill( reco2refDijetPtEtaPtEta, weight );
        
    } // if ( idRecoLead >= 0 && idRecoSubLead >= 0 )
    else {
        if ( fVerbose ) {
            std::cout << "Reco dijet found: [false]" << std::endl;
        }
    } // else
    
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processGenJets - end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processEvent(const Event* event) {

    //
    // Entry point for the analysis of the event
    //

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processEvent - start" << std::endl;
    }

    if ( !fHM ) return;

    //
    // Event quantities
    //

    // ptHat
    double ptHat = event->ptHat();
    // Vertex z position
    double vz = event->vz();
    // ptHat weight 
    double ptHatW = event->ptHatWeight();
    // Centrality
    double centrality = event->centrality();
    // Centrality weight
    double centW = event->centralityWeight();
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
    // Process gen jets
    //
    processGenJets(event, weight);

    //
    // Process reco jets
    //
    processRecoJets(event, weight);

    // Fill histograms with vertex z position
    fHM->hVz->Fill( vz, 1. );
    fHM->hVzCentWeighted->Fill( vz, centW );
    fHM->hVzPtHatWeighted->Fill( vz,  ptHatW );
    fHM->hVzWeighted->Fill( vz, weight );

    // Fill histograms with CMS definition of centrality in terms of HiBin
    double hiBin = event->hiBin();
    fHM->hHiBin->Fill( hiBin, 1. );
    fHM->hHiBinPtHatWeighted->Fill( hiBin, ptHatW );
    fHM->hHiBinWeighted->Fill( hiBin, weight );

    // Fill histograms with ptHat and ptHatW
    fHM->hPtHat->Fill( ptHat, 1. );
    fHM->hPtHatPtHatWeighted->Fill( ptHat, ptHatW );
    fHM->hPtHatCentWeighted->Fill( ptHat, centW );
    fHM->hPtHatWeighted->Fill( ptHat, weight );

    fHM->hPtHatWeight->Fill( ptHatW, 1. );
    fHM->hPtHatWeightWeighted->Fill( ptHatW, weight );

    // Fill histograms with collision centrality
    fHM->hCentrality->Fill( centrality, 1. );
    fHM->hCentralityPtHatWeighted->Fill( centrality, ptHatW );
    fHM->hCentralityWeighted->Fill( centrality, weight );

    // Fill histograms with vz, ptHat and centrality
    fHM->hVzPtHatCent->Fill( vz, ptHat, centrality, 1. );
    fHM->hVzPtHatCentPtHatWeighted->Fill( vz, ptHat, centrality, ptHatW );
    fHM->hVzPtHatCentWeighted->Fill( vz, ptHat, centrality, weight );

    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processEvent - end" << std::endl;
        std::cout << "====================================================" << std::endl;
    }
}

//________________
void JetESRAnalysis::finish() {
    // Save data and close files
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
