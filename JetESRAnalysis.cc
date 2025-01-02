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
    fCollisionSystem{2}, fIsPbGoingDir{false}, fPtHatRange{15., 30.},
    fIsLooseJetIdCut{false}, fVzWeight{nullptr} {
    /* Empty */
}

//________________
JetESRAnalysis::~JetESRAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void JetESRAnalysis::init() {
    // Initialize analysis

    // Print initial setup of the analysis
    print();
}

//________________
void JetESRAnalysis::print() {
    std::cout << "----------------------------------------\n";
    std::cout << "JetESRAnalysis initial parameters:\n";
    std::cout << "Use centrality weight       : " << fUseCentralityWeight << std::endl
              << "Histogram manager           : " << fHM << std::endl
              << "Is pPb                      : " << fIsPPb << std::endl
              << "Is Pb-going direction       : " << fIsPbGoingDir << std::endl
              << "eta shift                   : " << fEtaShift << std::endl
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl
              << "Leading jet pT              : " << fLeadJetPtLow << std::endl
              << "SubLeading jet pT           : " << fSubleadJetPtLow << std::endl
              << "Dijet phi cut               : " << fDijetPhiCut << std::endl
              << "Select jets in CM frame     : " << fSelectJetsInCMFrame << std::endl;
    std::cout << "----------------------------------------\n";
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
    double etaCut[2] {fSubleadJetEta[0], fSubleadJetEta[1]}; 
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
double JetESRAnalysis::deltaPhi(const double& phi1, const double &phi2) {
    double dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ) dphi -= TMath::TwoPi();
    if ( dphi < -TMath::Pi() ) dphi += TMath::TwoPi();
    return dphi;
}

//________________
bool JetESRAnalysis::isGoodDijet(const double& ptLead, const double& ptSubLead, const double& dphi) {
    bool isGood = ( ptLead > fLeadJetPtLow &&
                    ptSubLead > fSubLeadJetPtLow && 
                    dphi > fDijetPhiCut );
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::isGoodDijet " << isGood << " ";
        std::cout << Form("pTlead: %5.2f pTsub: %5.2f dphi: %4.2f\n", ptLead, ptSubLead, dphi);
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
bool DiJetAnalysis::isGoodJetId(const RecoJet* jet) {

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
void JetESRAnalysis::processGenJets(const Event* event, const double &weight) {
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
        fHM->hGenInclusiveJetPtEtaPhiFlavPtHatCent->Fill( jetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenInclusiveJetPtEtaPhiFlavPtHatCentWeighted->Fill( jetPtEtaPhiFlavPtHatCent, weight );

        if ( isGoodGenJet( (*genJetIter) ) ) {
            // Fill good jet histograms
            fHM->hGenInclusiveGoodJetPtEtaPhiFlavPtHatCent->Fill( jetPtEtaPhiFlavPtHatCent, 1. );
            fHM->hGenInclusiveGoodJetPtEtaPhiFlavPtHatCentWeighted->Fill( jetPtEtaPhiFlavPtHatCent, weight );
        } // if ( isGoodGenJet( (*genJetIter) ) )

        // Search for leading and subleading jets
        findLeadSubleadJets( pt, counter, ptLead, ptSubLead, idLead, idSubLead );

    } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )

    // Check if leading and subleading jets were found
    if ( idLead >= 0 && idSubLead >= 0 ) {
        if ( fVerbose ) {
            std::cout << "Gen dijet found: [true]" << std::endl;
        }

        // Get leading and subleading jets
        GenJet *leadJet = event->genJetCollection()->at( idLead );
        GenJet *subLeadJet = event->genJetCollection()->at( idSubLead );

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

        if ( !isGoodDijet( ptLead, ptSubLead, deltaPhi( phiLead, phiSubLead ) ) ) {
            break;
        }

        // Fill leading and subleading jet histograms
        double leadJetPtEtaPhiFlavPtHatCent[6] = { ptLead, etaLead, phiLead, (double)flavBLead, ptHat, (double)centrality };
        double subLeadJetPtEtaPhiFlavPtHatCent[6] = { ptSubLead, etaSubLead, phiSubLead, (double)flavBSubLead, ptHat, (double)centrality };

        fHM->hGenLeadJetPtEtaPhiFlavPtHatCent->Fill( leadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( leadJetPtEtaPhiFlavPtHatCent, weight );
        fHM->hGenSubLeadJetPtEtaPhiFlavPtHatCent->Fill( subLeadJetPtEtaPhiFlavPtHatCent, 1. );
        fHM->hGenSubLeadJetPtEtaPhiFlavPtHatCentWeighted->Fill( subLeadJetPtEtaPhiFlavPtHatCent, weight );

        // Calculate dijet quantities
        double dijetDphi = deltaPhi( phiLead, phiSubLead );
        double dijetEta = 0.5 * ( etaLead + etaSubLead );
        double dijetPt = 0.5 * ( ptLead + ptSubLead );

        fHM->hGenDijetPtEtaDphi->Fill( dijetPt, dijetEta, dijetDphi, 1. );
        fHM->hGenDijetPtEtaDphiWeighted->Fill( dijetPt, dijetEta, dijetDphi, weight );

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
void JetESRAnalysis::processRecoJets(const Event* event, const double &weight) {
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
    Int_t counter{0};
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
        double jetPtEtaPhiDummyIterPfNhfPfNefNumOfConstPfMufPfChfChargedMultPfCefNeutralMultPtHatCent[13] = 
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
        
        fHM->hNHF[dummyIter]->Fill( (*recoJetIter)->jtPfNHF(), weight );
        fHM->hNEmF[dummyIter]->Fill( (*recoJetIter)->jtPfNEF(), weight );
        fHM->hNumOfConst[dummyIter]->Fill( numberOfConstituents, weight );
        fHM->hMUF[dummyIter]->Fill( (*recoJetIter)->jtPfMUF(), weight );
        fHM->hCHF[dummyIter]->Fill( (*recoJetIter)->jtPfCHF(), weight );
        fHM->hChargedMult[dummyIter]->Fill( chargedMult, weight );
        fHM->hCEmF[dummyIter]->Fill( (*recoJetIter)->jtPfCEF(), weight );
        fHM->hNumOfNeutPart[dummyIter]->Fill( neutralMult, weight );

        


    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )
    
    if ( fVerbose ) {
        std::cout << "JetESRAnalysis::processGenJets - end" << std::endl;
    }
}

//________________
void JetESRAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    //std::cout << "JetESRAnalysis::processEvent" << std::endl;

    if ( !fHM ) return;

    Int_t leadJetIndex{-1}, currentIndex{0};
    Double_t leadJetPt{-1};
    RecoJetIterator recoJetIter;

    //
    // Searching for leading jet that also should have gen jet
    //

    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {
        Double_t pt = (*recoJetIter)->ptJECCorr();
        if ( pt > leadJetPt ) {
            leadJetIndex = currentIndex;
            leadJetPt = pt;
        }
        currentIndex++;
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    if ( leadJetIndex >= 0 ) {
        RecoJet *jet = event->recoJetCollection()->at( leadJetIndex );
        //std::cout << "Read leading jet info" << std::endl;
        if ( TMath::Abs( jet->eta() ) > 1.6 ||
             !jet->hasMatching() ) {
            //std::cout << "Bad event. Skip" << std::endl;
            return;
        }
    } //if ( leadJetIndex >= 0 )

    //
    // Event quantities
    //

    // ptHat, ptHatW, centrality, centW
    
    Double_t ptHat = event->ptHat();
    float ptHatW = event->ptHatWeight();
    Double_t centrality = event->centrality();
    Double_t centW = event->centralityWeight();
    //std::cout << "centrality weight: " << centW << std::endl;

    // Fill histograms with CMS definition of centrality in terms of HiBin
    fHM->hHiBin->Fill( event->hiBin() );
    fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW * centW );

    // Fill histograms with collision centrality
    fHM->hCentrality->Fill( centrality, centW );
    fHM->hCentralityWeighted->Fill( centrality, ptHatW * centW );

    // Fill histograms with vertex z position
    fHM->hVz->Fill( event->vz(),  centW );
    fHM->hVzWeighted->Fill( event->vz(), ptHatW * centW );
    Double_t vzPtHatCent[3] = { event->vz(), ptHat, centrality };
    fHM->hVzPtHatCent->Fill( vzPtHatCent, centW );
    fHM->hVzPtHatCentWeighted->Fill( vzPtHatCent, ptHatW * centW );
    
    // Fill histograms with ptHat and ptHatW
    fHM->hPtHat->Fill( ptHat, centW );
    fHM->hPtHatWeighted->Fill( ptHat, ptHatW * centW );
    fHM->hPtHatWeight->Fill( ptHatW, centW );

    // Number of overscaled jets
    fHM->hNBadJets[0]->Fill( event->numberOfOverscaledRecoJets() );
    if ( ptHat > 20 ) fHM->hNBadJets[1]->Fill( event->numberOfOverscaledRecoJets() );
    if ( ptHat > 40 ) fHM->hNBadJets[2]->Fill( event->numberOfOverscaledRecoJets() );
    if ( ptHat > 60 ) fHM->hNBadJets[3]->Fill( event->numberOfOverscaledRecoJets() );
    if ( ptHat > 80 ) fHM->hNBadJets[4]->Fill( event->numberOfOverscaledRecoJets() );

    //std::cout << "HiBin: " << event->hiBin() << " centrality: " << centrality << std::endl;


    //
    // Generated jet quantities
    //
    processGenJets(event, weight);

    //
    // Reco jets
    //

    // Counters for gen jets with pT cuts: >0, >20, >50, >80, >120 GeV
    Int_t nRecoJets[5] {0, 0, 0, 0, 0};
    Int_t nRefJets[5] {0, 0, 0, 0, 0};
    // Restart counter
    currentIndex = {0}; 

    // Loop over reconstructed particle flow jets
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        // For speed up purpose here
        // if ( !(*recoJetIter)->hasMatching() ) continue;

        // Reco jets
        Double_t eta = (*recoJetIter)->eta();
        Double_t phi = (*recoJetIter)->phi();
        Double_t ptRaw = (*recoJetIter)->pt();
        Double_t WTAeta = (*recoJetIter)->WTAEta();
        Double_t WTAphi = (*recoJetIter)->WTAPhi();
        Double_t dphi = TVector2::Phi_mpi_pi(phi - WTAphi);
        Double_t deltaR = TMath::Sqrt( (eta - WTAeta) * (eta - WTAeta) +
                                        dphi * dphi );

        Double_t recoJetRawPtEtaPhiCent[4]{ ptRaw, eta, phi, centrality };
        fHM->hRecoJetRawPtEtaPhiCent->Fill( recoJetRawPtEtaPhiCent, centW );

        // Corrected momentum of the reconstructed jet
        Double_t pt = (*recoJetIter)->ptJECCorr();
        nRecoJets[0]++;
        if ( pt > 20 ) nRecoJets[1]++;
        if ( pt > 50 ) nRecoJets[2]++;
        if ( pt > 80 ) nRecoJets[3]++;
        if ( pt > 120 ) nRecoJets[4]++;
        Double_t recoJetPtEtaPhiCent[4]{ pt, eta, phi, centrality };
        Double_t deltaRPtCent[3] = { deltaR, pt, centrality };
        fHM->hRecoJetDeltaRPtCent->Fill( deltaRPtCent, centW );

        Double_t flavB{-6};
        if ( !(*recoJetIter)->hasMatching() ) {
            Double_t tmp[4] {pt, -6., ptHat, centrality };
            fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( tmp, centW );
            fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( tmp, ptHatW * centW );
            fHM->hRecoUnmatchedJetPtFlavPtHatCent->Fill( tmp, centW);
            fHM->hRecoUnmatchedJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW);
        }

        // For the JES and other histograms matching is important
        if ( !(*recoJetIter)->hasMatching() ) continue;

        // Retrieve matched gen jet
        GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );

        Double_t genPt = matchedJet->pt();
        nRefJets[0]++;
        if ( genPt > 20 ) nRefJets[1]++;
        if ( genPt > 50 ) nRefJets[2]++;
        if ( genPt > 80 ) nRefJets[3]++;
        if ( genPt > 120 ) nRefJets[4]++;
        Double_t genEta = matchedJet->eta();
        Double_t genPhi = matchedJet->phi();
        switch ( matchedJet->flavorForB() )
        {
        case -99:
            flavB = -6;
            break;
        case -5:
            flavB = -5;
            break;
        case -4:
            flavB = -4;
            break;
        case -3:
            flavB = -3;
            break;
        case -2:
            flavB = -2;
            break;
        case -1:
            flavB = -1;
            break;
        case 0:
            flavB = 0;
            break;
        case 1:
            flavB = 1;
            break;
        case 2:
            flavB = 2;
            break;
        case 3:
            flavB = 3;
            break;
        case 4:
            flavB = 4;
            break;
        case 5:
            flavB = 5;
            break;
        case 21:
            flavB = 6;
            break;
        default:
            flavB = -6;
            break;
        }

        Double_t ptRawPtCorrPtGenCent[4] = { ptRaw, pt, genPt, centrality };
        fHM->hRecoJetRawPtCorrPtGenPtCent->Fill( ptRawPtCorrPtGenCent, centW );
        
        fHM->hRecoJetPtEtaPhiCent->Fill( recoJetPtEtaPhiCent, centW );
        fHM->hRecoJetPtEtaPhiCentWeighted->Fill( recoJetPtEtaPhiCent, ptHatW * centW);
        Double_t recoJetPtFlavPtHatCent[4] { pt, flavB, ptHat, centrality };
        fHM->hRecoJetPtFlavPtHatCent->Fill( recoJetPtFlavPtHatCent, centW );
        fHM->hRecoJetPtFlavPtHatCentWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW);
        fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( recoJetPtFlavPtHatCent, centW );
        fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW );

        // Fill the phi for lead jet (regardless of matching)
        if ( currentIndex == leadJetIndex ) {
            Double_t tmp[4] {pt, flavB, ptHat, centrality };
            fHM->hRecoLeadJetPtFlavPtHatCent->Fill( tmp, centW );
            fHM->hRecoLeadJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW );
        }

        //
        // Ref jets
        //

        Double_t refJetPtEtaPhiCent[4] { genPt, genEta, genPhi, centrality };
        Double_t refJetPtFlavPtHatCent[4] { genPt, flavB, ptHat, centrality };
        fHM->hRefJetPtEtaPhiCent->Fill( refJetPtEtaPhiCent, centW );
        fHM->hRefJetPtEtaPhiCentWeighted->Fill( refJetPtEtaPhiCent, ptHatW * centW );
        fHM->hRefJetPtFlavPtHatCent->Fill( refJetPtFlavPtHatCent, centW );
        fHM->hRefJetPtFlavPtHatCentWeighted->Fill( refJetPtFlavPtHatCent, ptHatW * centW);

        //
        // Jet Energy Scale
        //

        Double_t JESraw = ptRaw / genPt;
        Double_t JES = pt / genPt;
        Double_t jesPtEtaPhiCent[5] { JES, genPt, genEta, genPhi, centrality };
        Double_t jesRawPtEtaPhiCent[5] { JESraw, genPt, genEta, genPhi, centrality };
        Double_t jesPtFlavPtHatCent[5] { JES, genPt, flavB, ptHat, centrality };
        fHM->hJESPtEtaPhiCent->Fill( jesPtEtaPhiCent );
        fHM->hJESPtEtaPhiCentWeighted->Fill( jesPtEtaPhiCent, ptHatW );
        fHM->hJESRawPtFlavPtHatCent->Fill( jesRawPtEtaPhiCent );
        fHM->hJESRawPtFlavPtHatCentWeighted->Fill( jesRawPtEtaPhiCent, ptHatW );
        fHM->hJESPtFlavPtHatCent->Fill( jesPtFlavPtHatCent );
        fHM->hJESPtFlavPtHatCentWeighted->Fill( jesPtFlavPtHatCent, ptHatW );

        currentIndex++;
    } // for ( recoJetIter = event->recoJetCollection()->begin();

    //std::cout << "-------------------" << std::endl;
    for (Int_t i{0}; i<5; i++) {
        fHM->hNRecoJets[i]->Fill( nRecoJets[i] );
        fHM->hNGenJets[i]->Fill( nGenJets[i] );
        fHM->hNRefJets[i]->Fill( nRefJets[i] );
        // std::cout << Form("nRecoJets: %d nGenJets: %d nRefJets: %d\n", 
        //                   nRecoJets[i], nGenJets[i], nRefJets[i] );
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
