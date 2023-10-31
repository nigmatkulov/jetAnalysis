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

ClassImp(JetESRAnalysis)

//________________
JetESRAnalysis::JetESRAnalysis() : BaseAnalysis(), fDebug(kFALSE), fHM(nullptr) {
    /* Empty */
}

//________________
JetESRAnalysis::~JetESRAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void JetESRAnalysis::init() {
    // Initialize analysis
    //std::cout << "JetESRAnalysis::init" << std::endl;
}

//________________
void JetESRAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    //std::cout << "JetESRAnalysis::processEvent" << std::endl;

    if (fHM) {

        //
        // Event quantities
        //

        // ptHat weights
        float ptHatW = event->ptHatWeight();
        // ptHat
        Double_t ptHat = event->ptHat();

        fHM->hVz->Fill( event->vz() );
        fHM->hVzWeighted->Fill( event->vz(), ptHatW );
        fHM->hHiBin->Fill( event->hiBin() );
        fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW );
        fHM->hPtHat->Fill( ptHat );
        fHM->hPtHatWeighted->Fill( event->ptHat(), ptHatW );
        fHM->hPtHatWeight->Fill( ptHatW );

        // Collision centrality
        Double_t centrality = event->centrality();
        fHM->hCentrality->Fill( centrality );
        fHM->hCentralityWeighted->Fill( centrality, ptHatW );

        

        //std::cout << "HiBin: " << event->hiBin() << " centrality: " << centrality << std::endl;

        //
        // Gen jet quantities
        //

        fHM->hNGenJets->Fill( event->numberOfGenJets() );
        GenJetIterator genJetIter;
        for ( genJetIter = event->genJetCollection()->begin();
              genJetIter != event->genJetCollection()->end();
              genJetIter++ ) {

            Double_t pt = (*genJetIter)->pt();
            Double_t eta = (*genJetIter)->eta();
            Double_t phi = (*genJetIter)->phi();
            Double_t flavB{-6};
            switch ( (*genJetIter)->flavorForB() )
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

            Double_t genJetPtEtaPhiCent[4]{ pt, eta, phi, centrality };
            Double_t genJetPtFlavPtHatCent[4]{ pt, flavB, ptHat, centrality };

            fHM->hGenJetPtEtaPhiCent->Fill( genJetPtEtaPhiCent );
            fHM->hGenJetPtEtaPhiCentWeighted->Fill( genJetPtEtaPhiCent, ptHatW );
            fHM->hGenJetPtFlavPtHatCent->Fill( genJetPtFlavPtHatCent );
            fHM->hGenJetPtFlavPtHatCentWeighted->Fill( genJetPtFlavPtHatCent, ptHatW );
        } // for ( genJetIter = event->genJetCollection()->begin();


        //
        // Reco jet quantities
        //

        fHM->hNRecoJets->Fill( event->numberOfPFJets() );
        PartFlowJetIterator pfJetIter;
        Int_t nRefJets{0};
        for ( pfJetIter = event->pfJetCollection()->begin();
              pfJetIter != event->pfJetCollection()->end();
              pfJetIter++ ) {

            // For speed up purpose here
            // if ( !(*pfJetIter)->hasMatching() ) continue;

            // Reco jets
            Double_t eta = (*pfJetIter)->eta();
            Double_t phi = (*pfJetIter)->phi();
            Double_t ptRaw = (*pfJetIter)->pt();
            Double_t WTAeta = (*pfJetIter)->WTAEta();
            Double_t WTAphi = (*pfJetIter)->WTAPhi();
            Double_t deltaR = TMath::Sqrt( (eta - WTAeta) * (eta - WTAeta) +
                                           (phi - WTAphi) * (phi - WTAphi) );

            Double_t recoJetRawPtEtaPhiCent[4]{ ptRaw, eta, phi, centrality };
            fHM->hRecoJetRawPtEtaPhiCent->Fill( recoJetRawPtEtaPhiCent );

            // Corrected momentum of the reconstructed jet
            Double_t pt = (*pfJetIter)->ptJECCorr();
            Double_t recoJetPtEtaPhiCent[4]{ pt, eta, phi, centrality };
            Double_t deltaRPtCent[3] = { deltaR, pt, centrality };
            fHM->hRecoJetDeltaRPtCent->Fill( deltaRPtCent );

            Double_t flavB{-6};
            if ( !(*pfJetIter)->hasMatching() ) {
                Double_t tmp[4] {pt, -6., ptHat, centrality };
                fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( tmp );
                fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( tmp, ptHatW );
            }

            // For the JES and other histograms matching is important
            if ( !(*pfJetIter)->hasMatching() ) continue;

            nRefJets++;

            // Retrieve matched gen jet
            GenJet *matchedJet = event->genJetCollection()->at( (*pfJetIter)->genJetId() );

            Double_t genPt = matchedJet->pt();
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
            
            fHM->hRecoJetPtEtaPhiCent->Fill( recoJetPtEtaPhiCent );
            fHM->hRecoJetPtEtaPhiCentWeighted->Fill( recoJetPtEtaPhiCent, ptHatW );
            Double_t recoJetPtFlavPtHatCent[4] { pt, flavB, ptHat, centrality };
            fHM->hRecoJetPtFlavPtHatCent->Fill( recoJetPtFlavPtHatCent );
            fHM->hRecoJetPtFlavPtHatCentWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW );
            fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( recoJetPtFlavPtHatCent );
            fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW );

            //
            // Ref jets
            //

            Double_t refJetPtEtaPhiCent[4] { genPt, genEta, genPhi, centrality };
            Double_t refJetPtFlavPtHatCent[4] { genPt, flavB, ptHat, centrality };
            fHM->hRefJetPtEtaPhiCent->Fill( refJetPtEtaPhiCent );
            fHM->hRefJetPtEtaPhiCentWeighted->Fill( refJetPtEtaPhiCent, ptHatW );
            fHM->hRefJetPtFlavPtHatCent->Fill( refJetPtFlavPtHatCent );
            fHM->hRefJetPtFlavPtHatCentWeighted->Fill( refJetPtFlavPtHatCent, ptHatW );

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

        } // for ( pfJetIter = event->pfJetCollection()->begin();

        fHM->hNRefJets->Fill( nRefJets );
    } // if (fHM)
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