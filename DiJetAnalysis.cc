/**
 * @file DiJetAnalysis.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Dijet analysis
 * @version 0.1
 * @date 2024-01-09
 * 
 * @copyright Copyright (c) 2024
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
#include "DiJetAnalysis.h"

ClassImp(DiJetAnalysis)

//________________
DiJetAnalysis::DiJetAnalysis() : BaseAnalysis(), 
    fDebug{kFALSE}, fUseCentralityWeight{}, fHM{nullptr} {
    /* Empty */
}

//________________
DiJetAnalysis::~DiJetAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void DiJetAnalysis::init() {
    // Initialize analysis
    //std::cout << "DiJetAnalysis::init" << std::endl;
}

//________________
void DiJetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    //std::cout << "DiJetAnalysis::processEvent" << std::endl;

    if ( !fHM ) {
        std::cout << "[Warning] No histogram manager connected to the DiJetAnalysis\n";
    }

    // IMPORTANT!!! Skip events with more than two jets 
    if ( event->genJetCollection()->size() != 2 ) return;

    //
    // Event quantities
    //

    // ptHat weights
    float ptHatW = event->ptHatWeight();

    // ptHat
    Double_t ptHat = event->ptHat();
    Double_t centW = event->centralityWeight();
    //std::cout << "centrality weight: " << centW << std::endl;

    fHM->hHiBin->Fill( event->hiBin() );

    fHM->hVz->Fill( event->vz(),  centW );
    fHM->hVzWeighted->Fill( event->vz(), ptHatW * centW );

    fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW * centW );
    fHM->hPtHat->Fill( ptHat, centW );
    fHM->hPtHatWeighted->Fill( ptHat, ptHatW * centW );
    fHM->hPtHatWeight->Fill( ptHatW, centW );

    // Collision centrality
    Double_t centrality = event->centrality();
    fHM->hCentrality->Fill( centrality, centW );
    fHM->hCentralityWeighted->Fill( centrality, ptHatW * centW );

    Double_t vzPtHatCent[3] = { event->vz(), ptHat, centrality };
    fHM->hVzPtHatCent->Fill( vzPtHatCent, centW );
    fHM->hVzPtHatCentWeighted->Fill( vzPtHatCent, ptHatW * centW );

    //
    // Generated jets
    //

    Double_t ptLeadCut{30.};
    Double_t ptSubLeadCut{20.};
    Double_t phiDiJetCut{2. * TMath::Pi() / 3};

    Double_t ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
             phiLead{0.},  phiSubLead{0.};

    GenJetIterator genJetIter;
    for ( genJetIter = event->genJetCollection()->begin();
          genJetIter != event->genJetCollection()->end();
          genJetIter++ ) {

        Double_t pt = (*genJetIter)->pt();
        Double_t eta = (*genJetIter)->eta();
        Double_t phi = (*genJetIter)->phi();
        
        if ( pt > ptLead ) {
            ptSubLead = ptLead;
            etaSubLead = etaLead;
            phiSubLead = phiLead;
            ptLead = pt;
            etaLead = eta;
            phiLead = phi;
        }
        else if ( pt > ptSubLead ) {
            ptSubLead = pt;
            etaSubLead = eta;
            phiSubLead = phi;
        }
    } // for ( genJetIter = event->genJetCollection()->begin();

    
    /*
    std::cout << "ptLead: " << ptLead << " ptSubLead: " << ptSubLead
              << " etaLead: " << etaLead << " etaSubLead: " << etaSubLead 
              << " phiLead - phiSubLead: " << phiLead - phiSubLead << std::endl;
    */

    // Apply hardcoded cut
    if ( ptLead < 30. || ptSubLead < 20. ||
         etaLead < -3.465 || 2.535 < etaLead ||
         etaSubLead < -3.465 || 2.535 < etaSubLead ||
         /* TMath::Abs(etaLead) > 2.5 || TMath::Abs(etaSubLead) > 2.5 */
         TMath::Abs(phiLead - phiSubLead) < phiDiJetCut ) return;

    // std::cout << "Passed dijet cut" << std::endl;

    Double_t ptDiJet = 0.5 * (ptLead + ptSubLead);
    Double_t etaDiJet = 0.5 * (etaLead + etaSubLead) + fEtaShift;

    /*
    std::cout << "ptLead: " << ptLead << " ptSubLead: " << ptSubLead
              << " etaLead: " << etaLead << " etaSubLead: " << etaSubLead 
              << " phiLead - phiSubLead: " << phiLead - phiSubLead 
              << " ptDiJet: " << ptDiJet << " etaDiJet: " << etaDiJet << std::endl;
    */

    Double_t jetLeadPtEtaPhiSubLeadPtEtaPhi[7] { ptLead, etaLead, phiLead, ptSubLead, etaSubLead, phiSubLead, centrality };
    fHM->hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCent->Fill( jetLeadPtEtaPhiSubLeadPtEtaPhi );
    fHM->hGenLeadJetPtEtaPhiSubLeadPtEtaPhiCentWeighted->Fill( jetLeadPtEtaPhiSubLeadPtEtaPhi, ptHatW * centW );

    Double_t diJetPtEtaDeltaPhiCent[4] { ptDiJet, etaDiJet, phiLead - phiSubLead, centrality };
    fHM->hGenDiJetPtEtaDeltaPhiCent->Fill( diJetPtEtaDeltaPhiCent );
    fHM->hGenDiJetPtEtaDeltaPhiCentWeighted->Fill( diJetPtEtaDeltaPhiCent, ptHatW * centW );


    // Int_t leadJetIndex{-1}, currentIndex{0};
    // Double_t leadJetPt{-1};
    // PartFlowJetIterator pfJetIter;

    // // Loop to find leading jet
    // for ( pfJetIter = event->pfJetCollection()->begin();
    //         pfJetIter != event->pfJetCollection()->end();
    //         pfJetIter++ ) {
    //     Double_t pt = (*pfJetIter)->ptJECCorr();
    //     if ( pt > leadJetPt ) {
    //         leadJetIndex = currentIndex;
    //         leadJetPt = pt;
    //     }
    //     currentIndex++;
    // }

    // if ( leadJetIndex >= 0 ) {
    //     RecoJet *jet = event->pfJetCollection()->at( leadJetIndex );
    //     //std::cout << "Read leading jet info" << std::endl;
    //     if ( TMath::Abs( jet->eta() ) > 1.6 ||
    //             !jet->hasMatching() ) {
    //         //std::cout << "Bad event. Skip" << std::endl;
    //         return;
    //     }
    // }

    // //
    // // Event quantities
    // //

    // // ptHat weights
    // float ptHatW = event->ptHatWeight();
    // // ptHat
    // Double_t ptHat = event->ptHat();
    // Double_t centW = event->centralityWeight();
    // //std::cout << "centrality weight: " << centW << std::endl;

    // fHM->hHiBin->Fill( event->hiBin() );

    // fHM->hVz->Fill( event->vz(),  centW );
    // fHM->hVzWeighted->Fill( event->vz(), ptHatW * centW );

    // fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW * centW );
    // fHM->hPtHat->Fill( ptHat, centW );
    // fHM->hPtHatWeighted->Fill( ptHat, ptHatW * centW );
    // fHM->hPtHatWeight->Fill( ptHatW, centW );

    // // Collision centrality
    // Double_t centrality = event->centrality();
    // fHM->hCentrality->Fill( centrality, centW );
    // fHM->hCentralityWeighted->Fill( centrality, ptHatW * centW );

    // Double_t vzPtHatCent[3] = { event->vz(), ptHat, centrality };
    // fHM->hVzPtHatCent->Fill( vzPtHatCent, centW );
    // fHM->hVzPtHatCentWeighted->Fill( vzPtHatCent, ptHatW * centW );

    // fHM->hNBadJets[0]->Fill( event->numberOfOverscaledPFJets() );
    // if ( ptHat > 20 ) fHM->hNBadJets[1]->Fill( event->numberOfOverscaledPFJets() );
    // if ( ptHat > 40 ) fHM->hNBadJets[2]->Fill( event->numberOfOverscaledPFJets() );
    // if ( ptHat > 60 ) fHM->hNBadJets[3]->Fill( event->numberOfOverscaledPFJets() );
    // if ( ptHat > 80 ) fHM->hNBadJets[4]->Fill( event->numberOfOverscaledPFJets() );
    

    // //std::cout << "HiBin: " << event->hiBin() << " centrality: " << centrality << std::endl;

    // //
    // // Gen jet quantities
    // //

    // // Counters for gen jets with pT cuts: >0, >20, >50, >80, >120 GeV
    // Int_t nGenJets[5] {0, 0, 0, 0, 0};
    // GenJetIterator genJetIter;
    // for ( genJetIter = event->genJetCollection()->begin();
    //         genJetIter != event->genJetCollection()->end();
    //         genJetIter++ ) {

    //     Double_t pt = (*genJetIter)->pt();
    //     nGenJets[0]++;
    //     if ( pt > 20 ) nGenJets[1]++;
    //     if ( pt > 50 ) nGenJets[2]++;
    //     if ( pt > 80 ) nGenJets[3]++;
    //     if ( pt > 120 ) nGenJets[4]++;
    //     Double_t eta = (*genJetIter)->eta();
    //     Double_t phi = (*genJetIter)->phi();
    //     Double_t flavB{-6};
    //     switch ( (*genJetIter)->flavorForB() )
    //     {
    //     case -99:
    //         flavB = -6;
    //         break;
    //     case -5:
    //         flavB = -5;
    //         break;
    //     case -4:
    //         flavB = -4;
    //         break;
    //     case -3:
    //         flavB = -3;
    //         break;
    //     case -2:
    //         flavB = -2;
    //         break;
    //     case -1:
    //         flavB = -1;
    //         break;
    //     case 0:
    //         flavB = 0;
    //         break;
    //     case 1:
    //         flavB = 1;
    //         break;
    //     case 2:
    //         flavB = 2;
    //         break;
    //     case 3:
    //         flavB = 3;
    //         break;
    //     case 4:
    //         flavB = 4;
    //         break;
    //     case 5:
    //         flavB = 5;
    //         break;
    //     case 21:
    //         flavB = 6;
    //         break;
    //     default:
    //         flavB = -6;
    //         break;
    //     }

    //     Double_t genJetPtEtaPhiCent[4]{ pt, eta, phi, centrality };
    //     Double_t genJetPtFlavPtHatCent[4]{ pt, flavB, ptHat, centrality };

    //     fHM->hGenJetPtEtaPhiCent->Fill( genJetPtEtaPhiCent, centW);
    //     fHM->hGenJetPtEtaPhiCentWeighted->Fill( genJetPtEtaPhiCent, ptHatW * centW);
    //     fHM->hGenJetPtFlavPtHatCent->Fill( genJetPtFlavPtHatCent, centW);
    //     fHM->hGenJetPtFlavPtHatCentWeighted->Fill( genJetPtFlavPtHatCent, ptHatW * centW );
    // } // for ( genJetIter = event->genJetCollection()->begin();

    // //
    // // Reco jet quantities
    // //

    // // Counters for gen jets with pT cuts: >0, >20, >50, >80, >120 GeV
    // Int_t nRecoJets[5] {0, 0, 0, 0, 0};
    // Int_t nRefJets[5] {0, 0, 0, 0, 0};





    // currentIndex = {0}; // Restart counter

    // // Loop over reconstructed particle flow jets
    // for ( pfJetIter = event->pfJetCollection()->begin();
    //         pfJetIter != event->pfJetCollection()->end();
    //         pfJetIter++ ) {

    //     // For speed up purpose here
    //     // if ( !(*pfJetIter)->hasMatching() ) continue;

    //     // Reco jets
    //     Double_t eta = (*pfJetIter)->eta();
    //     Double_t phi = (*pfJetIter)->phi();
    //     Double_t ptRaw = (*pfJetIter)->pt();
    //     Double_t WTAeta = (*pfJetIter)->WTAEta();
    //     Double_t WTAphi = (*pfJetIter)->WTAPhi();
    //     Double_t dphi = TVector2::Phi_mpi_pi(phi - WTAphi);
    //     Double_t deltaR = TMath::Sqrt( (eta - WTAeta) * (eta - WTAeta) +
    //                                     dphi * dphi );

    //     Double_t recoJetRawPtEtaPhiCent[4]{ ptRaw, eta, phi, centrality };
    //     fHM->hRecoJetRawPtEtaPhiCent->Fill( recoJetRawPtEtaPhiCent, centW );

    //     // Corrected momentum of the reconstructed jet
    //     Double_t pt = (*pfJetIter)->ptJECCorr();
    //     nRecoJets[0]++;
    //     if ( pt > 20 ) nRecoJets[1]++;
    //     if ( pt > 50 ) nRecoJets[2]++;
    //     if ( pt > 80 ) nRecoJets[3]++;
    //     if ( pt > 120 ) nRecoJets[4]++;
    //     Double_t recoJetPtEtaPhiCent[4]{ pt, eta, phi, centrality };
    //     Double_t deltaRPtCent[3] = { deltaR, pt, centrality };
    //     fHM->hRecoJetDeltaRPtCent->Fill( deltaRPtCent, centW );

    //     Double_t flavB{-6};
    //     if ( !(*pfJetIter)->hasMatching() ) {
    //         Double_t tmp[4] {pt, -6., ptHat, centrality };
    //         fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( tmp, centW );
    //         fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( tmp, ptHatW * centW );
    //         fHM->hRecoUnmatchedJetPtFlavPtHatCent->Fill( tmp, centW);
    //         fHM->hRecoUnmatchedJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW);
    //     }

    //     // For the JES and other histograms matching is important
    //     if ( !(*pfJetIter)->hasMatching() ) continue;

    //     // Retrieve matched gen jet
    //     GenJet *matchedJet = event->genJetCollection()->at( (*pfJetIter)->genJetId() );

    //     Double_t genPt = matchedJet->pt();
    //     nRefJets[0]++;
    //     if ( genPt > 20 ) nRefJets[1]++;
    //     if ( genPt > 50 ) nRefJets[2]++;
    //     if ( genPt > 80 ) nRefJets[3]++;
    //     if ( genPt > 120 ) nRefJets[4]++;
    //     Double_t genEta = matchedJet->eta();
    //     Double_t genPhi = matchedJet->phi();
    //     switch ( matchedJet->flavorForB() )
    //     {
    //     case -99:
    //         flavB = -6;
    //         break;
    //     case -5:
    //         flavB = -5;
    //         break;
    //     case -4:
    //         flavB = -4;
    //         break;
    //     case -3:
    //         flavB = -3;
    //         break;
    //     case -2:
    //         flavB = -2;
    //         break;
    //     case -1:
    //         flavB = -1;
    //         break;
    //     case 0:
    //         flavB = 0;
    //         break;
    //     case 1:
    //         flavB = 1;
    //         break;
    //     case 2:
    //         flavB = 2;
    //         break;
    //     case 3:
    //         flavB = 3;
    //         break;
    //     case 4:
    //         flavB = 4;
    //         break;
    //     case 5:
    //         flavB = 5;
    //         break;
    //     case 21:
    //         flavB = 6;
    //         break;
    //     default:
    //         flavB = -6;
    //         break;
    //     }

    //     Double_t ptRawPtCorrPtGenCent[4] = { ptRaw, pt, genPt, centrality };
    //     fHM->hRecoJetRawPtCorrPtGenPtCent->Fill( ptRawPtCorrPtGenCent, centW );
        
    //     fHM->hRecoJetPtEtaPhiCent->Fill( recoJetPtEtaPhiCent, centW );
    //     fHM->hRecoJetPtEtaPhiCentWeighted->Fill( recoJetPtEtaPhiCent, ptHatW * centW);
    //     Double_t recoJetPtFlavPtHatCent[4] { pt, flavB, ptHat, centrality };
    //     fHM->hRecoJetPtFlavPtHatCent->Fill( recoJetPtFlavPtHatCent, centW );
    //     fHM->hRecoJetPtFlavPtHatCentWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW);
    //     fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( recoJetPtFlavPtHatCent, centW );
    //     fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW );

    //     // Fill the phi for lead jet (regardless of matching)
    //     if ( currentIndex == leadJetIndex ) {
    //         Double_t tmp[4] {pt, flavB, ptHat, centrality };
    //         fHM->hRecoLeadJetPtFlavPtHatCent->Fill( tmp, centW );
    //         fHM->hRecoLeadJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW );
    //     }

    //     //
    //     // Ref jets
    //     //

    //     Double_t refJetPtEtaPhiCent[4] { genPt, genEta, genPhi, centrality };
    //     Double_t refJetPtFlavPtHatCent[4] { genPt, flavB, ptHat, centrality };
    //     fHM->hRefJetPtEtaPhiCent->Fill( refJetPtEtaPhiCent, centW );
    //     fHM->hRefJetPtEtaPhiCentWeighted->Fill( refJetPtEtaPhiCent, ptHatW * centW );
    //     fHM->hRefJetPtFlavPtHatCent->Fill( refJetPtFlavPtHatCent, centW );
    //     fHM->hRefJetPtFlavPtHatCentWeighted->Fill( refJetPtFlavPtHatCent, ptHatW * centW);

    //     //
    //     // Jet Energy Scale
    //     //

    //     Double_t JESraw = ptRaw / genPt;
    //     Double_t JES = pt / genPt;
    //     Double_t jesPtEtaPhiCent[5] { JES, genPt, genEta, genPhi, centrality };
    //     Double_t jesRawPtEtaPhiCent[5] { JESraw, genPt, genEta, genPhi, centrality };
    //     Double_t jesPtFlavPtHatCent[5] { JES, genPt, flavB, ptHat, centrality };
    //     fHM->hJESPtEtaPhiCent->Fill( jesPtEtaPhiCent );
    //     fHM->hJESPtEtaPhiCentWeighted->Fill( jesPtEtaPhiCent, ptHatW );
    //     fHM->hJESRawPtFlavPtHatCent->Fill( jesRawPtEtaPhiCent );
    //     fHM->hJESRawPtFlavPtHatCentWeighted->Fill( jesRawPtEtaPhiCent, ptHatW );
    //     fHM->hJESPtFlavPtHatCent->Fill( jesPtFlavPtHatCent );
    //     fHM->hJESPtFlavPtHatCentWeighted->Fill( jesPtFlavPtHatCent, ptHatW );

    //     currentIndex++;
    // } // for ( pfJetIter = event->pfJetCollection()->begin();

    // //std::cout << "-------------------" << std::endl;
    // for (Int_t i{0}; i<5; i++) {
    //     fHM->hNRecoJets[i]->Fill( nRecoJets[i] );
    //     fHM->hNGenJets[i]->Fill( nGenJets[i] );
    //     fHM->hNRefJets[i]->Fill( nRefJets[i] );
    //     // std::cout << Form("nRecoJets: %d nGenJets: %d nRefJets: %d\n", 
    //     //                   nRecoJets[i], nGenJets[i], nRefJets[i] );
    // }
    //         case -99:
    //             flavB = -6;
    //             break;
    //         case -5:
    //             flavB = -5;
    //             break;
    //         case -4:
    //             flavB = -4;
    //             break;
    //         case -3:
    //             flavB = -3;
    //             break;
    //         case -2:
    //             flavB = -2;
    //             break;
    //         case -1:
    //             flavB = -1;
    //             break;
    //         case 0:
    //             flavB = 0;
    //             break;
    //         case 1:
    //             flavB = 1;
    //             break;
    //         case 2:
    //             flavB = 2;
    //             break;
    //         case 3:
    //             flavB = 3;
    //             break;
    //         case 4:
    //             flavB = 4;
    //             break;
    //         case 5:
    //             flavB = 5;
    //             break;
    //         case 21:
    //             flavB = 6;
    //             break;
    //         default:
    //             flavB = -6;
    //             break;
    //         }

    //         Double_t ptRawPtCorrPtGenCent[4] = { ptRaw, pt, genPt, centrality };
    //         fHM->hRecoJetRawPtCorrPtGenPtCent->Fill( ptRawPtCorrPtGenCent, centW );
            
    //         fHM->hRecoJetPtEtaPhiCent->Fill( recoJetPtEtaPhiCent, centW );
    //         fHM->hRecoJetPtEtaPhiCentWeighted->Fill( recoJetPtEtaPhiCent, ptHatW * centW);
    //         Double_t recoJetPtFlavPtHatCent[4] { pt, flavB, ptHat, centrality };
    //         fHM->hRecoJetPtFlavPtHatCent->Fill( recoJetPtFlavPtHatCent, centW );
    //         fHM->hRecoJetPtFlavPtHatCentWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW);
    //         fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( recoJetPtFlavPtHatCent, centW );
    //         fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW );

    //         // Fill the phi for lead jet (regardless of matching)
    //         if ( currentIndex == leadJetIndex ) {
    //             Double_t tmp[4] {pt, flavB, ptHat, centrality };
    //             fHM->hRecoLeadJetPtFlavPtHatCent->Fill( tmp, centW );
    //             fHM->hRecoLeadJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW );
    //         }

    //         //
    //         // Ref jets
    //         //

    //         Double_t refJetPtEtaPhiCent[4] { genPt, genEta, genPhi, centrality };
    //         Double_t refJetPtFlavPtHatCent[4] { genPt, flavB, ptHat, centrality };
    //         fHM->hRefJetPtEtaPhiCent->Fill( refJetPtEtaPhiCent, centW );
    //         fHM->hRefJetPtEtaPhiCentWeighted->Fill( refJetPtEtaPhiCent, ptHatW * centW );
    //         fHM->hRefJetPtFlavPtHatCent->Fill( refJetPtFlavPtHatCent, centW );
    //         fHM->hRefJetPtFlavPtHatCentWeighted->Fill( refJetPtFlavPtHatCent, ptHatW * centW);

    //         //
    //         // Jet Energy Scale
    //         //

    //         Double_t JESraw = ptRaw / genPt;
    //         Double_t JES = pt / genPt;
    //         Double_t jesPtEtaPhiCent[5] { JES, genPt, genEta, genPhi, centrality };
    //         Double_t jesRawPtEtaPhiCent[5] { JESraw, genPt, genEta, genPhi, centrality };
    //         Double_t jesPtFlavPtHatCent[5] { JES, genPt, flavB, ptHat, centrality };
    //         fHM->hJESPtEtaPhiCent->Fill( jesPtEtaPhiCent );
    //         fHM->hJESPtEtaPhiCentWeighted->Fill( jesPtEtaPhiCent, ptHatW );
    //         fHM->hJESRawPtFlavPtHatCent->Fill( jesRawPtEtaPhiCent );
    //         fHM->hJESRawPtFlavPtHatCentWeighted->Fill( jesRawPtEtaPhiCent, ptHatW );
    //         fHM->hJESPtFlavPtHatCent->Fill( jesPtFlavPtHatCent );
    //         fHM->hJESPtFlavPtHatCentWeighted->Fill( jesPtFlavPtHatCent, ptHatW );

    //         currentIndex++;
    //     } // for ( pfJetIter = event->pfJetCollection()->begin();

    //     //std::cout << "-------------------" << std::endl;
    //     for (Int_t i{0}; i<5; i++) {
    //         fHM->hNRecoJets[i]->Fill( nRecoJets[i] );
    //         fHM->hNGenJets[i]->Fill( nGenJets[i] );
    //         fHM->hNRefJets[i]->Fill( nRefJets[i] );
    //         // std::cout << Form("nRecoJets: %d nGenJets: %d nRefJets: %d\n", 
    //         //                   nRecoJets[i], nGenJets[i], nRefJets[i] );
    //     }
    //         case -99:
    //             flavB = -6;
    //             break;
    //         case -5:
    //             flavB = -5;
    //             break;
    //         case -4:
    //             flavB = -4;
    //             break;
    //         case -3:
    //             flavB = -3;
    //             break;
    //         case -2:
    //             flavB = -2;
    //             break;
    //         case -1:
    //             flavB = -1;
    //             break;
    //         case 0:
    //             flavB = 0;
    //             break;
    //         case 1:
    //             flavB = 1;
    //             break;
    //         case 2:
    //             flavB = 2;
    //             break;
    //         case 3:
    //             flavB = 3;
    //             break;
    //         case 4:
    //             flavB = 4;
    //             break;
    //         case 5:
    //             flavB = 5;
    //             break;
    //         case 21:
    //             flavB = 6;
    //             break;
    //         default:
    //             flavB = -6;
    //             break;
    //         }

    //         Double_t ptRawPtCorrPtGenCent[4] = { ptRaw, pt, genPt, centrality };
    //         fHM->hRecoJetRawPtCorrPtGenPtCent->Fill( ptRawPtCorrPtGenCent, centW );
            
    //         fHM->hRecoJetPtEtaPhiCent->Fill( recoJetPtEtaPhiCent, centW );
    //         fHM->hRecoJetPtEtaPhiCentWeighted->Fill( recoJetPtEtaPhiCent, ptHatW * centW);
    //         Double_t recoJetPtFlavPtHatCent[4] { pt, flavB, ptHat, centrality };
    //         fHM->hRecoJetPtFlavPtHatCent->Fill( recoJetPtFlavPtHatCent, centW );
    //         fHM->hRecoJetPtFlavPtHatCentWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW);
    //         fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( recoJetPtFlavPtHatCent, centW );
    //         fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW );

    //         // Fill the phi for lead jet (regardless of matching)
    //         if ( currentIndex == leadJetIndex ) {
    //             Double_t tmp[4] {pt, flavB, ptHat, centrality };
    //             fHM->hRecoLeadJetPtFlavPtHatCent->Fill( tmp, centW );
    //             fHM->hRecoLeadJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW );
    //         }

    //         //
    //         // Ref jets
    //         //

    //         Double_t refJetPtEtaPhiCent[4] { genPt, genEta, genPhi, centrality };
    //         Double_t refJetPtFlavPtHatCent[4] { genPt, flavB, ptHat, centrality };
    //         fHM->hRefJetPtEtaPhiCent->Fill( refJetPtEtaPhiCent, centW );
    //         fHM->hRefJetPtEtaPhiCentWeighted->Fill( refJetPtEtaPhiCent, ptHatW * centW );
    //         fHM->hRefJetPtFlavPtHatCent->Fill( refJetPtFlavPtHatCent, centW );
    //         fHM->hRefJetPtFlavPtHatCentWeighted->Fill( refJetPtFlavPtHatCent, ptHatW * centW);

    //         //
    //         // Jet Energy Scale
    //         //

    //         Double_t JESraw = ptRaw / genPt;
    //         Double_t JES = pt / genPt;
    //         Double_t jesPtEtaPhiCent[5] { JES, genPt, genEta, genPhi, centrality };
    //         Double_t jesRawPtEtaPhiCent[5] { JESraw, genPt, genEta, genPhi, centrality };
    //         Double_t jesPtFlavPtHatCent[5] { JES, genPt, flavB, ptHat, centrality };
    //         fHM->hJESPtEtaPhiCent->Fill( jesPtEtaPhiCent );
    //         fHM->hJESPtEtaPhiCentWeighted->Fill( jesPtEtaPhiCent, ptHatW );
    //         fHM->hJESRawPtFlavPtHatCent->Fill( jesRawPtEtaPhiCent );
    //         fHM->hJESRawPtFlavPtHatCentWeighted->Fill( jesRawPtEtaPhiCent, ptHatW );
    //         fHM->hJESPtFlavPtHatCent->Fill( jesPtFlavPtHatCent );
    //         fHM->hJESPtFlavPtHatCentWeighted->Fill( jesPtFlavPtHatCent, ptHatW );

    //         currentIndex++;
    //     } // for ( pfJetIter = event->pfJetCollection()->begin();

    //     //std::cout << "-------------------" << std::endl;
    //     for (Int_t i{0}; i<5; i++) {
    //         fHM->hNRecoJets[i]->Fill( nRecoJets[i] );
    //         fHM->hNGenJets[i]->Fill( nGenJets[i] );
    //         fHM->hNRefJets[i]->Fill( nRefJets[i] );
    //         // std::cout << Form("nRecoJets: %d nGenJets: %d nRefJets: %d\n", 
    //         //                   nRecoJets[i], nGenJets[i], nRefJets[i] );
    //     }
    //         case -99:
    //             flavB = -6;
    //             break;
    //         case -5:
    //             flavB = -5;
    //             break;
    //         case -4:
    //             flavB = -4;
    //             break;
    //         case -3:
    //             flavB = -3;
    //             break;
    //         case -2:
    //             flavB = -2;
    //             break;
    //         case -1:
    //             flavB = -1;
    //             break;
    //         case 0:
    //             flavB = 0;
    //             break;
    //         case 1:
    //             flavB = 1;
    //             break;
    //         case 2:
    //             flavB = 2;
    //             break;
    //         case 3:
    //             flavB = 3;
    //             break;
    //         case 4:
    //             flavB = 4;
    //             break;
    //         case 5:
    //             flavB = 5;
    //             break;
    //         case 21:
    //             flavB = 6;
    //             break;
    //         default:
    //             flavB = -6;
    //             break;
    //         }

    //         Double_t ptRawPtCorrPtGenCent[4] = { ptRaw, pt, genPt, centrality };
    //         fHM->hRecoJetRawPtCorrPtGenPtCent->Fill( ptRawPtCorrPtGenCent, centW );
            
    //         fHM->hRecoJetPtEtaPhiCent->Fill( recoJetPtEtaPhiCent, centW );
    //         fHM->hRecoJetPtEtaPhiCentWeighted->Fill( recoJetPtEtaPhiCent, ptHatW * centW);
    //         Double_t recoJetPtFlavPtHatCent[4] { pt, flavB, ptHat, centrality };
    //         fHM->hRecoJetPtFlavPtHatCent->Fill( recoJetPtFlavPtHatCent, centW );
    //         fHM->hRecoJetPtFlavPtHatCentWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW);
    //         fHM->hRecoJetPtFlavPtHatCentInclusive->Fill( recoJetPtFlavPtHatCent, centW );
    //         fHM->hRecoJetPtFlavPtHatCentInclusiveWeighted->Fill( recoJetPtFlavPtHatCent, ptHatW * centW );

    //         // Fill the phi for lead jet (regardless of matching)
    //         if ( currentIndex == leadJetIndex ) {
    //             Double_t tmp[4] {pt, flavB, ptHat, centrality };
    //             fHM->hRecoLeadJetPtFlavPtHatCent->Fill( tmp, centW );
    //             fHM->hRecoLeadJetPtFlavPtHatCentWeighted->Fill( tmp, ptHatW * centW );
    //         }

    //         //
    //         // Ref jets
    //         //

    //         Double_t refJetPtEtaPhiCent[4] { genPt, genEta, genPhi, centrality };
    //         Double_t refJetPtFlavPtHatCent[4] { genPt, flavB, ptHat, centrality };
    //         fHM->hRefJetPtEtaPhiCent->Fill( refJetPtEtaPhiCent, centW );
    //         fHM->hRefJetPtEtaPhiCentWeighted->Fill( refJetPtEtaPhiCent, ptHatW * centW );
    //         fHM->hRefJetPtFlavPtHatCent->Fill( refJetPtFlavPtHatCent, centW );
    //         fHM->hRefJetPtFlavPtHatCentWeighted->Fill( refJetPtFlavPtHatCent, ptHatW * centW);

    //         //
    //         // Jet Energy Scale
    //         //

    //         Double_t JESraw = ptRaw / genPt;
    //         Double_t JES = pt / genPt;
    //         Double_t jesPtEtaPhiCent[5] { JES, genPt, genEta, genPhi, centrality };
    //         Double_t jesRawPtEtaPhiCent[5] { JESraw, genPt, genEta, genPhi, centrality };
    //         Double_t jesPtFlavPtHatCent[5] { JES, genPt, flavB, ptHat, centrality };
    //         fHM->hJESPtEtaPhiCent->Fill( jesPtEtaPhiCent );
    //         fHM->hJESPtEtaPhiCentWeighted->Fill( jesPtEtaPhiCent, ptHatW );
    //         fHM->hJESRawPtFlavPtHatCent->Fill( jesRawPtEtaPhiCent );
    //         fHM->hJESRawPtFlavPtHatCentWeighted->Fill( jesRawPtEtaPhiCent, ptHatW );
    //         fHM->hJESPtFlavPtHatCent->Fill( jesPtFlavPtHatCent );
    //         fHM->hJESPtFlavPtHatCentWeighted->Fill( jesPtFlavPtHatCent, ptHatW );

    //         currentIndex++;
    //     } // for ( pfJetIter = event->pfJetCollection()->begin();

    //     //std::cout << "-------------------" << std::endl;
    //     for (Int_t i{0}; i<5; i++) {
    //         fHM->hNRecoJets[i]->Fill( nRecoJets[i] );
    //         fHM->hNGenJets[i]->Fill( nGenJets[i] );
    //         fHM->hNRefJets[i]->Fill( nRefJets[i] );
    //         // std::cout << Form("nRecoJets: %d nGenJets: %d nRefJets: %d\n", 
    //         //                   nRecoJets[i], nGenJets[i], nRefJets[i] );
    //     }
}

//________________
void DiJetAnalysis::finish() {
    // Save data and close files
    std::cout << "DiJetAnalysis::finish" << std::endl;
}

//________________
void DiJetAnalysis::report() {
    // Force to report everyone
}

//________________
TList* DiJetAnalysis::getOutputList() {
    TList *outputList = new TList();

    // Add list of settings for cuts

    return outputList;
}
