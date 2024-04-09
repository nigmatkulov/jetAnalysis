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
Double_t DiJetAnalysis::eventWeight(const Bool_t& isMc, const Bool_t& isPPb, 
                                    const Double_t& ptHat, const Double_t& vz) {
    Double_t weight{1.};
    Double_t genWeight{1.};
    Double_t vzWeight{1.};

    if (isMc && isPPb) {

        // Magic numbers are (cross section x Nevents generated)
		if (ptHat > 15.0 && ptHat <= 30.)       { genWeight = 1.0404701e-06 /* * 961104 */; }
		else if (ptHat > 30. && ptHat <= 50.)   { genWeight = 7.7966624e-08 /* * 952110 */; }
		else if (ptHat > 50. && ptHat <= 80.)   { genWeight = 1.0016052e-08 /* * 952554 */; }
		else if (ptHat > 80. && ptHat <= 120.)  { genWeight = 1.3018269e-09 /* * 996844 */; }
		else if (ptHat > 120.&& ptHat <= 170.)  { genWeight = 2.2648493e-10 /* * 964681 */; }
		else if (ptHat > 170. && ptHat <= 220.) { genWeight = 4.0879112e-11 /* * 999260 */; }
		else if (ptHat > 220. && ptHat <= 280.) { genWeight = 1.1898939e-11 /* * 964336 */; }
		else if (ptHat > 280. && ptHat <= 370.) { genWeight = 3.3364433e-12 /* * 995036 */; }
		else if (ptHat > 370. && ptHat <= 460.) { genWeight = 7.6612402e-13 /* * 958160 */; }
		else if (ptHat > 460. && ptHat <= 540.) { genWeight = 2.1341026e-13 /* * 981427 */; }
		else if (ptHat > 540.)                  { genWeight = 7.9191586e-14 /* * 1000000 */; }
		//evtgenWeight = (float) evtgenWeight/nevents;
		
		// Vz weighting
        TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
        VzWeightFunction->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);
        vzWeight = VzWeightFunction->Eval( vz );
        vzWeight = 1. / vzWeight;
    } // if (isMc && isPPb)

    weight = genWeight * vzWeight;

    return weight;
}

//________________
void DiJetAnalysis::processGenJets(const Event* event, Double_t ptHatW) {

    Double_t ptLeadCut{30.};
    Double_t ptSubLeadCut{20.};
    Double_t phiDiJetCut{2. * TMath::Pi() / 3};

    Double_t ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
             phiLead{0.},  phiSubLead{0.};

    GenJetIterator genJetIter;
    // Loop over generated jets
    for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

        Double_t pt = (*genJetIter)->pt();
        Double_t eta = (*genJetIter)->eta();
        Double_t phi = (*genJetIter)->phi();

        // Apply lab frame boost to CM for the pPb 
        if ( fIsPPb ) {
            if ( fIsPbGoingDir ) {
                eta += fEtaShift;
            }
            else {
                eta -= fEtaShift;
                eta = -eta;
            }
        } // if ( fIsPPb ) 
        
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
        // Fill inclusive jet pt
        fHM->hGenInclusiveJetPt->Fill(pt);
        fHM->hGenInclusiveJetPtEta->Fill(eta, pt);
    } // for ( genJetIter = event->genJetCollection()->begin();

    // Check the dijet selection on the MC level
    if ( !isGoodDijet(ptLead, ptSubLead, TMath::Abs(phiLead - phiSubLead)) ) continue;

    fHM->hGenPtLeadPtSublead->Fill(ptLead, ptSubLead);
    fHM->hGenEtaLeadEtaSublead->Fill(etaLead, etaSubLead);

    Double_t dijetPt = 0.5 * (ptLead + ptSubLead);
    Double_t dijetEta = 0.5 * (etaLead + etaSubLead);
    Double_t dijetDphi = TMath::Abs(phiLead - phiSubLead); // Should abs be used here

    Double_t genDijetLeadSublead[9] {dijetPt, dijetEta, dijetDphi, 
                                     ptLead, etaLead, phiLead, 
                                     ptSubLead, etaSubLead, phiSubLead };
    fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(genDijetLeadSublead);
    fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(genDijetLeadSublead, ptHatW);
}

//________________
void DiJetAnalysis::processRecoJets(const Event* event, Double_t ptHatW) {

    Double_t ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
             phiLead{0.},  phiSubLead{0.};

    // Loop over reconstructed jets
    PartFlowJetIterator pfJetIter;
    for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ ) {

        Double_t pt = (*pfJetIter)->ptJECCorr();
        Double_t eta = (*pfJetIter)->eta();
        Double_t phi = (*pfJetIter)->phi();

        // Apply lab frame boost to CM for the pPb 
        if ( fIsPPb ) {
            if ( fIsPbGoingDir ) {
                eta += fEtaShift;
            }
            else {
                eta -= fEtaShift;
                eta = -eta;
            }
        } // if ( fIsPPb ) 
        
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
        
        // Fill inclusive jet pt
        fHM->hRecoInclusiveJetPt->Fill(pt);
        fHM->hGenInclusiveJetPtEta->Fill(eta, pt);
    } // for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ )
}

//________________
Bool_t DiJetAnalysis::isGoodDijet(const Double_t& ptLead, const Double_t& ptSublead, const Double_t& dphi) {
    Bool_t isGood = ( ptLead > fLeadJetPtLow &&
                      ptSublead > fSubleadJetPtLow &&
                      dphi > fDijetPhiCut );
    return isGood;
}

//________________
void DiJetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    //std::cout << "DiJetAnalysis::processEvent" << std::endl;

    if ( !fHM ) {
        std::cout << "[Warning] No histogram manager connected to the DiJetAnalysis\n";
    }

    //
    // Event quantities
    //

    // ptHat
    Double_t ptHat = event->ptHat();
    Double_t vz = event->vz();
    // ptHat weight (a.k.a. event weight that includes the ptHat and vz)
    Double_t ptHatW{1.}
    // Check correct MC sample
    if ( fIsMc && fIsPPb ) {
        // Skip events with ptHat that is outside the ranged embedded
        if ( ptHat <= fPtHatRange[0] || ptHat > fPtHatRange[1] ) continue;
        ptHatW = eventWeight(fIsMc, fIsPPb, ptHat, vz);
    }

    // Check collision system and account for the lab frame rapidity
    if ( fIsPPb ) {
        if ( fIsPbGoingDir) {
            fEtaShift = {0.4654094531};
        }
        else {
            fEtaShift = {-0.4654094531};
        }
    }
    else {
        fEtaShift = {0};
    }

    Double_t centW = event->centralityWeight();
    centW = {1.}; // Do not apply weight for pPb
    //std::cout << "centrality weight: " << centW << std::endl;

    fHM->hHiBin->Fill( event->hiBin() );

    fHM->hVz->Fill( event->vz(),  centW );
    fHM->hVzWeighted->Fill( event->vz(), ptHatW * centW );

    fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW * centW );
    fHM->hPtHat->Fill( ptHat, centW );
    fHM->hPtHatWeighted->Fill( ptHat, ptHatW * centW );
    fHM->hPtHatWeight->Fill( ptHatW, centW );

    Double_t vzPtHat[2] = { vz, ptHat };
    fHM->hVzPtHat->Fill( vzPtHat, centW );
    fHM->hVzPtHatWeighted->Fill( vzPtHat, ptHatW * centW );

    if ( fIsMc ) {
        // Process and analyze gen jets
        processGenJets(event, ptHatW);
    }

    // Process and analyze reco jets
    processRecoJets(event, ptHatW);

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
