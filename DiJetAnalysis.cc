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
#include "TMath.h"

// C++ headers
#include <iostream>

// Jet analysis headers
#include "DiJetAnalysis.h"

ClassImp(DiJetAnalysis)

//________________
DiJetAnalysis::DiJetAnalysis() : BaseAnalysis(), 
    fDebug{kFALSE}, fUseCentralityWeight{}, fHM{nullptr},
    fEtaShift{0}, fIsMc{kFALSE}, fIsPPb{kTRUE},
    fLeadJetPtLow{50.}, fSubleadJetPtLow{30.},
    fDijetPhiCut{TMath::TwoPi() / 3},
    fIsPbGoingDir{kFALSE}, fVerbose{kFALSE},
    fNEventsInSample{1000000}, fUseEtaShiftAndSignSwap{kFALSE},
    fIsDijetFound{kFALSE}, fIsDijetJetIdFound{kFALSE},
    fEventCounter{0}, fCycleCounter{0} {
    fPtHatRange[0] = {15.};
    fPtHatRange[1] = {30.};
}

//________________
DiJetAnalysis::~DiJetAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void DiJetAnalysis::init() {
    // Initialize analysis
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::init" << std::endl;
        print();
    } 
}

//________________
void DiJetAnalysis::print() {
    std::cout << "----------------------------------------\n";
    std::cout << "DiJetAnalysis parameters:\n";
    std::cout << "Use centrality weight       : " << fUseCentralityWeight << std::endl
              << "Histogram manager           : " << fHM << std::endl
              << "Is MC                       : " << fIsMc << std::endl
              << "Is pPb                      : " << fIsPPb << std::endl
              << "Is Pb-going direction       : " << fIsPbGoingDir << std::endl
              << "eta shift                   : " << fEtaShift << std::endl
              << "Use eta shift and sign swap : " << fUseEtaShiftAndSignSwap << std::endl
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl
              << "Leading jet pT              : " << fLeadJetPtLow << std::endl
              << "SubLeading jet pT           : " << fSubleadJetPtLow << std::endl
              << "Dijet phi cut               : " << fDijetPhiCut << std::endl;
    std::cout << "----------------------------------------\n";
}

//________________
Double_t DiJetAnalysis::eventWeight(const Bool_t& isMc, const Bool_t& isPPb, 
                                    const Double_t& ptHat, const Double_t& vz) {
    Double_t weight{1.};
    Double_t genWeight{1.};
    Double_t vzWeight{1.};

    if (isMc && isPPb) {

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
        TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
        VzWeightFunction->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);
        vzWeight = VzWeightFunction->Eval( vz );
        vzWeight = 1. / vzWeight;
    } // if (isMc && isPPb)

    weight = genWeight * vzWeight;

    if ( fVerbose) {
        std::cout << "fNEventsInSample: " << fNEventsInSample << " genWeight: " 
                  << genWeight << " vzWeight: " << vzWeight 
                  << " weight: " << weight << std::endl;
    }

    return weight;
}

//________________
Double_t DiJetAnalysis::deltaPhi(const Double_t& phi1, const Double_t phi2) {
    Double_t dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ) dphi -= TMath::TwoPi();
    if ( dphi < -TMath::Pi() ) dphi += TMath::TwoPi();
    return dphi;
}

//________________
Bool_t DiJetAnalysis::isGoodGenJet(const GenJet* jet) {
    Bool_t goodJet{kFALSE};
    if ( jet->pt() > 20. && TMath::Abs( jet->eta() ) < 3.0 ) {
        goodJet = {kTRUE};
    }
    
    if ( fVerbose ) {
        std::cout << Form("Gen jet cut %s\n", goodJet ? "\t[passed]" : "\t[failed]" );
    }

    return goodJet;
}

//________________
Bool_t DiJetAnalysis::isGoodTrkMax(const RecoJet* jet) {
    Bool_t goodTrackMax = {kTRUE};
    Double_t rawPt = jet->rawPt();
    Double_t trackMaxPt = jet->trackMaxPt();
    if ( TMath::Abs( jet->eta() ) < 2.4 && 
         ( trackMaxPt/rawPt < 0.01 ||
           trackMaxPt/rawPt > 0.98) ) {
        goodTrackMax = {kFALSE};
    }

    if ( fVerbose ) {
        std::cout << "TrackMaxPt/rawPt: " << trackMaxPt/rawPt << ( (goodTrackMax) ? "[passed]" : "[failed]" ) 
                  << ( (trackMaxPt/rawPt < 0.01) ? " too low value " : "" ) << ( (trackMaxPt/rawPt > 0.98) ? " too large value " : "" )
                  << std::endl;
    }

    return goodTrackMax;
}

//_________________
Bool_t DiJetAnalysis::isGoodJetId(const RecoJet* jet) {

	Bool_t passJetId = {kFALSE};

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

    Bool_t passNHF{kFALSE};
    Bool_t passNEF{kFALSE};
    Bool_t passNumOfConstituents{kTRUE};
    Bool_t passMuonFrac{kTRUE};
    Bool_t passChargedFrac{kTRUE};
    Bool_t passChargedMult{kTRUE};
    Bool_t passChargedEmFrac{kTRUE};
    Bool_t passNeutralMult{kTRUE};
	
    // Check cuts depending on jet pseudorapdity
    if ( TMath::Abs( eta ) <= 2.7 ) {
        
        passNHF = ( nhf < neutFracCut ) ? kTRUE : kFALSE;
        passNEF = ( nhf < neutFracCut ) ? kTRUE : kFALSE;
        passNumOfConstituents = ( numberOfConstituents > 1 ) ? kTRUE : kFALSE;

        if ( !fIsLooseJetIdCut ) { 
            passMuonFrac = ( muf < 0.8 ) ? kTRUE : kFALSE; 
        } // if ( !fIsLooseJetIdCut )

        if( TMath::Abs( eta ) <= 2.4 ) {
            passChargedFrac = ( chf > 0 ) ? kTRUE : kFALSE;
            passChargedMult = ( chargedMult > 0 ) ? kTRUE : kFALSE;
            passChargedEmFrac = ( cef < chargedEmFracCut ) ? kTRUE : kFALSE;
        } // if( TMath::Abs( eta ) <= 2.4 )

    } // if ( TMath::Abs( eta ) <= 2.7 )
    else if ( TMath::Abs( eta ) <= 3.0) {

        passNEF = ( nef > 0.01 ) ? kTRUE : kFALSE;
        passNHF = ( nhf < 0.98 ) ? kTRUE : kFALSE;
        passNeutralMult = ( neutralMult > 2 ) ? kTRUE : kFALSE;

    } // else if ( TMath::Abs( eta ) <= 3.0)
    else  {
        passNEF = ( nef < 0.9 ) ? kTRUE : kFALSE;
        passNeutralMult = (neutralMult > 2 ) ? kTRUE : kFALSE; // CAUTION: JET MET says it should be >10
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
Bool_t DiJetAnalysis::isGoodRecoJet(const RecoJet* jet) {
    Bool_t goodJet{kFALSE};
    Bool_t goodKine{kFALSE};
    Bool_t hasMatching{kFALSE};

    if ( jet->ptJECCorr() > 20 && TMath::Abs( jet->eta() ) < 3.0 ) {
        goodKine = {kTRUE};
    }

    if ( fIsMc ) {
        if ( jet->hasMatching() ) {
            hasMatching = {kTRUE};
        }
    }
    else {
        hasMatching = {kTRUE};
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
void DiJetAnalysis::processGenJets(const Event* event, Double_t ptHatW) {

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processGenJets" << std::endl;
    }

    Double_t ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
             phiLead{0.},  phiSubLead{0.};
    Int_t idLead{-1}, idSubLead{-1};
    Bool_t isDijetFound{ kFALSE };

    GenJetIterator genJetIter;
    Int_t counter{0};

    // Loop over generated jets
    for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

        Double_t pt = (*genJetIter)->pt();
        Double_t eta = (*genJetIter)->eta();
        Double_t phi = (*genJetIter)->phi();

        if ( fVerbose ) {
            std::cout << "Gen jet #" << counter << " ";
            (*genJetIter)->print();
        }

        // Apply single-jet selection to gen jets
        //if ( !isGoodGenJet( *genJetIter ) ) continue;

        // Apply lab frame boost to CM for the pPb 
        if ( fUseEtaShiftAndSignSwap && fIsPPb ) {
            if ( fIsPbGoingDir ) {
                eta -= fEtaShift;
                eta = -eta;
            }
            else {
                eta += fEtaShift;
            }
        } // if ( fUseEtaShiftAndSignSwap && fIsPPb )
        
        if ( pt > ptLead ) {
            ptSubLead = ptLead;
            etaSubLead = etaLead;
            phiSubLead = phiLead;
            idSubLead = idLead;
            ptLead = pt;
            etaLead = eta;
            phiLead = phi;
            idLead = counter;
        }
        else if ( pt > ptSubLead ) {
            ptSubLead = pt;
            etaSubLead = eta;
            phiSubLead = phi;
            idSubLead = counter;
        }
        // Fill inclusive jet pt
        fHM->hGenInclusiveJetPt->Fill(pt, ptHatW);
        fHM->hGenInclusiveJetPtEta->Fill(eta, pt, ptHatW);

        if ( fVerbose ) {
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idLead: %d idSubLead: %d\n", 
                              ptLead, ptSubLead, idLead, idSubLead);
        }
        counter++;
    } // for ( genJetIter = event->genJetCollection()->begin();

    //
    // Check for gen dijet
    //
    if ( idLead>=0 && idSubLead>=0 ) {

        Bool_t goodLeadJet = isGoodGenJet( event->genJetCollection()->at( idLead ) );
        Bool_t goodSubLeadJet = isGoodGenJet( event->genJetCollection()->at( idSubLead ) );
        Bool_t goodDijet = isGoodDijet( ptLead, ptSubLead, TMath::Abs( deltaPhi(phiLead, phiSubLead) ) );
        isDijetFound = goodLeadJet && goodSubLeadJet && goodDijet;

        // Analyze gen dijets
        if ( isDijetFound ) {

            fHM->hGenPtLeadPtSublead->Fill(ptLead, ptSubLead, ptHatW);
            fHM->hGenEtaLeadEtaSublead->Fill(etaLead, etaSubLead, ptHatW);

            Double_t dijetPt = 0.5 * (ptLead + ptSubLead);
            Double_t dijetEta = 0.5 * (etaLead + etaSubLead);
            Double_t dijetDphi = deltaPhi(phiLead, phiSubLead);

            Double_t genDijetLeadSublead[9] {dijetPt, dijetEta, dijetDphi, 
                                             ptLead, etaLead, phiLead, 
                                             ptSubLead, etaSubLead, phiSubLead };
            fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(genDijetLeadSublead);
            fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(genDijetLeadSublead, ptHatW);
            fHM->hGenDijetEta->Fill(dijetEta, ptHatW);
            fHM->hGenDijetPtEtaDphi->Fill(dijetPt, dijetEta, dijetDphi, 1.);
            fHM->hGenDijetPtEtaDphiWeighted->Fill(dijetPt, dijetEta, dijetDphi, ptHatW);
        } // if ( isDijetFound )
    } // if ( idLead>=0 && idSubLead>=0 )

    if ( fVerbose ) {
        std::cout << "Gen dijet found: " << ( (isDijetFound) ? "[true]" : "[false]" ) << std::endl;
        std::cout << "Reporting from DiJetAnalysis::processGenJets - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRecoJets(const Event* event, Double_t ptHatW) {

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processRecoJets" << std::endl;
    }

    // Define variables
    Double_t ptRecoLead{-1.}, ptRecoSubLead{-1.},
             ptRawRecoLead{-1.}, ptRawRecoSubLead{-1.},
             etaRecoLead{0.}, etaRecoSubLead{0.},
             phiRecoLead{0.},  phiRecoSubLead{0.}, 
             ptRefLead{-1.}, ptRefSubLead{-1.},
             etaRefLead{0.}, etaRefSubLead{0.},
             phiRefLead{0.}, phiRefSubLead{0.},
             ptRecoLeadJetId{-1.}, ptRecoSubLeadJetId{-1.},
             etaRecoLeadJetId{0.}, etaRecoSubLeadJetId{0.},
             phiRecoLeadJetId{0.},  phiRecoSubLeadJetId{0.},
             ptRefLeadJetId{-1.}, ptRefSubLeadJetId{-1.},
             etaRefLeadJetId{0.}, etaRefSubLeadJetId{0.},
             phiRefLeadJetId{0.}, phiRefSubLeadJetId{0.};
    Int_t  idRecoLead{-1}, idRecoSubLead{-1}, idRecoLeadJetId{-1}, idRecoSubLeadJetId{-1};

    // Loop over reconstructed jets
    PartFlowJetIterator pfJetIter;
    Int_t counter{0};
    for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ ) {

        Double_t pt = (*pfJetIter)->ptJECCorr();
        Double_t eta = (*pfJetIter)->eta();
        Double_t phi = (*pfJetIter)->phi();
        Double_t ptRaw = (*pfJetIter)->pt();

        if ( fVerbose ) {
            std::cout << "Reco jet #" << counter << " ";
            (*pfJetIter)->print();
        }

        // JetId parameters
        int chargedMult = (*pfJetIter)->jtPfCHM() + (*pfJetIter)->jtPfCEM() + (*pfJetIter)->jtPfMUM();
        int neutralMult = (*pfJetIter)->jtPfNHM() + (*pfJetIter)->jtPfNEM();
        int numberOfConstituents = chargedMult + neutralMult;

        Int_t dummyIter{0};
        if ( TMath::Abs( eta ) <= 2.4 ) { dummyIter = {0}; }
        else if ( TMath::Abs( eta ) <= 2.7 ) { dummyIter = {1}; }
        else if ( TMath::Abs( eta ) <= 3.0 ) { dummyIter = {2}; }
        else { dummyIter = {3}; }

        fHM->hRecoInclusiveJetPt->Fill(pt, ptHatW);
        fHM->hRecoInclusiveAllJetPtVsEta->Fill(eta, pt, ptHatW);
        fHM->hRecoInclusiveJetPtVsEtaKineCut->Fill(eta, pt, ptHatW);

        // JetId histograms
        fHM->hNHF[dummyIter]->Fill( (*pfJetIter)->jtPfNHF(), ptHatW );
        fHM->hNEmF[dummyIter]->Fill( (*pfJetIter)->jtPfNEF(), ptHatW );
        fHM->hNumOfConst[dummyIter]->Fill( numberOfConstituents, ptHatW );
        fHM->hMUF[dummyIter]->Fill( (*pfJetIter)->jtPfMUF(), ptHatW );
        fHM->hCHF[dummyIter]->Fill( (*pfJetIter)->jtPfCHF(), ptHatW );
        fHM->hChargedMult[dummyIter]->Fill( chargedMult, ptHatW );
        fHM->hCEmF[dummyIter]->Fill( (*pfJetIter)->jtPfCEF(), ptHatW );
        fHM->hNumOfNeutPart[dummyIter]->Fill( neutralMult, ptHatW );

        // Local variables for the current jet analysis
        Bool_t passTrkMax{kFALSE};
        Bool_t passJetId{kFALSE};
        Bool_t hasMatching{kFALSE};
        GenJet *matchedJet{nullptr};
        Double_t genPt{-999};
        Double_t genEta{-999};
        Double_t genPhi{-999};
        Double_t JES {-999};
        Double_t dEta{-999};
        Double_t dPhi{-999};
        Double_t res[4] { JES, genPt, genEta, genPhi };

        // Check selection criteria
        passTrkMax = isGoodTrkMax( (*pfJetIter) );
        passJetId = isGoodJetId( (*pfJetIter) );

        // On MC check reco jet matching to gen
        if ( fIsMc ) {
            hasMatching = (*pfJetIter)->hasMatching();
            if ( hasMatching ) {
                matchedJet = event->genJetCollection()->at( (*pfJetIter)->genJetId() );
                genPt = matchedJet->pt();
                genEta = matchedJet->eta();
                genPhi = matchedJet->phi();

                JES = pt/genPt;
                dEta = eta - genEta;
                dPhi = phi - genPhi;
                res[0] = JES;
                res[1] = genPt; 
                res[2] = genEta;
                res[3] = genPhi;

                fHM->hRefInclusiveJetPt->Fill(genPt, ptHatW);
                fHM->hRefInclusiveJetPtEta->Fill(genEta, genPt, ptHatW);
                fHM->hRecoMatchedPtEta->Fill(eta, pt, ptHatW);

                Double_t correl[5] { pt, ptRaw, genPt, eta, genEta };
                fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, ptHatW);
                fHM->hJESInclusiveJetPtEtaPhi->Fill(res);
                fHM->hJESInclusiveJetPtEtaPhiWeighted->Fill(res, ptHatW);

                fHM->hRecoInclusiveMatchedJetPtVsEta->Fill(eta, pt, ptHatW);
                fHM->hRecoInclusiveMatchedJetPtVsEtaKineCut->Fill(eta, pt, ptHatW);
                fHM->hRecoInclusiveJetRefPtVsEtaKineCut->Fill(genEta, genPt, ptHatW);
                fHM->hRecoInclusiveJetJESPtEtaPhiKineCut->Fill(res, ptHatW);
                fHM->hRecoInclusiveJetDEtaPtEtaKineCut->Fill(dEta, genPt, genEta, ptHatW);
            } // if ( hasMatching )
            else {
                // Fill unmatched jets
                fHM->hRecoInclusiveUnmatchedJetPtVsEta->Fill(eta, pt, ptHatW);
                fHM->hRecoInclusiveUnmatchedJetPtVsEtaKineCut->Fill(eta, pt, ptHatW);
            } // else
        } // if ( fIsMc )


        //
        // For track max selection
        //
        if ( passTrkMax ) {

            //
            // Inclusive jet part
            //

            fHM->hRecoInclusiveJetPtVsEtaTrkMaxCut->Fill(eta, pt, ptHatW);

            // For MC only
            if ( fIsMc ) {
                if ( !hasMatching ) {
                    fHM->hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut->Fill(eta, pt, ptHatW);
                }
                else {
                    fHM->hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut->Fill(eta, pt, ptHatW);
                    fHM->hRecoInclusiveJetRefPtVsEtaTrkMaxCut->Fill(genEta, genPt, ptHatW);
                    fHM->hRecoInclusiveJetJESPtEtaPhiTrkMaxCut->Fill(res, ptHatW);
                    fHM->hRecoInclusiveJetDEtaPtEtaTrkMaxCut->Fill(dEta, genPt, genEta, ptHatW);
                } // else
            } // if ( fIsMc )

            //
            // Dijet part
            //

            // Check for leading and subleading jets
            if ( pt > ptRecoLead ) {
                ptRecoSubLead = ptRecoLead;
                etaRecoSubLead = etaRecoLead;
                phiRecoSubLead = phiRecoLead;
                ptRawRecoSubLead = ptRawRecoLead;
                idRecoSubLead = idRecoLead;
                ptRecoLead = pt;
                ptRawRecoLead = ptRaw;
                etaRecoLead = eta;
                phiRecoLead = phi;
                idRecoLead = counter;

                if ( fIsMc ) {
                    if ( hasMatching ) {
                        ptRefSubLead = ptRefLead;
                        etaRefSubLead = etaRefLead;
                        phiRefSubLead = phiRefLead;
                        ptRefLead = genPt;
                        etaRefLead = genEta;
                        phiRefLead = genPhi;
                    }
                    else {
                        ptRefSubLead = ptRefLead;
                        etaRefSubLead = etaRefLead;
                        phiRefSubLead = phiRefLead;
                        ptRefLead = {-999.};
                        etaRefLead = {-999.};
                        phiRefLead = {-999.};
                    }
                } // if ( fIsMc )
            } // if ( pt > ptRecoLead )
            else if ( pt > ptRecoSubLead ) {
                ptRecoSubLead = pt;
                ptRawRecoSubLead = ptRaw;
                etaRecoSubLead = eta;
                phiRecoSubLead = phi;
                idRecoSubLead = counter;

                if ( fIsMc ) {
                    if ( hasMatching ) {
                        ptRefSubLead = genPt;
                        etaRefSubLead = genEta;
                        phiRefSubLead = genPhi;
                    }
                    else {
                        ptRefSubLead = {-999.};
                        etaRefSubLead = {-999.};
                        phiRefSubLead = {-999.};
                    }
                } // if ( fIsMc )
            } // else if ( pt > ptRecoSubLead )
        } // if ( passTrkMax )

        //
        // For jetId selection
        //
        if ( passJetId ) {

            //
            // Inclusive jet part
            //

            fHM->hRecoInclusiveJetPtVsEtaJetIdCut->Fill(eta, pt, ptHatW);

            // For MC only
            if ( fIsMc ) { 
                // If not matched
                if ( !hasMatching ) {
                    fHM->hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut->Fill(eta, pt, ptHatW);
                }
                else { // If matched
                    fHM->hRecoInclusiveMatchedJetPtVsEtaJetIdCut->Fill(eta, pt, ptHatW);
                    fHM->hRecoInclusiveJetRefPtVsEtaJetIdCut->Fill(genEta, genPt, ptHatW);
                    fHM->hRecoInclusiveJetJESPtEtaPhiJetIdCut->Fill(res, ptHatW);
                    fHM->hRecoInclusiveJetDEtaPtEtaJetIdCut->Fill(dEta, genPt, genEta, ptHatW);
                }
            } // if ( fIsMc )

            //
            // Dijet part
            //

            if ( pt > ptRecoLeadJetId ) {
                ptRecoSubLeadJetId = ptRecoLeadJetId;
                etaRecoSubLeadJetId = etaRecoLeadJetId;
                phiRecoSubLeadJetId = phiRecoLeadJetId;
                idRecoSubLeadJetId = idRecoLeadJetId;
                ptRecoLeadJetId = pt;
                etaRecoLeadJetId = eta;
                phiRecoLeadJetId = phi;
                idRecoLeadJetId = counter;

                if ( fIsMc ) { 
                    if (hasMatching ) {
                        ptRefSubLeadJetId = ptRefLeadJetId;
                        etaRefSubLeadJetId = etaRefLeadJetId;
                        phiRefSubLeadJetId = phiRefLeadJetId;
                        ptRefLeadJetId = genPt;
                        etaRefLeadJetId = genEta;
                        phiRefLeadJetId = genPhi;
                    }
                    else {
                        ptRefSubLeadJetId = ptRefLeadJetId;
                        etaRefSubLeadJetId = etaRefLeadJetId;
                        phiRefSubLeadJetId = phiRefLeadJetId;
                        ptRefLeadJetId = {-999.};
                        etaRefLeadJetId = {-999.};
                        phiRefLeadJetId = {-999.};                       
                    }
                } // if ( fIsMc )
            } // if ( pt > ptRecoLeadJetId )
            else if ( pt > ptRecoSubLeadJetId ) {
                ptRecoSubLeadJetId = pt;
                etaRecoSubLeadJetId = eta;
                phiRecoSubLeadJetId = phi;
                idRecoSubLeadJetId = counter;

                if ( fIsMc ) {
                    if ( hasMatching ) {
                        ptRefSubLeadJetId = genPt;
                        etaRefSubLeadJetId = genEta;
                        phiRefSubLeadJetId = genPhi;
                    }
                    else {
                        ptRefLeadJetId = {-999.};
                        etaRefLeadJetId = {-999.};
                        phiRefLeadJetId = {-999.}; 
                    }
                } // if ( fIsMc )
            } // else if ( pt > ptRecoSubLeadJetId )
        }

        if ( fVerbose ) {
            std::cout << Form("TrkMax selection --> Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d\n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead);
            std::cout << Form("JetId selection  --> Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d\n", 
                              ptRecoLeadJetId, ptRecoSubLeadJetId, idRecoLeadJetId, idRecoSubLeadJetId);
        }

        // Increment counter
        counter++;
    } // for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ )


    //
    // TrkMax dijets
    //

    fIsDijetFound = {kFALSE};

    if ( idRecoLead>=0 && idRecoSubLead>=0 ) {

        fHM->hRecoLeadJetAllPtVsEta->Fill(etaRecoLead, ptRecoLead, ptHatW);
        fHM->hRecoSubLeadJetAllPtVsEta->Fill(etaRecoSubLead, ptRecoSubLead, ptHatW);   

        if ( fIsMc ) {
            // Check leading jet matching to gen
            if ( event->pfJetCollection()->at( idRecoLead )->hasMatching() ) {
                fHM->hRecoLeadJetMatchedPtVsEta->Fill(etaRecoLead, ptRecoLead, ptHatW);
            }
            else {
                fHM->hRecoLeadJetUnmatchedPtVsEta->Fill(etaRecoLead, ptRecoLead, ptHatW);
            }

            // Check subleading jet matching to gen
            if ( event->pfJetCollection()->at( idRecoSubLead )->hasMatching() ) {
                fHM->hRecoSubLeadJetMatchedPtVsEta->Fill(etaRecoSubLead, ptRecoSubLead, ptHatW);
            }
            else {
                fHM->hRecoSubLeadJetUnmatchedPtVsEta->Fill(etaRecoSubLead, ptRecoSubLead, ptHatW);
            }
        } // if ( fIsMc )

        Bool_t goodLeadJet = isGoodRecoJet( event->pfJetCollection()->at( idRecoLead ) );
        Bool_t goodSubLeadJet = isGoodRecoJet( event->pfJetCollection()->at( idRecoSubLead ) );
        Bool_t goodDijet = isGoodDijet( ptRecoLead, ptRecoSubLead, TMath::Abs( deltaPhi(phiRecoLead, phiRecoSubLead) ) );
        fIsDijetFound = goodLeadJet && goodSubLeadJet && goodDijet;

        // Analyze trkMax dijets
        if ( fIsDijetFound ) {

            // Dijet analysis
            Double_t dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
            Double_t dijetRecoEta = 0.5 * (etaRecoLead + etaRecoSubLead);
            Double_t dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);

            // Correlation between leading and subleading
            fHM->hRecoPtLeadPtSublead->Fill(ptRecoLead, ptRecoSubLead, ptHatW);
            fHM->hRecoEtaLeadEtaSublead->Fill(etaRecoLead, etaRecoSubLead, ptHatW);

            Double_t dijetRecoInfo[9] { dijetRecoPt, dijetRecoEta, dijetRecoDphi,
                                        ptRecoLead, etaRecoLead, phiRecoLead,
                                        ptRecoSubLead, etaRecoSubLead, phiRecoSubLead };
            fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(dijetRecoInfo);
            fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(dijetRecoInfo, ptHatW);
            fHM->hRecoDijetEta->Fill( dijetRecoEta, ptHatW );
            fHM->hRecoDijetPtEtaDphi->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, 1. );
            fHM->hRecoDijetPtEtaDphiWeighted->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, ptHatW );

            if ( fIsMc ) {

                Double_t dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
                Double_t dijetRefEta = 0.5 * (etaRefLead + etaRefSubLead);
                Double_t dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);            

                // Leading jet information
                Double_t correl[5] { ptRecoLead, ptRawRecoLead, ptRefLead, etaRecoLead, etaRefLead };
                fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, ptHatW);

                // Subleading jet information
                correl[0] = ptRecoSubLead;
                correl[1] = ptRawRecoSubLead;
                correl[2] = ptRefSubLead;
                correl[3] = etaRecoSubLead; 
                correl[4] = etaRefSubLead;
                fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, ptHatW);

                Double_t dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                                 ptRecoLead, etaRecoLead,
                                                 ptRecoSubLead, etaRecoSubLead,
                                                 dijetRefPt, dijetRefEta,
                                                 ptRefLead, etaRefLead,
                                                 ptRefSubLead, etaRefSubLead };

                fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Fill(dijetRecoUnfold);
                fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, ptHatW);
                fHM->hRefDijetEta->Fill( dijetRefEta, ptHatW );
                fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEta, dijetRefEta, ptHatW );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, 1.);
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, ptHatW);
                fHM->hRefDijetPtEtaDphi->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, 1. );
                fHM->hRefDijetPtEtaDphiWeighted->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW );
                        
            } // if ( fIsMc )

        } // if ( fIsDijetFound )
    } // if ( idRecoLead>=0 && idRecoSubLead>=0 )

    //
    // JetId dijets
    //

    fIsDijetJetIdFound = {kFALSE};

    if ( idRecoLeadJetId>=0 && idRecoSubLeadJetId>=0 ) {

        fHM->hRecoLeadJetAllPtVsEtaJetIdCut->Fill(etaRecoLeadJetId, ptRecoLeadJetId, ptHatW);
        fHM->hRecoSubLeadJetAllPtVsEtaJetIdCut->Fill(etaRecoSubLeadJetId, ptRecoSubLeadJetId, ptHatW);   

        if ( fIsMc ) {
            // Check leading jet matching to gen
            if ( event->pfJetCollection()->at( idRecoLeadJetId )->hasMatching() ) {
                fHM->hRecoLeadJetMatchedPtVsEtaJetIdCut->Fill(etaRecoLeadJetId, ptRecoLeadJetId, ptHatW);
            }
            else {
                fHM->hRecoLeadJetUnmatchedPtVsEtaJetIdCut->Fill(etaRecoLeadJetId, ptRecoLeadJetId, ptHatW);
            }

            // Check subleading jet matching to gen
            if ( event->pfJetCollection()->at( idRecoSubLeadJetId )->hasMatching() ) {
                fHM->hRecoSubLeadJetMatchedPtVsEtaJetIdCut->Fill(etaRecoSubLeadJetId, ptRecoSubLeadJetId, ptHatW);
            }
            else {
                fHM->hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut->Fill(etaRecoSubLeadJetId, ptRecoSubLeadJetId, ptHatW);
            }
        } // if ( fIsMc )
        
        Bool_t goodLeadJet = isGoodRecoJet( event->pfJetCollection()->at( idRecoLeadJetId ) );
        Bool_t goodSubLeadJet = isGoodRecoJet( event->pfJetCollection()->at( idRecoSubLeadJetId ) );
        Bool_t goodDijet = isGoodDijet( ptRecoLeadJetId, ptRecoSubLeadJetId, TMath::Abs( deltaPhi(phiRecoLeadJetId, phiRecoSubLeadJetId) ) );
        fIsDijetJetIdFound = goodLeadJet && goodSubLeadJet && goodDijet;

        // Analyze trkMax dijets
        if ( fIsDijetJetIdFound ) {

            // Dijet analysis
            Double_t dijetRecoPt = 0.5 * (ptRecoLeadJetId + ptRecoSubLeadJetId);
            Double_t dijetRecoEta = 0.5 * (etaRecoLeadJetId + etaRecoSubLeadJetId);
            Double_t dijetRecoDphi = deltaPhi(phiRecoLeadJetId, phiRecoSubLeadJetId);

            fHM->hRecoDijetPtEtaDphiJetId->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, ptHatW );

            if ( fIsMc ) {

                Double_t dijetRefPt = 0.5 * (ptRefLeadJetId + ptRefSubLeadJetId);
                Double_t dijetRefEta = 0.5 * (etaRefLeadJetId + etaRefSubLeadJetId);
                Double_t dijetRefDphi = deltaPhi(phiRefLeadJetId, phiRefSubLeadJetId);

                fHM->hRefDijetPtEtaDphiJetId->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, ptHatW);
            } // if ( fIsMc )

        } // if ( fIsDijetJetIdFound )
    } // if ( idRecoLeadJetId>=0 && idRecoSubLeadJetId>=0 )

    // Fill matching between the histograms
    if ( fIsDijetFound && fIsDijetJetIdFound ) {
        if ( (idRecoLead == idRecoLeadJetId) && (idRecoSubLead == idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(1., ptHatW);
        }
        else if ( (idRecoLead == idRecoLeadJetId) && (idRecoSubLead != idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(2., ptHatW);
        }
        else if ( (idRecoLead != idRecoLeadJetId) && (idRecoSubLead == idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(3., ptHatW);
        }
        else if ( (idRecoLead != idRecoLeadJetId) && (idRecoSubLead != idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(4., ptHatW);
        }
    }
    else if ( fIsDijetFound && !fIsDijetJetIdFound ) {
        fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(5., ptHatW);
    }
    else if ( !fIsDijetFound && fIsDijetJetIdFound ) {
        fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(6., ptHatW);   
    }
    else if ( !fIsDijetFound && !fIsDijetJetIdFound ) {
        fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(7., ptHatW);
    }

    if ( fVerbose ) {
        std::cout << "TrkMax dijet found: " << ( (fIsDijetFound) ? "[true]" : "[false]" ) << std::endl;
        std::cout << "JetId  dijet found: " << ( (fIsDijetJetIdFound) ? "[true]" : "[false]" ) << std::endl;
        std::cout << "Reporting from DiJetAnalysis::processRecoJets - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRefJets(const Event* event, Double_t ptHatW) {
    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processRefJets" << std::endl;
    }

    Double_t ptRecoLead{-1.}, ptRecoSubLead{-1.},
             ptRawRecoLead{-1.}, ptRawRecoSubLead{-1.},
             etaRecoLead{0.}, etaRecoSubLead{0.},
             phiRecoLead{0.},  phiRecoSubLead{0.}, 
             ptRefLead{-1.}, ptRefSubLead{-1.},
             etaRefLead{0.}, etaRefSubLead{0.},
             phiRefLead{0.}, phiRefSubLead{0.};
    Bool_t isDijetFound{kFALSE};
    Int_t  idRecoLead{-1}, idRecoSubLead{-1};

    // Loop over reconstructed jets
    PartFlowJetIterator pfJetIter;
    Int_t counter{0};
    for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ ) {

        if ( !(*pfJetIter)->hasMatching() ) continue;

        GenJet *matchedJet{nullptr};
        Double_t genPt{999.};
        Double_t genEta{-999.};
        Double_t genPhi{-999.};

        Double_t pt = (*pfJetIter)->ptJECCorr();
        Double_t eta = (*pfJetIter)->eta();
        Double_t phi = (*pfJetIter)->phi();
        Double_t ptRaw = (*pfJetIter)->pt();

        matchedJet = event->genJetCollection()->at( (*pfJetIter)->genJetId() );
        genPt = matchedJet->pt();
        genEta = matchedJet->eta();
        genPhi = matchedJet->phi();

        if ( fVerbose ) {
            std::cout << "Ref jet info for reco jet #" << counter;
            matchedJet->print();
            std::cout << "Reco jet #" << counter << " ";
            (*pfJetIter)->print();
        }

        // Apply lab frame boost to CM for the pPb 
        if ( fUseEtaShiftAndSignSwap && fIsPPb ) {
            if ( fIsMc ) { // For embedding: Pb goes to negative, p goes to positive
                if ( fIsPbGoingDir ) {
                    eta -= fEtaShift;
                    eta = -eta;
                    genEta -= fEtaShift;
                    genEta = -genEta;
                }
                else {
                    eta += fEtaShift;
                    genEta += fEtaShift;
                }
            }
            else { // For data: p goes to negative, Pb goes to positive
                if ( fIsPbGoingDir ) {
                    eta += fEtaShift;
                }
                else {
                    eta -= fEtaShift;
                    eta = -eta;
                }
            }

        } // if ( fIsPPb )
        
        if ( genPt > ptRefLead ) {
            ptRecoSubLead = ptRecoLead;
            etaRecoSubLead = etaRecoLead;
            phiRecoSubLead = phiRecoLead;
            ptRawRecoSubLead = ptRawRecoLead;
            idRecoSubLead = idRecoLead;
            ptRecoLead = pt;
            ptRawRecoLead = ptRaw;
            etaRecoLead = eta;
            phiRecoLead = phi;
            idRecoLead = counter;

            ptRefSubLead = ptRefLead;
            etaRefSubLead = etaRefLead;
            phiRefSubLead = phiRefLead;
            ptRefLead = genPt;
            etaRefLead = genEta;
            phiRefLead = genPhi;
        }
        else if ( genPt > ptRefSubLead ) {
            ptRecoSubLead = pt;
            ptRawRecoSubLead = ptRaw;
            etaRecoSubLead = eta;
            phiRecoSubLead = phi;
            idRecoSubLead = counter;

            ptRefSubLead = genPt;
            etaRefSubLead = genEta;
            phiRefSubLead = genPhi;
        }

        if ( fVerbose ) {
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d genId: %d \n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead, (*pfJetIter)->genJetId());
        }

        // Increment counter
        counter++;
    } // for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ )

    //
    // Check if leading and subleading jets were found
    //
    if (idRecoLead>=0 && idRecoSubLead>=0) {
        if ( fVerbose ) {
            std::cout << Form("Checking dijet Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d Lead genId: %d SubLead genId: %d\n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead, 
                              event->pfJetCollection()->at( idRecoLead )->genJetId(), 
                              event->pfJetCollection()->at( idRecoSubLead )->genJetId() );
        }
        Bool_t goodLeadJet{kFALSE};
        Bool_t goodSubLeadJet{kFALSE};
        Bool_t goodDijet{kFALSE};
        if ( event->pfJetCollection()->at( idRecoLead )->hasMatching() && 
             event->pfJetCollection()->at( idRecoSubLead )->hasMatching() ) {
            goodLeadJet = isGoodGenJet( event->genJetCollection()->at( event->pfJetCollection()->at( idRecoLead )->genJetId() ) );
            goodSubLeadJet = isGoodGenJet( event->genJetCollection()->at( event->pfJetCollection()->at( idRecoSubLead )->genJetId() ) );
            goodDijet = isGoodDijet( ptRefLead, ptRefSubLead, TMath::Abs( deltaPhi(phiRefLead, phiRefSubLead) ) );
        }
        isDijetFound = goodLeadJet && goodSubLeadJet && goodDijet;

        // Analyze trkMax dijets
        if ( isDijetFound ) {

            // Dijet analysis
            Double_t dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
            Double_t dijetRecoEta = 0.5 * (etaRecoLead + etaRecoSubLead);
            Double_t dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);

            Double_t dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
            Double_t dijetRefEta = 0.5 * (etaRefLead + etaRefSubLead);
            Double_t dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);

                // Dijet reco vs ref for unfolding
            Double_t dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                             ptRecoLead, etaRecoLead,
                                             ptRecoSubLead, etaRecoSubLead,
                                             dijetRefPt, dijetRefEta,
                                             ptRefLead, etaRefLead,
                                             ptRefSubLead, etaRefSubLead };    

            fHM->hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, ptHatW);
            fHM->hRefSelDijetPtEtaDphi->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, 1.);
            fHM->hRefSelDijetPtEtaDphiWeighted->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW);
            fHM->hRefSelDijetEta->Fill(dijetRefEta, ptHatW);
        } // if ( isDijetFound )
    } // if (idRecoLead>=0 && idRecoSubLead>=0)

    if ( fVerbose ) {
        std::cout << "Dijet found: " << isDijetFound << std::endl;
        std::cout << "Reporting from DiJetAnalysis::processRefJets - [DONE]" << std::endl;
    }
}

//________________
Bool_t DiJetAnalysis::isGoodDijet(const Double_t& ptLead, const Double_t& ptSublead, const Double_t& dphi) {
    Bool_t isGood = ( ptLead > fLeadJetPtLow &&
                      ptSublead > fSubleadJetPtLow && 
                      dphi > fDijetPhiCut );
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::isGoodDijet " << isGood << " ";
        std::cout << Form("pTlead: %5.2f pTsub: %5.2f dphi: %4.2f\n", ptLead, ptSublead, dphi);
    }
    return isGood;
}

//________________
void DiJetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    if ( fVerbose ) {
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "DiJetAnalysis::processEvent" << std::endl;
    }

    fEventCounter++;
    if ( fEventCounter >= 50000 ) {
        fCycleCounter++;
        std::cout << Form("DiJetAnalysis::processEvent [INFO] Events processed: %d Sample fraction: %3.2f%%", 
                          fCycleCounter * 50000, (Double_t)(fCycleCounter * 50000) / fNEventsInSample )
                  << std::endl;
        fEventCounter = {0};
    }

    if ( !fHM ) {
        std::cout << "[Warning] No histogram manager connected to the DiJetAnalysis\n";
    }

    // Must be flushed for each event !!!!
    fIsDijetFound = {kFALSE};
    fIsDijetJetIdFound = {kFALSE};

    //
    // Event quantities
    //

    // ptHat
    Double_t ptHat = event->ptHat();
    Double_t vz = event->vz();
    // ptHat weight (a.k.a. event weight that includes the ptHat and vz)
    Double_t ptHatW{1.};
    // Check correct MC sample
    if ( fIsMc && fIsPPb ) {
        // Skip events with ptHat that is outside the ranged embedded
        if ( ptHat <= fPtHatRange[0] || ptHat > fPtHatRange[1] ) {
            if ( fVerbose ) {
                std::cout << Form("[WARNING] Bad ptHat value: %4.1f < %4.1f <= %4.1f\n", fPtHatRange[0], ptHat, fPtHatRange[1]);
            }
            return;
        }

        // For MC we need to flip the direction of Pb-going to properly reweight distributions
        if ( fIsPbGoingDir ) {
            vz = -vz;
        }
        ptHatW = eventWeight(fIsMc, fIsPPb, ptHat, vz);
    }

    // Check collision system and account for the lab frame rapidity
    if ( fUseEtaShiftAndSignSwap && fIsPPb ) {
        if ( fIsMc ) { // For MC Pb is going to negative eta
            if ( fIsPbGoingDir) {
                fEtaShift = {-0.4654094531};
            }
            else { // For MC p is going to positive eta
                fEtaShift = {0.4654094531};
            }
        }
        else {
            if ( fIsPbGoingDir) { // For data Pb is going to positive eta
                fEtaShift = {0.4654094531};
            }
            else { // For data p is going to negative eta
                fEtaShift = {-0.4654094531};
            }
        }
    }
    else {
        fEtaShift = {0};
    }

    // Process and analyze reco jets
    processRecoJets(event, ptHatW);

    if ( fIsMc ) {
        // Process and analyze gen jets
        processGenJets(event, ptHatW);
        processRefJets(event, ptHatW);
    }

    Double_t centW = event->centralityWeight();
    centW = {1.}; // Do not apply weight for pPb
    //std::cout << "centrality weight: " << centW << std::endl;

    // For dijet analysis and reweighting purposes it is important 
    // to fill event histograms only when dijet is found
    if ( fIsDijetFound ) {
        fHM->hHiBin->Fill( event->hiBin() );

        fHM->hVz->Fill( event->vz(),  centW );
        fHM->hVzWeighted->Fill( event->vz(), ptHatW * centW );

        fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW * centW );

        fHM->hPtHatWeighted->Fill( ptHat, ptHatW * centW );
        fHM->hPtHatWeight->Fill( ptHatW, centW );

        Double_t vzPtHat[2] = { vz, ptHat };
        fHM->hVzPtHat->Fill( vzPtHat, centW );
        fHM->hVzPtHatWeighted->Fill( vzPtHat, ptHatW * centW );
    } // if ( fIsDijetFound ) 

    // Fill this outside to be sure that MC is calculated properly,
    // i.e. different ptHat samples are taken with the proper weight
    fHM->hPtHat->Fill( ptHat, centW );

    if ( fVerbose ) {
        std::cout << "Event quantities were read properly" << std::endl;
        //event->print();
        std::cout << "DiJetAnalysis::processEvent - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::finish() {
    // Save data and close files
    fCycleCounter++;
    std::cout << Form("DiJetAnalysis::processEvent [INFO] Total events processed: %d Sample fraction: %3.2f%%", 
                      (fCycleCounter * 50000) + fEventCounter, 
                      (Double_t)(fCycleCounter * 50000 + fEventCounter) / fNEventsInSample )
              << std::endl;
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
