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
    fLeadJetPtLow{30.}, fSubleadJetPtLow{20.},
    fDijetPhiCut{TMath::TwoPi() / 3},
    fIsPbGoingDir{kFALSE}, fVerbose{kFALSE},
    fNEventsInSample{1000000} {
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
    std::cout << "Use centrality weight   : " << fUseCentralityWeight << std::endl
              << "Histogram manager       : " << fHM << std::endl
              << "Is MC                   : " << fIsMc << std::endl
              << "Is pPb                  : " << fIsPPb << std::endl
              << "Is Pb-going direction   : " << fIsPbGoingDir << std::endl
              << "eta shift               : " << fEtaShift << std::endl
              << "ptHat range             : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl
              << "Leading jet pT          : " << fLeadJetPtLow << std::endl
              << "SubLeading jet pT       : " << fSubleadJetPtLow << std::endl
              << "Dijet phi cut           : " << fDijetPhiCut << std::endl;
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
void DiJetAnalysis::processGenJets(const Event* event, Double_t ptHatW) {

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processGenJets" << std::endl;
    }

    Double_t ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
             phiLead{0.},  phiSubLead{0.};
    Int_t idLead{-1}, idSubLead{-1};
    Bool_t isDijetFound{kFALSE};

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

        // Apply lab frame boost to CM for the pPb 
        if ( fIsPPb ) {
            if ( fIsPbGoingDir ) {
                eta -= fEtaShift;
                eta = -eta;
            }
            else {
                eta += fEtaShift;
            }
        } // if ( fIsPPb ) 
        
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

    // Check if two jets found
    if (idLead>=0 && idSubLead>=0) {
        isDijetFound = {kTRUE};
    }
    if ( !isDijetFound ) return;

    // Check the dijet selection on the MC level
    if ( !isGoodDijet(ptLead, ptSubLead, TMath::Abs( deltaPhi(phiLead, phiSubLead) ) ) ) return;

    fHM->hGenPtLeadPtSublead->Fill(ptLead, ptSubLead, ptHatW);
    fHM->hGenEtaLeadEtaSublead->Fill(etaLead, etaSubLead, ptHatW);

    Double_t dijetPt = 0.5 * (ptLead + ptSubLead);
    Double_t dijetEta = 0.5 * (etaLead + etaSubLead);
    Double_t dijetDphi = deltaPhi(phiLead, phiSubLead); // Should abs be used here

    Double_t genDijetLeadSublead[9] {dijetPt, dijetEta, dijetDphi, 
                                     ptLead, etaLead, phiLead, 
                                     ptSubLead, etaSubLead, phiSubLead };
    fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(genDijetLeadSublead);
    fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(genDijetLeadSublead, ptHatW);
    fHM->hGenDijetEta->Fill(dijetEta, ptHatW);
    fHM->hGenDijetPtEtaDphi->Fill(dijetPt, dijetEta, dijetDphi, ptHatW);

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processGenJets - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRecoJets(const Event* event, Double_t ptHatW) {

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processRecoJets" << std::endl;
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

        Double_t pt = (*pfJetIter)->ptJECCorr();
        Double_t eta = (*pfJetIter)->eta();
        Double_t phi = (*pfJetIter)->phi();
        Double_t ptRaw = (*pfJetIter)->pt();

        Double_t rawPt = (*pfJetIter)->rawPt();
        Double_t trackMaxPt = (*pfJetIter)->trackMaxPt();
        if ( TMath::Abs( eta ) < 2.4 && 
            ( trackMaxPt/rawPt < 0.01 ||
              trackMaxPt/rawPt > 0.98 ) ) {
            // Remove charge component at midrapidity
            continue;
        }

        if ( fVerbose ) {
            std::cout << "Reco jet #" << counter << " ";
            (*pfJetIter)->print();
        }
        
        GenJet *matchedJet{nullptr};
        Double_t genPt{999.};
        Double_t genEta{-999.};
        Double_t genPhi{-999.};

        // On MC will work with matching jets only
        if ( fIsMc ) {
            if ( !(*pfJetIter)->hasMatching() ) continue;
            matchedJet = event->genJetCollection()->at( (*pfJetIter)->genJetId() );
            genPt = matchedJet->pt();
            genEta = matchedJet->eta();
            genPhi = matchedJet->phi();
        }

        // Apply lab frame boost to CM for the pPb 
        if ( fIsPPb ) {
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

            ptRefSubLead = ptRefLead;
            etaRefSubLead = etaRefLead;
            phiRefSubLead = phiRefLead;
            ptRefLead = genPt;
            etaRefLead = genEta;
            phiRefLead = genPhi;
        }
        else if ( pt > ptRecoSubLead ) {
            ptRecoSubLead = pt;
            ptRawRecoSubLead = ptRaw;
            etaRecoSubLead = eta;
            phiRecoSubLead = phi;
            idRecoSubLead = counter;

            ptRefSubLead = genPt;
            etaRefSubLead = genEta;
            phiRefSubLead = genPhi;
        }
    
        // Fill inclusive jet information
        fHM->hRecoInclusiveJetPt->Fill(pt, ptHatW);
        fHM->hRecoMatchedPtEta->Fill(eta, pt, ptHatW);
        if ( fIsMc ) {
            fHM->hRefInclusiveJetPt->Fill(genPt, ptHatW);
            fHM->hRefInclusiveJetPtEta->Fill(genEta, genPt, ptHatW);

            Double_t correl[5] { pt, ptRaw, genPt, eta, genEta };
            fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
            fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl, ptHatW);

            Double_t res[5] { pt/genPt, genPt, eta, phi, event->ptHat() };
            fHM->hJESInclusiveJetPtEtaPhiPtHat->Fill(res);
            fHM->hJESInclusiveJetPtEtaPhiPtHatWeighted->Fill(res, ptHatW);
        }

        if ( fVerbose ) {
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d\n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead);
        }

        // Increment counter
        counter++;
    } // for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ )

    // Check if leading and subleading jets were found
    if (idRecoLead>=0 && idRecoSubLead>=0) {
        isDijetFound = kTRUE;
    }

    if ( fVerbose ) {
        std::cout << "Dijet found: " << isDijetFound << std::endl;
    }
    
    // Look only at events with dijets
    if ( !isDijetFound ) return;

    // Check the dijet selection on the MC level
    if ( !isGoodDijet(ptRecoLead, ptRecoSubLead, TMath::Abs( deltaPhi(phiRecoLead, phiRecoSubLead) ) ) ) return;

    if ( fIsMc ) {
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
    }

    // Correlation between leading and subleading
    fHM->hRecoPtLeadPtSublead->Fill(ptRecoLead, ptRecoSubLead, ptHatW);
    fHM->hRecoEtaLeadEtaSublead->Fill(etaRecoLead, etaRecoSubLead, ptHatW);

    // Dijet analysis
    Double_t dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
    Double_t dijetRecoEta = 0.5 * (etaRecoLead + etaRecoSubLead);
    Double_t dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead); // Should abs be used here?

    Double_t dijetRefPt{-999.};
    Double_t dijetRefEta{-999.};
    Double_t dijetRefDphi {-999.};

    if ( fIsMc ) {
        dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
        dijetRefEta = 0.5 * (etaRefLead + etaRefSubLead);
        dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);
    }

    Double_t dijetRecoInfo[9] { dijetRecoPt, dijetRecoEta, dijetRecoDphi,
                                ptRecoLead, etaRecoLead, phiRecoLead,
                                ptRecoSubLead, etaRecoSubLead, phiRecoSubLead };
    fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(dijetRecoInfo);
    fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(dijetRecoInfo, ptHatW);
    fHM->hRecoDijetEta->Fill( dijetRecoEta, ptHatW );
    fHM->hRecoDijetPtEtaDphi->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, ptHatW );

    // Dijet reco vs ref for unfolding
    Double_t dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                     ptRecoLead, etaRecoLead,
                                     ptRecoSubLead, etaRecoSubLead,
                                     dijetRefPt, dijetRefEta,
                                     ptRefLead, etaRefLead,
                                     ptRefSubLead, etaRefSubLead };
    if ( fIsMc ) {
        fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Fill(dijetRecoUnfold);
        fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, ptHatW);
        fHM->hRefDijetEta->Fill( dijetRefEta, ptHatW );
        fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEta, dijetRefEta, ptHatW );
        fHM->hRefDijetPtEtaDphi->Fill( dijetRefPt, dijetRefEta,  dijetRefDphi, ptHatW );
    }

    if ( fVerbose ) {
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
            std::cout << "Reco jet #" << counter << " ";
            (*pfJetIter)->print();
        }

        // Apply lab frame boost to CM for the pPb 
        if ( fIsPPb ) {
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
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d\n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead);
        }

        // Increment counter
        counter++;
    } // for ( pfJetIter = event->pfJetCollection()->begin(); pfJetIter != event->pfJetCollection()->end(); pfJetIter++ )

    // Check if leading and subleading jets were found
    if (idRecoLead>=0 && idRecoSubLead>=0) {
        isDijetFound = kTRUE;
    }

    if ( fVerbose ) {
        std::cout << "Dijet found: " << isDijetFound << std::endl;
    }
    
    // Look only at events with dijets
    if ( !isDijetFound ) return;

    // Check the dijet selection on the MC level
    if ( !isGoodDijet(ptRefLead, ptRefSubLead, TMath::Abs( deltaPhi(phiRefLead, phiRefSubLead) ) ) ) return;

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
    fHM->hRefSelDijetPtEtaDphi->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW);
    fHM->hRefSelDijetEta->Fill(dijetRefEta, ptHatW);

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processRefJets - [DONE]" << std::endl;
    }
}

//________________
Bool_t DiJetAnalysis::isGoodDijet(const Double_t& ptLead, const Double_t& ptSublead, const Double_t& dphi) {
    Bool_t isGood = ( ptLead > fLeadJetPtLow &&
                      ptSublead > fSubleadJetPtLow /* && dphi > fDijetPhiCut */ );
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
        std::cout << "DiJetAnalysis::processEvent" << std::endl;
    }

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
    if ( fIsPPb ) {
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

    if ( fVerbose ) {
        std::cout << "Event quantities were read properly" << std::endl;
        //event->print();
    }

    if ( fIsMc ) {
        // Process and analyze gen jets
        processGenJets(event, ptHatW);
        processRefJets(event, ptHatW);
    }

    // Process and analyze reco jets
    processRecoJets(event, ptHatW);

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processEvent - [DONE]" << std::endl;
    }
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
