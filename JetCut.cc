/**
 * @file JetCut.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Example of jet cut
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "JetCut.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

// C++ headers
#include <iostream>

ClassImp(JetCut)

//________________
JetCut::JetCut() : fRecoPt{0., 1e6}, fRecoConeR{1e6},
    fMustHaveGenMatching{kFALSE}, fRefPt{0., 1e6},
    fRefConeR{1e6}, fRefFlavorForB{-100000, 100000},
    fJetsPassed{0}, fJetsFailed{0} {
    /* Empty */
}

//________________
JetCut::~JetCut() {
    /* Empty */
}

//________________
void JetCut::report() {
    TString report = "\nReporting from JetCut";
    report += TString::Format( "reco pT         :\t %f - %f\n", fRecoPt[0], fRecoPt[1] );
    report += TString::Format( "reco R          :\t %f\n", fRecoConeR );
    report += TString::Format( "match to gen    :\t %d\n", fMustHaveGenMatching);
    report += TString::Format( "gen pT          :\t %f - %f\n", fRefPt[0], fRefPt[1] );
    report += TString::Format( "gen R           :\t %f\n", fRefConeR );
    report += TString::Format( "flavorForB      :\t %d - %d\n", fRefFlavorForB[0], fRefFlavorForB[1] );
    report += TString::Format( "Jets passed     :\t %lld\n", fJetsPassed );
    report += TString::Format( "Jets failed     :\t %lld\n", fJetsFailed );
    std::cout << report.Data() << std::endl;
}

//________________
Bool_t JetCut::pass(const Jet* jet) {
    if (fVerbose) {
        std::cout << "\n----- Jet cut -----\n";
    }

    Bool_t goodRecoPt = (fRecoPt[0] <= jet->recoJetPtJECCorr() &&
                         jet->recoJetPtJECCorr() <= fRecoPt[1]);
    if (fVerbose) {
        std::cout << Form("reco pT    : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fRecoPt[0], jet->recoJetPtJECCorr(), fRecoPt[1], ( goodRecoPt ) ? "true" : "false" );
    }

    Float_t recoR = TMath::Sqrt( jet->recoJetPhi() * jet->recoJetPhi() + 
                                 jet->recoJetEta() * jet->recoJetEta() );
    Bool_t goodRecoConeR = (recoR <= fRecoConeR);
    if (fVerbose) {
        std::cout << Form("reco cone R: %5.2f <= %5.2f \t %s \n",
                          recoR, fRecoConeR, ( goodRecoConeR ) ? "true" : "false" );
    }

    Bool_t goodMatching = kTRUE;
    if (fMustHaveGenMatching) {
        goodMatching = jet->hasMatching();
    }
    if (fVerbose) {
        std::cout << Form("has matching: \t %s \n",
                          ( goodMatching ) ? "true" : "false" );        
    }

    Bool_t goodRefPt = (fRefPt[0] <= jet->refJetPt() &&
                               jet->refJetPt() <= fRefPt[1]);
    if (fVerbose) {
        std::cout << Form("gen pT    : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fRefPt[0], jet->refJetPt(), fRefPt[1], ( goodRefPt ) ? "true" : "false" );
    }

    Float_t refR = TMath::Sqrt( jet->refJetPhi() * jet->refJetPhi() + 
                                jet->refJetEta() * jet->refJetEta() );
    Bool_t goodRefConeR = (refR <= fRefConeR);
    if (fVerbose) {
        std::cout << Form("ref cone R: %5.2f <= %5.2f \t %s \n",
                          refR, fRefConeR, ( goodRefConeR ) ? "true" : "false" );
    }

    Bool_t goodFlavorForB = ( fRefFlavorForB[0] <= jet->refFlavorForB() &&
                                    jet->refFlavorForB() <= fRefFlavorForB[1] );
    if (fVerbose) {
        std::cout << Form("ref flavorB: %d <= %d <= %d \t %s \n",
                          fRefFlavorForB[0], jet->refFlavorForB(), fRefFlavorForB[1], ( goodFlavorForB ) ? "true" : "false" );
    }

    Bool_t isGood = goodRecoPt && goodRecoConeR && goodMatching &&
                          goodRefPt && goodRefConeR && goodFlavorForB;

    if (fVerbose) {
        std::cout << Form("good jet     : \t %s \n", (isGood) ? "true" : "false");
    }

    return isGood;
}