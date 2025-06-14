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
#include "GenJet.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

// C++ headers
#include <iostream>

//________________
JetCut::JetCut() : fPt{0., 1e6}, fConeR{1e6},
    fMustHaveGenMatching{false}, fEta{-1e6, 1e6},
    fTrackMaxPtOverRawPt{0.01, 0.98},
	fVerbose{false}, fJetsPassed{0}, fJetsFailed{0} {
    /* Empty */
}

//________________
JetCut::~JetCut() {
    /* Empty */
}

//________________
void JetCut::report() {
    TString report = "\nReporting from JetCut";
    report += TString::Format( "pT              :\t %f - %f\n", fPt[0], fPt[1] );
    report += TString::Format( "cone R          :\t %f\n", fConeR );
    report += TString::Format( "match to gen    :\t %d\n", fMustHaveGenMatching);
    report += TString::Format( "eta             :\t %f - %f\n", fEta[0], fEta[1] );
    report += TString::Format( "Jets passed     :\t %ld\n", fJetsPassed );
    report += TString::Format( "Jets failed     :\t %ld\n", fJetsFailed );
    std::cout << report.Data() << std::endl;
}

//________________
bool JetCut::pass(const RecoJet* jet) {
    if (fVerbose) {
        std::cout << "\n----- Reco jet cut -----\n";
    }

    bool goodPt = (fPt[0] <= jet->ptJECCorr() &&
                     jet->ptJECCorr() <= fPt[1]);
    if (fVerbose) {
        std::cout << Form("pT : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fPt[0], jet->ptJECCorr(), fPt[1], ( goodPt ) ? "[true]" : "[false]" );
    }

    float recoR = TMath::Sqrt( jet->phi() * jet->phi() + jet->eta() * jet->eta() );
    bool goodConeR = (recoR <= fConeR);
    if (fVerbose) {
        std::cout << Form("cone R: %5.2f <= %5.2f \t %s \n",
                          recoR, fConeR, ( goodConeR ) ? "[true]" : "[false]" );
    }

    bool goodMatching {true};
    if (fMustHaveGenMatching) {
        goodMatching = jet->hasMatching();
    }

    if (fVerbose) {
        std::cout << Form("has matching: \t %s \n",
                          ( goodMatching ) ? "[true]" : "[false]" );
    }

    bool goodEta = ( fEta[0] <= jet->eta() && jet->eta() <= fEta[1] );
    if (fVerbose) {
        std::cout << Form("eta : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fEta[0], jet->eta(), fEta[1], ( goodEta ) ? "[true]" : "[false]" );
    }

    bool goodChargeComponent{true};
    // Float_t rawPt = jet->rawPt();
    // Float_t trackMaxPt = jet->trackMaxPt();
    // if ( TMath::Abs( jet->eta() ) < 2.4 && 
    //      ( trackMaxPt/rawPt < fTrackMaxPtOverRawPt[0] ||
    //        trackMaxPt/rawPt > fTrackMaxPtOverRawPt[1]) ) {
    //     goodChargeComponent = {false};
    // }

    // if (fVerbose) {
    //     std::cout << Form("rawPt: %5.2f trackMaxPt: %5.2f ptMax/rawPt: %3.2f lowCut: %3.2f highCut: %3.2f\n",
    //                       rawPt, trackMaxPt, trackMaxPt/rawPt, fTrackMaxPtOverRawPt[0], fTrackMaxPtOverRawPt[1] );
    // }

    // if ( goodMatching )
    // bool goodRefPt = (fRefPt[0] <= jet->refJetPt() &&
    //                            jet->refJetPt() <= fRefPt[1]);
    // if (fVerbose) {
    //     std::cout << Form("gen pT    : %5.2f <= %5.2f <= %5.2f \t %s \n",
    //                       fRefPt[0], jet->refJetPt(), fRefPt[1], ( goodRefPt ) ? "true" : "false" );
    // }

    // Float_t refR = TMath::Sqrt( jet->refJetPhi() * jet->refJetPhi() + 
    //                             jet->refJetEta() * jet->refJetEta() );
    // bool goodRefConeR = (refR <= fRefConeR);
    // if (fVerbose) {
    //     std::cout << Form("ref cone R: %5.2f <= %5.2f \t %s \n",
    //                       refR, fRefConeR, ( goodRefConeR ) ? "true" : "false" );
    // }

    // bool goodFlavorForB = ( fRefFlavorForB[0] <= jet->refFlavorForB() &&
    //                                 jet->refFlavorForB() <= fRefFlavorForB[1] );
    // if (fVerbose) {
    //     std::cout << Form("ref flavorB: %d <= %d <= %d \t %s \n",
    //                       fRefFlavorForB[0], jet->refFlavorForB(), fRefFlavorForB[1], ( goodFlavorForB ) ? "true" : "false" );
    // }


    bool isGood = goodPt && goodConeR && goodMatching && goodEta && goodChargeComponent;

    // bool isGood = goodRecoPt && goodRecoConeR && goodMatching &&
    //                       goodRefPt && goodRefConeR && goodFlavorForB;

    if (fVerbose) {
        std::cout << Form("good jet     : \t %s \n", (isGood) ? "[true]" : "[false]");
    }

    return isGood;
}

//________________
bool JetCut::pass(const GenJet* jet) {
    if (fVerbose) {
        std::cout << "\n----- Gen jet cut -----\n";
    }

    bool goodPt = (fPt[0] <= jet->pt() && jet->pt() <= fPt[1]);
    if (fVerbose) {
        std::cout << Form("pT : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fPt[0], jet->pt(), fPt[1], ( goodPt ) ? "[true]" : "[false]" );
    }
    Float_t genR = TMath::Sqrt( jet->phi() * jet->phi() + jet->eta() * jet->eta() );
    bool goodConeR = (genR <= fConeR);
    if (fVerbose) {
        std::cout << Form("cone R: %5.2f <= %5.2f \t %s \n",
                          genR, fConeR, ( goodConeR ) ? "[true]" : "[false]" );
    }
    bool goodEta = ( fEta[0] <= jet->eta() && jet->eta() <= fEta[1] );
    if (fVerbose) {
        std::cout << Form("eta : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fEta[0], jet->eta(), fEta[1], ( goodEta ) ? "[true]" : "[false]" );
    }

    return goodPt;
}
