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
#include <cmath>

//________________
JetCut::JetCut() : fPt{0., 1e6}, fConeR{1e6},
    fEtaLab{-1e6, 1e6}, fEtaCM{-1e6, 1e6},
    fSelectionMethod{0}, fLooseJetIdCut{true}, fVerbose{false} {
    /* Empty */
}

//________________
JetCut::~JetCut() {
    /* Empty */
}

//________________
JetCut& JetCut::operator=(const JetCut& other) {
    // Assignment operator
    if (this != &other) {
        fPt[0] = other.fPt[0];
        fPt[1] = other.fPt[1];
        fConeR = other.fConeR;
        fEtaLab[0] = other.fEtaLab[0];
        fEtaLab[1] = other.fEtaLab[1];
        fEtaCM[0] = other.fEtaCM[0];
        fEtaCM[1] = other.fEtaCM[1];
        fSelectionMethod = other.fSelectionMethod;
        fVerbose = other.fVerbose;
    }
    return *this;
}

//________________
bool JetCut::operator==(const JetCut& other) const {
    // Compare all members for equality
    return (fPt[0] == other.fPt[0] && fPt[1] == other.fPt[1] &&
            fConeR == other.fConeR && 
            fEtaLab[0] == other.fEtaLab[0] && fEtaLab[1] == other.fEtaLab[1] &&
            fEtaCM[0] == other.fEtaCM[0] && fEtaCM[1] == other.fEtaCM[1] &&
            fSelectionMethod == other.fSelectionMethod && fVerbose == other.fVerbose);
}

//________________
void JetCut::report() {
    TString report = "\nReporting from JetCut";
    report += TString::Format( "pT              :\t %f - %f\n", fPt[0], fPt[1] );
    report += TString::Format( "cone R          :\t %f\n", fConeR );
    report += TString::Format( "eta lab         :\t %f - %f\n", fEtaLab[0], fEtaLab[1] );
    report += TString::Format( "eta CM          :\t %f - %f\n", fEtaCM[0], fEtaCM[1] );
    report += TString::Format( "selection method:\t %d\n", fSelectionMethod );
    std::cout << report.Data() << std::endl;
}

//________________
bool JetCut::pass(const RecoJet* jet, bool isCM, bool isMC, bool requireMatching) {
    if (fVerbose) {

    }

    bool goodPt = (fPt[0] <= jet->ptJECCorr() && jet->ptJECCorr() <= fPt[1]);
    float recoR = sqrt( jet->phi() * jet->phi() + jet->eta() * jet->eta() );
    bool goodConeR = (recoR <= fConeR);
    bool goodMatching {true};
    if (isMC && requireMatching) {
        goodMatching = jet->hasMatching();
    }
    bool goodEta{false}; 
    if (isCM) {
        goodEta = ( fEtaCM[0] <= jet->eta() && jet->eta() <= fEtaCM[1] );
    } 
    else {
        goodEta = ( fEtaLab[0] <= jet->eta() && jet->eta() <= fEtaLab[1] );
    }

    bool goodSelection{false};
    if (fSelectionMethod == 0) {
        goodSelection = true; // no selection
    } 
    else if (fSelectionMethod == 1) {
        // trackMaxPt/RawPt selection
        goodSelection = jet->isGoodTrkMax();
    } 
    else if (fSelectionMethod == 2) {
        // jetId selection
        goodSelection = jet->isGoodJetId(fLooseJetIdCut);
    }

    bool isGood = goodPt && goodConeR && goodMatching && goodEta && goodSelection;

    if (fVerbose) {
        std::cout << "\n----- Reco jet cut -----\n";
        std::cout << Form("--> isCM: %d, isMC: %d, requireMatching: %s\n", isCM, isMC, (requireMatching) ? "[passed]" : "[failed]");
        std::cout << Form("--> pT : %5.2f <= %5.2f <= %5.2f \t %s \n", fPt[0], jet->ptJECCorr(), fPt[1], ( goodPt ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> cone R: %5.2f <= %5.2f \t %s \n", recoR, fConeR, ( goodConeR ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> has matching: \t %s \n", ( goodMatching ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> eta : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          (isCM) ? fEtaCM[0] : fEtaLab[0], jet->eta(),
                          (isCM) ? fEtaCM[1] : fEtaLab[1],
                          ( goodEta ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> selection method: %d \t %s \n", fSelectionMethod, 
                          ( goodSelection ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> good reco jet : \t %s \n", (isGood) ? "[passed]" : "[failed]");
    }

    return isGood;
}

//________________
bool JetCut::pass(const GenJet* jet, bool isCM) {
    if (fVerbose) {
        std::cout << "\n----- Gen jet cut -----\n";
        std::cout << Form("isCM: %d\n", isCM);
    }

    bool goodPt = (fPt[0] <= jet->pt() && jet->pt() <= fPt[1]);
    float genR = sqrt( jet->phi() * jet->phi() + jet->eta() * jet->eta() );
    bool goodConeR = (genR <= fConeR);

    bool goodEta{false}; 
    if (isCM) {
        goodEta = ( fEtaCM[0] <= jet->eta() && jet->eta() <= fEtaCM[1] );
    } 
    else {
        goodEta = ( fEtaLab[0] <= jet->eta() && jet->eta() <= fEtaLab[1] );
    }

    bool isGood = goodPt && goodConeR && goodEta;

    if (fVerbose) {
        std::cout << Form("--> pT : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          fPt[0], jet->pt(), fPt[1], ( goodPt ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> cone R: %5.2f <= %5.2f \t %s \n",
                          genR, fConeR, ( goodConeR ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> eta : %5.2f <= %5.2f <= %5.2f \t %s \n",
                          (isCM) ? fEtaCM[0] : fEtaLab[0], jet->eta(),
                          (isCM) ? fEtaCM[1] : fEtaLab[1],
                          ( goodEta ) ? "[passed]" : "[failed]" );
        std::cout << Form("--> good gen jet : \t %s \n", (isGood) ? "[passed]" : "[failed]");
    }

    return isGood;
}
