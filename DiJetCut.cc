// Jet analysis headers
#include "DiJetCut.h"

// ROOT headers
#include "TString.h"

// C++ headers
#include <iostream>
#include <cmath>

//________________
DiJetCut::DiJetCut() : fLeadJetPt{0.0f}, fSubLeadJetPt{0.0f},
                       fLeadJetEtaLab{-1000.f, 1000.f}, fSubLeadJetEtaLab{-1000.f, 1000.f},
                       fLeadJetEtaCM{-1000.f, 1000.f}, fSubLeadJetEtaCM{-1000.f, 1000.f},
                       fDijetDPhiCut{0.f}, fVerbose{false} {
    // Default constructor initializes all members to default values
}

//________________
DiJetCut::~DiJetCut() {
    // Destructor does nothing for now
}

//________________
DiJetCut& DiJetCut::operator=(const DiJetCut& other) {
    // Assignment operator
    if (this != &other) {
        fLeadJetPt = other.fLeadJetPt;
        fSubLeadJetPt = other.fSubLeadJetPt;
        fLeadJetEtaLab[0] = other.fLeadJetEtaLab[0];
        fLeadJetEtaLab[1] = other.fLeadJetEtaLab[1];
        fSubLeadJetEtaLab[0] = other.fSubLeadJetEtaLab[0];
        fSubLeadJetEtaLab[1] = other.fSubLeadJetEtaLab[1];
        fLeadJetEtaCM[0] = other.fLeadJetEtaCM[0];
        fLeadJetEtaCM[1] = other.fLeadJetEtaCM[1];
        fSubLeadJetEtaCM[0] = other.fSubLeadJetEtaCM[0];
        fSubLeadJetEtaCM[1] = other.fSubLeadJetEtaCM[1];
        fDijetDPhiCut = other.fDijetDPhiCut;
        fVerbose = other.fVerbose;
    }
    return *this;
}

//________________
bool DiJetCut::operator==(const DiJetCut& other) const {
    // Compare all members for equality
    return (fLeadJetPt == other.fLeadJetPt &&
            fSubLeadJetPt == other.fSubLeadJetPt &&
            fLeadJetEtaLab[0] == other.fLeadJetEtaLab[0] &&
            fLeadJetEtaLab[1] == other.fLeadJetEtaLab[1] &&
            fSubLeadJetEtaLab[0] == other.fSubLeadJetEtaLab[0] &&
            fSubLeadJetEtaLab[1] == other.fSubLeadJetEtaLab[1] &&
            fLeadJetEtaCM[0] == other.fLeadJetEtaCM[0] &&
            fLeadJetEtaCM[1] == other.fLeadJetEtaCM[1] &&
            fSubLeadJetEtaCM[0] == other.fSubLeadJetEtaCM[0] &&
            fSubLeadJetEtaCM[1] == other.fSubLeadJetEtaCM[1] &&
            fDijetDPhiCut == other.fDijetDPhiCut &&
            fVerbose == other.fVerbose);
}

//________________
void DiJetCut::report() {
    // Report the cut parameters
    std::cout << "DiJetCut parameters:" << std::endl;
    std::cout << "--> Minimum leading jet pT: " << fLeadJetPt << std::endl;
    std::cout << "--> Minimum subleading jet pT: " << fSubLeadJetPt << std::endl;
    std::cout << "--> Leading jet eta (lab frame): [" << fLeadJetEtaLab[0] << ", " << fLeadJetEtaLab[1] << "]" << std::endl;
    std::cout << "--> Subleading jet eta (lab frame): [" << fSubLeadJetEtaLab[0] << ", " << fSubLeadJetEtaLab[1] << "]" << std::endl;
    std::cout << "--> Leading jet eta (CM frame): [" << fLeadJetEtaCM[0] << ", " << fLeadJetEtaCM[1] << "]" << std::endl;
    std::cout << "--> Subleading jet eta (CM frame): [" << fSubLeadJetEtaCM[0] << ", " << fSubLeadJetEtaCM[1] << "]" << std::endl;
    std::cout << "--> Dijet dPhi cut: " << fDijetDPhiCut << std::endl;
}

//________________
bool DiJetCut::pass(const DiJet* dijet, bool isCM) {
    // Check if the dijet passes the cut criteria
    if (!dijet) {
        std::cerr << "DiJetCut::pass: dijet pointer is null!" << std::endl;
        return false;
    }

    float eta1 = ( isCM ) ? dijet->leadJetEtaCM() : dijet->leadJetEtaLab();
    float eta2 = ( isCM ) ? dijet->subLeadJetEtaCM() : dijet->subLeadJetEtaLab();
    float pt1 = dijet->leadJetPt();
    float pt2 = dijet->subLeadJetPt();
    float dphi = fabs( dijet->dPhi() );

    if ( fVerbose ) {
        std::cout << "DiJetCut::pass: Checking dijet with parameters in " << (isCM ? "CM" : "lab") << " frame:" << std::endl;
        std::cout << Form("--> Lead jet pT: %.2f > %.2f  %s\n", pt1, fLeadJetPt, ((pt1 > fLeadJetPt) ? "[passed]" : "[failed]"));
        std::cout << Form("--> Sublead jet pT: %.2f > %.2f  %s\n", pt2, fSubLeadJetPt, (pt2 > fSubLeadJetPt) ? "[passed]" : "[failed]");
        std::cout << Form("--> Lead jet eta: %.2f in [%f, %f]  %s\n", eta1, fLeadJetEtaLab[0], fLeadJetEtaLab[1], 
                          ((eta1 >= fLeadJetEtaLab[0] && eta1 <= fLeadJetEtaLab[1]) ? "[passed]" : "[failed]"));
        std::cout << Form("--> Sublead jet eta: %.2f in [%f, %f]  %s\n", eta2, fSubLeadJetEtaLab[0], fSubLeadJetEtaLab[1], 
                          ((eta2 >=  fSubLeadJetEtaLab[0] && eta2 <= fSubLeadJetEtaLab[1]) ? "[passed]" : "[failed]"));
        std::cout << Form("--> Dijet dPhi: %.2f > %.2f  %s\n", dphi, fDijetDPhiCut,
                          ((dphi > fDijetDPhiCut) ? "[passed]" : "[failed]"));
    }

    if (pt1 < fLeadJetPt) return false;
    if (pt2 < fSubLeadJetPt) return false;
    if ( isCM ) {
        if (eta1 < fLeadJetEtaCM[0] || eta1 > fLeadJetEtaCM[1]) return false;
        if (eta2 < fSubLeadJetEtaCM[0] || eta2 > fSubLeadJetEtaCM[1]) return false;
    } 
    else {
        if (eta1 < fLeadJetEtaLab[0] || eta1 > fLeadJetEtaLab[1]) return false;
        if (eta2 < fSubLeadJetEtaLab[0] || eta2 > fSubLeadJetEtaLab[1]) return false;
    }
    if (dphi < fDijetDPhiCut) return false;

    return true;
}
