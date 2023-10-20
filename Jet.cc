/**
 * @file Jet.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet class description
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "Jet.h"

// ROOT headers
#include "TString.h"

// C++ headers
#include <iostream>

ClassImp(Jet)

//________________
Jet::Jet() : TObject(), fRecoPt{-999.f}, 
    fRecoEta{-999.f}, fRecoPhi{-999.f}, 
    fRecoPtJECCorr{-999.f}, fRecoWTAeta{-999.f}, fRecoWTAphi{-999.f},
    fRecoJetPtWeight{-999.f}, fRecoJetPtSmearingWeight{-999.f},
    fRefPt{-999.f}, fRefEta{-999.f}, fRefPhi{-999.f}, 
    fRefFlavor{-999}, fRefFlavorForB{-99}, fRefPtWeight{-999.f} {
    // Empty
}

//________________
Jet::~Jet() {
    // Empty
}

//________________
void Jet::print() {
    std::cout << Form("---------------------------------\n")
              << Form("Reconstructed jet params:\n")
              << Form("pT: %5.2f  eta: %3.2f  phi: %3.2f  pTcorr: %5.2f  WTAeta: %3.2f  WTAphi: %3.2f\n",
                      fRecoPt, fRecoEta, fRecoPhi, fRecoPtJECCorr, fRecoWTAeta, fRecoWTAphi)
              << Form("Matched jet params:\n")
              << Form("pT: %5.2f  eta: %3.2f  phi: %3.2f  flavor: %d  flavorForB: %d\n",
                      fRefPt, fRefEta, fRefPhi, refFlavor(), refFlavorForB() )
              << Form("---------------------------------\n");
}
