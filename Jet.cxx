// Jet analysis headers
#include "Jet.h"

// ROOT headers
#include "TString.h"

// C++ headers
#include <iostream>

ClassImp(Jet)

//________________
Jet::Jet() : TObject(), fRecoPt{-999.}, 
    fRecoEta{-999.}, fRecoPhi{-999.}, fRefPt{-999.},
    fRefEta{-999.}, fRefPhi{-999.}, 
    fRefFlavor{-999}, fRefFlavorForB{-99}, 
    fGenPt{-999.}, fGenEta{-999}, fGenPhi{-999.}, 
    fDebug{kFALSE} {
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
              << Form("pT: %5.2f  eta: %3.2f  phi: %3.2f\n",
                      fRecoPt, fRecoEta, fRecoPhi)
              << Form("Matched jet params:\n")
              << Form("pT: %5.2f  eta: %3.2f  phi: %3.2f  flavor: %d  flavorForB: %d\n",
                      fRefPt, fRefEta, fRefPhi, refFlavor(), refFlavorForB() )
              << Form("Generated jet params:\n")
              << Form("pT: %5.2f  eta: %3.2f  phi: %3.2f\n",
                      fGenPt, fGenEta, fGenPhi)
              << Form("---------------------------------\n");
}
