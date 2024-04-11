// JetAnalysis headers
#include "BaseJet.h"
#include "TString.h"

// C++ headers
#include <iostream>

ClassImp(BaseJet)

//________________
BaseJet::BaseJet() : TObject(), fPt{0}, fEta{0}, fPhi{0}, 
    fWTAEta{0}, fWTAPhi{0} {
    /* Empty */
}

//________________
void BaseJet::print() {
    std::cout << Form("pT: %5.2f eta: %3.2f phi: %3.2f\n", fPt, fEta, fPhi);
}