// JetAnalysis headers
#include "BaseJet.h"

// C++ headers
#include <cmath>

//________________
BaseJet::BaseJet() : TObject(), fPt{0}, fEta{0}, fPhi{0}, 
    fWTAEta{0}, fWTAPhi{0} {
    /* Empty */
}

//________________
BaseJet& BaseJet::operator=(const BaseJet& other) {
    // Assignment operator
    if (this != &other) {
        fPt = other.fPt;
        fEta = other.fEta;
        fPhi = other.fPhi;
        fWTAEta = other.fWTAEta;
        fWTAPhi = other.fWTAPhi;
    }
    return *this;
}

//________________
bool BaseJet::operator==(const BaseJet& other) const {
    return ( fabs(fPt - other.fPt) < 1e-6 && 
             fabs(fEta - other.fEta) < 1e-6 &&
             fabs(fPhi - other.fPhi) < 1e-6 );
}
