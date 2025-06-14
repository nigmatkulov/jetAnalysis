// JetAnalysis headers
#include "BaseJet.h"

//________________
BaseJet::BaseJet() : TObject(), fPt{0}, fEta{0}, fPhi{0}, 
    fWTAEta{0}, fWTAPhi{0} {
    /* Empty */
}

//________________
bool BaseJet::operator==(const BaseJet& other) const {
    return ( fabs(fPt - other.fPt) < 1e-6 && 
             fabs(fEta - other.fEta) < 1e-6 &&
             fabs(fPhi - other.fPhi) < 1e-6 );
}
