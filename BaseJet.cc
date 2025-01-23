// JetAnalysis headers
#include "BaseJet.h"

//________________
BaseJet::BaseJet() : TObject(), fPt{0}, fEta{0}, fPhi{0}, 
    fWTAEta{0}, fWTAPhi{0} {
    /* Empty */
}