// JetAnalysis headers
#include "BaseJet.h"

ClassImp(BaseJet)

//________________
BaseJet::BaseJet() : TObject(), fPt{0}, fEta{0}, fPhi{0}, 
    fWTAEta{0}, fWTAPhi{0} {
    /* Empty */
}