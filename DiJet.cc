// JetAnalysis headers
#include "DiJet.h"

// ROOT headers
#include "TString.h"

// C++ headers
#include <iostream>


//________________
DiJet::DiJet(float leadJetPt, float leadJetEtaLab, float leadJetEtaCM, float leadJetPhi,
             float subLeadJetPt, float subLeadJetEtaLab, float subLeadJetEtaCM, float subLeadJetPhi) : 
    fLeadJetPt(leadJetPt), fLeadJetEtaLab(leadJetEtaLab), 
    fLeadJetEtaCM(leadJetEtaCM), fLeadJetPhi(leadJetPhi),
    fSubLeadJetPt(subLeadJetPt), fSubLeadJetEtaLab(subLeadJetEtaLab), 
    fSubLeadJetEtaCM(subLeadJetEtaCM), fSubLeadJetPhi(subLeadJetPhi) {
    // Empty constructor body
}

//________________
float DiJet::deltaPhi(const float& phi1, const float &phi2) {
    float dphi = phi1 - phi2;
    if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
    return dphi;
}

//________________
void DiJet::print() {
    std::cout << "DiJet parameters:" << std::endl;
    std::cout << Form("--> Leading jet: pTave: %.2f etaLab: %.2f etaCM: %.2f phi: %.2f",
                 fLeadJetPt, fLeadJetEtaLab, fLeadJetEtaCM, fLeadJetPhi) << std::endl;
    std::cout << Form("--> Subleading jet: pTave: %.2f etaLab: %.2f etaCM: %.2f phi: %.2f",
                 fSubLeadJetPt, fSubLeadJetEtaLab, fSubLeadJetEtaCM, fSubLeadJetPhi) << std::endl;
    std::cout << Form("--> pTave: %.2f etaLab: %.2f etaCM: %.2f dEtaCM: %.2f dphi: %.2f",
                 this->ptAve(), this->etaLab(), this->etaCM(), this->dEta(), this->dPhi()) << std::endl;
}
