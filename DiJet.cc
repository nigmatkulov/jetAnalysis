// JetAnalysis headers
#include "DiJet.h"

// ROOT headers
#include "TString.h"
#include "TMath.h"

// C++ headers
#include <iostream>
#include <cmath>

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
DiJet& DiJet::operator=(const DiJet& other) {
    // Assignment operator
    if (this != &other) {
        fLeadJetPt = other.fLeadJetPt;
        fLeadJetEtaLab = other.fLeadJetEtaLab;
        fLeadJetEtaCM = other.fLeadJetEtaCM;
        fLeadJetPhi = other.fLeadJetPhi;
        fSubLeadJetPt = other.fSubLeadJetPt;
        fSubLeadJetEtaLab = other.fSubLeadJetEtaLab;
        fSubLeadJetEtaCM = other.fSubLeadJetEtaCM;
        fSubLeadJetPhi = other.fSubLeadJetPhi;
    }
    return *this;
}

//________________
bool DiJet::operator==(const DiJet& other) const {
    // Compare all members for equality
    return (fLeadJetPt == other.fLeadJetPt &&
            fLeadJetEtaLab == other.fLeadJetEtaLab &&
            fLeadJetEtaCM == other.fLeadJetEtaCM &&
            fLeadJetPhi == other.fLeadJetPhi &&
            fSubLeadJetPt == other.fSubLeadJetPt &&
            fSubLeadJetEtaLab == other.fSubLeadJetEtaLab &&
            fSubLeadJetEtaCM == other.fSubLeadJetEtaCM &&
            fSubLeadJetPhi == other.fSubLeadJetPhi);
}

//________________
float DiJet::deltaPhi(const float& phi1, const float &phi2) {
    float dphi = phi1 - phi2;
    if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
    return dphi;
}

//________________
float wrapTo0to2Pi(const float &angle) {
    float wrapped = std::fmod(angle, TMath::TwoPi());
    return wrapped < 0 ? wrapped + TMath::TwoPi() : wrapped;
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

//________________
void DiJet::cleanParameters() {
    // Reset all parameters to zero
    fLeadJetPt = 0.0f;
    fLeadJetEtaLab = 0.0f;
    fLeadJetEtaCM = 0.0f;
    fLeadJetPhi = 0.0f;

    fSubLeadJetPt = 0.0f;
    fSubLeadJetEtaLab = 0.0f;
    fSubLeadJetEtaCM = 0.0f;
    fSubLeadJetPhi = 0.0f;
}

//________________
float DiJet::phi() const {
    // Calculate the azimuthal angle of the dijet
    float phiLead = fLeadJetPhi;
    float phiSubLead = fSubLeadJetPhi;

    float pxLead = fLeadJetPt * TMath::Cos(phiLead);
    float pyLead = fLeadJetPt * TMath::Sin(phiLead);
    float pxSubLead = fSubLeadJetPt * TMath::Cos(phiSubLead);
    float pySubLead = fSubLeadJetPt * TMath::Sin(phiSubLead);

    float pxDijet = pxLead + pxSubLead;
    float pyDijet = pyLead + pySubLead;

    return TMath::ATan2(pyDijet, pxDijet);
}
