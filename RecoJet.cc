/**
 * @file RecoJet.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Class describes reconstructed jet parameters
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "RecoJet.h"

// ROOT headers
#include "TString.h"
#include "TMath.h"

// C++ headers
#include <iostream>

//________________
RecoJet::RecoJet() : BaseJet{}, fPtJECCorr{0}, fGenJetId{-99}, fTrackPtMax{0},
    fJtPfNHF{0}, fJtPfNEF{0}, fJtPfCHF{0}, fJtPfMUF{0}, fJtPfCEF{0}, 
    fJtPfCHM{0}, fJtPfCEM{0}, fJtPfNHM{0}, fJtPfNEM{0}, fJtPfMUM{0} {
    /* Empty */
}

//________________
void RecoJet::print() {
    std::cout << Form("--  Reconstructed jet info  --\n")
              << Form("pT: %5.2f  eta: %3.2f  phi: %3.2f  pTcorr: %5.2f  WTAeta: %3.2f  WTAphi: %3.2f  GenJetId: %d NHF: %3.2f  NEF: %3.2f  CHF: %3.2f  MUF: %3.2f  CEF: %3.2f  CHM: %d  CEM: %d  NHM: %d  NEM: %d  MUM: %d \n",
                      this->pt(), this->eta(), this->phi(), fPtJECCorr, this->WTAEta(), this->WTAPhi(), this->genJetId(),
                      fJtPfNHF, fJtPfNEF, fJtPfCHF, fJtPfMUF, fJtPfCEF, (Int_t)fJtPfCHM, (Int_t)fJtPfCEM, (Int_t)fJtPfNHM, (Int_t)fJtPfNEM, (Int_t)fJtPfMUM);
}

//________________
bool RecoJet::isGoodTrkMax() const {

    return ( TMath::Abs(this->eta()) < 2.4 && 
             ( ((fTrackPtMax / this->pt()) < 0.01) || ((fTrackPtMax / this->pt()) > 0.98) ) ) ? false : true;
}

//________________
bool RecoJet::isGoodJetId(const bool& useLooseJetIdCut) const {

    bool passJetId = {false};

    int chm = this->jtPfCHM();
    int cem = this->jtPfCEM();
    int mum = this->jtPfMUM();
    int nhm = this->jtPfNHM();
    int nem = this->jtPfNEM();

    float chf = this->jtPfCHF();
    float cef = this->jtPfCEF();
    float nhf = this->jtPfNHF();
    float nef = this->jtPfNEF();
    float muf = this->jtPfMUF();

    int chargedMult = chm + cem + mum;
    int neutralMult = nhm + nem;
    int numberOfConstituents = chargedMult + neutralMult;

    float eta = this->eta();

    float chargedEmFracCut{1.}, neutFracCut{1.};
    if ( !useLooseJetIdCut ) {
        chargedEmFracCut = {0.9};
        neutFracCut = {0.9};
    }
    else {
        chargedEmFracCut = {0.99};
        neutFracCut = {0.99};
    }

    bool passNHF{false};
    bool passNEF{false};
    bool passNumOfConstituents{true};
    bool passMuonFrac{true};
    bool passChargedFrac{true};
    bool passChargedMult{true};
    bool passChargedEmFrac{true};
    bool passNeutralMult{true};
	
    // Check cuts depending on jet pseudorapdity
    if ( TMath::Abs( eta ) <= 2.7 ) {
        
        passNHF = ( nhf < neutFracCut ) ? true : false;
        passNEF = ( nhf < neutFracCut ) ? true : false;
        passNumOfConstituents = ( numberOfConstituents > 1 ) ? true : false;

        if ( !useLooseJetIdCut ) { 
            passMuonFrac = ( muf < 0.8 ) ? true : false; 
        } // if ( !useLooseJetIdCut )

        if( TMath::Abs( eta ) <= 2.4 ) {
            passChargedFrac = ( chf > 0 ) ? true : false;
            passChargedMult = ( chargedMult > 0 ) ? true : false;
            passChargedEmFrac = ( cef < chargedEmFracCut ) ? true : false;
        } // if( TMath::Abs( eta ) <= 2.4 )

    } // if ( TMath::Abs( eta ) <= 2.7 )
    else if ( TMath::Abs( eta ) <= 3.0) {

        passNEF = ( nef > 0.01 ) ? true : false;
        passNHF = ( nhf < 0.98 ) ? true : false;
        passNeutralMult = ( neutralMult > 2 ) ? true : false;

    } // else if ( TMath::Abs( eta ) <= 3.0)
    else  {
        passNEF = ( nef < 0.9 ) ? true : false;
        passNeutralMult = (neutralMult > 10 ) ? true : false; // CAUTION: JET MET says it should be >10
    } // else 

    passJetId = passNHF && passNEF && passNumOfConstituents && passMuonFrac && 
                passChargedFrac && passChargedMult && passChargedEmFrac && passNeutralMult;

    // if ( fVerbose ) {
    //     std::cout << "JetId selection results: " << ( (passJetId) ? "[passed]" : "[failed]" ) << " Reasons ";
    //     std::cout << Form("passNHF: %d \tpassNEF: %d \tpassNumConst: %d \tpassMuonFrac: %d \tpassChFrac: %d \tpassChMult: %d \tpassChEmFrac: %d \tpassNeutMult: %d\n", 
    //                       passNHF, passNEF, passNumOfConstituents, passMuonFrac, passChargedFrac, 
    //                       passChargedMult , passChargedEmFrac , passNeutralMult);
    // }
		
	return passJetId;
}
