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

// C++ headers
#include <iostream>

ClassImp(RecoJet)

//________________
RecoJet::RecoJet() : BaseJet{}, fPtJECCorr{0}, fGenJetId{-99}, 
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
