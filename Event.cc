/**
 * @file Event.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Internal event structure of the framework
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "Event.h"

// ROOT headers
#include "TString.h"

// C++ headers
#include <iostream>

//________________
Event::Event() : TObject(), fRunId{0}, fEventId{0}, fLumi{0},
                 fVx{0}, fVy{0}, fVz{0}, fHiBin{-1}, fCentralityWeight{1.}, 
                 fPtHat{-1}, fPtHatWeight{-1}, 
                 fNBadPFJets{0},  fNBadCaloJets{0}, fMult{0},
                 fGenJetsCollectionIsFilled{kFALSE} {
    fPFJetCollection = new PartFlowJetCollection{};
    fCaloJetCollection = new CaloJetCollection{};
    fGenJetCollection = new GenJetCollection{};
    fTrackCollection = new TrackCollection{};
    fGenTrackCollection = new GenTrackCollection{};
    fTrigAndSkim = new TriggerAndSkim{};
}

//________________
Event::Event(const UInt_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
             const Float_t& vx, const Float_t& vy, const Float_t& vz, 
             const Int_t& hiBin, const Float_t& centW, const Float_t& ptHat, 
             const Float_t& w, const Int_t& nBadPFJets, 
             const Int_t& nBadCaloJets, const Int_t& mult) : TObject(),
    fRunId{runId}, fEventId{eventId}, fLumi{lumi}, 
    fVx{vx}, fVy{vy}, fVz{vz},
    fHiBin{(Short_t)hiBin}, fCentralityWeight{centW}, fPtHat{ptHat}, fPtHatWeight{w}, 
    fNBadPFJets{(UChar_t)nBadPFJets}, fNBadCaloJets{(UChar_t)nBadCaloJets},
    fMult{(UShort_t)mult}, fGenJetsCollectionIsFilled{kFALSE} {
    
    // Create new collections 
    fPFJetCollection = new PartFlowJetCollection{};
    fCaloJetCollection = new CaloJetCollection{};
    fGenJetCollection = new GenJetCollection{};
    fTrackCollection = new TrackCollection{};
    fGenTrackCollection = new GenTrackCollection{};
    fTrigAndSkim = new TriggerAndSkim{};
}

//________________
Event::~Event() {
    // Clean collection of particle jets
    for (PartFlowJetIterator iter=fPFJetCollection->begin();
         iter!=fPFJetCollection->end(); iter++) {
        delete *iter;
    }
    // Clean collection of calorimeter jets
    for (CaloJetIterator iter=fCaloJetCollection->begin();
         iter!=fCaloJetCollection->end(); iter++) {
        delete *iter;
    }
    // Clean collection of generated jets
    for (GenJetIterator iter=fGenJetCollection->begin();
         iter!=fGenJetCollection->end(); iter++) {
        delete *iter;
    }
    // Clean track collection
    for (TrackIterator iter=fTrackCollection->begin();
         iter!=fTrackCollection->end(); iter++) {
        delete *iter;
    }
    // Clean MC track collection
    for (GenTrackIterator iter=fGenTrackCollection->begin();
         iter!=fGenTrackCollection->end(); iter++) {
        delete *iter;
    }
    // Clear trigger and skim instance
    if (fTrigAndSkim) delete fTrigAndSkim;
}

//________________
void Event::print() {
    std::cout << Form("-------------------------------------\n")
              << Form("runId: %d  eventId: %llu  lumi: %d  vx: %5.2f  vy: %5.2f  vz: %5.2f\n", fRunId, fEventId, fLumi, fVx, fVy, fVz)
              << Form("hiBin: %d  ptHat: %3.2f  ptHatWeight: %4.2f \n", hiBin(), fPtHat, fPtHatWeight)
              << Form("-------------------------------------\n");
}