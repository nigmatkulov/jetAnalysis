// Jet analysis headers
#include "Event.h"

// ROOT headers
#include "TString.h"

// C++ headers
#include <iostream>

//________________
Event::Event() : TObject(), fRunId{0}, fEventId{0}, fLumi{0},
                 fVz{0}, fHiBin{-1}, fPtHat{-1}, fWeight{-1}, 
                 fJetTriggerBit{-1} {
    fPFJetCollection = new PartFlowJetCollection();
    fCaloJetCollection = new CaloJetCollection();
    fTrackCollection = new TrackCollection();
    fGenTrackCollection = new GenTrackCollection();
}

//________________
Event::Event(const UInt_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
             const Float_t& vz, const Int_t& hiBin, const Float_t& ptHat, 
             const Float_t& w, const Int_t& bit) : TObject(),
    fRunId{runId}, fEventId{eventId}, fLumi{lumi}, fVz{vz},
    fHiBin{(Short_t)hiBin}, fPtHat{ptHat}, fWeight{w}, fJetTriggerBit{bit} {
    
    // Create new collections 
    fPFJetCollection = new PartFlowJetCollection();
    fCaloJetCollection = new CaloJetCollection();
    fTrackCollection = new TrackCollection();
    fGenTrackCollection = new GenTrackCollection();
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
}

//________________
void Event::print() {
    std::cout << Form("-------------------------------------\n")
              << Form("runId: %d  eventId: %llu  lumi: %d  vz: %5.2f\n", fRunId, fEventId, fLumi, fVz)
              << Form("hiBin: %d  ptHat: %3.2f  eventWeight: %4.2f  jetTriggerBit: %d\n",
                      hiBin(), fPtHat, fWeight, fJetTriggerBit)
              << Form("-------------------------------------\n");
}