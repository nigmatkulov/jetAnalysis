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
    // Empty 
}

//________________
Event::Event(const UInt_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
              const Float_t& vz, const Int_t& hiBin, const Float_t& ptHat, 
              const Float_t& w, const Float_t& bit) : Event() {
    fRunId = runId;
    fEventId = eventId;
    fLumi = lumi;
    fVz = vz;
    fHiBin = (Short_t)hiBin;
    fPtHat = ptHat;
    fWeight = w;
    fJetTriggerBit = bit;
}

//________________
Event::~Event() {
    // Empty
}

//________________
void Event::print() {
    std::cout << Form("-------------------------------------\n")
              << Form("runId: %d  eventId: %d  lumi: %d  vz: %5.2f\n", fRunId, fEventId, fLumi, fVz)
              << Form("hiBin: %d  ptHat: %3.2f  eventWeight: %4.2f  jetTriggerBit: %d\n",
                      hiBin(), fPtHat, fWeight, fJetTriggerBit)
              << Form("-------------------------------------\n");

}