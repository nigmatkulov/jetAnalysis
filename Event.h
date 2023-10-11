#ifndef EVENT_H
#define EVENT_H

// ROOT headers
#include "TObject.h"
#include "Collections.h"

//________________
class Event : public TObject {
  public:
    /// @brief Default constructor
    Event();
    /// @brief Parametrized constructor
    Event(const UInt_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
          const Float_t& vz, const Int_t& hiBin, const Float_t& ptHat, 
          const Float_t& w, const Int_t& bit);
    /// @brief Destructor
    virtual ~Event();

    //
    // Setters
    //

    /// @brief Set run index
    void setRunId(const UInt_t& id)      { fRunId = id; }
    /// @brief Set event index
    void setEventId(const ULong64_t& id) { fEventId = id; }
    /// @brief Set luminosity
    void setLumi(const UInt_t& lumi)     { fLumi = lumi; }
    /// @brief Set vertex z position
    void setVz(const Float_t& vz)        { fVz = vz; }
    /// @brief Set centrality bin (0-200)
    void setHiBin(const Int_t& hiBin)    { fHiBin = hiBin; }
    /// @brief Set ptHat
    void setPtHat(const Float_t& ptHat)  { fPtHat = ptHat; }
    /// @brief Set event weight
    void setPtHatWeight(const Float_t& w) { fPtHatWeight = w; }
    /// @brief Set jet trigger bit
    void setJetTriggerBit(const Int_t& bit) { fJetTriggerBit = bit; }
    /// @brief  Print event information
    void print();

    //
    // Getters
    //

    /// @brief Return run index 
    UInt_t runId() const      { return fRunId; }
    /// @brief Return event index 
    ULong64_t eventId() const { return fEventId; }
    /// @brief Return luminosity 
    UInt_t lumi() const       { return fLumi; }
    /// @brief Return vertex z 
    Float_t vz() const        { return fVz; }
    /// @brief Return hiBin bin 
    Int_t hiBin() const       { return (Int_t)fHiBin; }
    /// @brief Return centrality bin
    Double_t centrality() const  { return 100. - Double_t(200 - fHiBin) * 0.5; }
    /// @brief Return ptHat 
    Float_t ptHat() const     { return fPtHat; }
    /// @brief Return event weight 
    Float_t ptHatWeight() const    { return fPtHatWeight; }
    /// @brief Return jet trigger bit 
    Int_t jetTriggerBit() const { return fJetTriggerBit; }

    /// @brief Return pointer to a collection of tracks
    TrackCollection *trackCollection() const { return fTrackCollection; }
    /// @brief Return pointer to a collection of MC tracks 
    GenTrackCollection *genTrackCollection() const { return fGenTrackCollection; }
    /// @brief Return pointer to a collection of particle flow jets 
    PartFlowJetCollection *pfJetCollection() const { return fPFJetCollection; }
    /// @brief Return pointer to a collection of calorimeter jets 
    CaloJetCollection *caloJetCollection() const { return fCaloJetCollection; }

  private:
    /// @brief Run index
    UInt_t    fRunId;
    /// @brief Event index
    ULong64_t fEventId;
    /// @brief Luminosity value
    UInt_t    fLumi;
    /// @brief Vertex z position
    Float_t   fVz;
    /// @brief Centrality bin
    Short_t   fHiBin;
    /// @brief pthat sclaing
    Float_t   fPtHat;
    /// @brief Event weight scaling
    Float_t   fPtHatWeight;

    /// @brief Trigger scheme
    Int_t     fJetTriggerBit;

    /// @brief Particle flow jet collection
    PartFlowJetCollection *fPFJetCollection;
    /// @brief Calorimeter jet collection
    CaloJetCollection *fCaloJetCollection;
    /// @brief Track collection
    TrackCollection *fTrackCollection;
    /// @brief MC track collection
    GenTrackCollection *fGenTrackCollection;

    ClassDef(Event, 1)
};

#endif // #define EVENT_H