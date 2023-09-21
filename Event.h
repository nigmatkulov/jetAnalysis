#ifndef EVENT_H
#define EVENT_H

// ROOT headers
#include "TObject.h"

//________________
class Event : public TObject {
  public:
    /// @brief Default constructor
    Event();
    Event(const UInt_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
          const Float_t& vz, const Int_t& hiBin, const Float_t& ptHat, 
          const Float_t& w, const Float_t& bit);
    /// @brief destructor
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
    void setWeight(const Float_t& w)     { fWeight = w; }
    /// @brief Set jet trigger bit
    void seetJetTriggerBit(const Int_t& bit) { fJetTriggerBit = bit; }
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
    /// @brief Return centrality bin 
    Int_t hiBin() const       { return (Int_t)fHiBin; }
    /// @brief Return ptHat 
    Float_t ptHat() const     { return fPtHat; }
    /// @brief Return event weight 
    Float_t weight() const    { return fWeight; }
    /// @brief Return jet trigger bit 
    Int_t jetTriggerBit() const { return fJetTriggerBit; }

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
    Short_t     fHiBin;
    /// @brief pthat sclaing
    Float_t   fPtHat;
    /// @brief Event weight scaling
    Float_t   fWeight;
    /// @brief Trigger scheme
    Int_t     fJetTriggerBit;

    ClassDef(Event, 1)
};

#endif // #define EVENT_H