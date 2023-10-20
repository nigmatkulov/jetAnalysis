/**
 * @file Event.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Internal event structure of the framework
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef Event_h
#define Event_h

// ROOT headers
#include "TObject.h"
#include "Collections.h"
#include "TriggerAndSkim.h"

//________________
class Event : public TObject {
  public:
    /// @brief Default constructor
    Event();
    /// @brief Parametrized constructor
    Event(const UInt_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
          const Float_t& vx, const Float_t& vy, const Float_t& vz, 
          const Int_t& hiBin, const Float_t& ptHat, 
          const Float_t& w, const Int_t& nBadPFJets, 
          const Int_t& nBadCaloJets, const Int_t& mult);
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
    /// @brief Set vertex x position
    void setVx(const Float_t& vx)        { fVx = vx; }
    /// @brief Set vertex y position
    void setVy(const Float_t& vy)        { fVy = vy; }
    /// @brief Set vertex z position
    void setVz(const Float_t& vz)        { fVz = vz; }
    /// @brief Set vertex x, y, and z coordinates
    void setVertex(const Float_t& x, const Float_t& y, const Float_t &z) { fVx = x; fVy = y; fVz = z; }
    /// @brief Set centrality bin (0-200)
    void setHiBin(const Int_t& hiBin)    { fHiBin = hiBin; }
    /// @brief Set ptHat
    void setPtHat(const Float_t& ptHat)  { fPtHat = ptHat; }
    /// @brief Set event weight
    void setPtHatWeight(const Float_t& w) { fPtHatWeight = w; }
    /// @brief Set number of particle flow jets with pT > pThat
    void setNumberOfOverscaledPFJets(const Int_t& n) { fNBadPFJets = (UChar_t)n; }
    /// @brief Set number of calorimeter jets with pT > pThat
    void setNumberOfOverscaledCaloJets(const Int_t& n) { fNBadCaloJets = (UChar_t)n; }
    /// @brief Set reference multiplicity (CMS way)
    void setMultiplicity(const Int_t& mult) { fMult = (UShort_t)mult; }
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
    /// @brief Return number of particle flow jets in even with pT > pThat 
    Int_t numberOfOverscaledPFJets() const { return (Int_t)fNBadPFJets; }
    /// @brief Return number of calorimeter jets in even with pT > pThat 
    Int_t numberOfOverscaledCaloJets() const { return (Int_t)fNBadCaloJets; }
    /// @brief Return reference multiplicity (CMS way) 
    Int_t multiplicity() const { return (Int_t)fMult; }

    /// @brief Return pointer to a collection of tracks
    TrackCollection *trackCollection() const { return fTrackCollection; }
    /// @brief Return pointer to a collection of MC tracks 
    GenTrackCollection *genTrackCollection() const { return fGenTrackCollection; }
    /// @brief Return pointer to a collection of particle flow jets 
    PartFlowJetCollection *pfJetCollection() const { return fPFJetCollection; }
    /// @brief Return pointer to a collection of calorimeter jets 
    CaloJetCollection *caloJetCollection() const { return fCaloJetCollection; }
    /// @brief Return pointer to a trigger and skimming information 
    TriggerAndSkim *trigAndSkim() const { return fTrigAndSkim; }

  private:
    /// @brief Run index
    UInt_t    fRunId;
    /// @brief Event index
    ULong64_t fEventId;
    /// @brief Luminosity value
    UInt_t    fLumi;
    /// @brief Vertex x position
    Float_t   fVx;
    /// @brief Vertex y position
    Float_t   fVy;
    /// @brief Vertex z position
    Float_t   fVz;
    /// @brief Centrality bin
    Short_t   fHiBin;
    /// @brief pthat sclaing
    Float_t   fPtHat;
    /// @brief Event weight scaling
    Float_t   fPtHatWeight;
    /// @brief Number of particle flow jets with pT > pThat
    UChar_t   fNBadPFJets;
    /// @brief Number of calorimeter jets with pT > pThat
    UChar_t   fNBadCaloJets;
    /// @brief Reference charged track multiplicity (CMS way)
    UShort_t  fMult;

    /// @brief Particle flow jet collection
    PartFlowJetCollection *fPFJetCollection;
    /// @brief Calorimeter jet collection
    CaloJetCollection *fCaloJetCollection;
    /// @brief Track collection
    TrackCollection *fTrackCollection;
    /// @brief MC track collection
    GenTrackCollection *fGenTrackCollection;

    /// @brief Trigger and skimming information
    TriggerAndSkim *fTrigAndSkim;

    ClassDef(Event, 1)
};

#endif // #define Event_h