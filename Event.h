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
    Event(const Int_t& runId, const ULong64_t& eventId, const UInt_t& lumi, 
          const Float_t& vx, const Float_t& vy, const Float_t& vz, 
          const Int_t& hiBin, const Float_t& centW, const Float_t& ptHat, 
          const Float_t& w, const Int_t& nBadRecoJets, const Int_t& mult);
    /// @brief Destructor
    virtual ~Event();

    //
    // Setters
    //

    /// @brief Set run index
    void setRunId(const Int_t& id)       { fRunId = id; }
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
    /// @brief Set centrality weight
    void setCentralityWeight(const Float_t& w) { fCentralityWeight = w; }
    /// @brief Set ptHat
    void setPtHat(const Float_t& ptHat)  { fPtHat = ptHat; }
    /// @brief Set event weight
    void setPtHatWeight(const Float_t& w) { fPtHatWeight = w; }
    /// @brief Set number of particle flow jets with pT > pThat
    void setNumberOfOverscaledRecoJets(const Int_t& n) { fNBadRecoJets = (UChar_t)n; }
    /// @brief Set reference multiplicity (CMS way)
    void setMultiplicity(const Int_t& mult) { fMult = (UShort_t)mult; }
    /// @brief Set the flag that generated jet collection is filled to true
    void setGenJetCollectionIsFilled()      { fGenJetsCollectionIsFilled = kTRUE; }
    /// @brief  Print event information
    void print();

    //
    // Getters
    //

    /// @brief Return run index 
    Int_t runId() const       { return fRunId; }
    /// @brief Return event index 
    ULong64_t eventId() const { return fEventId; }
    /// @brief Return luminosity 
    UInt_t lumi() const       { return fLumi; }
    /// @brief Return vertex z 
    Float_t vz() const        { return fVz; }
    /// @brief Return hiBin bin 
    Int_t hiBin() const       { return (Int_t)fHiBin; }
    /// @brief Return centrality bin
    Double_t centrality() const  { return (fHiBin < 0) ? -5 : 100. - Double_t(200 - fHiBin) * 0.5; }
    /// @brief Return centrality weight 
    Double_t centralityWeight() const { return (Double_t)fCentralityWeight; }
    /// @brief Return ptHat 
    Float_t ptHat() const     { return fPtHat; }
    /// @brief Return event weight 
    Float_t ptHatWeight() const    { return fPtHatWeight; }
    /// @brief Return number of particle flow jets in even with pT > pThat 
    Int_t numberOfOverscaledRecoJets() const { return (Int_t)fNBadRecoJets; }
    /// @brief Return reference multiplicity (CMS way) 
    Int_t multiplicity() const { return (Int_t)fMult; }
    /// @brief Return if generated jet collection filled 
    Bool_t isGenJetCollectionFilled() const { return fGenJetsCollectionIsFilled; }
    /// @brief Return number of particle flow jets 
    UInt_t numberOfRecoJets() const   { return this->recoJetCollection()->size(); }
    /// @brief Return number of generated jets 
    UInt_t numberOfGenJets() const  { return this->genJetCollection()->size(); }

    /// @brief Return pointer to a collection of tracks
    TrackCollection *trackCollection() const { return fTrackCollection; }
    /// @brief Return pointer to a collection of MC tracks 
    GenTrackCollection *genTrackCollection() const { return fGenTrackCollection; }
    /// @brief Return pointer to a collection of particle flow jets 
    RecoJetCollection *recoJetCollection() const { return fRecoJetCollection; }
    /// @brief Return pointer to a collection of generated jets 
    GenJetCollection *genJetCollection() const { return fGenJetCollection; }
    /// @brief Return pointer to a trigger and skimming information 
    TriggerAndSkim *trigAndSkim() const { return fTrigAndSkim; }

  private:
    /// @brief Run index
    Int_t    fRunId;
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
    /// @brief Centrality weight
    Float_t   fCentralityWeight;
    /// @brief pthat sclaing
    Float_t   fPtHat;
    /// @brief Event weight scaling
    Float_t   fPtHatWeight;
    /// @brief Number of particle flow jets with pT > pThat
    UChar_t   fNBadRecoJets;
    /// @brief Reference charged track multiplicity (CMS way)
    UShort_t  fMult;
    /// @brief Check if collection of generated jets is filled
    Bool_t    fGenJetsCollectionIsFilled;

    /// @brief Particle flow jet collection
    RecoJetCollection *fRecoJetCollection;
    /// @brief Generated jet collection
    GenJetCollection *fGenJetCollection;

    /// @brief Track collection
    TrackCollection *fTrackCollection;
    /// @brief MC track collection
    GenTrackCollection *fGenTrackCollection;


    /// @brief Trigger and skimming information
    TriggerAndSkim *fTrigAndSkim;

    ClassDef(Event, 1)
};

#endif // #define Event_h