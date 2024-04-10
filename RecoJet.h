/**
 * @file RecoJet.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Class describes reconstructed jet parameters
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef RecoJet_h
#define RecoJet_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

// JetAnalysis headers
#include "BaseJet.h"

//________________
class RecoJet : public BaseJet {
  public:
    /// @brief Default constructor
    RecoJet();
    /// @brief Destructor
    virtual ~RecoJet() { /* empty */ }

    //
    // Setters
    //

    /// @brief Set reconstructed jet JEC-corrected pt
    void setPtJECCorr(const Float_t& pt) { fPtJECCorr = {pt}; }
    /// @brief Set index of the matched GenJet
    void setGenJetId(const Int_t& id)    { fGenJetId = (Char_t)id; }
    /// @brief Set transverse momentum of tracks in jet
    void setRawPt(const Float_t& pt)     { fRawPt = pt; }
    /// @brief Set momentum of the track with the highest pt in the jet
    void setTrackMaxPt(const Float_t& pt){ fTrackPtMax = pt; }
    /// @brief Print parameters of the given jet
    void print();

    //
    // Getters
    //

    /// @brief Transverse momentum after JEC
    Float_t ptJECCorr() const { return fPtJECCorr; }
    /// @brief Return reconstructed jet parameters
    TVector3 vecJECCorr() const 
    { TVector3 v; v.SetPtEtaPhi(fPtJECCorr, this->eta(), this->phi()); return v; }
    /// @brief Check if reconstructed jet has matched MC jet 
    Bool_t hasMatching() const  
    { return (fGenJetId < 0) ? kFALSE : kTRUE; }
    /// @brief Index of the matched GenJet
    Int_t genJetId() const  { return (Int_t)fGenJetId; }
    /// @brief Return sum of tracks pT in the jet
    Float_t rawPt() const      { return fRawPt; }
    /// @brief Return transverse momentum of the track with highest pT in the jet
    Float_t trackMaxPt() const { return fTrackPtMax; }

  private:

    /// @brief Transverse momentum after JEC
    Float_t fPtJECCorr;
    /// @brief Index of the matched Monte Carlo jet (-99 if not matched)
    Char_t   fGenJetId;
    /// @brief Raw pT of tracks in the jet
    Float_t fRawPt;
    /// @brief Track in the jet with the highest pT
    Float_t fTrackPtMax;
    
    ClassDef(RecoJet, 2)
};

#endif // #define RecoJet_h
