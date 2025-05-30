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
    void setGenJetId(const Int_t& id)    { fGenJetId = (Short_t)id; }
    /// @brief Set transverse momentum of tracks in jet
    void setRawPt(const Float_t& pt)     { fRawPt = pt; }
    /// @brief Set momentum of the track with the highest pt in the jet
    void setTrackMaxPt(const Float_t& pt){ fTrackPtMax = pt; }
    /// @brief Set neutral hadron fraction
    void setJtPfNHF(const Float_t& x) { fJtPfNHF = x; }
    /// @brief Set neutral ElectroMagnetic (EM) fraction
    void setJtPfNEF(const Float_t& x) { fJtPfNEF = x; }
    /// @brief Set charged hadron fraction
    void setJtPfCHF(const Float_t& x) { fJtPfCHF = x; }
    /// @brief Set muon fraction
    void setJtPfMUF(const Float_t& x) { fJtPfMUF = x; }
    /// @brief Set charged ElectroMagnetic fraction
    void setJtPfCEF(const Float_t& x) { fJtPfCEF = x; }
    /// @brief Set charged hadron multiplicity
    void setJtPfCHM(const Int_t& x) { fJtPfCHM = (UChar_t)x; }
    /// @brief Set charged EM multiplicity
    void setJtPfCEM(const Int_t& x) { fJtPfCEM = (UChar_t)x; }
    /// @brief Set neutral hadron multiplicity
    void setJtPfNHM(const Int_t& x) { fJtPfNHM = (UChar_t)x; }
    /// @brief Set neutral EM multiplicity
    void setJtPfNEM(const Int_t& x) { fJtPfNEM = (UChar_t)x; }
    /// @brief Set muon multiplicity
    void setJtPfMUM(const Int_t& x) { fJtPfMUM = (UChar_t)x; }


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

    /// @brief Return neutral hadron fraction
    Float_t jtPfNHF() const { return fJtPfNHF; }
    /// @brief Return neutral ElectroMagnetic (EM) fraction
    Float_t jtPfNEF() const { return fJtPfNEF; }
    /// @brief Return charged hadron fraction
    Float_t jtPfCHF() const { return fJtPfCHF; }
    /// @brief Return muon fraction
    Float_t jtPfMUF() const { return fJtPfMUF; }
    /// @brief Return charged ElectroMagnetic fraction
    Float_t jtPfCEF() const { return fJtPfCEF; }
    /// @brief Return charged hadron multiplicity
    Int_t jtPfCHM() const { return (Int_t)fJtPfCHM; }
    /// @brief Return charged EM multiplicity
    Int_t jtPfCEM() const { return (Int_t)fJtPfCEM; }
    /// @brief Return neutral hadron multiplicity
    Int_t jtPfNHM() const { return (Int_t)fJtPfNHM; }
    /// @brief Return neutral EM multiplicity
    Int_t jtPfNEM() const { return (Int_t)fJtPfNEM; }
    /// @brief Return muon multiplicity
    Int_t jtPfMUM() const { return (Int_t)fJtPfMUM; }

  private:

    /// @brief Transverse momentum after JEC
    Float_t fPtJECCorr;
    /// @brief Index of the matched Monte Carlo jet (-99 if not matched)
    Short_t   fGenJetId;
    /// @brief Raw pT of tracks in the jet
    Float_t fRawPt;
    /// @brief Track in the jet with the highest pT
    Float_t fTrackPtMax;
    /// @brief Neutral hadron fraction
    Float_t fJtPfNHF;
    /// @brief Neutral ElectroMagnetic (EM) fraction
    Float_t fJtPfNEF;
    /// @brief Charged hadron fraction
    Float_t fJtPfCHF;
    /// @brief Muon fraction
    Float_t fJtPfMUF;
    /// @brief Charged ElectroMagnetic fraction
    Float_t fJtPfCEF;
    /// @brief Charged hadron multiplicity
    UChar_t fJtPfCHM;
    /// @brief Charged EM multiplicity
    UChar_t fJtPfCEM;
    /// @brief Neutral hadron multiplicity
    UChar_t fJtPfNHM;
    /// @brief Neutral EM multiplicity
    UChar_t fJtPfNEM;
    /// @brief Muon multiplicity
    UChar_t fJtPfMUM;
    
    ClassDef(RecoJet, 4)
};

#endif // #define RecoJet_h
