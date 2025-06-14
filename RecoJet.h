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

    /// @brief Equality operator to compare two RecoJet objects
    bool operator==(const RecoJet& other) const;

    //
    // Setters
    //

    /// @brief Set reconstructed jet JEC-corrected pt
    void setPtJECCorr(const float& pt) { fPtJECCorr = pt; }
    /// @brief Set index of the matched GenJet
    void setGenJetId(const int& id)    { fGenJetId = (Short_t)id; }
    /// @brief Set transverse momentum of tracks in jet
    void setRawPt(const float& pt)     { this->setPt(pt); }
    /// @brief Set momentum of the track with the highest pt in the jet
    void setTrackMaxPt(const float& pt){ fTrackPtMax = pt; }
    /// @brief Set neutral hadron fraction
    void setJtPfNHF(const float& x) { fJtPfNHF = x; }
    /// @brief Set neutral ElectroMagnetic (EM) fraction
    void setJtPfNEF(const float& x) { fJtPfNEF = x; }
    /// @brief Set charged hadron fraction
    void setJtPfCHF(const float& x) { fJtPfCHF = x; }
    /// @brief Set muon fraction
    void setJtPfMUF(const float& x) { fJtPfMUF = x; }
    /// @brief Set charged ElectroMagnetic fraction
    void setJtPfCEF(const float& x) { fJtPfCEF = x; }
    /// @brief Set charged hadron multiplicity
    void setJtPfCHM(const int& x) { fJtPfCHM = (UChar_t)x; }
    /// @brief Set charged EM multiplicity
    void setJtPfCEM(const int& x) { fJtPfCEM = (UChar_t)x; }
    /// @brief Set neutral hadron multiplicity
    void setJtPfNHM(const int& x) { fJtPfNHM = (UChar_t)x; }
    /// @brief Set neutral EM multiplicity
    void setJtPfNEM(const int& x) { fJtPfNEM = (UChar_t)x; }
    /// @brief Set muon multiplicity
    void setJtPfMUM(const int& x) { fJtPfMUM = (UChar_t)x; }


    /// @brief Print parameters of the given jet
    void print();

    //
    // Getters
    //

    /// @brief Transverse momentum after JEC
    float ptJECCorr() const { return fPtJECCorr; }
    /// @brief Return reconstructed jet parameters
    TVector3 vecJECCorr() const 
    { TVector3 v; v.SetPtEtaPhi(fPtJECCorr, this->eta(), this->phi()); return v; }
    /// @brief Check if reconstructed jet has matched MC jet 
    bool hasMatching() const  
    { return (fGenJetId < 0) ? kFALSE : kTRUE; }
    /// @brief Index of the matched GenJet
    int genJetId() const  { return (int)fGenJetId; }
    /// @brief Return sum of tracks pT in the jet
    float rawPt() const      { return this->pt(); }
    /// @brief Return transverse momentum of the track with highest pT in the jet
    float trackMaxPt() const { return fTrackPtMax; }

    /// @brief Return neutral hadron fraction
    float jtPfNHF() const { return fJtPfNHF; }
    /// @brief Return neutral ElectroMagnetic (EM) fraction
    float jtPfNEF() const { return fJtPfNEF; }
    /// @brief Return charged hadron fraction
    float jtPfCHF() const { return fJtPfCHF; }
    /// @brief Return muon fraction
    float jtPfMUF() const { return fJtPfMUF; }
    /// @brief Return charged ElectroMagnetic fraction
    float jtPfCEF() const { return fJtPfCEF; }
    /// @brief Return charged hadron multiplicity
    int jtPfCHM() const { return (int)fJtPfCHM; }
    /// @brief Return charged EM multiplicity
    int jtPfCEM() const { return (int)fJtPfCEM; }
    /// @brief Return neutral hadron multiplicity
    int jtPfNHM() const { return (int)fJtPfNHM; }
    /// @brief Return neutral EM multiplicity
    int jtPfNEM() const { return (int)fJtPfNEM; }
    /// @brief Return muon multiplicity
    int jtPfMUM() const { return (int)fJtPfMUM; }

    /// @brief Check if jet passes trackMaxPt/jetRawPt cut in various eta ranges
    bool isGoodTrkMax() const;
    /// @brief Check if jet passes jetId selection criteria for various eta ranges
    bool isGoodJetId(const bool& useLooseJetIdCut = false) const;

  private:

    // Note: rawPt is stored in BaseJet class, so no need to redefine it here

    /// @brief Jet-energy-corrected transverse momentum
    Float_t fPtJECCorr;
    /// @brief Index of the matched Monte Carlo jet (-99 if not matched)
    Short_t   fGenJetId;
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
    
    ClassDef(RecoJet, 5)
};

#endif // #define RecoJet_h
