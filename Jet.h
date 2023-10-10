#ifndef JET_H
#define JET_H

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//________________
class Jet : public TObject {
  public:
    /// @brief Default constructor
    Jet();
    /// @brief Destructor
    virtual ~Jet();

    //
    // Setters
    //

    /// @brief Set reconstructed jet pt 
    void setRecoJetPt(const Float_t& pt)     { fRecoPt = pt; }
    /// @brief Set reconstructed jet eta
    void setRecoJetEta(const Float_t& eta)   { fRecoEta = eta; }
    /// @brief Set reconstructed jet phi
    void setRecoJetPhi(const Float_t& phi)   { fRecoPhi = phi; }
    /// @brief Set reconstructed jet JES-corrected pt
    void setRecoJetPtJESCorr(const Float_t& pt) { fRecoPtJESCorr = pt; }
    /// @brief Set WTA eta axis
    void setRecoJetWTAeta(const Float_t& eta) { fRecoWTAeta = eta; }
    /// @brief Set WTA phi axis
    void setRecoJetWTAphi(const Float_t& phi) {fRecoWTAphi = phi; }
    /// @brief Set jet weight (for MC). Jet pT-smearing is not applied
    void setRecoJetPtWeight(const Float_t& w) { fRecoJetPtWeight = w; }
    /// @brief Set jet pT smearing weight (for MC)
    void setRecoJetPtSmearingWeight(const Float_t& w) { fRecoJetPtSmearingWeight = w; }

    /// @brief Set jet-matched generated jet transverse momentum
    void setRefJetPt(const Float_t& pt)   { fRefPt = pt; }
    /// @brief Set jet-matched generated jet eta
    void setRefJetEta(const Float_t& eta) { fRefEta = eta; }
    /// @brief Set jet-matched generated jet phi
    void setRefJetPhi(const Float_t& phi) { fRefPhi = phi; }
    /// @brief Set flavor for generated jet
    void setRefFlavor(const Int_t& flav)
    { fRefFlavor = (Short_t)flav; }
    /// @brief Set flavor for B
    void setRefFlavorForB(const Int_t& flav) 
    { fRefFlavorForB = (Char_t)flav; }
    /// @brief Set jet-matched pT weight
    void setRefJetPtWeight(const Float_t& w) { fRefPtWeight = w; }
    /// @brief Print parameters of the given jet
    void print();

    //
    // Getters
    //

    /// @brief Reconstructed jet transverse momentum 
    Float_t recoJetPt() const  { return fRecoPt; }
    /// @brief Reconstructed jet pseudorapidity 
    Float_t recoJetEta() const { return fRecoEta; }
    /// @brief Reconstructed jet phi 
    Float_t recoJetPhi() const { return fRecoPhi; }
    /// @brief Reconstructed jet JES-corrected pt 
    Float_t recoJetPtJESCorr() const { return fRecoPtJESCorr; }
    /// @brief Reconstructed jet WTA eta axis 
    Float_t recoJetWTAeta() const { return fRecoWTAeta; }
    /// @brief Reconstructed jet WTA phi axis
    Float_t recoJetWTAphi() const { return fRecoWTAphi; }
    /// @brief Return reconstructed jet parameters
    TVector3 recoJetVec() const 
    { TVector3 v; v.SetPtEtaPhi(fRecoPt, fRecoEta, fRecoPhi); return v; }
    /// @brief Return reconstructed jet parameters
    TVector3 recoJetVecJESCorr() const 
    { TVector3 v; v.SetPtEtaPhi(fRecoPtJESCorr, fRecoEta, fRecoPhi); return v; }
    /// @brief Return jet pT weight
    Float_t recoJetPtWeight() const { return fRecoJetPtWeight; }
    /// @brief Return jet pT smearing weight 
    Float_t recoJetPtSmearingWeight() const { return fRecoJetPtSmearingWeight; }
    /// @brief Return jet weight (pTweight x pTsmearWeight) 
    Float_t recoJetWeight() const { return recoJetPtWeight() * recoJetPtSmearingWeight(); }
    /// Check if reconstructed jet has matched generated jet
    Bool_t hasMatching() const  
    { return (fRefPt < -998.f) ? kFALSE : kTRUE; }

    /// @brief Generated jet pt that matched reconstructed one
    Float_t refJetPt() const  { return fRefPt; }
    /// @brief Generated jet eta that matched reconstructed one
    Float_t refJetEta() const { return fRefEta; }
    /// @brief Generated jet phi that matched reconstructed one
    Float_t refJetPhi() const { return fRefPhi; }
    /// @brief Return vector of generated jet that matched reconstructed one
    TVector3 refJetVec() const  
    { TVector3 v; v.SetPtEtaPhi(fRefPt, fRefEta, fRefPhi); return v; }
    /// @brief Return flavor of the generated jet that was matched with the reconstructed one
    Int_t refFlavor() const      { return (Int_t)fRefFlavor; }
    /// @brief Return flavor from B
    Int_t refFlavorForB() const { return (Int_t)fRefFlavorForB; }
    /// @brief Return MC jet pT weight
    Float_t refJetPtWeight() const { return fRefPtWeight; }

  private:

    /// @brief Jet transverse momentum (-999. for non-existing reco jet)
    Float_t fRecoPt;
    /// @brief Jet pseudorapidity (-999. for non-existing reco jet)
    Float_t fRecoEta;
    /// @brief Jet azimuthal angle (-999. for non-existing reco jet)
    Float_t fRecoPhi;
    /// @brief Jet transverse momentum after energy correction
    Float_t fRecoPtJESCorr;
    /// @brief Jet WTA eta axis
    Float_t fRecoWTAeta;
    /// @brief Jet WTA phi axis
    Float_t fRecoWTAphi;
    /// @brief Jet weight (for MC). Without smearing weight
    Double_t fRecoJetPtWeight;
    /// @brief Jet smearing pT weight
    Double_t fRecoJetPtSmearingWeight;

    /// @brief Reco-matched MC pT (-999. for non-existing MC jet)
    Float_t fRefPt;
    /// @brief Reco-matched MC eta (-999. for non-existing MC jet)
    Float_t fRefEta;
    /// @brief Reco-matched MC phi (-999. for non-existing MC jet)
    Float_t fRefPhi;
    /// @brief Reco-matched MC flavor
    Short_t fRefFlavor;
    /// @brief Reco-matched MC flavor (-5 - antib, -4 - antic, -3 - antis, -2 - antiu, -1 - antid, 0 - unknown, 1 - d, 2 - u, 3 - s, 4 - c, 5 - b, 21 - gluon, -99 - for non-existing MC jet)
    Char_t  fRefFlavorForB;
    /// @brief Reco-matched MC weight
    Float_t fRefPtWeight;
    
    ClassDef(Jet, 1)
};

#endif // #define JET_H
