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

    /// @brief Set reconstructed jet parameters
    void setRecoJet(const Float_t& pt = -999.f, const Float_t& eta = -999.f, const Float_t& phi = -999.f) 
    { fRecoPt = pt; fRecoEta = eta; fRecoPhi = phi; }
    /// @brief Set parameters of generated jet parameters that matched reconstructed one
    void setRefJet(const Float_t& pt = -999.f, const Float_t& eta = -999.f, const Float_t& phi = -999.f, const Int_t& flav = -999, const Int_t& flavB = -99) 
    { fRefPt = pt; fRefEta = eta; fRefPhi = phi; fRefFlavor = flav; fRefFlavorForB = flavB; }
    void setGenJet(const Float_t& pt = -999.f, const Float_t& eta = -999.f, const Float_t& phi = -999.f) 
    { fGenPt = pt; fGenEta = eta; fGenPt = phi; }
    /// @brief Set flavor for generated jet
    void setRefFlavor(const Int_t& flav)
    { fRefFlavor = (Short_t)flav; }
    /// @brief Set flavor for B
    void setRefFlavorForB(const Int_t& flav) 
    { fRefFlavorForB = (Char_t)flav; }
    /// @brief Set debug
    void setDebug(Bool_t debug)  { fDebug = debug; }
    /// @brief Print parameters of the given jet
    void print();

    //
    // Getters
    //

    /// @brief Return reconstructed jet parameters
    TVector3 recoJetVec() const 
    { TVector3 v; v.SetPtEtaPhi(fRecoPt, fRecoEta, fRecoPhi); return v; }
    /// Check if reconstructed jet has matched generated jet
    Bool_t hasMatching() const  
    { return (fRecoPt < 998.f) ? kFALSE : kTRUE; }
    /// @brief Return vector of generated jet that matched reconstructed one
    TVector3 refJetVec() const  
    { TVector3 v; v.SetPtEtaPhi(fRefPt, fRefEta, fRefPhi); return v; }
    /// @brief Return vector of generated jet
    TVector3 genJetVec() const
    { TVector3 v; v.SetPtEtaPhi(fGenPt, fGenEta, fGenPhi); return v; }

    /// @brief Return flavor of the generated jet that was matched with the reconstructed one
    Int_t refFlavor() const      { return (Int_t)fRefFlavor; }
    /// @brief Return flavor from B
    Int_t refFlavorForB() const { return (Int_t)fRefFlavorForB; }

  private:

    /// @brief Jet transverse momentum (-999. for non-existing reco jet)
    Float_t fRecoPt;
    /// @brief Jet pseudorapidity (-999. for non-existing reco jet)
    Float_t fRecoEta;
    /// @brief Jet azimuthal angle (-999. for non-existing reco jet)
    Float_t fRecoPhi;

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

    /// @brief enerated pT (-999 for non-MC jet)
    Float_t fGenPt;
    /// @brief Generated pT (-999 for non-MC jet)
    Float_t fGenEta;
    /// @brief Generated pT (-999 for non-MC jet)
    Float_t fGenPhi;

    /// Debug flag (default kFALSE)
    Bool_t fDebug;
    
    ClassDef(Jet, 1)
};

#endif // #define JET_H
