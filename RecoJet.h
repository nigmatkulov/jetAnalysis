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

  private:

    /// @brief Transverse momentum after JEC
    Float_t fPtJECCorr;
    /// @brief Index of the matched Monte Carlo jet (-99 if not matched)
    Char_t   fGenJetId;
    
    ClassDef(RecoJet, 1)
};

#endif // #define RecoJet_h
