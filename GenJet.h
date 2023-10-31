/**
 * @file GenJet.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Class describes jet from the MС level
 * @version 0.1
 * @date 2023-10-23
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef GenJet_h
#define GenJet_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

// JetAnalysis headers
#include "BaseJet.h"

//________________
class GenJet : public BaseJet {
  public:
    /// @brief Constructor
    GenJet();
    /// @brief Destructor
    virtual ~GenJet() { /* Empty*/ }

    //
    // Setters
    //

    /// @brief Set flavor for generated jet
    void setFlavor(const Int_t& flav)  { fFlavor = {(Short_t)flav}; }
    /// @brief Set flavor for B
    void setFlavorForB(const Int_t& flav) { fFlavorForB = {(Char_t)flav}; }
    /// @brief Set jet-matched pT weight
    void setPtWeight(const Float_t& w)    { fPtWeight = {w}; }
    /// @brief Print parameters of the given jet
    void print();

    //
    // Getters
    //

    /// @brief Flavor (-999 for non-matched to RecoJet)
    Int_t   flavor() const     { return (Int_t)fFlavor; }
    /// @brief Flavor (-5 - antib, -4 - antic, -3 - antis, -2 - antiu, -1 - antid, 0 - unknown, 1 - d, 2 - u, 3 - s, 4 - c, 5 - b, 21 - gluon, -99 for non-matched to RecoJet)
    Int_t   flavorForB() const { return (Int_t)fFlavorForB; }

  private:

    /// @brief Flavor (-990 for non-matched to RecoJet)
    Short_t fFlavor; 
    /// @brief Flavor for B (-5 - antib, -4 - antic, -3 - antis, -2 - antiu, 
    /// -1 - antid, 0 - unknown, 1 - d, 2 - u, 3 - s, 4 - c, 5 - b, 
    /// 21 - gluon, -99 for non-matched to RecoJet)
    Char_t  fFlavorForB;
    /// @brief Transverse momentum weight
    Float_t fPtWeight;

    ClassDef(GenJet, 1)
};

#endif