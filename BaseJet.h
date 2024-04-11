/**
 * @file BaseJet.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Base class for jet description
 * @version 0.1
 * @date 2023-10-23
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef BaseJet_h
#define BaseJet_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"

//________________
class BaseJet : public TObject {
  public:
    /// @brief Constructor
    BaseJet();
    /// @brief Destructor
    virtual ~BaseJet() { /* Empty */ }
    /// @brief Print jet information
    virtual print();

    //
    // Setters
    //

    /// @brief Set jet-matched generated jet transverse momentum
    void setPt(const Float_t& pt)      { fPt = {pt}; }
    /// @brief Set jet-matched generated jet eta
    void setEta(const Float_t& eta)    { fEta = {eta}; }
    /// @brief Set jet-matched generated jet phi
    void setPhi(const Float_t& phi)    { fPhi = {phi}; }
    /// @brief Set jet-matched generated jet WTA eta
    void setWTAEta(const Float_t& eta) { fWTAEta = {eta}; }
    /// @brief Set jet-matched generated jet WTA phi
    void setWTAPhi(const Float_t& phi) { fWTAPhi = {phi}; }

    //
    // Getters
    //

    /// @brief Transverse momentum 
    Float_t pt() const         { return fPt; }
    /// @brief Pseudorapidity 
    Float_t eta() const        { return fEta; }
    /// @brief Azimuthal angle 
    Float_t phi() const        { return fPhi; }
    /// @brief Return reconstructed jet parameters
    TVector3 vec() const 
    { TVector3 v; v.SetPtEtaPhi(fPt, fEta, fPhi); return v; }
    /// @brief Pseudorapidity of the WTA axis 
    Float_t WTAEta() const     { return fWTAEta; }
    /// @brief Azimuthal angle of the WTA axis
    Float_t WTAPhi() const     { return fWTAPhi; }

  private:

    /// @brief  Transverse momentum
    Float_t fPt;
    /// @brief Pseudorapidity
    Float_t fEta;
    /// @brief Azimuthal angle
    Float_t fPhi;
    /// @brief Pseudorapidity of the WTA axis
    Float_t fWTAEta;
    /// @brief Azimuthal angle of the WTA axis
    Float_t fWTAPhi;

  ClassDef(BaseJet, 1)
};

#endif