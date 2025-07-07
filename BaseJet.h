/**
 * @file BaseJet.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Base class for jet description
 * @version 1.0
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
    /// @brief Copy constructor
    BaseJet(const BaseJet& other) : TObject(other),
        fPt(other.fPt), fEta(other.fEta), fPhi(other.fPhi),
        fWTAEta(other.fWTAEta), fWTAPhi(other.fWTAPhi) {}

    /// @brief Assignment operator
    BaseJet& operator=(const BaseJet& other);
    /// @brief Equality operator to compare two jets
    bool operator==(const BaseJet& other) const;
    /// @brief Inequality operator to compare two jets
    bool operator!=(const BaseJet& other) const { return !(*this == other); }

    //
    // Setters
    //

    /// @brief Set jet unique identifier
    void setId(const unsigned int& id) { fId = static_cast<UInt_t>(id); }
    /// @brief Set jet unique identifier
    void setId(const int& id) { fId = static_cast<UInt_t>(id); }
    /// @brief Set jet-matched generated jet transverse momentum
    void setPt(const float& pt)      { fPt = {pt}; }
    /// @brief Set jet-matched generated jet transverse momentum
    void setPt(const double& pt)     { fPt = static_cast<Float_t>(pt); }
    /// @brief Set jet-matched generated jet eta
    void setEta(const float& eta)    { fEta = {eta}; }
    /// @brief Set jet-matched generated jet eta
    void setEta(const double& eta)   { fEta = static_cast<Float_t>(eta); }
    /// @brief Set jet-matched generated jet phi
    void setPhi(const float& phi)    { fPhi = {phi}; }
    /// @brief Set jet-matched generated jet phi
    void setPhi(const double& phi)   { fPhi = static_cast<Float_t>(phi); }
    /// @brief Set jet-matched generated jet WTA eta
    void setWTAEta(const float& eta) { fWTAEta = {eta}; }
    /// @brief Set jet-matched generated jet WTA eta
    void setWTAEta(const double& eta) { fWTAEta = static_cast<Float_t>(eta); }
    /// @brief Set jet-matched generated jet WTA phi
    void setWTAPhi(const float& phi) { fWTAPhi = {phi}; }
    /// @brief Set jet-matched generated jet WTA phi
    void setWTAPhi(const double& phi) { fWTAPhi = static_cast<Float_t>(phi); }

    //
    // Getters
    //

    /// @brief Get jet unique identifier
    unsigned int id() const { return static_cast<unsigned int>(fId); }
    /// @brief Transverse momentum 
    float pt() const         { return static_cast<float>(fPt); }
    /// @brief Pseudorapidity 
    float eta() const        { return static_cast<float>(fEta); }
    /// @brief Azimuthal angle 
    float phi() const        { return static_cast<float>(fPhi); }
    /// @brief Return reconstructed jet parameters
    TVector3 vec() const 
    { TVector3 v; v.SetPtEtaPhi(fPt, fEta, fPhi); return v; }
    /// @brief Pseudorapidity of the WTA axis 
    float WTAEta() const     { return static_cast<float>(fWTAEta); }
    /// @brief Azimuthal angle of the WTA axis
    float WTAPhi() const     { return static_cast<float>(fWTAPhi); }

  private:

    /// @brief Jet unique identifier
    UInt_t fId; 
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
