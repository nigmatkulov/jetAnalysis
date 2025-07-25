/** @file DiJet.h
 * 
 * @brief Header file for the DiJet class, which represents a pair of jets in high-energy physics.
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @version 1.0
 * @date 2025-06-05
 * 
 * @copyright Copyright (c) 2025
 */

#ifndef DiJet_h
#define DiJet_h

// ROOT headers
#include "TObject.h"

//________________
class DiJet : public TObject {
  public:
    /// @brief Default constructor
    DiJet() : TObject(), fLeadJetPt(0), fLeadJetEtaLab(0), fLeadJetEtaCM(0), fLeadJetPhi(0),
              fSubLeadJetPt(0), fSubLeadJetEtaLab(0), fSubLeadJetEtaCM(0), fSubLeadJetPhi(0) {}

    /// @brief Constructor with parameters
    DiJet(float leadJetPt, float leadJetEtaLab, float leadJetEtaCM, float leadJetPhi,
          float subLeadJetPt, float subLeadJetEtaLab, float subLeadJetEtaCM, float subLeadJetPhi);

    /// @brief Copy constructor
    DiJet(const DiJet& other) : TObject(other),
        fLeadJetPt(other.fLeadJetPt), fLeadJetEtaLab(other.fLeadJetEtaLab),
        fLeadJetEtaCM(other.fLeadJetEtaCM), fLeadJetPhi(other.fLeadJetPhi),
        fSubLeadJetPt(other.fSubLeadJetPt),
        fSubLeadJetEtaLab(other.fSubLeadJetEtaLab), fSubLeadJetEtaCM(other.fSubLeadJetEtaCM),
        fSubLeadJetPhi(other.fSubLeadJetPhi) {};

    /// @brief Assignment operator
    DiJet& operator=(const DiJet& other);
    /// @brief Compare two DiJet objects
    bool operator==(const DiJet& other) const;
    /// @brief Compare two DiJet objects
    bool operator!=(const DiJet& other) const { return !(*this == other); }

    /// @brief Destructor
    virtual ~DiJet() { /* Empty*/ }

    /// @brief Print parameters of the given jet
    void print();

    //
    // Setters
    //

    /// @brief Set the transverse momentum of the leading jet
    void setLeadJetPt(const float& pt) { fLeadJetPt = pt; }
    /// @brief Set the pseudorapidity of the leading jet in the lab frame
    void setLeadJetEtaLab(const float& eta) { fLeadJetEtaLab = eta; }
    /// @brief Set the pseudorapidity of the leading jet in the center of mass frame
    void setLeadJetEtaCM(const float& eta) { fLeadJetEtaCM = eta; }
    /// @brief Set the azimuthal angle of the leading jet
    void setLeadJetPhi(const float& phi) { fLeadJetPhi = phi; }

    /// @brief Get the transverse momentum of the leading jet
    float leadJetPt() const { return fLeadJetPt; }
    /// @brief Get the pseudorapidity of the leading jet in the lab frame
    float leadJetEtaLab() const { return fLeadJetEtaLab; }
    /// @brief Get the pseudorapidity of the leading jet in the center of mass frame
    float leadJetEtaCM() const { return fLeadJetEtaCM; }
    /// @brief Get the azimuthal angle of the leading jet
    float leadJetPhi() const { return fLeadJetPhi; }

    /// @brief Set the subleading jet transverse momentum
    void setSubLeadJetPt(const float& pt) { fSubLeadJetPt = pt; }
    /// @brief Set the pseudorapidity of the subleading jet in the lab frame
    void setSubLeadJetEtaLab(const float& eta) { fSubLeadJetEtaLab = eta; }
    /// @brief Set the pseudorapidity of the subleading jet in the center of mass frame
    void setSubLeadJetEtaCM(const float& eta) { fSubLeadJetEtaCM = eta; }
    /// @brief Set the azimuthal angle of the subleading jet
    void setSubLeadJetPhi(const float& phi) { fSubLeadJetPhi = phi; }

    /// @brief Get the subleading jet transverse momentum
    float subLeadJetPt() const { return fSubLeadJetPt; }
    /// @brief Get the pseudorapidity of the subleading jet in the lab frame
    float subLeadJetEtaLab() const { return fSubLeadJetEtaLab; }
    /// @brief Get the pseudorapidity of the subleading jet in the center of mass frame
    float subLeadJetEtaCM() const { return fSubLeadJetEtaCM; }
    /// @brief Get the azimuthal angle of the subleading jet
    float subLeadJetPhi() const { return fSubLeadJetPhi; }

    /// @brief Get the average transverse momentum of the dijet
    float ptAve() const { return 0.5 * (fLeadJetPt + fSubLeadJetPt); }
    /// @brief Get the total transverse momentum of the dijet
    float pt() const { return fLeadJetPt + fSubLeadJetPt; }
    /// @brief Get the pseudorapidity of the dijet in the lab frame
    float etaLab() const { return 0.5 * (fLeadJetEtaLab + fSubLeadJetEtaLab); }
    /// @brief Get the pseudorapidity of the dijet in the center of mass frame
    float etaCM() const { return 0.5 * (fLeadJetEtaCM + fSubLeadJetEtaCM); }
    /// @brief Get the opening azimuthal angle of the dijet
    float dPhi() const { return deltaPhi(fLeadJetPhi, fSubLeadJetPhi); }
    /// @brief Get the azimuthal angle of the dijet
    float phi() const;
    /// @brief Get the relative pseudorapidity between two jets
    float dEta() const { return 0.5 * (fLeadJetEtaCM - fSubLeadJetEtaCM); }
    /// @brief Alias for dEta() in the center of mass frame
    float dEtaCM() const { return dEta(); }

    /// @brief Get relative azimuthal angle between two jets 
    static float deltaPhi(const float& phi1, const float &phi2);
    /// @brief Wrap angle to the range [0, 2π) 
    static float wrapTo0to2Pi(const float &angle);

    /// @brief Clean parameters of the dijet
    void cleanParameters();

  private:

    /// Leading jet pT
    Float_t fLeadJetPt;
    /// Leading jet pseudorapidity in the lab frame
    Float_t fLeadJetEtaLab;
    /// Leading jet pseudorapidity in the center of mass frame
    Float_t fLeadJetEtaCM;
    /// Leading jet azimuthal angle
    Float_t fLeadJetPhi;

    /// Subleading jet pT
    Float_t fSubLeadJetPt;
    /// Subleading jet pseudorapidity in the lab frame
    Float_t fSubLeadJetEtaLab;
    /// Subleading jet pseudorapidity in the center of mass frame
    Float_t fSubLeadJetEtaCM;
    /// Subleading jet azimuthal angle
    Float_t fSubLeadJetPhi;

    ClassDef(DiJet, 1)
};

#endif // DiJet_h

