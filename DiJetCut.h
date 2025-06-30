/**
 * @file DiJetCut.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Class for dijet selection
 * @version 0.1
 * @date 2025-06-29
 * 
 * @copyright Copyright (c) 2025
 */

#ifndef DiJetCut_h
#define DiJetCut_h

// Jet analysis headers
#include "DiJet.h"

//________________
class DiJetCut {
  public:
    /// @brief Default constructor
    DiJetCut();
    /// @brief Destructor
    virtual ~DiJetCut();
    /// @brief Copy constructor (allow compiler to generate it)
    DiJetCut(const DiJetCut& other) = default;
    /// @brief Assignment operator
    DiJetCut& operator=(const DiJetCut& other);
    /// @brief Compare two DiJetCut objects
    bool operator==(const DiJetCut& other) const;
    /// @brief Compare two DiJetCut objects
    bool operator!=(const DiJetCut& other) const { return !(*this == other); }

    /// @brief Report cut parameters
    void report();

    //
    // Setters
    //

    /// @brief Check if dijet passes the cut
    virtual bool pass(const DiJet* dijet, bool isCM = false);

    /// @brief Set minimum pT for the leading jet
    void setLeadJetPtMinimum(const double& pTmin) { fLeadJetPt = static_cast<float>(pTmin); }
    /// @brief Set minimum pT for the leading jet
    void setLeadJetPtMinimum(const float& pTmin) { fLeadJetPt = pTmin; }
    /// @brief Set minimum pT for the subleading jet
    void setSubLeadJetPtMinimum(const double& pTmin) { fSubLeadJetPt = static_cast<float>(pTmin); }
    /// @brief Set minimum pT for the subleading jet
    void setSubLeadJetPtMinimum(const float& pTmin) { fSubLeadJetPt = pTmin; }
    /// @brief Set pseudorapidity range for the leading jet in the lab frame
    void setLeadJetEtaLab(const double& min, const double& max) {
        fLeadJetEtaLab[0] = static_cast<float>(min); fLeadJetEtaLab[1] = static_cast<float>(max);
    }
    /// @brief Set pseudorapidity range for the leading jet in the lab frame
    void setLeadJetEtaLab(const float& min, const float& max) { fLeadJetEtaLab[0] = min; fLeadJetEtaLab[1] = max;}
    /// @brief Set pseudorapidity range for the subleading jet in the lab frame
    void setSubLeadJetEtaLab(const double& min, const double& max) {
        fSubLeadJetEtaLab[0] = static_cast<float>(min); fSubLeadJetEtaLab[1] = static_cast<float>(max);
    }
    /// @brief Set pseudorapidity range for the subleading jet in the lab frame
    void setSubLeadJetEtaLab(const float& min, const float& max) { fSubLeadJetEtaLab[0] = min; fSubLeadJetEtaLab[1] = max; }
    /// @brief Set pseudorapidity range for the leading jet in the center-of-mass frame
    void setLeadJetEtaCM(const double& min, const double& max) {
        fLeadJetEtaCM[0] = static_cast<float>(min); fLeadJetEtaCM[1] = static_cast<float>(max);
    }
    /// @brief Set pseudorapidity range for the leading jet in the center-of-mass frame
    void setLeadJetEtaCM(const float& min, const float& max) { fLeadJetEtaCM[0] = min; fLeadJetEtaCM[1] = max; }
    /// @brief Set pseudorapidity range for the subleading jet in the center-of-mass frame
    void setSubLeadJetEtaCM(const double& min, const double& max) {
        fSubLeadJetEtaCM[0] = static_cast<float>(min); fSubLeadJetEtaCM[1] = static_cast<float>(max);
    }
    /// @brief Set pseudorapidity range for the subleading jet in the center-of-mass frame
    void setSubLeadJetEtaCM(const float& min, const float& max) { fSubLeadJetEtaCM[0] = min; fSubLeadJetEtaCM[1] = max; }   
    /// @brief Set the azimuthal angle between leading and subleading jets
    void setDijetDPhi(const double& dphi) { fDijetDPhiCut = static_cast<float>(dphi); }
    /// @brief Set the azimuthal angle between leading and subleading jets
    void setDijetDPhi(const float& dphi) { fDijetDPhiCut = dphi; }
    /// @brief Set verbose mode
    void setVerbose() { fVerbose = true; }

  private:

    /// @brief Minimum leading jet pT
    float fLeadJetPt;
    /// @brief Minimum subleading jet pT
    float fSubLeadJetPt;
    /// @brief Leading jet pseudorapidity range in the lab frame [min, max]
    float fLeadJetEtaLab[2]; 
    /// @brief Subleading jet pseudorapidity range in the lab frame [min, max]
    float fSubLeadJetEtaLab[2];
    /// @brief Leading jet pseudorapidity range in the center-of-mass frame [min, max]
    float fLeadJetEtaCM[2];
    /// @brief Subleading jet pseudorapidity range in the center-of-mass frame [min, max]
    float fSubLeadJetEtaCM[2];
    /// @brief Azimuthal angle between leading and subleading jets
    float fDijetDPhiCut;
    /// @brief Flag for verbose output
    bool fVerbose = false;
};

#endif // #ifndef DiJetCut_h
