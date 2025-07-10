/**
 * @file JetCut.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Example of jet cut
 * @version 1.1
 * @date 2025-06-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef JetCut_h
#define JetCut_h

// Jet analysis headers
#include "RecoJet.h"
#include "GenJet.h"

// ROOT headers
#include "Rtypes.h"

// C++ headers
#include <limits>
#include <iostream>

//________________
class JetCut {
  public:
    /// @brief Construtor
    JetCut();
    /// @brief Destructor
    virtual ~JetCut();
    /// @brief Assignment operator
    JetCut& operator=(const JetCut& other);
    /// @brief Copy constructor (allow compiler to generate it)
    JetCut(const JetCut& other) = default;
    /// @brief Compare two JetCut objects
    bool operator==(const JetCut& other) const;
    /// @brief Compare two JetCut objects
    bool operator!=(const JetCut& other) const { return !(*this == other); }

    //
    // Setters
    //

    /// @brief Jet momentum cut
    void setPt(const double& lo, const double& hi) { fPt[0] = static_cast<float>(lo); fPt[1] = static_cast<float>(hi); }
    /// @brief Jet momentum cut
    void setPt(const float& lo, const float& hi)   { fPt[0] = lo; fPt[1] = hi; }
    /// @brief Jet R of cone radius maximum 
    void setConeR(const double& max) { fConeR = static_cast<float>(max); }
    /// @brief Jet R of cone radius maximum 
    void setConeR(const float& max)   { fConeR = max; }
    /// @brief Preudorapidity of the reconstructed jet in the laboratory frame
    void setEtaLab(const double& lo, const double& hi) { fEtaLab[0] = static_cast<float>(lo); fEtaLab[1] = static_cast<float>(hi); }
    /// @brief Preudorapidity of the reconstructed jet in the laboratory frame
    void setEtaLab(const float& lo, const float& hi)   { fEtaLab[0] = lo; fEtaLab[1] = hi; }
    /// @brief Preudorapidity of the reconstructed jet in the center-of-mass frame
    void setEtaCM(const double& lo, const double& hi) { fEtaCM[0] = static_cast<float>(lo); fEtaCM[1] = static_cast<float>(hi); }
    /// @brief Preudorapidity of the reconstructed jet in the center-of-mass frame
    void setEtaCM(const float& lo, const float& hi)   { fEtaCM[0] = lo; fEtaCM[1] = hi; }
    /// @brief Set selection method: 0 - no selection (default), 1 - trkMaxPt/RawPt, 2 - jetId
    void setSelectionMethod(const int& method) {
        if (method < 0 || method > 2) {
            std::cout << "JetCut::setSelectionMethod: method must be 0, 1 or 2, but got " << method << ". Setting to default (0)." << std::endl;
            fSelectionMethod = 0;
        } else {
            fSelectionMethod = method;
        }
    }

    /// @brief Set loose or tight jetId cut (default: true = loose)
    void setLooseJetIdCut(const bool& loose = true) {
        fLooseJetIdCut = loose; // true = loose, false = tight
    }

    /// @brief Set verbose mode
    void setVerbose() { fVerbose = true; }

    /// @brief Report cut limits and passed/failed statistics
    void report();
    /// @brief Check if jet passes the cut 
    virtual bool pass(const RecoJet* jet, bool isCM, bool isMC, bool requireMatching);
    /// @brief Check if jet passes the cut 
    virtual bool pass(const GenJet* jet, bool isCM);


  private:
    
    /// @brief Jet pT
    float fPt[2];
    /// @brief Jet cone radius
    float fConeR;
    /// @brief Pseudorapidity of the jet in the laboratory frame
    float fEtaLab[2];
    /// @brief Pseudorapidity of the jet in the center-of-mass frame
    float fEtaCM[2];
    /// @brief Selection method: 0 - no selection (default), 1 - trkMaxPt/RawPt, 2 - jetId
    int   fSelectionMethod;
    /// @brief Loose or tight jetId cut (default: true = loose)
    bool fLooseJetIdCut; 
    /// @brief Print status for each jet
    bool   fVerbose; 
};

#endif // #define JetCut_h
