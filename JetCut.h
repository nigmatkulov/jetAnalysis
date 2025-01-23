/**
 * @file JetCut.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Example of jet cut
 * @version 0.1
 * @date 2023-10-19
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

//________________
class JetCut {
  public:
    /// @brief Construtor
    JetCut();
    /// @brief Destructor
    virtual ~JetCut();

    //
    // Setters
    //

    /// @brief Jet momentum cut
    void setPt(const double& lo, const double& hi) { fPt[0]=lo; fPt[1]=hi; }
    /// @brief Jet R of cone radius maximum 
    void setConeR(const double& max) { fConeR=max; }
    /// @brief Require reconstructed jet to have a matching to generated jet
    void setMustHaveGenMathing() { fMustHaveGenMatching=true; }
    /// @brief Preudorapidity of the reconstructed jet
    void setEta(const double& lo, const double& hi) { fEta[0]=lo; fEta[1]=hi; }
    /// @brief Cut jets with values below low and above high at midrapidity
    void setTrackMaxPtOverRawPt(const double& lo, const double& hi) { fTrackMaxPtOverRawPt[0]=lo; fTrackMaxPtOverRawPt[1]=hi; }
    /// @brief Report cut limits and passed/failed statistics
    void report();
    /// @brief Check if jet passes the cut 
    virtual bool pass(const RecoJet* jet);
    /// @brief Check if jet passes the cut 
    virtual bool pass(const GenJet* jet);
    /// @brief Set verbose mode
    void setVerbose() { fVerbose = true; }

  private:
    
    /// @brief Jet pT
    double fPt[2];
    /// @brief Jet cone radius
    double fConeR;
    /// @brief Jet must have generated jet matching
    bool   fMustHaveGenMatching;
    /// @brief Pseudorapidity of the jet
    double fEta[2];
    /// @brief Cut on charged particle fraction inside jet (at midrapidity cut jets with val<[0] and val>[1])
    double fTrackMaxPtOverRawPt[2];
    /// @brief Print status for each jet
    bool   fVerbose; 
    /// @brief Number of jets passed cut
    long int fJetsPassed;
    /// @brief Number of jet failed cut
    long int fJetsFailed;

    ClassDef(JetCut, 0)
};

#endif // #define JetCut_h