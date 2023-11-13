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
    void setPt(const Double_t& lo, const Double_t& hi) { fPt[0]=lo; fPt[1]=hi; }
    /// @brief Jet R of cone radius maximum 
    void setConeR(const Double_t& max) { fConeR=max; }
    /// @brief Require reconstructed jet to have a matching to generated jet
    void setMustHaveGenMathing() { fMustHaveGenMatching=kTRUE; }
    /// @brief Preudorapidity of the reconstructed jet
    void setEta(const Double_t& lo, const Double_t& hi) { fEta[0]=lo; fEta[1]=hi; }
    /// @brief Report cut limits and passed/failed statistics
    void report();
    /// @brief Check if jet passes the cut 
    virtual Bool_t pass(const RecoJet* jet);
    /// @brief Check if jet passes the cut 
    virtual Bool_t pass(const GenJet* jet);
    /// @brief Set verbose mode
    void setVerbose() { fVerbose = kTRUE; }

  private:
    
    /// @brief Jet pT
    Double_t fPt[2];
    /// @brief Jet cone radius
    Double_t fConeR;
    /// @brief Jet must have generated jet matching
    Bool_t   fMustHaveGenMatching;
    /// @brief Pseudorapidity of the jet
    Double_t fEta[2];
    /// @brief Print status for each jet
    Bool_t   fVerbose; 
    /// @brief Number of jets passed cut
    Long64_t fJetsPassed;
    /// @brief Number of jet failed cut
    Long64_t fJetsFailed;

    ClassDef(JetCut, 0)
};

#endif // #define JetCut_h