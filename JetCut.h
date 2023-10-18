#ifndef JetCut_h
#define JetCut_h

// Jet analysis headers
#include "Jet.h"

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

    /// @brief Reconstructed jet momentum cut
    void setRecoPt(const Float_t& lo, const Float_t& hi) { fRecoPt[0]=lo; fRecoPt[1]=hi; }
    /// @brief Reconstructed jet R of cone radius maximum 
    void setRecoConeR(const Float_t& max) { fRecoConeR=max; }
    /// @brief Require reconstructed jet to have a matching to generated jet
    void setMustHaveGenMathing() { fMustHaveGenMatching=kTRUE; }
    /// @brief Generated jet momentum cut
    void setRefPt(const Float_t& lo, const Float_t& hi) { fRefPt[0]=lo; fRefPt[1]=hi; }
    /// @brief Generated jet R of cone radius maximum 
    void setRefConeR(const Float_t& max) { fRefConeR=max; }
    /// @brief Generated jet flavor for b
    void setRefFlavorForB(const Int_t& lo, const Int_t& hi) { fRefFlavorForB[0]=lo; fRefFlavorForB[1]=hi; }
    /// @brief Report cut limits and passed/failed statistics
    void report();
    /// @brief Check if jet passes the cut 
    virtual Bool_t pass(const Jet* jet);
    /// @brief Set verbose mode
    void setVerbose() { fVerbose = kTRUE; }

  private:
    
    /// @brief Reconstructed jet pT
    Float_t fRecoPt[2];
    /// @brief Reconstructed jet cone radius
    Float_t fRecoConeR;
    /// @brief Reonctructed jet must have generated jet matching
    Bool_t  fMustHaveGenMatching;
    /// @brief Generated jet pT
    Float_t fRefPt[2];
    /// @brief Generated jet cone radius
    Float_t fRefConeR;
    /// @brief Generated jet flavor for b
    Int_t   fRefFlavorForB[2];
    /// @brief Print status for each jet
    Bool_t  fVerbose; 

    /// @brief Number of jets passed cut
    Long64_t fJetsPassed;
    /// @brief Number of jet failed cut
    Long64_t fJetsFailed;

    ClassDef(JetCut, 0)
};

#endif // #define JetCut_h