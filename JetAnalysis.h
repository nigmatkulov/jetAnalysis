#ifndef JETANALYSIS_H
#define JETANALYSIS_H

// Load ROOT libraries
#include "TObject.h"
#include "TString.h"
#include "Rtypes.h"
#include "TChain.h"

// Jet analysis headers
#include "BaseAnalysis.h"
#include "Event.h"

//________________
class JetAnalysis : public BaseAnalysis {
  public:
    /// @brief Default constructor
    JetAnalysis();
    /// @brief Destructor
    virtual ~JetAnalysis();

    /// @brief Initialize variables and functions
    void init();
    /// @brief Process event
    void processEvent(const Event* ev);
    /// @brief Finish analysis
    void finish();

    /// @brief Returns reports of all cuts applied and correlation functions being done
    virtual void report();
    /// @brief Return a TList of objects to be written as output
    virtual TList* getOutputList();

    /// @brief Set debug information
    void setDebug(const Bool_t debug) { fDebug = debug; }

  private:

    /// @brief Pring debug information
    Bool_t fDebug;

  ClassDef(JetAnalysis, 0)
};

#endif // #define JETANALYSIS_H