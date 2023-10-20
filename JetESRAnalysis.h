/**
 * @file JetESRAnalysis.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Jet energy scale and resolution analysis
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef JetESRAnalysis_h
#define JetESRAnalysis_h

// Load ROOT libraries
#include "TObject.h"
#include "TString.h"
#include "Rtypes.h"
#include "TChain.h"

// Jet analysis headers
#include "BaseAnalysis.h"
#include "BasicHistoManager.h"
#include "Event.h"

//________________
class JetESRAnalysis : public BaseAnalysis {
  public:
    /// @brief Default constructor
    JetESRAnalysis();
    /// @brief Destructor
    virtual ~JetESRAnalysis();

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
    /// @brief Add histogram manager to the analysis
    void addHistoManager(BasicHistoManager *hm) { fHM = hm; }

  private:

    /// @brief Pring debug information
    Bool_t fDebug;
    /// @brief Histogram manager
    BasicHistoManager *fHM;

  ClassDef(JetESRAnalysis, 0)
};

#endif // #define JetESRAnalysis_h