/**
 * @file Manager.h
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Program manager - main entry point to the framework
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef Manager_h
#define Manager_h

// Jet analysis headers
#include "BaseAnalysis.h"
#include "BaseReader.h"
#include "Collections.h"

// ROOT headers
#include "TObject.h"
#include "Rtypes.h"

//________________
class Manager {
  public:
    /// @brief Constructor
    Manager();
    /// @brief Destructor
    virtual ~Manager();

    /// @brief Initialize all objects used
    void init();
    /// @brief Run analysis over events
    void performAnalysis();
    /// @brief Finish
    void finish();

    /// @brief Report details
    void report();

    /// @brief Add analysis to the collection
    void addAnalysis(BaseAnalysis* ana) { fAnalysisCollection->push_back(ana); }
    /// @brief Set event reader
    void setEventReader(BaseReader* reader) { fEventReader = reader; }

  private:
    /// @brief Pointer to analysis collection
    AnalysisCollection *fAnalysisCollection;
    /// @brief Poiter to event reader
    BaseReader *fEventReader;
    /// Number of events in input
    Long64_t fEventsInChain;

  ClassDef(Manager, 0)
};

#endif // #define ANALYSISMANAGER_H