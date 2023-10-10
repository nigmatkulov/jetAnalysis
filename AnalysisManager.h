#ifndef ANALYSISMANAGER_H
#define ANALYSISMANAGER_H

// Jet analysis headers
#include "BaseAnalysis.h"
#include "BaseReader.h"
#include "Collections.h"

// ROOT headers
#include "TObject.h"
#include "Rtypes.h"

//________________
class AnalysisManager {
  public:
    /// @brief Constructor
    AnalysisManager();
    /// @brief Destructor
    virtual ~AnalysisManager();

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

  ClassDef(AnalysisManager, 0)
};

#endif // #define ANALYSISMANAGER_H