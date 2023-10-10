#ifndef BASEANALYSIS_H
#define BASEANALYSIS_H

// ROOT headers
#include "TObject.h"
#include "TList.h"

// Forward declaration
class Event;

//_________________
class BaseAnalysis {

  public:
    /// @brief Default constructor
    BaseAnalysis()          { /* noop */ }
    /// @brief Destructor
    virtual ~BaseAnalysis() { /* noop */ }

    /// @brief Initialize analysis
    virtual void init() = 0; ///<

    /// @brief Returns reports of all cuts
    virtual void report() = 0; //!<

    /// @brief Obtain number of objects to be written as an output
    virtual TList *getOutputList() = 0; ///<

    /// @brief Event processing
    virtual void processEvent(const Event *) = 0; ///<

    /// @brief Finish analysis
    virtual void finish() = 0; ///<

    ClassDef(BaseAnalysis, 0)
};

#endif // #define BASEANALYSIS_H
