#ifndef BASEREADER_H
#define BASEREADER_H

// ROOT headers
#include "TObject.h"
#include "Rtypes.h"

// JetAnalysis headers
#include "Event.h"

// C++ headers
#include <iostream>

//_________________
class BaseReader {
  public:
    /// @brief Default constructor
    BaseReader();
    /// @brief Destructor
    virtual ~BaseReader()  { /* empty */ }

    /// @brief Read and return a pointer to the Event from the given event
    /// @return Instance of Event class
    virtual Event* returnEvent() = 0;

    /// @brief Initialize event reader
    /// @return Base return function
    virtual Int_t init()  { std::cout << "BaseReader::init() - Do nothing\n"; return 0; }

    /// @brief Finish reading
    virtual void finish() { /* empty*/ }

    /// @brief Return reader status
    /// @return 0 - good, else - bad
    Int_t status() const { return fReaderStatus; }

    /// @brief Report reader status including cuts
    virtual void report();

    virtual Long64_t nEventsTotal() const { return 0; }

  protected:
    /// @brief Reader status. 0 - good, 1 - error, 2 - EOF
    Int_t fReaderStatus;

    ClassDef(BaseReader, 0)
};

#endif // #define BASEREADER_H