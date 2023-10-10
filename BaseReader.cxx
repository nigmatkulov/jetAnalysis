// JetAnalysis headers
#include "BaseReader.h"

// ROOT headers
#include "TString.h"

ClassImp(BaseReader)

//_________________
BaseReader::BaseReader() : fReaderStatus{0} {
    /* empty */
}

//_________________
void BaseReader::report() {
    std::cout << Form("\nReporting from the BaseReader class\n");

    // Need to call event, track, jet and other reportes here
    // if (mEventCut) {
    //     temp += mEventCut->report();
    // }
    // else {
    //     temp += "NONE";
    // }
    std::cout << "\n" << std::endl;
}