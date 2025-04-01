/**
 * @file BaseReader.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Base class for event readers
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// JetAnalysis headers
#include "BaseReader.h"

// ROOT headers
#include "TString.h"

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