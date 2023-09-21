// C++ headers
#include <iostream>

// Jet analysis headers
#include "JetAnalysis.h"

//________________
/// @brief The prorgram that launches the physics analysis
/// @param argc Number of arguments
/// @param argv Argument list
/// @return 0 in case of OKAY
int main(int argc, char const *argv[]) {
    /* code */

    JetAnalysis *ana = new JetAnalysis();
    ana->init();
    ana->processData();
    ana->finish();
    
    return 0;
}
