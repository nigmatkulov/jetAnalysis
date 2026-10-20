#include "JetCorrector.h"
#include "TSystem.h"
#include "test_jec_lookup.C"

int main() {
    try {
        test_jec_lookup();
    } catch (const std::exception &error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
    return 0;
}
