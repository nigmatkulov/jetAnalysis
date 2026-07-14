#include "AnalysisContract.h"

#include <cassert>
#include <cmath>

int main() {
    using analysis_contract::passesDijetSelection;
    assert(passesDijetSelection(50.f, 40.f, 0.f, 0.f, 2.2f, 1.9f));
    assert(!passesDijetSelection(49.99f, 40.f, 0.f, 0.f, 2.2f, 1.9f));
    assert(!passesDijetSelection(50.f, 39.99f, 0.f, 0.f, 2.2f, 1.9f));
    assert(passesDijetSelection(50.f, 40.f, 1.9f, 0.f, 2.2f, 1.9f));
    assert(passesDijetSelection(50.f, 40.f, 0.f, 0.f,
                                analysis_contract::kTwoPiOverThree, 1.9f));
    return 0;
}
