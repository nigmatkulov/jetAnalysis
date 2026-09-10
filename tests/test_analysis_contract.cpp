#include "AnalysisContract.h"

#undef NDEBUG // Keep contract checks active in Release builds.
#include <cassert>
#include <cmath>
#include <initializer_list>

int main() {
    using analysis_contract::passesDijetSelection;
    assert(passesDijetSelection(50.f, 40.f, 0.f, 0.f, 2.2f, 1.9f));
    assert(!passesDijetSelection(49.99f, 40.f, 0.f, 0.f, 2.2f, 1.9f));
    assert(!passesDijetSelection(50.f, 39.99f, 0.f, 0.f, 2.2f, 1.9f));
    assert(passesDijetSelection(50.f, 40.f, 1.9f, 0.f, 2.2f, 1.9f));
    assert(passesDijetSelection(50.f, 40.f, 0.f, 0.f,
                                analysis_contract::kTwoPiOverThree, 1.9f));
    // Exercise both data/MC orientations, nominal and alternative shifts.
    for (bool isMc : {false, true}) {
        for (bool isPbGoing : {false, true}) {
            for (float eta : {-2.4f, -1.9f, -0.2f, 0.f, 0.2f, 1.9f, 2.4f}) {
                const bool flip = isMc ? isPbGoing : !isPbGoing;
                const float expectedLab = flip ? -eta : eta;
                assert(analysis_contract::etaLab(eta, isPbGoing, isMc) == expectedLab);
                assert(analysis_contract::etaLab(eta, isPbGoing, isMc, true) == eta);
                for (float shift : {0.f, 0.460f, analysis_contract::kNominalEtaShift, 0.490f}) {
                    const float expectedCM = flip ? -(eta + shift) : eta - shift;
                    assert(analysis_contract::etaCM(eta, shift, isPbGoing, isMc) == expectedCM);
                    assert(analysis_contract::etaCM(eta, shift, isPbGoing, isMc, false) == expectedCM);
                    assert(analysis_contract::etaCM(eta, shift, isPbGoing, isMc, true) == eta);
                }
            }
        }
    }
    // A boundary dijet enters pp acceptance but fails the legacy pPb boost.
    const float ppEta = analysis_contract::etaCM(-1.8f, 0.465f, false, true, true);
    const float pPbEta = analysis_contract::etaCM(-1.8f, 0.465f, false, true);
    assert(passesDijetSelection(60.f, 45.f, ppEta, ppEta, 3.f, 1.9f));
    assert(!passesDijetSelection(60.f, 45.f, pPbEta, pPbEta, 3.f, 1.9f));
    return 0;
}
