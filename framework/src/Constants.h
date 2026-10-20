#pragma once
#include <array>
#include <vector>
#include <RtypesCore.h>
namespace jetanalysis {
inline const char* RED="\033[1;31m";
inline const char* GREEN="\033[1;32m";
inline const char* RESET="\033[0m";
// Eta shifts for pPb 8.16 TeV collisions, used for etaCM calculation
const int nEtaShifts = 13;
static constexpr std::array<float, nEtaShifts> etaShift{0.460, 0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470, 0.475, 0.480, 0.485, 0.490 };
const int nEtaCuts = 8;
static constexpr std::array<float, nEtaCuts> etaCuts{1.5, 1.6, 1.7, 1.8, 1.9, 2.4, 2.5, 3.0};
// Run IDs used for the data-taking-period comparison histograms.
const int nDiagnosticRuns = 5;
const int nDiagnosticRunHistograms = nDiagnosticRuns + 1;
static constexpr std::array<UInt_t, nDiagnosticRuns> diagnosticRunIds{
    285480, 285505, 285517, 285832, 285993
};
static constexpr std::array<float, nDiagnosticRuns> diagnosticRunPileup{
    0.04f, 0.25f, 0.1f, 0.004f, 0.2f
};

// One event-level split is reused by truth, reco, response, miss and fake fills.
#ifdef DO_MIXED_PRIOR
// Trial alternative prior for unfolding regularization. The target is an
// equal mixture of two Gaussians, expressed as a ratio to the configured
// nominal Gaussian centered at zero in dijet etaCM.
constexpr double nominalPriorMean = 0.;
constexpr double nominalPriorSigma = 0.9;
constexpr double leftPriorMean = -0.2;
constexpr double leftPriorSigma = 1.2 * nominalPriorSigma;
constexpr double rightPriorMean = 0.5;
constexpr double rightPriorSigma = 0.8 * nominalPriorSigma;
#endif

namespace unfolding_diagnostics {
    enum PairCategory : int {
        kMatchedDirect = 1,
        kMatchedSwapped,
        kGenPassRecoFail,
        kRecoPassGenFail,
        kSelectedPairMismatch,
        kMissingOneRecoReference,
        kMissingBothRecoReferences,
        kNPairCategories
    };
    constexpr int nPairCategories = kNPairCategories - 1;
} // namespace unfolding_diagnostics
inline const std::vector<float> etaBinEdges = {-3.60f, -3.46f, -3.31f, -3.17f, -3.02f, -2.88f, -2.74f, -2.59f, -2.45f, -2.30f, -2.16f, -2.02f, -1.87f, -1.73f, -1.58f, -1.44f, -1.30f, -1.15f, -1.01f, -0.86f, -0.72f, -0.58f, -0.43f, -0.29f, -0.14f, 0.00f, 0.14f, 0.29f, 0.43f, 0.58f, 0.72f, 0.86f, 1.01f, 1.15f, 1.30f, 1.44f, 1.58f, 1.73f, 1.87f, 2.02f, 2.16f, 2.30f, 2.45f, 2.59f, 2.74f, 2.88f, 3.02f, 3.17f, 3.31f, 3.46f, 3.60f};

}
