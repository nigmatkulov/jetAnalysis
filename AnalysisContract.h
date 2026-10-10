#ifndef ANALYSIS_CONTRACT_H
#define ANALYSIS_CONTRACT_H

#include <cmath>

// Shared producer/test contract. Momenta are in GeV; angles are in radians.
// These constants define the fiducial selection, not histogram display ranges.
namespace analysis_contract {

constexpr float kNominalEtaShift = 0.465f;
constexpr float kNominalJetEtaMax = 1.9f;
constexpr float kLeadingJetPtMin = 50.f;
constexpr float kSubleadingJetPtMin = 40.f;
constexpr float kTwoPiOverThree = 2.09439510239319549f;
constexpr float kDijetMatchDeltaR = 0.2f;

// pp retains the stored detector orientation; pPb uses the legacy beam convention.
// First orient the lab axis, then subtract shift: etaCM = etaLab - shift.
// The different data/MC signs encode the forest coordinate convention; the
// direction flag alone is not sufficient to interpret the stored eta sign.
// This is a shifted pseudorapidity observable (massless approximation), not
// an exact Lorentz boost of a massive jet four-vector or true rapidity y.
inline float etaCM(float eta, float shift, bool isPbGoing, bool isMc,
                   bool isPp = false) {
    if (isPp) return eta;
    if (isMc) return isPbGoing ? -1.f * (eta + shift) : eta - shift;
    return isPbGoing ? eta - shift : -1.f * (eta + shift);
}

inline float etaLab(float eta, bool isPbGoing, bool isMc, bool isPp = false) {
    if (isPp) return eta;
    if (isMc) return isPbGoing ? -eta : eta;
    return isPbGoing ? eta : -eta;
}

// Select the globally leading two jets first, then require BOTH to be inside
// etaMax. This differs from finding the two leading jets inside acceptance.
// Boundaries are inclusive; deltaPhi must already be wrapped into [-pi, pi].
inline bool passesDijetSelection(float leadingPt, float subleadingPt,
                                 float leadingEtaCM, float subleadingEtaCM,
                                 float deltaPhi, float etaMax) {
    return leadingPt >= kLeadingJetPtMin &&
           subleadingPt >= kSubleadingJetPtMin &&
           std::abs(leadingEtaCM) <= etaMax &&
           std::abs(subleadingEtaCM) <= etaMax &&
           std::abs(deltaPhi) >= kTwoPiOverThree;
}

} // namespace analysis_contract

#endif
