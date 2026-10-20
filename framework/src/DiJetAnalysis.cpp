#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TString.h"
#include "THnSparse.h"

#include "JetCorrector.h"
#include "JetUncertainty.h"
#include "AnalysisContract.h"

#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <cstring>


#include "jetanalysis/DiJetAnalysis.h"
#include "jetanalysis/AnalysisManager.h"
#include <TDirectory.h>
#include "jetanalysis/Output.h"
#include "Histograms.h"
#include "StageTimings.h"
namespace jetanalysis {
using analysis_contract::etaCM;
using analysis_contract::etaLab;
struct JERSmearingFuncVsEta{
    double fPtLow;
    double fPtHigh;
    std::unique_ptr<TF1> fJERSmearFuncVsEta;
};

enum class JERSmearingMode : int {
    Hybrid = 0,
    Stochastic = 1,
    TruthToPseudoReco = 2
};

//________________
class JERSmearingHelper {
    public:
        /// Constructor
        JERSmearingHelper() {
            fJERSmearFuncVsPt = std::make_unique<TF1>("fJERSmearFuncVsPt", "sqrt( [0]*[0] + [1]*[1]/x )", 30., 800.);
            fJERSmearFuncVsPt->SetParameters(0.04360, 0.81815); // Default parameters (Pb-going direction for pPb 8.16 TeV collisions)
        }

        /// Destructor
        ~JERSmearingHelper() {
            if (fJERSmearFuncVsPt) fJERSmearFuncVsPt.reset();
        }

        /// Set resolution smearing function parameters
        void setJERSmearingFuncVsPtParams(const double &param0, const double &param1) {
            if (!fJERSmearFuncVsPt) throw std::runtime_error("Invalid TF1 passed to JERSmearingHelper");
            fJERSmearFuncVsPt->SetParameters(param0, param1);
        }

        /// Add a JER smearing function for a given pt range
        void addJERSmearingFuncVsEta(const double &ptLow, const double &ptHigh, std::unique_ptr<TF1> fJERSmearFunc) {
            if (!fJERSmearFunc) throw std::runtime_error("Invalid TF1 passed to JERSmearingHelper");
            JERSmearingFuncVsEta jerFunc;
            jerFunc.fPtLow = ptLow;
            jerFunc.fPtHigh = ptHigh;
            jerFunc.fJERSmearFuncVsEta = std::move(fJERSmearFunc);
            fJERSmearFuncVsEtaVec.push_back(std::move(jerFunc));
        }

        /// Return the JER value for given pt and eta
        double evalJERScalingFactorForPtEta(const double &pt, const double &eta) const {
            if ( std::abs(eta) < 0.8 || std::abs(eta) > 2.6 ) {
                return 1.0; // No EXTRA eta multiplier; the base pT resolution still applies
            }
            else {
                for (const auto &jerFunc : fJERSmearFuncVsEtaVec) {
                    if (pt >= jerFunc.fPtLow && pt < jerFunc.fPtHigh) {
                        return jerFunc.fJERSmearFuncVsEta->Eval(eta);
                    }
                }                
            }
            return 1.0; // Return 1 if no matching pt range found
        }
        
        /// Return JER scaling factor for the default, up and down variations. Default is -99 (return 1.0)
        double JERScaleFactorSystematics(const float &eta, const int &variation = -99) const {
            // Implementation for getting JER scaling factors
            if (variation == 0) {
                // Default JER scaling factor
                for (size_t i = 0; i < fJerEtaLow.size(); ++i) {
                    if (eta >= fJerEtaLow[i] && eta < fJerEtaHi[i]) {
                        return fJerDef[i];
                    }
                }
            } 
            else if (variation == -1) {
                // Down variation
                for (size_t i = 0; i < fJerEtaLow.size(); ++i) {
                    if (eta >= fJerEtaLow[i] && eta < fJerEtaHi[i]) {
                        return fJerLow[i];
                    }
                }
            } 
            else if (variation == 1) {
                // Up variation
                for (size_t i = 0; i < fJerEtaLow.size(); ++i) {
                    if (eta >= fJerEtaLow[i] && eta < fJerEtaHi[i]) {
                        return fJerHi[i];
                    }
                }
            }
            else {
                return 1.0; // No variation specified, return 1
            }
            return 1.0; // No scale factor is defined outside the configured eta bins
        }

        // Fractional MC resolution sigma(pT)/pT, distinct from the data/MC SF.
        // Clamp the evaluation point, not the jet momentum, to the fit domain.
        double resolution(const float pt, const float eta, const bool useEtaDependence) const {
            if ( !fJERSmearFuncVsPt ) {
                throw std::runtime_error("JERSmearingHelper: JER smearing function not set.");
            }
            double value{0.};
            if ( pt <= 30.) {
                value = fJERSmearFuncVsPt->Eval(31.);
            }
            else if ( pt >= 800.) {
                value = fJERSmearFuncVsPt->Eval(799.);
            }
            else {
                value = fJERSmearFuncVsPt->Eval(pt);
            }
            if (useEtaDependence) value *= evalJERScalingFactorForPtEta(pt, eta);
            return value;
        }

        /// Hash-based standard-normal draw. The entropy follows the CMS JER
        /// recommendation: jet kinematics, event rho, and event identity. The
        /// round distinguishes optional pseudo-experiment replicas. Because the
        /// JER variation is deliberately absent: all branches share one draw.
        /// The caller freezes the initial jet index before any response sort
        /// and combines it with the replica index to form round.
        double normalDeviate(const float jetPt, const float jetEta, const float rho,
                             const UInt_t runId, const UInt_t lumiId,
                             const ULong64_t eventId, const int round) const {
            auto floatBits = [](const float value) {
                std::uint32_t bits{0};
                std::memcpy(&bits, &value, sizeof(bits));
                return static_cast<std::uint64_t>(bits);
            };
            auto mix = [](std::uint64_t value) {
                value += 0x9e3779b97f4a7c15ULL;
                value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
                value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
                return value ^ (value >> 31);
            };
            std::uint64_t state = mix(eventId);
            state = mix(state ^ (static_cast<std::uint64_t>(runId) << 32) ^ lumiId);
            state = mix(state ^ floatBits(jetPt));
            state = mix(state ^ (floatBits(jetEta) << 1));
            state = mix(state ^ (floatBits(rho) << 2));
            state = mix(state ^ static_cast<std::uint64_t>(round));
            const double u1 = (static_cast<double>((mix(state) >> 11) + 1ULL)) /
                              (static_cast<double>(1ULL << 53) + 1.0);
            const double u2 = (static_cast<double>((mix(state ^ 0xd1b54a32d192ed03ULL) >> 11) + 1ULL)) /
                              (static_cast<double>(1ULL << 53) + 1.0);
            return std::sqrt(-2.0 * std::log(u1)) * std::cos(TMath::TwoPi() * u2);
        }

        float smearRecoMomentum(const float recoPt, const float refPt, const float eta,
                                const double scaleFactor, const bool useEtaDependence,
                                const JERSmearingMode mode, const bool hasGoodGenMatch,
                                const double normal) const {
            // Hybrid rescales an existing reco-truth residual. Stochastic adds
            // only the extra width needed beyond MC; SF<1 cannot sharpen an
            // unmatched jet, so its additional stochastic width is zero.
            // TruthToPseudoReco instead models the FULL resolution around truth.
            const double jer = resolution(recoPt, eta, useEtaDependence);
            double smearedPt{recoPt};
            if (mode == JERSmearingMode::Hybrid && hasGoodGenMatch) {
                smearedPt = refPt + scaleFactor * (recoPt - refPt);
            }
            else if (mode != JERSmearingMode::TruthToPseudoReco) {
                const double extraWidth = jer * std::sqrt(std::max(scaleFactor * scaleFactor - 1.0, 0.0));
                smearedPt = recoPt * (1.0 + normal * extraWidth);
            }
            else {
                smearedPt = refPt * (1.0 + normal * scaleFactor * jer);
            }
            if (!std::isfinite(smearedPt) || smearedPt <= 0.) {
                return mode == JERSmearingMode::TruthToPseudoReco ? refPt : recoPt;
            }
            return static_cast<float>(smearedPt);
        }

        float smearTruthMomentum(const float refPt, const float eta, const double scaleFactor,
                                 const bool useEtaDependence, const double normal) const {
            const double smearedPt = refPt * (1.0 + normal * scaleFactor * resolution(refPt, eta, useEtaDependence));
            if (!std::isfinite(smearedPt) || smearedPt <= 0.) return refPt;
            return static_cast<float>(smearedPt);
        }

    private:
        // JER as a function of pt for -0.8 < eta < 0.8
        std::unique_ptr<TF1> fJERSmearFuncVsPt;
        // Vector to hold JER smearing functions for different pt and eta ranges
        std::vector<JERSmearingFuncVsEta> fJERSmearFuncVsEtaVec;
        // Data/MC JER scale factors from Summer16_25nsV1b_MC_SF_AK4PF.txt:
        // https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Summer16_25nsV1b_MC/Summer16_25nsV1b_MC_SF_AK4PF.txt
        std::vector<float> fJerEtaLow{-5.191, -3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139};
        std::vector<float> fJerEtaHi{-3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139, 5.191};
        std::vector<float> fJerDef{1.1922, 1.1869, 1.7788, 1.3418, 1.2963, 1.1512, 1.1426, 1.1000, 1.1278, 1.1609, 1.1464, 1.1948, 1.1595, 1.1595, 1.1948, 1.1464, 1.1609, 1.1278, 1.1000, 1.1426, 1.1512, 1.2963, 1.3418, 1.7788, 1.1869, 1.1922};
        std::vector<float> fJerLow{1.0434, 1.0626, 1.578, 1.1327, 1.0592, 1.0372, 1.0212, 0.9921, 1.0292, 1.0584, 1.0832, 1.1296, 1.095, 1.095, 1.1296, 1.0832, 1.0584, 1.0292, 0.9921, 1.0212, 1.0372, 1.0592, 1.1327, 1.578, 1.0626, 1.0434};
        std::vector<float> fJerHi{1.341, 1.3112, 1.9796, 1.5509, 1.5334, 1.2652, 1.264, 1.2079, 1.2264, 1.2634, 1.2096, 1.26, 1.224, 1.224, 1.26, 1.2096, 1.2634, 1.2264, 1.2079, 1.264, 1.2652, 1.5334, 1.5509, 1.9796, 1.3112, 1.341};
};


namespace {
struct Physics : EventInfo {
    void usage() { throw std::runtime_error("Invalid analysis configuration"); }
    StageTimings timings;
    const Event* currentEvent=nullptr;
    Histograms hs;
    std::unique_ptr<TF1> fVzWeight;
    std::unique_ptr<JERSmearingHelper> fJerSmearingHelper;
    std::unique_ptr<JetCorrector> jecOwner;
    std::unique_ptr<JetUncertainty> jeuOwner;
    bool isPp=false,isMc=false,isPbGoing=false,isPythia=false,useEtaDependentJER=false;
    int triggerId=0,jetSelectionMethod=2,vertexFilterSelection=0;
    JERSmearingMode jerSmearingMode=JERSmearingMode::Hybrid;
    std::uint64_t inputRecoCount=0,inputGenCount=0;
    int nEventsProcessed=0,nGoodEvents=0;
    bool isGenOverweight=false,isRecoOverweight=false;
    double weight=1.;
    std::vector<GenJet> genJets;
    std::vector<RecoJet> recoJets;
    int printEvery=100000;
    Long64_t nextPrintAt=100000,nEntries=0;
    std::chrono::steady_clock::time_point startTime;
    std::mt19937 sampleSplitGenerator{0x5eed1234u};
    std::bernoulli_distribution fillControlledSample{0.60};
    std::string output;
    std::vector<std::string> exclusions;
    bool isFillControlledSample=true;
    float ptHatRange[2]{0.f,8160.f};
void cacheRecoJetSmearing(RecoJet &jet, const JERSmearingHelper &jer,
                         JERSmearingMode mode, bool useEtaDependentJER,
                         double normal) {
    const bool hasReference = jet.refPt > 0.;
    auto momentum = [&](double sf, bool etaDependent) -> float {
        if (!hasReference && mode == JERSmearingMode::TruthToPseudoReco) return -1.f;
        const double dEta = hasReference ? jet.recoEta - jet.refEta : 0.;
        const double dPhi = hasReference ? std::remainder(
            static_cast<double>(jet.recoPhi - jet.refPhi), TMath::TwoPi()) : 0.;
        const double resolution = jer.resolution(jet.recoPt, jet.recoEta, etaDependent);
        const bool goodMatch = hasReference && std::hypot(dEta, dPhi) < 0.2 &&
            std::abs(jet.recoPt - jet.refPt) < 3.0 * resolution * jet.recoPt;
        return jer.smearRecoMomentum(jet.recoPt, hasReference ? jet.refPt : jet.recoPt,
                                     jet.recoEta, sf, etaDependent, mode, goodMatch, normal);
    };
    jet.refPtSmeared = momentum(1., false);
    jet.refPtSmearedEta = momentum(1., true);
    for (int scale = 0; scale < 9; ++scale) {
        double sf = scale == 1 ? 2. : scale == 2 ? 0.9f : 1.;
        if (scale >= 3) sf = jer.JERScaleFactorSystematics(
            jet.recoEta, (scale == 3 || scale == 6) ? -1 :
                         (scale == 4 || scale == 7) ? 0 : 1);
        jet.jerMomenta[scale] = momentum(sf, useEtaDependentJER || scale >= 6);
    }
    jet.refPtSmearedJerDown = jet.jerMomenta[3];
    jet.refPtSmearedJerDef = jet.jerMomenta[4];
    jet.refPtSmearedJerUp = jet.jerMomenta[5];
}

//________________
/// Find the two highest-pT jets without sorting the full collection.
/// Fall back to the original copy-and-sort behavior when equal or non-finite
/// pT values make the ordering ambiguous.
template <typename Jet, typename GetPt>
std::pair<const Jet*, const Jet*> selectLeadingTwo(const std::vector<Jet> &jets,
                                                   GetPt getPt,
                                                   std::vector<Jet> &sortFallback) {
    const Jet *leading = &jets[0];
    const Jet *subleading = &jets[1];
    if (getPt(*subleading) > getPt(*leading)) {
        std::swap(leading, subleading);
    }

    bool requiresSortFallback = !std::isfinite(getPt(*leading)) || !std::isfinite(getPt(*subleading));
    for (size_t iJet = 2; iJet < jets.size(); ++iJet) {
        const Jet *jet = &jets[iJet];
        const auto pt = getPt(*jet);
        if (!std::isfinite(pt)) requiresSortFallback = true;

        if (pt > getPt(*leading)) {
            subleading = leading;
            leading = jet;
        }
        else if (pt > getPt(*subleading)) {
            subleading = jet;
        }
    }

    const auto leadingPt = getPt(*leading);
    const auto subleadingPt = getPt(*subleading);
    for (const auto &jet : jets) {
        const Jet *jetAddress = &jet;
        if ((jetAddress != leading && getPt(jet) == leadingPt) ||
            (jetAddress != subleading && getPt(jet) == subleadingPt)) {
            requiresSortFallback = true;
            break;
        }
    }

    if (requiresSortFallback) {
        sortFallback.assign(jets.begin(), jets.end());
        std::sort(sortFallback.begin(), sortFallback.end(), [&](const Jet &a, const Jet &b) {
            return getPt(a) > getPt(b);
        });
        return {&sortFallback[0], &sortFallback[1]};
    }

    return {leading, subleading};
}
double eventWeight(const float &ptHat, const float &vz, TF1 &fVzWeight, Long64_t eventsInSample,
                   bool isPp = false, Float_t storedWeight = 1.f) {
    if (isPp) {
        if (!std::isfinite(storedWeight)) {
            throw std::runtime_error("Non-finite hiEvtAnalyzer/HiTree/weight in pp event");
        }
        return storedWeight;
    }

    // P04: eventsInSample is the entries in THIS invocation, before cuts.
    // Splitting one pT-hat sample into separately normalized outputs and adding
    // them changes its total normalization. Use a campaign-level convention
    // before interpreting merged pPb yields as a cross section.
    double genWeight = 1.0;
    double vzWeight = 1./ fVzWeight.Eval( vz );
    // vzWeight = 1.;
    // Half-open lookup intervals [low, high), consistent with isGoodEvent.
    // Legacy per-event cross-section coefficients times reference generated
    // counts. Their provenance/units require campaign validation; this code
    // does not obtain a cross section from the input forest metadata.
    if (ptHat >= 15.0 && ptHat < 30.)       { genWeight = 1.0404701e-06 * 961104 ; }
    else if (ptHat >= 30. && ptHat < 50.)   { genWeight = 7.7966624e-08 * 952110 ; }
    else if (ptHat >= 50. && ptHat < 80.)   { genWeight = 1.0016052e-08 * 952554 ; }
    else if (ptHat >= 80. && ptHat < 120.)  { genWeight = 1.3018269e-09 * 996844 ; }
    else if (ptHat >= 120. && ptHat < 170.)  { genWeight = 2.2648493e-10 * 964681 ; }
    else if (ptHat >= 170. && ptHat < 220.) { genWeight = 4.0879112e-11 * 999260 ; }
    else if (ptHat >= 220. && ptHat < 280.) { genWeight = 1.1898939e-11 * 964336 ; }
    else if (ptHat >= 280. && ptHat < 370.) { genWeight = 3.3364433e-12 * 995036 ; }
    else if (ptHat >= 370. && ptHat < 460.) { genWeight = 7.6612402e-13 * 958160 ; }
    else if (ptHat >= 460. && ptHat < 540.) { genWeight = 2.1341026e-13 * 981427 ; }
    else if (ptHat >= 540.)                  { genWeight = 7.9191586e-14 * 1000000; }
    genWeight /= eventsInSample;

    return genWeight * vzWeight;
}

//________________
/// Check if a reconstructed jet passes the selection criteria
/// @param recoJet: the reconstructed jet object
/// @param jetSelectionMethod: the jet selection method (0 - no selection, 1 - trackMax, 2 - jetIdTightLept, 3 - jetIdLoose)
/// @return true if the jet is good, false otherwise
bool isGoodRecoJet(const RecoJet &recoJet, const int &jetSelectionMethod) {
    // jetSelectionMethod: 0 - no selection, 1 - trackMax, 2 - jetIdTightLept, 3 - jetIdLoose

    bool debug = {false};
    if (debug) {
        std::cout << Form("\nJet parameters: ptRaw = %.2f, ptCorr = %.2f, eta = %.2f, phi = %.2f, NHF = %.3f, NEF = %.3f, CHF = %.3f, MUF = %.3f, CEF = %.3f, CHM = %d, PfCEM = %d, PfNHM = %d, PfNEM = %d, PfMUM = %d, trackMaxPt = %.2f, refPt = %.2f, refEta = %.2f jeuUp = %.2f, jeuDown = %.2f, refPtSmeared = %.2f\n",
                        recoJet.recoPtRaw, recoJet.recoPt, recoJet.recoEta, recoJet.recoPhi,
                        recoJet.recoPfNHF, recoJet.recoPfNEF, recoJet.recoPfCHF,
                        recoJet.recoPfMUF, recoJet.recoPfCEF, recoJet.recoPfCHM,
                        recoJet.recoPfCEM, recoJet.recoPfNHM, recoJet.recoPfNEM,
                        recoJet.recoPfMUM, recoJet.recoTrackMaxPt,
                        recoJet.refPt, recoJet.refEta, 
                        recoJet.recoPtJeuUp, recoJet.recoPtJeuDown, 
                        recoJet.refPtSmeared);
        std::cout << Form("Jet selection method: %d\n", jetSelectionMethod);
    }

    bool isGood = false;
    if (jetSelectionMethod == 0) {
        isGood = true; // No selection, all jets are good
    }
    else if (jetSelectionMethod == 1) {
        isGood = true;
        if ( std::abs(recoJet.recoEta) < 2.4 ) {
            if ( recoJet.recoTrackMaxPt / recoJet.recoPtRaw < 0.01 ) {
                isGood = false;
            }
            if ( recoJet.recoTrackMaxPt / recoJet.recoPtRaw > 0.98 ) {
                isGood = false;
            }
        }
    }
    else if (jetSelectionMethod == 2 || jetSelectionMethod == 3) {
        isGood = true;
        int chargedMult = recoJet.recoPfCHM + recoJet.recoPfCEM + recoJet.recoPfMUM;
        int neutralMult = recoJet.recoPfNHM + recoJet.recoPfNEM;
        int numberOfConstituents = chargedMult + neutralMult;

        float chargedEmFracCut{1.}, neutFracCut{1.};
        if ( jetSelectionMethod == 2 ) { 
            chargedEmFracCut = {0.9f};
            neutFracCut = {0.9f};
        }
        else {
            chargedEmFracCut = {0.99f};
            neutFracCut = {0.99f};
        }

        bool passNHF{false};
        bool passNEF{false};
        bool passNumOfConstituents{true};
        bool passMuonFrac{true};
        bool passChargedFrac{true};
        bool passChargedMult{true};
        bool passChargedEmFrac{true};
        bool passNeutralMult{true};

        if (debug) {
            std::cout << Form("chargedMult: %d, neutralMult: %d\n", chargedMult, neutralMult);
        }

        // Check cuts depending on jet pseudorapdity
        if ( std::abs( recoJet.recoEta ) <= 2.7 ) {
            
            passNHF = ( recoJet.recoPfNHF < neutFracCut ) ? true : false;
            passNEF = ( recoJet.recoPfNEF < neutFracCut ) ? true : false;
            passNumOfConstituents = ( numberOfConstituents > 1 ) ? true : false;

            if (debug) {
                std::cout << Form("|eta|<=2.7 NHF < neutHadronFrac: %4.3f < %4.3f -> %s\n", 
                                  recoJet.recoPfNHF, neutFracCut, 
                                  (passNHF)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                std::cout << Form("|eta|<=2.7 NEF < neutEMFrac: %4.3f < %4.3f -> %s\n", 
                                  recoJet.recoPfNEF, neutFracCut, 
                                  (passNEF)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                std::cout << Form("|eta|<=2.7 numOfConstituents > 1: %d > 1 -> %s\n", 
                                  numberOfConstituents, 
                                  (passNumOfConstituents)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
            }

            if ( jetSelectionMethod == 2 ) { 
                passMuonFrac = ( recoJet.recoPfMUF < 0.8 ) ? true : false;
                if (debug) {
                    std::cout << Form("tightJetId |eta|<=2.7 MUF < 0.8: %4.3f < 0.8 -> %s\n", recoJet.recoPfMUF, (passMuonFrac)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                }
            } // if ( jetSelectionMethod == 2 )

            if( std::abs( recoJet.recoEta ) <= 2.4 ) {
                passChargedFrac = ( recoJet.recoPfCHF > 0 ) ? true : false;
                passChargedMult = ( chargedMult > 0 ) ? true : false;
                passChargedEmFrac = ( recoJet.recoPfCEF < chargedEmFracCut ) ? true : false;
                if (debug) {
                    std::cout << Form("|eta|<=2.4 CHF > 0: %4.3f > 0 -> %s\n", 
                                      recoJet.recoPfCHF, 
                                      (passChargedFrac)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                    std::cout << Form("|eta|<=2.4 chargedMult > 0: %d > 0 -> %s\n", 
                                      chargedMult, 
                                      (passChargedMult)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                    std::cout << Form("|eta|<=2.4 CEF < chargedEmFracCut: %4.3f < %4.3f -> %s\n", 
                                      recoJet.recoPfCEF, chargedEmFracCut, 
                                      (passChargedEmFrac)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                }
            } // if( std::abs( recoJetEta[jetIndex] ) <= 2.4 )

            

        } // if ( std::abs( recoJetEta[jetIndex] ) <= 2.7 )
        else if ( std::abs( recoJet.recoEta ) <= 3.0) {

            passNEF = ( recoJet.recoPfNEF > 0.01 ) ? true : false;
            passNHF = ( recoJet.recoPfNHF < 0.98 ) ? true : false;
            passNeutralMult = ( neutralMult > 2 ) ? true : false;

            if (debug) {
                std::cout << Form("|eta|<=3.0 NEF > 0.01: %4.3f > 0.01 -> %s\n", 
                                  recoJet.recoPfNEF, 
                                  (passNEF)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                std::cout << Form("|eta|<=3.0 NHF < 0.98: %4.3f < 0.98 -> %s\n", 
                                  recoJet.recoPfNHF, 
                                  (passNHF)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                std::cout << Form("|eta|<=3.0 neutralMult > 2: %d > 2 -> %s\n", 
                                  neutralMult, 
                                  (passNeutralMult)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
            }

        } // else if ( std::abs( recoJetEta[jetIndex] ) <= 3.0)
        else  {
            passNEF = ( recoJet.recoPfNEF < 0.9 ) ? true : false;
            passNeutralMult = (neutralMult > 10 ) ? true : false; // CAUTION: JET MET says it should be >10
            if (debug) {
                std::cout << Form("|eta|>3.0 NEF < 0.9: %4.3f < 0.9 -> %s\n", 
                                  recoJet.recoPfNEF, 
                                  (passNEF)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
                std::cout << Form("|eta|>3.0 neutralMult > 10: %d > 10 -> %s\n", 
                                  neutralMult, 
                                  (passNeutralMult)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
            }
        } // else 

        isGood = passNHF && passNEF && passNumOfConstituents && passMuonFrac && passChargedFrac && 
                 passChargedMult && passChargedEmFrac && passNeutralMult;

        if (debug) {
            std::cout << Form("Jet selection method: %d, isGood: %s\n", 
                              jetSelectionMethod, 
                              (isGood)?Form("%sPASS%s",GREEN,RESET):Form("%sFAIL%s",RED,RESET));
        }

        // std::cout << Form("Jet %d: passNHF: %d, passNEF: %d, passNumOfConstituents: %d, passMuonFrac: %d, passChargedFrac: %d, passChargedMult: %d, passChargedEmFrac: %d, passNeutralMult: %d \t [isGood: %d]\n",
        //                   jetIndex, passNHF, passNEF, passNumOfConstituents, passMuonFrac, passChargedFrac, 
        //                   passChargedMult, passChargedEmFrac, passNeutralMult, isGood);
    } // else if (jetSelectionMethod == 2 || jetSelectionMethod == 3)
    else {
        std::cerr << RED << Form("Invalid jet selection method: %d. No selection will be applied.", jetSelectionMethod) << RESET << std::endl;
        isGood = true;
    }

    return isGood;
}

//________________
/// Wrap the delta phi value to be within the range [-pi, pi]
/// @param phi: the delta phi value to be wrapped
void wrapDeltaPhi(float &phi) {
    while (phi > TMath::Pi()) phi -= TMath::TwoPi();
    while (phi <= -TMath::Pi()) phi += TMath::TwoPi();
}

//________________
/// Calculate the delta phi between two angles, wrapping the result to be within [-pi, pi]
/// @param phi1: the first angle
/// @param phi2: the second angle
/// @return the wrapped delta phi value
float deltaPhi(const float &phi1, const float &phi2) {
    float dPhi = phi1 - phi2;
    wrapDeltaPhi(dPhi);
    return dPhi;
}

//________________
// Create the histograms for the analysis
/// @param hs: the histograms object
/// @param isMc: true if the data is from Monte Carlo, false otherwise
void setPtHatRange(const int &ptHatSample) {
    std::cout << "Setting pT hat range for sample " << ptHatSample;
    if (ptHatSample == 15)       { ptHatRange[0] = 15.; ptHatRange[1] = 30.; }
    else if (ptHatSample == 30)  { ptHatRange[0] = 30.; ptHatRange[1] = 50.; }
    else if (ptHatSample == 50)  { ptHatRange[0] = 50.; ptHatRange[1] = 80.; }
    else if (ptHatSample == 80)  { ptHatRange[0] = 80.; ptHatRange[1] = 120.; }
    else if (ptHatSample == 120) { ptHatRange[0] = 120.; ptHatRange[1] = 170.; }
    else if (ptHatSample == 170) { ptHatRange[0] = 170.; ptHatRange[1] = 220.; }
    else if (ptHatSample == 220) { ptHatRange[0] = 220.; ptHatRange[1] = 280.; }
    else if (ptHatSample == 280) { ptHatRange[0] = 280.; ptHatRange[1] = 370.; }
    else if (ptHatSample == 370) { ptHatRange[0] = 370.; ptHatRange[1] = 460.; }
    else if (ptHatSample == 460) { ptHatRange[0] = 460.; ptHatRange[1] = 540.; }
    else if (ptHatSample == 540) { ptHatRange[0] = 540.; ptHatRange[1] = 8160.; }
    std::cout << Form(": [%.1f, %.1f)", ptHatRange[0], ptHatRange[1]) << std::endl;
}

//________________
/// Check if event passes selection criteria
/// @param isPbGoing: true if Pb is going in the positive direction, false otherwise
/// @param isMc: true if the event is from Monte Carlo, false otherwise
/// @param triggerId: the trigger ID (0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100)
/// @return true if the event is good, false otherwise
bool passesVertexFilter(const int &vertexFilterSelection) {
    if (vertexFilterSelection == 0) return pVertexFilterCutdz1p0 == 1;
    if (vertexFilterSelection == 1) return pVertexFilterCutGplus == 1;
    if (vertexFilterSelection == 2) return pVertexFilterCutVtx1 == 1;
    return false;
}

bool isGoodEvent(const bool &isPbGoing, const bool &isMc, const int &triggerId,
                 const int &vertexFilterSelection, bool isPp = false) {

    // MC bypasses HLT, but retains pPb skim cuts; pp uses only vz and HLT.
    // triggerId=0 means no explicit HLT requirement, not a checked MB trigger.
    // pPb sample cuts and weight lookup both use [low, high), assigning each
    // shared pT-hat boundary to the sample beginning at that value.
    bool goodTrigger{false};
    if (isMc) {
        if (!isPp && (ptHat < ptHatRange[0] || ptHat >= ptHatRange[1])) {
            return false;
        }
        goodTrigger = {true};
    }
    else {
        if (triggerId == 0) {
            goodTrigger = {true};
        }
        else if ( triggerId == 1 && HLT_PAAK4PFJet60_Eta5p1_v4 == 1 ) {
            goodTrigger = {true};
        }
        else if ( triggerId == 2 && HLT_PAAK4PFJet80_Eta5p1_v3 == 1 ) {
            goodTrigger = {true};
        }
        else if ( triggerId == 3 && HLT_PAAK4PFJet100_Eta5p1_v3 == 1 ) {
            goodTrigger = {true};
        }
    }

    if (isPp) return std::abs(vz) <= 15. && goodTrigger;

    return ( ( std::abs(vz) <= 15. ) && 
             ( pBeamScrapingFilter == 1 ) && 
             ( pPAprimaryVertexFilter == 1 ) && 
             ( HBHENoiseFilterResultRun2Loose == 1 ) && 
             ( phfCoincFilter == 1 ) && 
             passesVertexFilter(vertexFilterSelection) &&
             goodTrigger );
}

void fillOverweightHistograms(std::vector<GenJet> &genJets, std::vector<RecoJet> &recoJets, const double &weight, 
                              Histograms &hs, const bool &isPythia, bool &isGenOverweight, bool &isRecoOverweight) {
    
    // Tail veto acts on the whole event before configured reco jet ID and
    // before gen/reco fiducial dijet cuts. It is an MC sample-cleaning model,
    // not a detector efficiency; all downstream physics histograms inherit it.
    isGenOverweight = {false};
    isRecoOverweight = {false};

    double c0{0.}, c1{0.}, c2{0.}, c3{0.}, c4{0.};
    double d0{0.}, d1{0.}, d2{0.}, d3{0.}, d4{0.};

    float genDijetPtAve{0.};
    if ( genJets.size() >= 2) {
        genDijetPtAve = 0.5 * (genJets[0].pt + genJets[1].pt);
        hs.hGenDijetPtAveOverPtHatVsPtHat->Fill(ptHat, genDijetPtAve / ptHat, weight);
        hs.hGenLeadJetPtOverPtHatVsPtHat->Fill(ptHat, genJets[0].pt / ptHat, weight);
        if (isPythia) { // >0.5% overweight
            c0 = 1.1397348; c1 = 0.12215761; c2 = 0.0025413349; c3 = 0.1221576; c4 = 0.0025413349;
        }
        else { // >0.5% overweight
            c0 = 1.1333168; c1 = 0.12438097; c2 = 0.0023870616; c3 = 0.12438097; c4 = 0.0023870616;
        }
        
        if (genDijetPtAve / ptHat > (c0 + c1 * std::exp(-c2 * ptHat) + c3 * std::exp(-c4 * ptHat))) {
            isGenOverweight = {true};
        }
    }

    float recoDijetPtAve{0.};
    if ( recoJets.size() >= 1 ) {

        if (isPythia) { // >0.5% overweight
            c0 = 1.1486378; c1 = 0.11949062; c2 = 0.0025387593; c3 = 0.11949062; c4 = 0.0025387593;
            d0 = 1.2134125; d1 = 0.54101789; d2 = 0.0013676456; d3 = 0.32960156; d4 = 0.02549354;
        }
        else { // >0.5% overweight
            c0 = 1.135187; c1 = 1.7571184; c2 = 0.093903768; c3 = 0.25480253; c4 = 0.0022177174;
            d0 = 1.296563; d1 = 0.51835045; d2 = 0.0020646687; d3 = 0.99496387; d4 = 0.05675552;
        }

        hs.hRecoLeadJetPtOverPtHatVsPtHat->Fill(ptHat, recoJets[0].recoPt / ptHat, weight);
        if (recoJets[0].recoPt / ptHat > (d0 + d1 * std::exp(-d2 * ptHat) + d3 * std::exp(-d4 * ptHat))) {
            isRecoOverweight = {true};
        }

        if (recoJets.size() >= 2) {
            recoDijetPtAve = 0.5 * (recoJets[0].recoPt + recoJets[1].recoPt);
            hs.hRecoDijetPtAveOverPtHatVsPtHat->Fill(ptHat, recoDijetPtAve / ptHat, weight);
            if (recoDijetPtAve / ptHat > (c0 + c1 * std::exp(-c2 * ptHat) + c3 * std::exp(-c4 * ptHat))) {
                isRecoOverweight = {true};
            }
        }
    }

    if (!isGenOverweight && !isRecoOverweight) {
        if (genJets.size() >= 2) {
            hs.hGenDijetPtAveOverPtHatVsPtHatPass->Fill(ptHat, genDijetPtAve / ptHat, weight);
            hs.hGenLeadJetPtOverPtHatVsPtHatPass->Fill(ptHat, genJets[0].pt / ptHat, weight);
        }
        if (!recoJets.empty()) {
            hs.hRecoLeadJetPtOverPtHatVsPtHatPass->Fill(ptHat, recoJets[0].recoPt / ptHat, weight);
            if (recoJets.size() >= 2) {
                hs.hRecoDijetPtAveOverPtHatVsPtHatPass->Fill(ptHat, recoDijetPtAve / ptHat, weight);
            }
            
        }   
    }
}

//________________
void calculateScaleSmearing(float &ptDef, const float &eta) {
    float scaleFactor = 1.f;
    // if ( 30.f < ptDef && ptDef < 50.f ) {
    //     // Find bin for the given eta
    //     // scaleFactor = jesValues_pt_30_50(iBin)
    // }
    // else if ( 50.f <= ptDef && ptDef < 80.f ) {
    //     // scaleFactor = 1.01f; // 1% up
    // }
    // else if ( 80.f <= ptDef && ptDef < 120.f ) {
    //     // scaleFactor = 1.005f; // 0.5% up
    // }
    ptDef *= scaleFactor;
}

//________________
void processGenJets(const bool isPp, const bool &isPbGoing, const bool &isMc, const double &weight,
                    std::vector<GenJet> &genJets, Histograms &hs, const double &ptHat) {

    float genJetEtaLabFlipped{0.};
    float genJetEtaCM{0.};
    float pt{0.};
    
    // Loop over gen jets and fill histograms
    for (const auto &genJet : genJets) {
        pt = genJet.pt;
        genJetEtaLabFlipped = etaLab(genJet.eta, isPbGoing, isMc, isPp);
        genJetEtaCM = etaCM(genJet.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);

        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted->Fill(pt, genJet.eta, 1.);
        hs.hGenInclusiveJetPtEtaLabUnflipped->Fill(pt, genJet.eta, weight);
        hs.hGenInclusiveJetPtEtaLab->Fill(pt, genJetEtaLabFlipped, weight);
        hs.hGenInclusiveJetPtEtaCM->Fill(pt, genJetEtaCM, weight);

        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            genJetEtaCM = etaCM(genJet.eta, etaShift[iShift], isPbGoing, isMc, isPp);
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Fill( genJetEtaCM, 1.); // Unweighted
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Fill(genJetEtaCM, weight);
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        hs.hGenInclusiveJetPtEtaStdBins->Fill(pt, genJet.eta, weight);
        hs.hGenInclusiveJetPtPtHat->Fill(pt, ptHat, weight);
    } // for (const auto &genJet : genJets)

    // Let's do a little cleanup of the genJets vector
    // for (auto it = genJets.begin(); it != genJets.end(); ) {
    //     if ( it->pt < 20.f || std::abs(it->eta) > 3.0 ) {
    //         it = genJets.erase(it);  // returns iterator to next element
    //     } 
    //     else {
    //         ++it;
    //     }
    // } // for (auto it = genJets.begin(); it != genJets.end(); ) {
}

//________________
void processGenDijets(const bool isPp, const bool &isPbGoing, const bool &isMc, const double &weight,
                      std::vector<GenJet> &genJets, Histograms &hs,
                      JERSmearingHelper *fJERSmearingHelper) {

    // Must be at least 2 jets to form a dijet system
    if (genJets.size() < 2) return;

    // Perform initial sorting of gen jets by pT in descending order
    std::sort(genJets.begin(), genJets.end(), [](const GenJet &a, const GenJet &b) { return a.pt > b.pt; });

    const auto &leadingJet = genJets[0];
    const auto &subleadingJet = genJets[1];
    if (leadingJet.pt >= 50.f && subleadingJet.pt >= 40.) {
        float dphi = deltaPhi(leadingJet.phi, subleadingJet.phi);
        if (std::abs(dphi) >= TMath::TwoPi() / 3.) {
            float dijetPtAve = 0.5f * (leadingJet.pt + subleadingJet.pt);
            float leadEtaLabUnflipped = leadingJet.eta;
            float leadEtaLabFlipped = etaLab(leadingJet.eta, isPbGoing, isMc, isPp);
            float leadEtaCM = etaCM(leadingJet.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
            float subleadEtaLabUnflipped = subleadingJet.eta;
            float subleadEtaLabFlipped = etaLab(subleadingJet.eta, isPbGoing, isMc, isPp);
            float subleadEtaCM = etaCM(subleadingJet.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);

            if (std::abs(leadEtaCM) <= 1.9 && std::abs(subleadEtaCM) <= 1.9) {
                float dijetEtaLabUnflipped = 0.5f * (leadEtaLabUnflipped + subleadEtaLabUnflipped);
                float dijetEtaLabFlipped = 0.5f * (leadEtaLabFlipped + subleadEtaLabFlipped);
                float dijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

                hs.hGenDijetPtEtaLabUnflipped->Fill(dijetPtAve, dijetEtaLabUnflipped, weight);
                hs.hGenDijetPtEtaLab->Fill(dijetPtAve, dijetEtaLabFlipped, weight);
                hs.hGenDijetPtEtaCM->Fill(dijetPtAve, dijetEtaCM, weight);
            } // if (std::abs(leadEtaCM) <= 1.9 && std::abs(subleadEtaCM) <= 1.9)

            // Verify forward/backward dijet etaCM bias due to the eta shift
            for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
                float leadingEtaCMShifted = etaCM(leadingJet.eta, etaShift[iShift], isPbGoing, isMc, isPp);
                float subleadingEtaCMShifted = etaCM(subleadingJet.eta, etaShift[iShift], isPbGoing, isMc, isPp);
                float dijetEtaCMShifted = 0.5 * (leadingEtaCMShifted + subleadingEtaCMShifted);
                
                if (std::abs(leadingEtaCMShifted) > 1.9 || std::abs(subleadingEtaCMShifted) > 1.9) continue;

                // if (dijetPtAve < 60. || dijetPtAve > 80.) continue;
                hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Fill(dijetEtaCMShifted);
                hs.hGenDijetEtaCMShifted[iShift]->Fill(dijetEtaCMShifted, weight);
            } // for (int iShift = 0; iShift < nEtaShifts; ++iShift)
        } // if (std::abs(dphi) >= TMath::TwoPi() / 3.)
    } // if (leadingJet.pt >= 50.f && subleadingJet.pt >= 40.)

    // Re-select the leading pair for each pT accessor so nominal and smeared
    // generator distributions use the same thresholds and angular selection.
    std::vector<GenJet> genJetSortFallback;
    genJetSortFallback.reserve(genJets.size());
    auto fillGenDijetSystematic = [&](auto getPt,
                                      std::unique_ptr<TH2D> *genFullArr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(genJets, getPt, genJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return;
        float systDphi = deltaPhi(leadingJetSyst.phi, subleadingJetSyst.phi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaCM = etaCM(leadingJetSyst.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        float subleadEtaCM = etaCM(subleadingJetSyst.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

#ifdef DO_MIXED_PRIOR
        auto gaussianRatio = [](double eta, double targetMean, double targetSigma) {
            const double nominalPull = (eta - nominalPriorMean) / nominalPriorSigma;
            const double targetPull = (eta - targetMean) / targetSigma;
            return (nominalPriorSigma / targetSigma) *
                   std::exp(0.5 * (nominalPull * nominalPull - targetPull * targetPull));
        };
        const double mixedPriorWeight = 0.5 * gaussianRatio(
                                                systDijetEtaCM, leftPriorMean, leftPriorSigma) +
                                        0.5 * gaussianRatio(
                                                systDijetEtaCM, rightPriorMean, rightPriorSigma);
#endif

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;

            genFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
            if (genFullArr == hs.hGenDijetPtEtaCMArr) {
                // LO massless two-jet estimators, not generator-record Bjorken
                // x: higher-order radiation, jet masses and eta vs y enter as
                // approximations. For pp, the historical xPb name means the
                // negative-z proton estimator; there is no Pb beam.
                // Both supported collision systems use sqrt(sNN) = 8.16 TeV.
                const double xp = (leadingJetSyst.pt * std::exp(leadEtaCM) +
                                   subleadingJetSyst.pt * std::exp(subleadEtaCM)) / 8160.;
                const double xPb = (leadingJetSyst.pt * std::exp(-leadEtaCM) +
                                    subleadingJetSyst.pt * std::exp(-subleadEtaCM)) / 8160.;
                const double valuesPb[3] = {systDijetPtAve, systDijetEtaCM, xPb};
                const double valuesP[3] = {systDijetPtAve, systDijetEtaCM, xp};
                const double valuesRatio[3] = {systDijetPtAve, systDijetEtaCM, xPb / xp};
                hs.hGenDijetPtEtaCMXPbArr[iCut]->Fill(valuesPb, weight);
                hs.hGenDijetPtEtaCMXpArr[iCut]->Fill(valuesP, weight);
                hs.hGenDijetPtEtaCMXPbOverXpArr[iCut]->Fill(valuesRatio, weight);
                auto *sampleArr = isFillControlledSample
                    ? hs.hGenDijetPtEtaCMControlledSampleArr
                    : hs.hGenDijetPtEtaCMStudiedSampleArr;
                sampleArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
#ifdef DO_MIXED_PRIOR
                hs.hGenDijetPtEtaCMMixedPriorArr[iCut]->Fill(
                    systDijetPtAve, systDijetEtaCM, weight * mixedPriorWeight);
#endif
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };

    // Fill dijet histograms for systematic variations in lab frame (etaLab)
    auto fillGenDijetLabSystematic = [&](auto getPt,
                                         std::unique_ptr<TH2D> *genFullArr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(genJets, getPt, genJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        // pT cuts for systematic variations (can be adjusted as needed)
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return;
        float systDphi = deltaPhi(leadingJetSyst.phi, subleadingJetSyst.phi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaLab = etaLab(leadingJetSyst.eta, isPbGoing, isMc, isPp);
        float subleadEtaLab = etaLab(subleadingJetSyst.eta, isPbGoing, isMc, isPp);
        float systDijetEtaLab = 0.5f * (leadEtaLab + subleadEtaLab);

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            if (std::abs(leadEtaLab) > etaCuts[iCut] || std::abs(subleadEtaLab) > etaCuts[iCut]) continue;
            genFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaLab, weight);
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };
    
    //
    // Gen level dijet histograms (no smearing, no extra scale factors)
    //
    fillGenDijetSystematic(
        [](const GenJet &jet) { return jet.pt; },
        hs.hGenDijetPtEtaCMArr
    );

    fillGenDijetLabSystematic(
        [](const GenJet &jet) { return jet.pt; },
        hs.hGenDijetPtEtaLabArr
    );

    //
    // Smear gen jets and fill smeared-gen histograms for several extraScale factors.
    //

    // Number of smearing iterations (runs) per scale.
    // Diagnostic replicas carry the full event weight (no 1/nSmearRuns).
    // Dijet spectra receive 20 fills/event, inclusive diagnostics only 10.
    // Compare normalized shapes or account explicitly for these multiplicities;
    // Sumw2 alone does not encode correlations between same-event replicas.
    int nSmearRuns = 20; // adjust as needed

    // Scale index map: 0: 1.0f (Default), 1: 2.0f (Double), 2: 0.9f (Ninety), 3: 1.0f (Extra/eta-dependent)
    const std::array<float, 4> extraScales = {1.0f, 2.0f, 0.9f, 1.0f};

    for (int iScale = 0; iScale < static_cast<int>(extraScales.size()); ++iScale) {
        float extraScale = extraScales[iScale];
        // Choose the target histogram arrays based on the scale index
        std::unique_ptr<TH2D> *fullArr{nullptr};
        float genJetEtaLabFlipped = {-999.f};

        if (iScale == 0) {  // x1.0f (Default)
            fullArr = hs.hGenDijetDefPtEtaCMArr;
        }
        else if (iScale == 1) {  // x2.0f (Double)
            fullArr = hs.hGenDijetDefDoublePtEtaCMArr;
        }
        else if (iScale == 2) {  // x0.9f (Ninety)
            fullArr = hs.hGenDijetDefNinetyPtEtaCMArr;
        }
        else if (iScale == 3) {  // x1.0f (Extra/eta-dependent)
            fullArr = hs.hGenDijetDefExtraPtEtaCMArr;
        }

        // If any of the target histogram arrays weren't created, skip this scale
        if (!fullArr) continue;

        for (int iRun = 0; iRun < nSmearRuns; ++iRun) {
            // Apply smearing to each gen jet and store in smearedPt
            for (size_t iJet = 0; iJet < genJets.size(); ++iJet) {
                auto &jet = genJets[iJet];
                jet.smearedPt = jet.pt;
                genJetEtaLabFlipped = etaLab(jet.eta, isPbGoing, isMc, isPp);

                // std::cout << "gen jet pt: " << jet.pt << ", eta: " << jet.eta << ", smeared pt: " << jet.smearedPt;

                const double normal = fJERSmearingHelper->normalDeviate(
                    jet.pt, jet.eta, eventRho, run, lumi, evt,
                    iRun * static_cast<int>(genJets.size()) + static_cast<int>(iJet));
                jet.smearedPt = fJERSmearingHelper->smearTruthMomentum(
                    jet.pt, jet.eta, extraScale, iScale == 3, normal);

                // std::cout << ", smeared pt after smearing: " << jet.smearedPt << std::endl;

                // Fill histograms for the first 10 smearing runs to avoid overfilling
                if (iRun > 9) continue;

                if (iScale == 0 && hs.hGenInclusiveJetJESDefPtEta) {  // x1.0f (Default)
                    hs.hGenInclusiveJetJESDefPtEta->Fill(jet.smearedPt / jet.pt, jet.pt, jet.eta, weight);
                    hs.hGenInclusiveJetDefPtEtaFlipped->Fill(genJetEtaLabFlipped, jet.smearedPt, weight);
                }
                else if (iScale == 1 && hs.hGenInclusiveJetJESDefDoublePtEta) {  // x2.0f (Double)
                    hs.hGenInclusiveJetJESDefDoublePtEta->Fill(jet.smearedPt / jet.pt, jet.pt, jet.eta, weight);
                    hs.hGenInclusiveJetDefDoublePtEtaFlipped->Fill(genJetEtaLabFlipped, jet.smearedPt, weight);
                }
                else if (iScale == 2 && hs.hGenInclusiveJetJESDefNinetyPtEta) {  // x0.9f (Ninety)
                    hs.hGenInclusiveJetJESDefNinetyPtEta->Fill(jet.smearedPt / jet.pt, jet.pt, jet.eta, weight);
                    hs.hGenInclusiveJetDefNinetyPtEtaFlipped->Fill(genJetEtaLabFlipped, jet.smearedPt, weight);
                }
                else if (iScale == 3 && hs.hGenInclusiveJetJESDefExtraPtEta) {  // x1.0f (Extra/eta-dependent)
                    hs.hGenInclusiveJetJESDefExtraPtEta->Fill(jet.smearedPt / jet.pt, jet.pt, jet.eta, weight);
                    hs.hGenInclusiveJetDefExtraPtEtaFlipped->Fill(genJetEtaLabFlipped, jet.smearedPt, weight);
                }
            } // for (auto &jet : genJets)

            // Fill histograms using smeared pT accessor
            fillGenDijetSystematic(
                [](const GenJet &jet) { return jet.smearedPt; },
                fullArr
            );
        } // for runs
    } // for iScale
}

//________________
void prepareUnfolding(const bool isPp, const bool &isPbGoing, const bool &isMc, const float &weight,
                      std::vector<GenJet> &genJets, std::vector<RecoJet> &recoJets, 
                      std::unique_ptr<THnSparseF> *responseArr,
                      std::unique_ptr<TH2D> *missArr,
                      std::unique_ptr<TH2D> *fakeArr,
                      std::unique_ptr<TH1D> *classificationArr,
                      std::unique_ptr<THnSparseF> *controlledResponseArr = nullptr,
                      std::unique_ptr<TH2D> *controlledMissArr = nullptr,
                      std::unique_ptr<TH2D> *controlledFakeArr = nullptr,
                      std::unique_ptr<THnSparseF> *studiedResponseArr = nullptr,
                      std::unique_ptr<TH2D> *studiedMissArr = nullptr,
                      std::unique_ptr<TH2D> *studiedFakeArr = nullptr,
                      int jerVariation = -1,
                      UnfoldingPairHistograms *pairHistograms = nullptr,
                      UnfoldingPairHistograms *controlledPairHistograms = nullptr,
                      UnfoldingPairHistograms *studiedPairHistograms = nullptr,
                      ExtendedUnfoldingHistograms *extendedHistograms = nullptr,
                      ExtendedUnfoldingHistograms *controlledExtendedHistograms = nullptr,
                      ExtendedUnfoldingHistograms *studiedExtendedHistograms = nullptr) {

    if (!isMc) return; // Only prepare unfolding for MC
    
    // Lambda function to compute the squared delta R between two jets
    auto deltaR2 = [&](float eta1, float phi1, float eta2, float phi2) {
        float dEta = eta1 - eta2;
        float dPhi = deltaPhi(phi1, phi2);
        return dEta * dEta + dPhi * dPhi;
    };

    // Lambda function to check if a pair of jets passes the dijet selection criteria
    auto goodPair = [](float leadPt, float subleadPt, float leadEtaCM, float subleadEtaCM, float dphi, float etaCut) {
        return analysis_contract::passesDijetSelection(
            leadPt, subleadPt, leadEtaCM, subleadEtaCM, dphi, etaCut);
    };

    // Sort gen jets by pT in descending order to identify leading and subleading jets
    std::sort(genJets.begin(), genJets.end(), [&](const GenJet &a, const GenJet &b) {
            return a.pt > b.pt;
    });

    // Cached momenta and stable draw identities travel with jets through sorting.
    // The default accessor preserves direct helper calls with refPtSmeared.
    auto measuredPt = [jerVariation](const RecoJet &jet) {
        return jerVariation < 0 ? jet.refPtSmeared : jet.jerMomenta[jerVariation];
    };
    // Sort reco jets by the caller-prepared smeared pT so the selected measured
    // dijet matches the response, miss, fake, and classification outputs for
    // the requested JER variation.
    std::sort(recoJets.begin(), recoJets.end(), [&](const RecoJet &a, const RecoJet &b) {
        return measuredPt(a) > measuredPt(b);
    });

    //
    // Process gen-level (truth) leading and subleading jets and dijet properties
    //

    GenJet genLeadingJet{};
    GenJet genSubleadingJet{};

    float genLeadingJetPt = {-99.f};
    float genLeadingJetEtaCM = {-99.f};
    float genLeadingJetPhi = {0.f};

    float genSubleadingJetPt = {-99.f};
    float genSubleadingJetEtaCM = {-99.f};
    float genSubleadingJetPhi = {0.f};

    float genDijetPtAve = {-99.f};
    float genDijetEtaCM = {-99.f};
    float genDijetDeltaPhi = {0.f};

    bool has2genJets = (genJets.size() >= 2);

    // Truth-level (gen) leading and subleading jets and dijet properties
    // If both leading and subleading gen jets pass the selection criteria, 
    // compute their properties and call truth distribution
    if (has2genJets) {
        genLeadingJet = genJets[0];
        genSubleadingJet = genJets[1];

        genLeadingJetPt = genLeadingJet.pt;
        genLeadingJetEtaCM = etaCM(genLeadingJet.eta,
                                   analysis_contract::kNominalEtaShift,
                                   isPbGoing, isMc, isPp);
        genLeadingJetPhi = genLeadingJet.phi;

        genSubleadingJetPt = genSubleadingJet.pt;
        genSubleadingJetEtaCM = etaCM(genSubleadingJet.eta,
                                      analysis_contract::kNominalEtaShift,
                                      isPbGoing, isMc, isPp);
        genSubleadingJetPhi = genSubleadingJet.phi;

        genDijetPtAve = 0.5 * (genLeadingJetPt + genSubleadingJetPt);
        genDijetEtaCM = 0.5 * (genLeadingJetEtaCM + genSubleadingJetEtaCM);
        genDijetDeltaPhi = deltaPhi(genLeadingJetPhi, genSubleadingJetPhi);
    } // if (has2genJets)

    //
    // Process reco-level (reconstructed) leading and subleading jets and dijet properties
    //
    RecoJet recoLeadingJet;
    RecoJet recoSubleadingJet;

    bool directMatching = {false};
    bool swapMatching = {false};
    bool goodMatching = {false};
    bool hasRecoReferenceMatches = {false};
    int nRecoReferenceMatches = {0};

    bool has2recoJets = (recoJets.size() >= 2);
    float recoLeadingJetPt = {-99.f};
    float recoLeadingJetEtaCM = {-99.f};
    float recoLeadingJetPhi = {0.f};

    float recoSubleadingJetPt = {-99.f};
    float recoSubleadingJetEtaCM = {-99.f};
    float recoSubleadingJetPhi = {0.f};

    float recoDijetPtAve = {-99.f};
    float recoDijetEtaCM = {-99.f};
    float recoDijetDeltaPhi = {0.f};

    // If both leading and subleading reco jets pass the selection criteria, 
    // compute their properties and call reco (measured) distribution
    if (has2recoJets) {
        // Get leading and subleading reco jets
        recoLeadingJet = recoJets[0];
        recoSubleadingJet = recoJets[1];

        recoLeadingJetPt = measuredPt(recoLeadingJet);
        recoLeadingJetEtaCM = etaCM(recoLeadingJet.recoEta,
                                    analysis_contract::kNominalEtaShift,
                                    isPbGoing, isMc, isPp);
        recoLeadingJetPhi = recoLeadingJet.recoPhi;

        recoSubleadingJetPt = measuredPt(recoSubleadingJet);
        recoSubleadingJetEtaCM = etaCM(recoSubleadingJet.recoEta,
                                       analysis_contract::kNominalEtaShift,
                                       isPbGoing, isMc, isPp);
        recoSubleadingJetPhi = recoSubleadingJet.recoPhi;

        recoDijetPtAve = 0.5 * (recoLeadingJetPt + recoSubleadingJetPt);
        recoDijetEtaCM = 0.5 * (recoLeadingJetEtaCM + recoSubleadingJetEtaCM);
        recoDijetDeltaPhi = deltaPhi(recoLeadingJetPhi, recoSubleadingJetPhi);

        nRecoReferenceMatches = static_cast<int>(recoLeadingJet.refPt > 0.f) +
                                static_cast<int>(recoSubleadingJet.refPt > 0.f);
        hasRecoReferenceMatches = nRecoReferenceMatches == 2;

        // Define the deltaR cut for matching
        if (has2genJets && hasRecoReferenceMatches) {
            constexpr float deltaRCut2 = analysis_contract::kDijetMatchDeltaR *
                                         analysis_contract::kDijetMatchDeltaR;
            directMatching = (deltaR2(recoLeadingJet.recoEta, recoLeadingJetPhi, genLeadingJet.eta, genLeadingJetPhi) < deltaRCut2) &&
                             (deltaR2(recoSubleadingJet.recoEta, recoSubleadingJetPhi, genSubleadingJet.eta, genSubleadingJetPhi) < deltaRCut2);
            swapMatching = (deltaR2(recoLeadingJet.recoEta, recoLeadingJetPhi, genSubleadingJet.eta, genSubleadingJetPhi) < deltaRCut2) &&
                           (deltaR2(recoSubleadingJet.recoEta, recoSubleadingJetPhi, genLeadingJet.eta, genLeadingJetPhi) < deltaRCut2);
            goodMatching = directMatching || swapMatching;
        } // if (has2genJets && hasRecoReferenceMatches)
    } // if (has2recoJets)

    // Existence of an injective two-jet correspondence. Requiring distinct
    // indices prevents one jet inside both cones from counting as two jets.
    auto distinctCounterparts = [](const auto &jets, auto matchMask) {
        int firstA = -1, firstB = -1;
        int countA = 0, countB = 0;
        for (size_t i = 0; i < jets.size(); ++i) {
            const int mask = matchMask(jets[i]);
            if (mask & 1) { if (firstA < 0) firstA = static_cast<int>(i); ++countA; }
            if (mask & 2) { if (firstB < 0) firstB = static_cast<int>(i); ++countB; }
        }
        return countA > 0 && countB > 0 &&
               (firstA != firstB || countA > 1 || countB > 1);
    };
    constexpr float matchCut2 = analysis_contract::kDijetMatchDeltaR *
                                analysis_contract::kDijetMatchDeltaR;
    const bool truthRecoFound = has2genJets && distinctCounterparts(recoJets, [&](const RecoJet &jet) {
        if (jet.refPt <= 0.f) return 0;
        return (deltaR2(jet.recoEta, jet.recoPhi, genLeadingJet.eta, genLeadingJet.phi) < matchCut2 ? 1 : 0) |
               (deltaR2(jet.recoEta, jet.recoPhi, genSubleadingJet.eta, genSubleadingJet.phi) < matchCut2 ? 2 : 0);
    });
    const bool recoGenFound = has2recoJets && hasRecoReferenceMatches &&
        distinctCounterparts(genJets, [&](const GenJet &jet) {
            return (deltaR2(jet.eta, jet.phi, recoLeadingJet.recoEta, recoLeadingJet.recoPhi) < matchCut2 ? 1 : 0) |
                   (deltaR2(jet.eta, jet.phi, recoSubleadingJet.recoEta, recoSubleadingJet.recoPhi) < matchCut2 ? 2 : 0);
        });

    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

        // Explicitly list the conditions for each category to avoid ambiguity 
        // and ensure clarity in the unfolding diagnostics

        // Gen (truth) dijet passes selection
        const bool genPairPass = has2genJets && goodPair(
            genLeadingJetPt, genSubleadingJetPt,
            genLeadingJetEtaCM, genSubleadingJetEtaCM,
            genDijetDeltaPhi, etaCuts[iCut]);
        // Reco (measured) dijet passes selection
        const bool recoPairPass = has2recoJets && goodPair(
            recoLeadingJetPt, recoSubleadingJetPt,
            recoLeadingJetEtaCM, recoSubleadingJetEtaCM,
            recoDijetDeltaPhi, etaCuts[iCut]);

        // Event-level response, deliberately independent of jet matching.
        // Truth includes every available global gen-leading pair, not only
        // pairs that produced selected reco; otherwise inefficiency is lost.
        auto fillExtended = [&](ExtendedUnfoldingHistograms *extended) {
            if (!extended) return;
            const double truthPoint[3] = {genDijetPtAve, genDijetEtaCM, genPairPass ? 2. : 1.};
            if (has2genJets) {
                extended->truth[iCut]->Fill(truthPoint, weight);
                if (recoPairPass) {
                    const double responsePoint[5] = {genDijetPtAve, genDijetEtaCM,
                        recoDijetPtAve, recoDijetEtaCM, truthPoint[2]};
                    extended->response[iCut]->Fill(responsePoint, weight);
                    extended->matchedTruth[iCut]->Fill(truthPoint, weight);
                    extended->matchedMeasured[iCut]->Fill(recoDijetPtAve, recoDijetEtaCM, weight);
                } else {
                    extended->miss[iCut]->Fill(truthPoint, weight);
                }
            }
            if (recoPairPass) {
                extended->measured[iCut]->Fill(recoDijetPtAve, recoDijetEtaCM, weight);
                if (!has2genJets) extended->fake[iCut]->Fill(recoDijetPtAve, recoDijetEtaCM, weight);
            }
        };
        fillExtended(extendedHistograms);
        fillExtended(isFillControlledSample ? controlledExtendedHistograms : studiedExtendedHistograms);

        // Check if both gen and reco dijets pass selection
        const bool bothPairsPass = genPairPass && recoPairPass;
        // Check if the reco dijet matches the gen dijet (either direct or swapped)
        const bool matchedPair = bothPairsPass && goodMatching;
        // Exclusive subdivisions of the existing miss/fake definitions.
        // "Found" includes counterparts displaced from the leading pair.
        auto fillPairDiagnostics = [&](UnfoldingPairHistograms *pairs) {
            if (!pairs) return;
            if (matchedPair) {
                pairs->hGenDijetPtEtaCMMatchedRecoSelected[iCut]->Fill(genDijetPtAve, genDijetEtaCM, weight);
                pairs->hRecoDijetPtEtaCMMatchedGenSelected[iCut]->Fill(recoDijetPtAve, recoDijetEtaCM, weight);
            } else {
                if (genPairPass) {
                    auto &histogram = truthRecoFound
                        ? pairs->hGenDijetPtEtaCMSelectedRecoFound[iCut]
                        : pairs->hGenDijetPtEtaCMSelectedRecoNotFound[iCut];
                    histogram->Fill(genDijetPtAve, genDijetEtaCM, weight);
                }
                if (recoPairPass) {
                    auto &histogram = recoGenFound
                        ? pairs->hRecoDijetPtEtaCMMatchedGenFound[iCut]
                        : pairs->hRecoDijetPtEtaCMSelectedGenNotFound[iCut];
                    histogram->Fill(recoDijetPtAve, recoDijetEtaCM, weight);
                }
            }
        };
        fillPairDiagnostics(pairHistograms);
        fillPairDiagnostics(isFillControlledSample ? controlledPairHistograms : studiedPairHistograms);
        // Check if the composition of the reco dijet has changed compared to the gen dijet
        const bool pairCompositionChanged = bothPairsPass &&
                                            hasRecoReferenceMatches &&
                                            !goodMatching;
        // Check if one or both reco jets are missing a reference match
        const bool missingOneRecoReference = bothPairsPass &&
                                             nRecoReferenceMatches == 1;
        // Check if both reco jets are missing reference matches
        const bool missingBothRecoReferences = bothPairsPass &&
                                               nRecoReferenceMatches == 0;

        // Determine the category for unfolding diagnostics
        int pairCategory = 0;
        if (matchedPair) {
            pairCategory = directMatching
                ? unfolding_diagnostics::kMatchedDirect
                : unfolding_diagnostics::kMatchedSwapped;
        }
        else if (pairCompositionChanged) {
            pairCategory = unfolding_diagnostics::kSelectedPairMismatch;
        }
        else if (missingOneRecoReference) {
            pairCategory = unfolding_diagnostics::kMissingOneRecoReference;
        }
        else if (missingBothRecoReferences) {
            pairCategory = unfolding_diagnostics::kMissingBothRecoReferences;
        }
        else if (genPairPass) {
            pairCategory = unfolding_diagnostics::kGenPassRecoFail;
        }
        else if (recoPairPass) {
            pairCategory = unfolding_diagnostics::kRecoPassGenFail;
        }

        if (pairCategory != 0 && classificationArr) {
            classificationArr[iCut]->Fill(pairCategory, weight);
        }

        // Axis contract: [gen pTave, gen eta, reco pTave, reco eta]. Both
        // selected pairs must match directly or with swapped leading roles.
        // A selected but mismatched pair contributes BOTH a truth miss and a
        // reco fake, while classification gets exactly one category entry.
        // Fill the response matrix for matched pairs, and fill the appropriate 
        // histograms for gen-only or reco-only pairs
        if (matchedPair) {
            const double response[4] = {
                genDijetPtAve, genDijetEtaCM, recoDijetPtAve, recoDijetEtaCM
            };
            responseArr[iCut]->Fill(response, weight);
            if (controlledResponseArr && studiedResponseArr) {
                auto *sampleResponseArr = isFillControlledSample
                    ? controlledResponseArr
                    : studiedResponseArr;
                sampleResponseArr[iCut]->Fill(response, weight);
            }
        }
        else if (genPairPass) {
            missArr[iCut]->Fill(
                genDijetPtAve, genDijetEtaCM, weight);
            if (controlledMissArr && studiedMissArr) {
                auto *sampleMissArr = isFillControlledSample
                    ? controlledMissArr
                    : studiedMissArr;
                sampleMissArr[iCut]->Fill(genDijetPtAve, genDijetEtaCM, weight);
            }
        }

        if (recoPairPass && !matchedPair) {
            fakeArr[iCut]->Fill(
                recoDijetPtAve, recoDijetEtaCM, weight);
            if (controlledFakeArr && studiedFakeArr) {
                auto *sampleFakeArr = isFillControlledSample
                    ? controlledFakeArr
                    : studiedFakeArr;
                sampleFakeArr[iCut]->Fill(recoDijetPtAve, recoDijetEtaCM, weight);
            }
        }

    } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
}


//________________
void processRecoJets(const bool isPp, const bool &isPbGoing, const bool &isMc, const double &weight,
                     std::vector<RecoJet> &recoJets, Histograms &hs,
                     JetUncertainty &jeu,
                     const int &jetSelectionMethod, const double &ptHat) {
    
    float recoJetEtaLabFlipped{0.};
    float recoJetEtaCM{0.};
    float refJetEtaLabFlipped{0.};
    float refJetEtaCM{0.};
    bool  hasMatchingGenJet{false};

    // Loop over reconstructed jets
    for (auto it = recoJets.begin(); it != recoJets.end(); ) {

        // Reference to the current reco jet
        auto &recoJet = *it;

        hasMatchingGenJet = (recoJet.refPt > 0.f);
        recoJetEtaLabFlipped = etaLab(recoJet.recoEta, isPbGoing, isMc, isPp);
        recoJetEtaCM = etaCM(recoJet.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        if (hasMatchingGenJet) {
            refJetEtaLabFlipped = etaLab(recoJet.refEta, isPbGoing, isMc, isPp);
            refJetEtaCM = etaCM(recoJet.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        }

        // Acceptance for reco jets without any selection applied
        if ( isGoodRecoJet(recoJet, 0) ) {
            hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.recoEta, weight);
            hs.hRecoInclusiveJetNoSelPtEtaLab->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            hs.hRecoInclusiveJetNoSelPtEtaCM->Fill(recoJet.recoPt, recoJetEtaCM, weight);
            if (hasMatchingGenJet) {
                hs.hRefInclusiveJetNoSelPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.refEta, weight);
                hs.hRefInclusiveJetNoSelPtEtaLab->Fill(recoJet.recoPt, refJetEtaLabFlipped, weight);
                hs.hRefInclusiveJetNoSelPtEtaCM->Fill(recoJet.recoPt, refJetEtaCM, weight);
                hs.hRefSelInclusiveJetNoSelPtEtaLabUnflipped->Fill(recoJet.refPt, recoJet.refEta, weight);
                hs.hRefSelInclusiveJetNoSelPtEtaLab->Fill(recoJet.refPt, refJetEtaLabFlipped, weight);
                hs.hRefSelInclusiveJetNoSelPtEtaCM->Fill(recoJet.refPt, refJetEtaCM, weight);
            }
        }
        // Acceptance for reco jets with trkMax selection applied
        if ( isGoodRecoJet(recoJet, 1) ) {
            hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.recoEta, weight);
            hs.hRecoInclusiveJetTrkMaxPtEtaLab->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            hs.hRecoInclusiveJetTrkMaxPtEtaCM->Fill(recoJet.recoPt, recoJetEtaCM, weight);
            if (hasMatchingGenJet) {
                hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.refEta, weight);
                hs.hRefInclusiveJetTrkMaxPtEtaLab->Fill(recoJet.recoPt, refJetEtaLabFlipped, weight);
                hs.hRefInclusiveJetTrkMaxPtEtaCM->Fill(recoJet.recoPt, refJetEtaCM, weight);
                hs.hRefSelInclusiveJetTrkMaxPtEtaLabUnflipped->Fill(recoJet.refPt, recoJet.refEta, weight);
                hs.hRefSelInclusiveJetTrkMaxPtEtaLab->Fill(recoJet.refPt, refJetEtaLabFlipped, weight);
                hs.hRefSelInclusiveJetTrkMaxPtEtaCM->Fill(recoJet.refPt, refJetEtaCM, weight);
            }
        }
        // This inclusive family always uses method 2, independently of the
        // runtime method used to prune the collection for dijets below.
        // Acceptance for reco jets with jetId selection applied
        if ( isGoodRecoJet(recoJet, 2) ) {
            hs.hRecoInclusiveJetPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.recoEta,  weight);
            hs.hRecoInclusiveJetPtEtaLab->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            hs.hRecoInclusiveJetPtEtaCM->Fill(recoJet.recoPt, recoJetEtaCM, weight);
            if (hasMatchingGenJet) {
                hs.hRefInclusiveJetPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.refEta, weight);
                hs.hRefInclusiveJetPtEtaLab->Fill(recoJet.recoPt, refJetEtaLabFlipped, weight);
                hs.hRefInclusiveJetPtEtaCM->Fill(recoJet.recoPt, refJetEtaCM, weight);
                hs.hRefSelInclusiveJetPtEtaLabUnflipped->Fill(recoJet.refPt, recoJet.refEta, weight);
                hs.hRefSelInclusiveJetPtEtaLab->Fill(recoJet.refPt, refJetEtaLabFlipped, weight);
                hs.hRefSelInclusiveJetPtEtaCM->Fill(recoJet.refPt, refJetEtaCM, weight);
            }   
        }

        // Remove reco jets that do not pass the selection criteria based on the specified method
        if ( !isGoodRecoJet(recoJet, jetSelectionMethod) ) {
            it = recoJets.erase(it);  // returns iterator to next element
            continue;
        } 
        else {
            ++it;
        }

        // Fill after the configured jet selection so the raw and corrected
        // histograms contain the same selected jet population.
        hs.hRecoInclusiveJetRawPtEtaLabUnflipped->Fill(recoJet.recoPtRaw, recoJet.recoEta, weight);

        // std::cout << Form("Jet: eta = %.2f, phi = %.2f, rawPt = %.1f, recoPt = %.1f", 
        //                   recoJet.recoEta, recoJet.recoPhi, recoJet.recoRawPt, recoJet.recoPt) << std::endl;
        hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted->Fill(recoJet.recoPt, recoJet.recoEta, 1.);
        hs.hRecoInclusiveJetPtEtaStdBins->Fill(recoJet.recoPt, recoJet.recoEta, weight);

        // for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
        //     recoJetEtaCM = etaCM(recoJet.recoEta, etaShift[iShift], isPbGoing, isMc, isPp);
        //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Fill( recoJetEtaCM, 1.); // Unweighted
        //     hs.hRecoInclusiveJetEtaCMShifted[iShift]->Fill(recoJetEtaCM, weight);
        // } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        // Fill histograms for MC jets with and without matching gen jets
        if (isMc) {
            if (hasMatchingGenJet) {
                // Fill histograms for reco jets with matching gen jets
                hs.hRecoInclusiveJetPtPtHat->Fill(recoJet.recoPt, ptHat, weight);
                // Fill histograms for reco jets with matching gen jets (efficiency calculation)
                hs.hRecoInclusiveJetPtEtaLabMatched->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
                // Fill histogram for JES and JER extraction
                hs.hRecoInclusiveJetJESPtEta->Fill(recoJet.recoPt / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);

                // Fill reference jet histograms
                hs.hRefInclusiveJetPtEtaStdBins->Fill(recoJet.refPt, recoJet.refEta, weight);
            } // if (hasMatchingGenJet)
            else {
                // Fill histograms for reco jets without matching gen jets (fake rate calculation)
                hs.hRecoInclusiveJetPtEtaLabUnmatched->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            }
        } // if (isMc )

        // Perform JEU calculations
        if ( !isMc ) {
            jeu.SetJetPT( recoJet.recoPt );
            jeu.SetJetEta( recoJet.recoEta );
            jeu.SetJetPhi( recoJet.recoPhi );
            recoJet.recoPtJeuUp = recoJet.recoPt * (1. + jeu.GetUncertainty().first);
            recoJet.recoPtJeuDown = recoJet.recoPt * (1. - jeu.GetUncertainty().second);
        } // if ( !isMc )
    } // for (auto it = recoJets.begin(); it != recoJets.end();

}

//________________
void processRecoDijets(const bool isPp, const bool &isPbGoing, const bool &isMc, const float &weight,
                       std::vector<RecoJet> &recoJets, std::vector<GenJet> &genJets,
                       Histograms &hs,
                       JERSmearingHelper *fJERSmearingHelper,
                       const JERSmearingMode jerSmearingMode,
                       const bool useEtaDependentJER) {

    // std::cout << "processRecoDijets input parameters:" << std::endl;
    // std::cout << "  isPbGoing: " << (isPbGoing ? "true" : "false") << std::endl;
    // std::cout << "  isMc: " << (isMc ? "true" : "false") << std::endl;
    // std::cout << "  weight: " << weight << std::endl;
    // std::cout << "  recoJets.size(): " << recoJets.size() << std::endl;
    // std::cout << "  genJets.size(): " << genJets.size() << std::endl;
    // std::cout << "  hs: " << &hs << std::endl;
    // std::cout << "  fJERSmearingHelper: " << fJERSmearingHelper << std::endl;

    // With fewer than two selected reco jets, no momentum variation can form
    // a measured pair. Still account for a selected truth pair as a miss in
    // BOTH response families, once per event (not once per smearing replica).
    // prepareUnfolding guards the absent pair before inspecting its momenta.
    if (recoJets.size() < 2) {
        if (isMc) {
            prepareUnfolding(isPp, isPbGoing, isMc, weight, genJets, recoJets,
                             hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr,
                             hs.hGenDijetPtEtaCMMissArr,
                             hs.hRecoDijetPtEtaCMFakeArr,
                             hs.hUnfoldingPairClassificationArr,
                             hs.hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr,
                             hs.hGenDijetPtEtaCMMissControlledSampleArr,
                             hs.hRecoDijetPtEtaCMFakeControlledSampleArr,
                             hs.hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr,
                             hs.hGenDijetPtEtaCMMissStudiedSampleArr,
                             hs.hRecoDijetPtEtaCMFakeStudiedSampleArr, -1,
                             &hs.unfoldingPairs, &hs.unfoldingPairsControlled, &hs.unfoldingPairsStudied,
                             &hs.unfoldingExtended, &hs.unfoldingExtendedControlled, &hs.unfoldingExtendedStudied);
            prepareUnfolding(isPp, isPbGoing, isMc, weight, genJets, recoJets,
                             hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr,
                             hs.hGenDijetPtEtaCMMissJerDefExtraArr,
                             hs.hRecoDijetPtEtaCMFakeJerDefExtraArr,
                             hs.hUnfoldingPairClassificationJerDefExtraArr,
                             nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, -1,
                             &hs.unfoldingPairsJerDefExtra, nullptr, nullptr,
                             &hs.unfoldingExtendedJerDefExtra);
        }
        return;
    }

    std::sort(recoJets.begin(), recoJets.end(), [](const RecoJet &a, const RecoJet &b) { return a.recoPt > b.recoPt; });

    const auto& leadingJet = recoJets[0];
    const auto& subleadingJet = recoJets[1];
    if (leadingJet.recoPt >= 50.f && subleadingJet.recoPt >= 40.f) {
        float dphi = deltaPhi(leadingJet.recoPhi, subleadingJet.recoPhi);
        if (std::abs(dphi) >= TMath::TwoPi() / 3.) {
            float recoDijetPtAve = 0.5f * (leadingJet.recoPt + subleadingJet.recoPt);
            float recoLeadEtaLabUnflipped = leadingJet.recoEta;
            float recoLeadEtaLab = etaLab(leadingJet.recoEta, isPbGoing, isMc, isPp);
            float recoLeadEtaCM = etaCM(leadingJet.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
            float recoSubleadEtaLabUnflipped = subleadingJet.recoEta;
            float recoSubleadEtaLab = etaLab(subleadingJet.recoEta, isPbGoing, isMc, isPp);
            float recoSubleadEtaCM = etaCM(subleadingJet.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
            if (std::abs(recoLeadEtaCM) <= 1.9 && std::abs(recoSubleadEtaCM) <= 1.9) {
                float dijetEtaLabUnflipped = 0.5f * (recoLeadEtaLabUnflipped + recoSubleadEtaLabUnflipped);
                float dijetEtaLab = 0.5f * (recoLeadEtaLab + recoSubleadEtaLab);
                float dijetEtaCM = 0.5f * (recoLeadEtaCM + recoSubleadEtaCM);

                hs.hRecoDijetPtEtaLabUnflipped->Fill(recoDijetPtAve, dijetEtaLabUnflipped, weight);
                hs.hRecoDijetPtEtaLab->Fill(recoDijetPtAve, dijetEtaLab, weight);
                hs.hRecoDijetPtEtaCM->Fill(recoDijetPtAve, dijetEtaCM, weight);

                if (!isMc) {
                    hs.hRecoDijetPtEtaCMRunArr[nDiagnosticRuns]->Fill(
                        recoDijetPtAve, dijetEtaCM, weight);
                    for (int iRun{0}; iRun < nDiagnosticRuns; ++iRun) {
                        if (run == diagnosticRunIds[iRun]) {
                            hs.hRecoDijetPtEtaCMRunArr[iRun]->Fill(
                                recoDijetPtAve, dijetEtaCM, weight);
                            break;
                        }
                    }
                }

                if (leadingJet.refPt > 0. && subleadingJet.refPt > 0.) {
                    float refLeadEtaLabUnflipped = leadingJet.refEta;
                    float refLeadEtaLab = etaLab(leadingJet.refEta, isPbGoing, isMc, isPp);
                    float refLeadEtaCM = etaCM(leadingJet.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
                    float refSubleadEtaLabUnflipped = subleadingJet.refEta;
                    float refSubleadEtaLab = etaLab(subleadingJet.refEta, isPbGoing, isMc, isPp);
                    float refSubleadEtaCM = etaCM(subleadingJet.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
                    
                    float refDijetEtaLabUnflipped = 0.5f * (refLeadEtaLabUnflipped + refSubleadEtaLabUnflipped);
                    float refDijetEtaLab = 0.5f * (refLeadEtaLab + refSubleadEtaLab);
                    float refDijetEtaCM = 0.5f * (refLeadEtaCM + refSubleadEtaCM);

                    hs.hRefDijetPtEtaLabUnflipped->Fill(recoDijetPtAve, refDijetEtaLabUnflipped, weight);
                    hs.hRefDijetPtEtaLab->Fill(recoDijetPtAve, refDijetEtaLab, weight);
                    hs.hRefDijetPtEtaCM->Fill(recoDijetPtAve, refDijetEtaCM, weight);
                } // if (leadingJet.refPt > 0. && subleadingJet.refPt > 0.)
            } // if (std::abs(recoLeadEtaCM) <= 1.9 && std::abs(recoSubleadEtaCM) <= 1.9)
        } // if (std::abs(dphi) >= TMath::TwoPi() / 3.)
    } // if (leadingJet.recoPt >= 50.f && subleadingJet.recoPt >= 40.f)

    // fillRecoDijetSystematic:
    // - `getPt`: callable taking `const RecoJet&` and returning the pT value
    //   used for ordering and dijet pT calculations (float-like).
    // - Finds the two highest-pT jets without sorting the full collection.
    // - Applies pT (lead >= 50, sublead >= 40) and |Δφ| >= 2π/3 cuts, computes
    //   the dijet pT average using `getPt`, and fills the provided per-eta-cut
    //   full 2D histograms; analysis folds hemispheres from these.
    // - Matched-reference arrays retain the selected reconstructed pTave on x
    //   and use reference dijet eta on y. Pointing arrays instead apply each
    //   eta acceptance to the reference jets.
    // - Passing different accessors (e.g. `recoPt`, `recoPtJerUp`) lets this
    //   lambda fill nominal and systematic variations without duplicating code.
    // Re-selection allows pT rank changes and threshold migrations, including
    // events that fail nominal pT but pass a variation. Acceptance is imposed
    // AFTER choosing the pair; no replacement by lower-ranked in-range jets.
    std::vector<RecoJet> recoJetSortFallback;
    recoJetSortFallback.reserve(recoJets.size());
    auto fillRecoDijetSystematic = [&](auto getPt,
                                       std::unique_ptr<TH2D> *recoFullArr,
                                       std::unique_ptr<TH2D> *refFullArr = nullptr,
                                       std::unique_ptr<TH2D> *recoFullArrPointing = nullptr,
                                       std::unique_ptr<TH2D> *refFullArrPointing = nullptr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(recoJets, getPt, recoJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        // pT cuts for systematic variations (can be adjusted as needed)
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return; 
        float systDphi = deltaPhi(leadingJetSyst.recoPhi, subleadingJetSyst.recoPhi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaCM = etaCM(leadingJetSyst.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        float subleadEtaCM = etaCM(subleadingJetSyst.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

        bool  leadHasRef = (leadingJetSyst.refPt > 0.);
        bool  subleadHasRef = (subleadingJetSyst.refPt > 0.);
        bool  hasRefPair = (leadHasRef && subleadHasRef);
        float refLeadEtaCM = 0.;
        float refSubleadEtaCM = 0.;
        float refDijetEtaCM = 0.;
        const float refDijetPtAve = 0.5f * (leadingJetSyst.refPt + subleadingJetSyst.refPt);

        if (hasRefPair) {
            refLeadEtaCM = etaCM(leadingJetSyst.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
            refSubleadEtaCM = etaCM(subleadingJetSyst.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
            refDijetEtaCM = 0.5f * (refLeadEtaCM + refSubleadEtaCM);
        }

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            // Standard dijet processing
            if (std::abs(leadEtaCM) <= etaCuts[iCut] && std::abs(subleadEtaCM) <= etaCuts[iCut]) {
                // Fill reco dijet histograms
                recoFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);

                // Systematic calls reuse this helper; split only the nominal MC sample.
                if (isMc && recoFullArr == hs.hRecoDijetPtEtaCMArr) {
                    auto *sampleArr = isFillControlledSample
                        ? hs.hRecoDijetPtEtaCMControlledSampleArr
                        : hs.hRecoDijetPtEtaCMStudiedSampleArr;
                    sampleArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
                }
                if (isMc && recoFullArr == hs.hRecoDijetPtEtaCMArrJerDefUnfold) {
                    auto *sampleArr = isFillControlledSample
                        ? hs.hRecoDijetPtEtaCMJerDefUnfoldControlledSampleArr
                        : hs.hRecoDijetPtEtaCMJerDefUnfoldStudiedSampleArr;
                    sampleArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
                }

                // Fill ref dijet histograms (if a ref pair is available)
                if (refFullArr && hasRefPair) {
                    refFullArr[iCut]->Fill(systDijetPtAve, refDijetEtaCM, weight);
                    if (refFullArr == hs.hRefDijetPtEtaCMArr) {
                        auto *sampleArr = isFillControlledSample
                            ? hs.hRefDijetPtEtaCMControlledSampleArr
                            : hs.hRefDijetPtEtaCMStudiedSampleArr;
                        sampleArr[iCut]->Fill(systDijetPtAve, refDijetEtaCM, weight);
                    }
                }

            } // if (std::abs(leadEtaCM) < etaCuts[iCut] && std::abs(subleadEtaCM) < etaCuts[iCut])

            // Pointing resolution histograms for reference jets (if available)
            if (hasRefPair && recoFullArrPointing && refFullArrPointing) {
                // Reco and legacy Ref use reco pTave; the new Ref uses matched ref pTave.
                // Apply each acceptance in the
                // coordinate system whose dijet eta is being filled.
                if (std::abs(leadEtaCM) <= etaCuts[iCut] && std::abs(subleadEtaCM) <= etaCuts[iCut]) {
                    recoFullArrPointing[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
                }
                if (std::abs(refLeadEtaCM) <= etaCuts[iCut] && std::abs(refSubleadEtaCM) <= etaCuts[iCut]) {
                    refFullArrPointing[iCut]->Fill(systDijetPtAve, refDijetEtaCM, weight);
                    hs.hRefDijetPtEtaCMPointingResRefPtAveArr[iCut]->Fill(refDijetPtAve, refDijetEtaCM, weight);
                }
            } // if (hasRefPair && recoFullArrPointing && refFullArrPointing)
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };

    auto fillRecoDijetLabSystematic = [&](auto getPt,
                                       std::unique_ptr<TH2D> *recoFullArr,
                                       std::unique_ptr<TH2D> *refFullArr = nullptr,
                                       std::unique_ptr<TH2D> *recoFullArrPointing = nullptr,
                                       std::unique_ptr<TH2D> *refFullArrPointing = nullptr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(recoJets, getPt, recoJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return;
        float systDphi = deltaPhi(leadingJetSyst.recoPhi, subleadingJetSyst.recoPhi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaLab = etaLab(leadingJetSyst.recoEta, isPbGoing, isMc, isPp);
        float subleadEtaLab = etaLab(subleadingJetSyst.recoEta, isPbGoing, isMc, isPp);
        float systDijetEtaLab = 0.5f * (leadEtaLab + subleadEtaLab);

        bool  leadHasRef = (leadingJetSyst.refPt > 0.);
        bool  subleadHasRef = (subleadingJetSyst.refPt > 0.);
        bool  hasRefPair = (leadHasRef && subleadHasRef);
        float refLeadEtaLab = 0.;
        float refSubleadEtaLab = 0.;
        float refDijetEtaLab = 0.;
        const float refDijetPtAve = 0.5f * (leadingJetSyst.refPt + subleadingJetSyst.refPt);

        if (hasRefPair) {
            refLeadEtaLab = etaLab(leadingJetSyst.refEta, isPbGoing, isMc, isPp);
            refSubleadEtaLab = etaLab(subleadingJetSyst.refEta, isPbGoing, isMc, isPp);
            refDijetEtaLab = 0.5f * (refLeadEtaLab + refSubleadEtaLab);
        }

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            if (std::abs(leadEtaLab) <= etaCuts[iCut] && std::abs(subleadEtaLab) <= etaCuts[iCut]) {
                // Fill reco dijet histograms
                recoFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaLab, weight);

                // Fill ref dijet histograms (if a ref pair is available)
                if (refFullArr && hasRefPair) {
                    refFullArr[iCut]->Fill(systDijetPtAve, refDijetEtaLab, weight);
                }
            } // if (std::abs(leadEtaLab) < etaCuts[iCut] && std::abs(subleadEtaLab) < etaCuts[iCut])

            if (hasRefPair && recoFullArrPointing && refFullArrPointing) {
                if (std::abs(leadEtaLab) <= etaCuts[iCut] && std::abs(subleadEtaLab) <= etaCuts[iCut]) {
                    recoFullArrPointing[iCut]->Fill(systDijetPtAve, systDijetEtaLab, weight);
                }
                if (std::abs(refLeadEtaLab) <= etaCuts[iCut] && std::abs(refSubleadEtaLab) <= etaCuts[iCut]) {
                    refFullArrPointing[iCut]->Fill(systDijetPtAve, refDijetEtaLab, weight);
                    hs.hRefDijetPtEtaLabPointingResRefPtAveArr[iCut]->Fill(refDijetPtAve, refDijetEtaLab, weight);
                }
            } // if (hasRefPair && recoFullArrPointing && refFullArrPointing)

        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    }; // auto fillRecoDijetLabSystematic

    // Fill nominal dijet histograms
    if (isMc) {
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPt; },
                                hs.hRecoDijetPtEtaCMArr,
                                hs.hRefDijetPtEtaCMArr,
                                hs.hRecoDijetPtEtaCMPointingResArr,
                                hs.hRefDijetPtEtaCMPointingResArr);

        fillRecoDijetLabSystematic([](const RecoJet &jet) { return jet.recoPt; },
                                    hs.hRecoDijetPtEtaLabArr,
                                    hs.hRefDijetPtEtaLabArr,
                                    hs.hRecoDijetPtEtaLabPointingResArr,
                                    hs.hRefDijetPtEtaLabPointingResArr);
    }
    else {
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPt; },
                                hs.hRecoDijetPtEtaCMArr);

        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJeuUp; },
                                hs.hRecoDijetPtEtaCMArrJeuUp);

        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJeuDown; },
                                hs.hRecoDijetPtEtaCMArrJeuDown);
        fillRecoDijetLabSystematic([](const RecoJet &jet) { return jet.recoPt; },
                                    hs.hRecoDijetPtEtaLabArr);
    }

    // Fill JER systematic dijet histograms if requested
    if (isMc) {

        // Number of iterations for smearing
        const int nSmearRounds = jerSmearingMode == JERSmearingMode::TruthToPseudoReco ? 10 : 1;
        // Freeze identity before the first response sort, including later replicas.
        for (size_t i = 0; i < recoJets.size(); ++i)
            recoJets[i].smearingIndex = static_cast<int>(i);

        for (int iIter = 0; iIter < nSmearRounds; ++iIter) {
            for (auto &jet : recoJets) {
                const double normal = fJERSmearingHelper->normalDeviate(
                    jet.recoPt, jet.recoEta, eventRho, run, lumi, evt,
                    iIter * static_cast<int>(recoJets.size()) + jet.smearingIndex);
                cacheRecoJetSmearing(jet, *fJERSmearingHelper, jerSmearingMode,
                                    useEtaDependentJER, normal);
            }
            for (int iScale = 0; iScale < 9; ++iScale) {
                for (auto &recoJet : recoJets) {
                    const bool hasReference = recoJet.refPt > 0.;
                    // Retain at most the first ten diagnostic replicas; the
                    // official Hybrid/Stochastic modes execute only iIter=0.
                    if (iIter > 9) continue; 

                    if (!hasReference) continue;
                    if (iScale == 0) { // JER def no SF (x1.0)
                        hs.hRecoInclusiveJetJESDefNoSFPtEta->Fill(recoJet.jerMomenta[iScale] / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 1) { // JER def double SF (x2.0)
                        hs.hRecoInclusiveJetJESDefDoublePtEta->Fill(recoJet.jerMomenta[iScale] / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 2) { // JER def ninety SF (x0.9)
                        hs.hRecoInclusiveJetJESDefNinetyPtEta->Fill(recoJet.jerMomenta[iScale] / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 3) { // JER down
                        hs.hRecoInclusiveJetJESDownPtEta->Fill(recoJet.jerMomenta[iScale] / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }   
                    if (iScale == 4) { // JER def
                        hs.hRecoInclusiveJetJESDefPtEta->Fill(recoJet.jerMomenta[iScale] / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 5) { // JER up
                        hs.hRecoInclusiveJetJESUpPtEta->Fill(recoJet.jerMomenta[iScale] / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                } // for (const auto &recoJet : recoJets)

                if (iScale == 0) { // JER def no SF (x1.0)
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDefNoSF);
                }
                else if (iScale == 1) { // JER def double SF (x2.0)
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDefDoubleSF);
                }
                else if (iScale == 2) { // JER def ninety SF (x0.9)
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDefNinetySF);
                }
                else if (iScale == 3) { // JER down
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDown,
                                            hs.hRefDijetPtEtaCMArrJerDown);
                }
                else if (iScale == 4) { // JER def
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDef,
                                            hs.hRefDijetPtEtaCMArrJerDef);

                    // For the first iteration, fill the nominal histograms with the default 
                    // JER smearing for unfolding purposes
                    if (iIter == 0) {
                        fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                                hs.hRecoDijetPtEtaCMArrJerDefUnfold);
                        prepareUnfolding(isPp, isPbGoing, isMc, weight, genJets, recoJets,
                                         hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr,
                                         hs.hGenDijetPtEtaCMMissArr,
                                         hs.hRecoDijetPtEtaCMFakeArr,
                                         hs.hUnfoldingPairClassificationArr,
                                         hs.hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr,
                                         hs.hGenDijetPtEtaCMMissControlledSampleArr,
                                         hs.hRecoDijetPtEtaCMFakeControlledSampleArr,
                                         hs.hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr,
                                         hs.hGenDijetPtEtaCMMissStudiedSampleArr,
                                         hs.hRecoDijetPtEtaCMFakeStudiedSampleArr, iScale,
                                         &hs.unfoldingPairs, &hs.unfoldingPairsControlled, &hs.unfoldingPairsStudied,
                                         &hs.unfoldingExtended, &hs.unfoldingExtendedControlled, &hs.unfoldingExtendedStudied);
                    }
                }
                else if (iScale == 5) { // JER up
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerUp,
                                            hs.hRefDijetPtEtaCMArrJerUp);
                }
                else if (iScale == 6) { // JER down with eta-dependent scaling
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDownExtra);
                }
                else if (iScale == 7) { // JER default with eta-dependent scaling
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerDefExtra);
                    if (iIter == 0) {
                        fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                                hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold);
                        prepareUnfolding(isPp, isPbGoing, isMc, weight, genJets, recoJets,
                                         hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr,
                                         hs.hGenDijetPtEtaCMMissJerDefExtraArr,
                                         hs.hRecoDijetPtEtaCMFakeJerDefExtraArr,
                                         hs.hUnfoldingPairClassificationJerDefExtraArr,
                                         nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, iScale,
                                         &hs.unfoldingPairsJerDefExtra, nullptr, nullptr,
                                         &hs.unfoldingExtendedJerDefExtra);
                    }
                }
                else if (iScale == 8) { // JER up with eta-dependent scaling
                    fillRecoDijetSystematic([iScale](const RecoJet &jet) { return jet.jerMomenta[iScale]; },
                                            hs.hRecoDijetPtEtaCMArrJerUpExtra);
                }
            } // for (int iScale = 0; iScale < 9; ++iScale)
        } // for (int iIter = 0; iIter < nSmearRounds; ++iIter)

    } // if (isMc)
}

//________________
void processRefDijets(const bool isPp, const bool &isPbGoing, const bool &isMc, const float &weight,
                       std::vector<RecoJet> &recoJets, Histograms &hs) {

    // std::cout << "processRefDijets input parameters:" << std::endl;
    // std::cout << "  isPbGoing: " << (isPbGoing ? "true" : "false") << std::endl;
    // std::cout << "  isMc: " << (isMc ? "true" : "false") << std::endl;
    // std::cout << "  weight: " << weight << std::endl;
    // std::cout << "  recoJets.size(): " << recoJets.size() << std::endl;

    // Must be at least 2 jets to form a dijet system
    if (recoJets.size() < 2) return;

    // fillRefDijetSystematic:
    // - `getPt`: callable taking `const RecoJet&` and returning the pT value
    //   used for ordering and dijet pT calculations (float-like).
    // - Finds the two highest-pT jets without sorting the full collection.
    // - Applies pT (lead >= 50, sublead >= 40) and |Δφ| >= 2π/3 cuts, computes
    //   the dijet pT average using `getPt`, and fills the provided per-eta-cut
    //   corresponding reference-level 2D histograms.
    std::vector<RecoJet> refJetSortFallback;
    refJetSortFallback.reserve(recoJets.size());
    auto fillRefDijetSystematic = [&](auto getPt,
                                       std::unique_ptr<TH2D> *refFullArr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(recoJets, getPt, refJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        // pT cuts for systematic variations (can be adjusted as needed)
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return; 
        float systDphi = deltaPhi(leadingJetSyst.refPhi, subleadingJetSyst.refPhi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaCM = etaCM(leadingJetSyst.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        float subleadEtaCM = etaCM(subleadingJetSyst.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc, isPp);
        float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;
            refFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };

    // RefSel is the matched-reference collection BEFORE configured reco ID
    // pruning (see processEvents order). It is not an independent inclusive
    // generator sample and should not be interpreted as a truth efficiency.
    // Fill nominal dijet histograms
    fillRefDijetSystematic([](const RecoJet &jet) { return jet.refPt; },
                           hs.hRefSelDijetPtEtaCMArr);
}

void loadGenJets(std::vector<GenJet>& jets) {
    jets.assign(currentEvent->genJets.begin(), currentEvent->genJets.end());
    std::sort(jets.begin(), jets.end(), [](const GenJet& a,const GenJet& b){return a.pt>b.pt;});
}
void loadRecoJets(std::vector<RecoJet>& jets, JetCorrector& jec, const bool& isMc) {
    jets.clear();
    for (const auto& raw : currentEvent->recoJets) {
        jec.SetJetPT(raw.recoPtRaw); jec.SetJetEta(raw.recoEta); jec.SetJetPhi(raw.recoPhi);
        const float pt=jec.GetCorrectedPT();
        RecoJet jet=raw;
        jet.recoPt=pt;
        if (!isMc) { jet.refPt=-1.f; jet.refEta=0.f; jet.refPhi=0.f; }
        jets.push_back(jet);
    }
    std::sort(jets.begin(), jets.end(), [](const RecoJet& a,const RecoJet& b){return a.recoPt>b.recoPt;});
}
void initializeSmearing() {
    std::cout << "Start event processing..." << std::endl;

    // Vz weight function for MC, derived from data/MC ratio of vz distribution in minimum bias events
    fVzWeight = std::make_unique<TF1>("fVzWeight", "pol8", -15.1, 15.1);
    fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);

    // JER smearing helper for calculating the smearing factor based on eta and pT
    fJerSmearingHelper = std::make_unique<JERSmearingHelper>();

    auto fJERSmearVsEtaFunc_pt_30_50 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_30_50", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_50_80 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_50_80", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_80_100 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_80_100", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_100_120 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_100_120", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_120_150 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_120_150", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_150_180 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_150_180", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_180_220 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_180_220", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_220_280 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_220_280", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_280_350 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_280_350", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_350_450 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_350_450", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_450_540 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_450_540", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);
    auto fJERSmearVsEtaFunc_pt_540_8160 = std::make_unique<TF1>("fJERSmearVsEtaFunc_pt_540_8160", "gaus(0) + gaus(3) + pol2(6)", -2.8, 2.8);

    if (isMc) {
        if (isPythia) {
            if (isPbGoing) {
                fJerSmearingHelper->setJERSmearingFuncVsPtParams(0.04360, 0.81815); // -0.8 < eta < 0.8, 30 < pt (GeV) < 800
                fJERSmearVsEtaFunc_pt_30_50->SetParameters(0.297573, -1.45445, 0.280701, 0.322475, 1.47685, 0.291071, 0.989036, -0.00955405, 0.025316);
                fJERSmearVsEtaFunc_pt_50_80->SetParameters(0.293634, -1.43467, 0.277413, 0.345306, 1.44334, 0.310660, 0.994119, -0.0106212, 0.00165775);
                fJERSmearVsEtaFunc_pt_80_100->SetParameters(0.293292, -1.41581, 0.269992, 0.322632, 1.41280, 0.311503, 0.998526, -0.00943762, -0.0117420);
                fJERSmearVsEtaFunc_pt_100_120->SetParameters(0.292156, -1.40265, 0.219616, 0.321592, 1.38074, 0.272997, 1.00353, -0.0102294, -0.0148459);
                fJERSmearVsEtaFunc_pt_120_150->SetParameters(0.277961, -1.38076, 0.218392, 0.303384, 1.40581, 0.273054, 1.00547, -0.0142556, -0.0246123);
                fJERSmearVsEtaFunc_pt_150_180->SetParameters(0.305976, -1.37269, 0.169487, 0.290438, 1.36259, 0.283952, 1.00294, -0.0202212, -0.0234219);
                fJERSmearVsEtaFunc_pt_180_220->SetParameters(0.249333, -1.34427, 0.202333, 0.291205, 1.34517, 0.237504, 1.00372, -0.0115382, -0.0216865);
                fJERSmearVsEtaFunc_pt_220_280->SetParameters(0.261671, -1.34257, 0.191754, 0.280618, 1.33863, 0.274849, 1.00387, -0.0154493, -0.0247799);
                fJERSmearVsEtaFunc_pt_280_350->SetParameters(0.276569, -1.33522, 0.168794, 0.313315, 1.34124, 0.271635, 1.00615, -0.0219329, -0.0321017);
                fJERSmearVsEtaFunc_pt_350_450->SetParameters(0.288649, -1.33524, 0.165996, 0.360752, 1.33645, 0.289183, 1.00648, -0.0420297, -0.0420253);
                fJERSmearVsEtaFunc_pt_450_540->SetParameters(0.330524, -1.33153, 0.159125, 0.403670, 1.33649, 0.288587, 1.00779, -0.0517481, -0.0513578);
                fJERSmearVsEtaFunc_pt_540_8160->SetParameters(0.380368, -1.33655, 0.142485, 0.537708, 1.42214, 0.413657, 0.997332, -0.0970631, -0.0743566);
            }
            else {
                fJerSmearingHelper->setJERSmearingFuncVsPtParams(0.04305, 0.82470); // -0.8 < eta < 0.8, 30 < pt (GeV) < 800
                fJERSmearVsEtaFunc_pt_30_50->SetParameters(0.310346, -1.49539, 0.343705, 0.332612, 1.41243, 0.300797, 0.986439, 0.00256937, 0.0210452);
                fJERSmearVsEtaFunc_pt_50_80->SetParameters(0.313000, -1.46533, 0.280744, 0.336817, 1.43451, 0.299505, 0.997038, -0.00414605, 0.00216964);
                fJERSmearVsEtaFunc_pt_80_100->SetParameters(0.282013, -1.42527, 0.262737, 0.329079, 1.40888, 0.276463, 1.00096, 0.00160294, -0.0166477);
                fJERSmearVsEtaFunc_pt_100_120->SetParameters(0.290421, -1.39473, 0.262791, 0.306412, 1.39970, 0.255765, 0.997020, -0.000126060, -0.0166709);
                fJERSmearVsEtaFunc_pt_120_150->SetParameters(0.269288, -1.37681, 0.219038, 0.281362, 1.37189, 0.273782, 1.00072, -0.00553830, -0.0192238);
                fJERSmearVsEtaFunc_pt_150_180->SetParameters(0.276496, -1.33329, 0.280118, 0.289559, 1.36278, 0.249674, 0.994538, 0.00157502, -0.0228272);
                fJERSmearVsEtaFunc_pt_180_220->SetParameters(0.277700, -1.35258, 0.199943, 0.266430, 1.36892, 0.244938, 1.00112, 0.00169535, -0.0244297);
                fJERSmearVsEtaFunc_pt_220_280->SetParameters(0.290909, -1.33599, 0.207813, 0.303474, 1.32557, 0.240638, 1.00184, 0.00191870, -0.0218494);
                fJERSmearVsEtaFunc_pt_280_350->SetParameters(0.284738, -1.31639, 0.209728, 0.282979, 1.32785, 0.216577, 1.00362, 0.00343919, -0.0278720);
                fJERSmearVsEtaFunc_pt_350_450->SetParameters(0.295510, -1.31097, 0.264653, 0.295145, 1.33895, 0.237917, 1.00112, 0.00920355, -0.0342905);
                fJERSmearVsEtaFunc_pt_450_540->SetParameters(0.302981, -1.30795, 0.223667, 0.332526, 1.32390, 0.228734, 1.00695, 0.00227566, -0.0339503);
                fJERSmearVsEtaFunc_pt_540_8160->SetParameters(0.301793, -1.30626, 0.242770, 0.364985, 1.31322, 0.239480, 1.00411, -0.0106817, -0.0309593);
            }
        }
        else { // Embdedding
            if (isPbGoing) {
                fJerSmearingHelper->setJERSmearingFuncVsPtParams(0.04360, 0.81815); // -0.8 < eta < 0.8, 30 < pt (GeV) < 800
                fJERSmearVsEtaFunc_pt_30_50->SetParameters(0.251023, -1.44091, 0.247814, 0.259187, 1.43814, 0.284808, 0.980647, -0.00616568, 0.0127762);
                fJERSmearVsEtaFunc_pt_50_80->SetParameters(2.80975e-01, -1.43118e+00, 2.62786e-01, 3.11693e-01, 1.43402e+00, 2.98206e-01, 9.97888e-01, -9.32355e-03, -4.82988e-03);
                fJERSmearVsEtaFunc_pt_80_100->SetParameters(2.94796e-01, -1.41661e+00, 2.34356e-01, 3.06224e-01, 1.43206e+00, 2.88265e-01, 1.00413e+00, -8.87238e-03, -1.53309e-02);
                fJERSmearVsEtaFunc_pt_100_120->SetParameters(0.268538, -1.38565, 0.228907, 0.322785, 1.38592, 0.26244, 1.00201, -0.0128027, -0.021899);
                fJERSmearVsEtaFunc_pt_120_150->SetParameters(0.270962, -1.35472, 0.222158, 0.295618, 1.40682, 0.296922, 1.00146, -0.0130908, -0.023406);
                fJERSmearVsEtaFunc_pt_150_180->SetParameters(0.291945, -1.37874, 0.170634, 0.274105, 1.39584, 0.282491, 1.00488, -0.0178596, -0.0266459);
                fJERSmearVsEtaFunc_pt_180_220->SetParameters(0.248293, -1.33017, 0.206415, 0.279631, 1.35826, 0.274813, 1.00113, -0.0154812, -0.0243555);
                fJERSmearVsEtaFunc_pt_220_280->SetParameters(0.268718, -1.33038, 0.19311, 0.298086, 1.33168, 0.282377, 0.998327, -0.0210509, -0.0276671);
                fJERSmearVsEtaFunc_pt_280_350->SetParameters(0.274991, -1.33119, 0.176131, 0.327383, 1.32507, 0.231403, 1.00344, -0.0174497, -0.0290836);
                fJERSmearVsEtaFunc_pt_350_450->SetParameters(0.293324, -1.34356, 0.166471, 0.328376, 1.34978, 0.288753, 1.00615, -0.0339521, -0.0397907);
                fJERSmearVsEtaFunc_pt_450_540->SetParameters(0.324936, -1.32021, 0.164281, 0.369252, 1.34049, 0.286628, 1.00614, -0.0389711, -0.0465924);
                fJERSmearVsEtaFunc_pt_540_8160->SetParameters(0.37107, -1.34014, 0.148853, 0.495634, 1.37641,  0.350713, 1.00368, -0.0782385, -0.0654875);
            }
            else {
                fJerSmearingHelper->setJERSmearingFuncVsPtParams(0.04181, 0.8614); // -0.8 < eta < 0.8, 30 < pt (GeV) < 800
                fJERSmearVsEtaFunc_pt_30_50->SetParameters(0.236966, -1.51017, 0.335154, 0.305743, 1.44193, 0.292118, 0.986551, 6.73397e-05, 0.0118094);
                fJERSmearVsEtaFunc_pt_50_80->SetParameters(0.296463, -1.46524, 0.303889, 0.332084, 1.45171, 0.295073, 0.99803, -0.00414073, -0.00244323);
                fJERSmearVsEtaFunc_pt_80_100->SetParameters(0.281362, -1.436, 0.269096, 0.277846, 1.40117, 0.297503, 0.999944, 0.000818982, -0.0206659);
                fJERSmearVsEtaFunc_pt_100_120->SetParameters(0.285419, -1.40383, 0.20919, 0.293156, 1.38149, 0.261957, 1.00248, -0.00378822, -0.0210735);
                fJERSmearVsEtaFunc_pt_120_150->SetParameters(0.286209, -1.35991, 0.261485, 0.277341, 1.34677, 0.276339, 0.999703, 0.00635947, -0.0228129);
                fJERSmearVsEtaFunc_pt_150_180->SetParameters(0.264685, -1.38645, 0.207738, 0.265336, 1.37243, 0.246389, 1.00305, -0.00340848, -0.0242265);
                fJERSmearVsEtaFunc_pt_180_220->SetParameters(0.28098, -1.32577, 0.240958, 0.292022, 1.34173, 0.246079, 0.995414, 0.00611453, -0.0235741);
                fJERSmearVsEtaFunc_pt_220_280->SetParameters(0.269833, -1.37708, 0.364279, 0.257267, 1.33253, 0.216335, 0.993059, 0.0303286, -0.0349511);
                fJERSmearVsEtaFunc_pt_280_350->SetParameters(0.297607, -1.31883, 0.207127, 0.294736, 1.34136, 0.212788, 1.00491, 0.00283823, -0.0309229);
                fJERSmearVsEtaFunc_pt_350_450->SetParameters(0.291975, -1.32442, 0.247271, 0.317328, 1.32188, 0.221408, 1.00603, 0.0079703, -0.0378725);
                fJERSmearVsEtaFunc_pt_450_540->SetParameters(3.06494e-01, -1.27420, 2.57302e-01, 3.37387e-01, 1.31932e+00, 2.29305e-01, 1.00217e+00, 6.98099e-03, -3.63918e-02);
                fJERSmearVsEtaFunc_pt_540_8160->SetParameters(2.65766e-01, -1.30667e+00, 2.73191e-01, 3.68171e-01, 1.31464e+00, 2.50415e-01, 1.00047e+00, -1.67119e-02, -2.60016e-02);
            }
        }

        fJerSmearingHelper->addJERSmearingFuncVsEta(30., 50., std::move(fJERSmearVsEtaFunc_pt_30_50));
        fJerSmearingHelper->addJERSmearingFuncVsEta(50., 80., std::move(fJERSmearVsEtaFunc_pt_50_80));
        fJerSmearingHelper->addJERSmearingFuncVsEta(80., 100., std::move(fJERSmearVsEtaFunc_pt_80_100));
        fJerSmearingHelper->addJERSmearingFuncVsEta(100., 120., std::move(fJERSmearVsEtaFunc_pt_100_120));
        fJerSmearingHelper->addJERSmearingFuncVsEta(120., 150., std::move(fJERSmearVsEtaFunc_pt_120_150));
        fJerSmearingHelper->addJERSmearingFuncVsEta(150., 180., std::move(fJERSmearVsEtaFunc_pt_150_180));
        fJerSmearingHelper->addJERSmearingFuncVsEta(180., 220., std::move(fJERSmearVsEtaFunc_pt_180_220));
        fJerSmearingHelper->addJERSmearingFuncVsEta(220., 280., std::move(fJERSmearVsEtaFunc_pt_220_280));
        fJerSmearingHelper->addJERSmearingFuncVsEta(280., 350., std::move(fJERSmearVsEtaFunc_pt_280_350));
        fJerSmearingHelper->addJERSmearingFuncVsEta(350., 450., std::move(fJERSmearVsEtaFunc_pt_350_450));
        fJerSmearingHelper->addJERSmearingFuncVsEta(450., 540., std::move(fJERSmearVsEtaFunc_pt_450_540));
        fJerSmearingHelper->addJERSmearingFuncVsEta(540., 8160., std::move(fJERSmearVsEtaFunc_pt_540_8160));
    }

    // auto fJERSmearVsEtaFunc = std::make_unique<TF1>("fJERSmearVsEtaFunc", "gaus(0) + gaus(3) + pol2(6)", -2.6, 2.6);
    // fJERSmearVsEtaFunc->SetParameters(0.280798, -1.42826, 0.240604, 0.350327, 1.41888, 0.260303, 1.00702, -0.0255149, -0.032875);

    startTime=std::chrono::steady_clock::now();
}
void processEvent(const Event& event) {
    currentEvent=&event;
    auto& jec=*jecOwner;
    auto& jeu=*jeuOwner;
        if (isMc) {
            isFillControlledSample = fillControlledSample(sampleSplitGenerator);
        }

        if (!event.hasRecoJets || (isMc && !event.hasGenJets))
            throw std::runtime_error("DiJetAnalysis requires loaded reco jets and MC gen jets");
        if (static_cast<int>(event.sample) != (isMc ? (isPythia ? 2 : 1) : 0))
            throw std::runtime_error("Reader sample type does not match analysis configuration");
        static_cast<EventInfo&>(*this)=event.info;
        nEventsProcessed++;
        inputRecoCount += event.recoJets.size(); inputGenCount += event.genJets.size();
        Long64_t entriesPassed = nEventsProcessed;
        if ((entriesPassed >= nextPrintAt) || (entriesPassed == nEntries)) {
            auto now = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = now - startTime;
            double elapsedSeconds = elapsed.count();
            double etaSeconds = 0.0;
            if (entriesPassed > 0 && entriesPassed < nEntries) {
                etaSeconds = elapsedSeconds * static_cast<double>(nEntries - entriesPassed) / static_cast<double>(entriesPassed);
            }

            Long64_t elapsedTotalSeconds = static_cast<Long64_t>(std::llround(elapsedSeconds));
            Long64_t elapsedMinutes = elapsedTotalSeconds / 60;
            Long64_t elapsedRemainingSeconds = elapsedTotalSeconds % 60;

            Long64_t etaTotalSeconds = static_cast<Long64_t>(std::llround(etaSeconds));
            Long64_t etaMinutes = etaTotalSeconds / 60;
            Long64_t etaRemainingSeconds = etaTotalSeconds % 60;

            std::cout << Form("Processed %lld / %lld entries | elapsed %lldm:%llds | ETA %lldm:%llds",
                              entriesPassed, nEntries,
                              elapsedMinutes, elapsedRemainingSeconds,
                              etaMinutes, etaRemainingSeconds) << std::endl;

            nextPrintAt += printEvery;
        }

        // P05: unlike etaLab/etaCM this legacy vz transform has no pp guard.
        // pp uses an even |vz| cut and stored weights, so selection is unchanged,
        // but its signed vertex diagnostics still depend on isPbGoing.
        // Flip vz to properly account for the p-going and Pb-going directions
        // in data and MC. Definition is proton beam going in +z direction
        if (isMc) {
            if (isPbGoing) vz = -vz; 
        }
        else {
            if (!isPbGoing) vz = -vz; 
        }

        hs.hVzRaw->Fill(vz, 1.);

        // std::cout << "\n========================================\n";
        // std::cout << Form("Event: %d ptHat = %.1f GeV, vz = %.2f cm", (nEventsProcessed-1), ptHat, vz) << std::endl;
        // std::cout << Form("%.1f <= ptHat <= %.1f ? %s", ptHatRange[0], ptHatRange[1], 
        //                  (ptHat >= ptHatRange[0] && ptHat <= ptHatRange[1]) ? "true" : "false") << std::endl;
        // std::cout << Form("|vz| <= 15 cm ? %s", (std::abs(vz) <= 15.) ? "true" : "false") << std::endl;

        // Check the event is satisfies basic event selection
        if (!isGoodEvent(isPbGoing, isMc, triggerId, vertexFilterSelection, isPp)) return;

        if ( isMc ) {
            weight = eventWeight(ptHat, vz, *fVzWeight, nEntries, isPp, inputEventWeight);
        }

        // Load reconstructed jets
        timings.measure(StageTimings::RecoLoad, [&] { loadRecoJets(recoJets, jec, isMc); });

        if (isMc) {
            // Load generated jets
            timings.measure(StageTimings::GenLoad, [&] { loadGenJets( genJets ); });
            // Fill histograms for overweight events and skip them
            fillOverweightHistograms(genJets, recoJets, weight, hs, isPythia, isGenOverweight, isRecoOverweight);
            if (isGenOverweight) return;
            if (isRecoOverweight) return;
        } // if ( isMc )

        nGoodEvents++;

        hs.hPtHatUnweighted->Fill(ptHat, 1.);
        hs.hPtHat->Fill(ptHat, weight);
        hs.hVzUnweighted->Fill(vz, 1.);
        hs.hVz->Fill(vz, weight);

        //
        // Gen level processing
        //
        if (isMc) {
            timings.measure(StageTimings::GenJets, [&] { processGenJets(isPp, isPbGoing, isMc, weight, genJets, hs, ptHat); });
            timings.measure(StageTimings::GenDijets, [&] { processGenDijets(isPp, isPbGoing, isMc, weight, genJets, hs, fJerSmearingHelper.get()); });
            timings.measure(StageTimings::RefDijets, [&] { processRefDijets(isPp, isPbGoing, isMc, weight, recoJets, hs); });
        }

        //
        // Reco level processing
        //
        timings.measure(StageTimings::RecoJets, [&] { processRecoJets(isPp, isPbGoing, isMc, weight, recoJets, hs, jeu, jetSelectionMethod, ptHat); });
        timings.measure(StageTimings::RecoDijets, [&] { processRecoDijets(isPp, isPbGoing, isMc, weight, recoJets, genJets, hs,
                          fJerSmearingHelper.get(), jerSmearingMode, useEtaDependentJER); });
}
void reportCounts() {
    std::cout << "Input objects: reco=" << inputRecoCount << " gen=" << inputGenCount << std::endl;
    std::cout << Form("Processed %d events, of which %d passed selection.", nEventsProcessed, nGoodEvents) << std::endl;
}
void initialize(const AnalysisConfig& config, const RunContext& context) {
    timings.enabled=config.profileStages;
    const char* input=config.reader.input.c_str(); const char* output=config.output.c_str();
    const char* collisionSystem=config.reader.pp ? "pp" : "pPb";
    const char* path2auxFiles=config.auxDirectory.c_str();
    const int mcType=config.reader.mcType, isPbGoingDir=config.isPbGoingDir,
        ptHatSample=config.ptHatSample, triggerId=config.triggerId,
        recoJetSelMethod=config.recoJetSelMethod, vertexFilterSelection=config.vertexFilterSelection,
        jerSmearingMode=config.jerSmearingMode, useEtaDependentJER=config.useEtaDependentJER;
    if (collisionSystem == nullptr ||
        (std::string(collisionSystem) != "pPb" && std::string(collisionSystem) != "pp")) {
        std::cerr << RED << "Invalid collisionSystem. Must be pPb or pp." << RESET << std::endl;
        usage();
    }
    const bool isPp = std::string(collisionSystem) == "pp";

    // Validate input parameters
    if (mcType < 0 || mcType > 2) {
        std::cerr << RED << "Invalid mcType parameter. Must be 0 (data), 1 (embedding), or 2 (pythia). Terminating." << RESET << std::endl;
        usage();
    }
    if (isPbGoingDir < 0 || isPbGoingDir > 1) {
        std::cerr << RED << "Invalid isPbGoingDir parameter. Must be 0 (p-going) or 1 (Pb-going). Terminating." << RESET << std::endl;
        usage();
    }
    if (triggerId < 0 || triggerId > 3) {
        std::cerr << RED << "Invalid triggerId parameter. Must be 0 (no trigger/MB), 1 (jet60), 2 (jet80), or 3 (jet100). Terminating." << RESET << std::endl;
        usage();
    }
    if (recoJetSelMethod < 0 || recoJetSelMethod > 3) {
        std::cerr << RED << "Invalid recoJetSelMethod parameter. Must be 0 (no selection), 1 (trkMaxPt/RawPt), 2 (jetIdTightLeptVeto), or 3 (jetIdLoose). Terminating." << RESET << std::endl;
        usage();
    }
    if (vertexFilterSelection < 0 || vertexFilterSelection > 2) {
        std::cerr << RED << "Invalid vertexFilterSelection parameter. Must be 0 (dz1p0), 1 (Gplus), or 2 (Vtx1). Terminating." << RESET << std::endl;
        usage();
    }
    if (jerSmearingMode < 0 || jerSmearingMode > 2) {
        std::cerr << RED << "Invalid jerSmearingMode. Must be 0 (hybrid), 1 (stochastic), or 2 (truth-to-pseudo-reco)." << RESET << std::endl;
        usage();
    }
    if (useEtaDependentJER < 0 || useEtaDependentJER > 1) {
        std::cerr << RED << "Invalid useEtaDependentJER. Must be 0 or 1." << RESET << std::endl;
        usage();
    }

    bool isMc = (mcType != 0);
    bool isPythia = ((mcType != 0) && (mcType == 2));
    // Select by sample type, not the Boolean MC flag. Standalone PYTHIA
    // (including pp MC) uses the direction-specific pPb Unembedded payload.
    bool isEmbedding = (mcType == 1);
    bool isPbGoing = (isPbGoingDir == 1);

    // Define the standard eta bins for L2L3 corrections
    TString generator = isPythia ? "pythia" : "embedding";
    TString tag = isPythia ? "unembedded" : "embedded";
    TString directionStr = isPbGoing ? "Pbgoing" : "pgoing";

    std::vector< std::string > fJECFileNames;
    std::string fJEUFileName;

    // The correction family is selected by MC/direction, not collisionSystem:
    // pp currently reuses these pPb calibrations (P05). A pp physics comparison
    // requires establishing that this is the intended calibration convention.
    if (isMc) {
        if (isPbGoingDir) {
            if (isEmbedding) {
                // PYTHIA+EPOS (embedding)
                fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
            }
            else {
                // PYTHIA
                fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_Pbgoing_Unembedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
            }
        }
        else {
            if (isEmbedding) {
                // PYTHIA+EPOS (embedding)
                fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
            }
            else {
                // PYTHIA
                fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
            }
        }
    }
    else {
        if (isPbGoingDir) { // Remember to flip to p-going for data
            fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_pgoing_Embedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
        }
        else {
            fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
        }
        // Old correction
        // JECFileDataName = "Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF.txt";
        // New correction (with updated residuals)
        fJECFileNames.push_back( Form("%s/JEC/Summer16_07Aug2017GH_V11_DATA_L2L3Residual_AK4PF.txt", path2auxFiles) );
        
        // JEU (Jet Energy Uncertainty)
        fJEUFileName = Form("%s/JEC/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF.txt", path2auxFiles);
    } // else


    if (config.reader.jetTree != "ak4PFJetAnalyzer/t")
        throw std::runtime_error("DiJetAnalysis currently supports AK4PF only; other collections require approved collection-specific JEC/JER calibrations. Use --reader-only.");
    if (!config.jecFiles.empty()) fJECFileNames=config.jecFiles;
    if (!config.jeuFile.empty()) fJEUFileName=config.jeuFile;
    // Set pT hat range for MC samples
    if (isMc && !isPp) setPtHatRange(ptHatSample);

    // Setup input filename
    TString fname;
    if (input == nullptr || input[0] == '\0') {
        fname = Form("~/cernbox/ana/pPb8160/%s/%s/forest/HiForestSkim_pPb_MC_pthat%d_%s_%s.root", 
                     generator.Data(), directionStr.Data(), ptHatSample, directionStr.Data(), tag.Data());
    }
    else {
        fname = input;
    }
 
    // Setup output filename
    TString oFileName;
    if (output == nullptr || output[0] == '\0') {
        oFileName = Form("macro/eta_shift/%s_%s_ptHat%d.root", generator.Data(), directionStr.Data(), ptHatSample);
    }
    else {
        oFileName = output;
    }

    // Print code input parameters here
    std::cout << "Collision system: " << collisionSystem << std::endl;
    std::cout << "Code input parameters:" << std::endl;
    std::cout << "  generator: " << GREEN << generator.Data() << RESET << std::endl;
    std::cout << "  direction: " << GREEN << directionStr.Data() << RESET << std::endl;
    std::cout << "  tag: " << GREEN << tag.Data() << RESET << std::endl;
    std::cout << "  ptHatSample: " << GREEN << ptHatSample << RESET << std::endl;
    std::cout << "  isMc: " << GREEN << (isMc?"true":"false") << RESET << std::endl;
    std::cout << "  isPbGoing: " << GREEN << (isPbGoing?"true":"false") << RESET << std::endl;
    std::cout << "  mcType: " << GREEN << (mcType?"MC":"Data") << RESET << std::endl;
    std::cout << "  isPythia: " << GREEN << (isPythia?"true":"false") << RESET << std::endl;
    std::cout << "  input override: " << GREEN << (input?input:"(null)") << RESET << std::endl;
    std::cout << "  output override: " << GREEN << (output?output:"(null)") << RESET << std::endl;
    std::cout << "  triggerId: " << GREEN << triggerId << RESET << std::endl;
    std::cout << "  recoJetSelMethod: " << GREEN << recoJetSelMethod << RESET << std::endl;
    const std::array<const char *, 3> jerModeNames{"hybrid", "stochastic", "truth-to-pseudo-reco"};
    std::cout << "  JER smearing mode: " << GREEN << jerModeNames[jerSmearingMode] << RESET << std::endl;
    std::cout << "  JER eta dependence: " << GREEN << (useEtaDependentJER ? "enabled" : "disabled") << RESET << std::endl;
    std::cout << "  JEC files:" << std::endl;
    for (const auto &fn : fJECFileNames) std::cout << "    " << GREEN << fn.c_str() << RESET << std::endl;

    // TFile *f = TFile::Open( fname );

    { TDirectory::TContext detached(nullptr); createHistograms(hs, isMc); }
    jecOwner=std::make_unique<JetCorrector>(fJECFileNames);
    jeuOwner=std::make_unique<JetUncertainty>(fJEUFileName);
    this->isPp=isPp; this->isMc=isMc; this->isPbGoing=isPbGoing; this->isPythia=isPythia;
    this->triggerId=triggerId; this->jetSelectionMethod=recoJetSelMethod;
    this->vertexFilterSelection=vertexFilterSelection;
    this->jerSmearingMode=static_cast<JERSmearingMode>(jerSmearingMode);
    this->useEtaDependentJER=useEtaDependentJER!=0;
    nEntries=context.entries;
    this->output=oFileName.Data(); exclusions=config.excludeOutputPrefixes;
    initializeSmearing();
}
void finish() {
    reportCounts();
    timings.measure(StageTimings::Output, [&] {
        jetanalysis::writeOutput(output,hs,isMc,static_cast<int>(jerSmearingMode),useEtaDependentJER,
            OutputSelection{exclusions});
    });
    timings.report();
}

};
} // namespace
struct DiJetAnalysis::Impl {
    Physics physics;
    bool initialized=false,finished=false;
};
DiJetAnalysis::DiJetAnalysis(AnalysisConfig config):config_(std::move(config)),impl_(std::make_unique<Impl>()) {}
DiJetAnalysis::~DiJetAnalysis()=default;
void DiJetAnalysis::initialize(const RunContext& context) {
    if(impl_->initialized) throw std::logic_error("DiJetAnalysis is single-use");
    impl_->initialized=true;
    impl_->physics.initialize(config_,context);
}
void DiJetAnalysis::process(const Event& event) {
    if(!impl_->initialized || impl_->finished) throw std::logic_error("DiJetAnalysis is not active");
    impl_->physics.processEvent(event);
}
void DiJetAnalysis::finish() {
    if(!impl_->initialized || impl_->finished) throw std::logic_error("DiJetAnalysis is not active");
    impl_->finished=true;
    impl_->physics.finish();
}
void DiJetAnalysis::run(EventReader& reader) {runTask(reader,*this,config_.profileStages);}
} // namespace jetanalysis
