#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"
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
#include <utility>

const char* RED    = "\033[1;31m";
const char* GREEN  = "\033[1;32m";
const char* RESET  = "\033[0m";

//________________
void usage() {
    std::cout << GREEN << "Usage: ./processForestSimple  [options]\n"
                << "Options:\n"
                << "  <input_file_or_list>  Input ROOT file or file list (default: empty)\n"
                << "  <output_file>         Output ROOT file name (default: empty)\n"
                << "  <0|1|2>               MC type (default 2): 0 - data, 1 - embedding, 2 - pythia\n"
                << "  <0|1>                 Direction (default 1): 0 - p-going, 1 - Pb-going\n"
                << "  <value>               pT hat sample for MC (default: 30)\n"
                << "  <value>               Trigger ID (default 0): 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100\n"
                << "  <value>               Reco jet selection method (default 2): 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetIdTightLeptVeto, 3 - jetIdLoose\n"
                << RESET;
    exit(EXIT_FAILURE);
}

// Eta shifts for pPb 8.16 TeV collisions, used for etaCM calculation
const int nEtaShifts = 13;
static constexpr std::array<float, nEtaShifts> etaShift{0.460, 0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470, 0.475, 0.480, 0.485, 0.490 };
const int nEtaCuts = 8;
static constexpr std::array<float, nEtaCuts> etaCuts{1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.5, 3.0};

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

std::vector<float> etaBinEdges = {-3.60f, -3.46f, -3.31f, -3.17f, -3.02f, -2.88f, -2.74f, -2.59f, -2.45f, -2.30f, -2.16f, -2.02f, -1.87f, -1.73f, -1.58f, -1.44f, -1.30f, -1.15f, -1.01f, -0.86f, -0.72f, -0.58f, -0.43f, -0.29f, -0.14f, 0.00f, 0.14f, 0.29f, 0.43f, 0.58f, 0.72f, 0.86f, 1.01f, 1.15f, 1.30f, 1.44f, 1.58f, 1.73f, 1.87f, 2.02f, 2.16f, 2.30f, 2.45f, 2.59f, 2.74f, 2.88f, 3.02f, 3.17f, 3.31f, 3.46f, 3.60f};

// Range for the selection
float ptHatRange[2] { 0.f, 8160.f};

// Event info variables
float ptHat;
float vz;

// Skim/event filter variables
int pBeamScrapingFilter;
int pPAprimaryVertexFilter;
int HBHENoiseFilterResultRun2Loose;
int phfCoincFilter;
int pVertexFilterCutdz1p0;
int pVertexFilterCutGplus;
int pVertexFilterCutVtx1;

// Trigger variables
int HLT_PAAK4PFJet60_Eta5p1_v4;
int HLT_PAAK4PFJet80_Eta5p1_v3;
int HLT_PAAK4PFJet100_Eta5p1_v3;

// Maximum number of jets to store info for (adjust as needed)
const int maxJets = 1500;

// Gen jet info
int   nGenJets;           // number of generated jets
float genJetPt[maxJets];  // generated jet pT
float genJetEta[maxJets]; // generated jet eta
float genJetPhi[maxJets]; // generated jet phi

// Reco and ref jet info
int   nRecoJets;                  // number of reconstructed jets

float recoJetPtRaw[maxJets];      // reconstructed jet pT without JEC
float recoJetTrackMaxPt[maxJets]; // reconstructed track pT max in the jet
float recoJetEta[maxJets];        // reconstructed jet eta
float recoJetPhi[maxJets];        // reconstructed jet phi
float recoJetPfNHF[maxJets];      // reconstructed jet neutral hadron fraction
float recoJetPfNEF[maxJets];      // reconstructed jet neutral electromagnetic fraction
float recoJetPfCHF[maxJets];      // reconstructed jet charged hadron fraction
float recoJetPfMUF[maxJets];      // reconstructed jet muon fraction
float recoJetPfCEF[maxJets];      // reconstructed jet charged electromagnetic fraction
int   recoJetPfCHM[maxJets];      // reconstructed jet charged muon fraction
int   recoJetPfCEM[maxJets];      // reconstructed jet charged electromagnetic fraction
int   recoJetPfNHM[maxJets];      // reconstructed jet neutral hadron fraction
int   recoJetPfNEM[maxJets];      // reconstructed jet neutral electromagnetic fraction
int   recoJetPfMUM[maxJets];      // reconstructed jet muon fraction

float refJetPt[maxJets];   // reference jet pT
float refJetEta[maxJets];  // reference jet eta
float refJetPhi[maxJets];  // reference jet phi

//________________
struct JERSmearingFuncVsEta{
    double fPtLow;
    double fPtHigh;
    std::unique_ptr<TF1> fJERSmearFuncVsEta;
};

//________________
class JERSmearingHelper {
    public:
        /// Constructor
        JERSmearingHelper() {
            fJERSmearFuncVsPt = std::make_unique<TF1>("fJERSmearFuncVsPt", "sqrt( [0]*[0] + [1]*[1]/x )", 30., 800.);
            fJERSmearFuncVsPt->SetParameters(0.04360, 0.81815); // Default parameters (Pb-going direction for pPb 8.16 TeV collisions)
            fRndm = std::make_unique<TRandom3>(0); // Seed with 0 for random seed
        }

        /// Destructor
        ~JERSmearingHelper() {
            if (fJERSmearFuncVsPt) fJERSmearFuncVsPt.reset();
            if (fRndm) fRndm.reset();
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
                return 1.0; // No smearing for |eta| < 0.8 or |eta| > 2.6
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

        /// @brief Returns the JER smeared pT for a given
        /// @param ptOriginal The original pT of the jet
        /// @param recoEta The reconstructed eta of the jet
        /// @param ptSmeared The smeared pT of the jet (output)
        /// @param variation The variation for JER scaling factor:
        ///                  -99: default JER (SF = 1.0), 
        ///                  -1: JER down with SF, 
        ///                   0: default JER with SF, 
        ///                   1: JER up with SF, 
        ///                   2: use eta-dependent JER scaling factor
        ///                   3: use provided scaling factor
        ///                  19: JER down with SF x eta-dependent scaling factor
        ///                  20: default JER with SF x eta-dependent scaling factor
        ///                  21: JER up with SF x eta-dependent scaling factor
        /// @param scaleFactor The JER scaling factor to be used if variation is set to 2
        void smearMomentum(const float &ptOriginal, const float &eta, float &ptSmeared, const int &variation = -99, const float &scaleFactor = 1.0) const {
    
            double sigmaSmear{0.};
            double evalValue{0.};
            double scaleFactorToUse{1.0};

            if ( !fJERSmearFuncVsPt ) {
                throw std::runtime_error("JERSmearingHelper: JER smearing function not set.");
            }
            if ( !fRndm ) {
                throw std::runtime_error("JERSmearingHelper: Random number generator not set.");
            }

            if (variation == -99) { // Default JER (no scaling)
                scaleFactorToUse = 1.0; 
            }
            else if (variation == -1) { // JER down with SF
                scaleFactorToUse = JERScaleFactorSystematics(eta, -1); 
            }
            else if (variation == 0) {  // Default JER with SF
                scaleFactorToUse = JERScaleFactorSystematics(eta, 0); 
            }
            else if (variation == 1) { // JER up with SF
                scaleFactorToUse = JERScaleFactorSystematics(eta, 1);
            }
            else if (variation == 2) { // Use eta-dependent JER scaling factor
                scaleFactorToUse = evalJERScalingFactorForPtEta(ptOriginal, eta);
            }
            else if (variation == 3) { // Use provided scaling factor
                scaleFactorToUse = scaleFactor; 
            }
            else if (variation == 19) { // JER down with SF x eta-dependent scaling factor
                scaleFactorToUse = JERScaleFactorSystematics(eta, -1) * evalJERScalingFactorForPtEta(ptOriginal, eta);
            }
            else if (variation == 20) { // Default JER with SF x eta-dependent scaling factor
                scaleFactorToUse = JERScaleFactorSystematics(eta, 0) * evalJERScalingFactorForPtEta(ptOriginal, eta);
            }
            else if (variation == 21) { // JER up with SF x eta-dependent scaling factor
                scaleFactorToUse = JERScaleFactorSystematics(eta, 1) * evalJERScalingFactorForPtEta(ptOriginal, eta);
            }
            else {
                throw std::invalid_argument("Invalid variation value for JER smearing.");
            }


            if ( ptOriginal <= 30.) {
                evalValue = fJERSmearFuncVsPt->Eval( 31. );
            }
            else if ( ptOriginal >= 800 ) {
                evalValue = fJERSmearFuncVsPt->Eval( 799. );
            }
            else {
                evalValue = fJERSmearFuncVsPt->Eval( ptOriginal );
            }
            sigmaSmear = scaleFactorToUse * evalValue;

            double smearFactor = 1. + fRndm->Gaus( 0., sigmaSmear );
            while ( smearFactor < 0. ) {
                smearFactor = 1. + fRndm->Gaus( 0., sigmaSmear );
            }
            ptSmeared = ptOriginal * smearFactor;
        }

    private:
        // JER as a function of pt for -0.8 < eta < 0.8
        std::unique_ptr<TF1> fJERSmearFuncVsPt;
        // Random number generator for smearing
        std::unique_ptr<TRandom3> fRndm;
        // Vector to hold JER smearing functions for different pt and eta ranges
        std::vector<JERSmearingFuncVsEta> fJERSmearFuncVsEtaVec;
        // JER scaling factors for different eta ranges
        std::vector<float> fJerEtaLow{-5.191, -3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139};
        std::vector<float> fJerEtaHi{-3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139, 5.191};
        std::vector<float> fJerDef{1.1922, 1.1869, 1.7788, 1.3418, 1.2963, 1.1512, 1.1426, 1.1000, 1.1278, 1.1609, 1.1464, 1.1948, 1.15958, 1.15958, 1.1948, 1.1464, 1.1609, 1.1278, 1.1000, 1.1426, 1.1512, 1.2963, 1.3418, 1.7788, 1.1869, 1.1922};
        std::vector<float> fJerLow{1.0434, 1.0626, 1.578, 1.1327, 1.0592, 1.0372, 1.0212, 0.9921, 1.0292, 1.0584, 1.0832, 1.1296, 1.095, 1.095, 1.1296, 1.0832, 1.0584, 1.0292, 0.9921, 1.0212, 1.0372, 1.0592, 1.1327, 1.578, 1.0626, 1.0434};
        std::vector<float> fJerHi{1.341, 1.3112, 1.9796, 1.5509, 1.5334, 1.2652, 1.264, 1.2079, 1.2264, 1.2634, 1.2096, 1.26, 1.224, 1.224, 1.26, 1.2096, 1.2634, 1.2264, 1.2079, 1.264, 1.2652, 1.5334, 1.5509, 1.9796, 1.3112, 1.341};
};

//________________
/// Function to calculate the resolution smearing for a given reference pT and reconstructed eta
/// @param refPt The reference pT of the jet
/// @param recoEta The reconstructed eta of the jet
/// @param ptReco The reconstructed pT of the jet
/// @param scaleFactor The scaling factor for the JER
/// @param fJERSmearFunc The JER smearing function
/// @param fRndm The random number generator
void calculateResolutionSmearing(const float &refPt, const float &recoEta, float &ptReco, const float &scaleFactor,
                                 TF1 *fJERSmearFunc, TRandom3 *fRndm) {
    
    double sigmaSmear{0.};
    double evalValue{0.};
    if ( refPt <= 30.) {
        evalValue = fJERSmearFunc->Eval( 31. );
    }
    else if ( refPt >= 800 ) {
        evalValue = fJERSmearFunc->Eval( 799. );
    }
    else {
        evalValue = fJERSmearFunc->Eval( refPt );
    }
    sigmaSmear = scaleFactor * evalValue;

    double smearFactor = 1. + fRndm->Gaus( 0., sigmaSmear );
    while ( smearFactor < 0. ) {
        smearFactor = 1. + fRndm->Gaus( 0., sigmaSmear );
    }
    ptReco *= smearFactor;
}

//________________
struct GenJet {
    float pt;
    float eta;
    float phi;
    float smearedPt;
};

//________________
struct RecoJet {
    float recoPtRaw;
    float recoPt;
    float recoEta;
    float recoPhi;
    float recoPfNHF;
    float recoPfNEF;
    float recoPfCHF;
    float recoPfMUF;
    float recoPfCEF;
    int  recoPfCHM;
    int  recoPfCEM;
    int  recoPfNHM;
    int  recoPfNEM;
    int  recoPfMUM;
    float recoTrackMaxPt;
    float refPt;
    float refEta;
    float refPhi;

    float recoPtJeuUp;
    float recoPtJeuDown;

    float refPtSmeared; // reference jet pT after smearing
};

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

//________________
/// Define the histograms for the analysis
struct Histograms {
    //
    // Event level histograms
    //
    std::unique_ptr<TH1D> hVzRaw;
    std::unique_ptr<TH1D> hPtHatUnweighted;
    std::unique_ptr<TH1D> hPtHat;
    std::unique_ptr<TH1D> hVzUnweighted;
    std::unique_ptr<TH1D> hVz;

    std::unique_ptr<TH2D> hGenDijetPtAveOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hGenDijetPtAveOverPtHatVsPtHatPass;
    std::unique_ptr<TH2D> hGenLeadJetPtOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hGenLeadJetPtOverPtHatVsPtHatPass;
    std::unique_ptr<TH2D> hRecoDijetPtAveOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hRecoDijetPtAveOverPtHatVsPtHatPass;
    std::unique_ptr<TH2D> hRecoLeadJetPtOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hRecoLeadJetPtOverPtHatVsPtHatPass;

    //
    // Gen level histograms
    //

    // Gen jet histograms
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaLabUnflippedUnweighted;
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaLab;
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaCM;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hGenInclusiveJetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaStdBins;
    std::unique_ptr<TH2D> hGenInclusiveJetPtPtHat;

    std::unique_ptr<TH3D> hGenInclusiveJetJESDefPtEta;
    std::unique_ptr<TH3D> hGenInclusiveJetJESDefExtraPtEta;
    std::unique_ptr<TH3D> hGenInclusiveJetJESDefDoublePtEta;
    std::unique_ptr<TH3D> hGenInclusiveJetJESDefNinetyPtEta;

    std::unique_ptr<TH2D> hGenInclusiveJetDefPtEtaFlipped;
    std::unique_ptr<TH2D> hGenInclusiveJetDefExtraPtEtaFlipped;
    std::unique_ptr<TH2D> hGenInclusiveJetDefDoublePtEtaFlipped;
    std::unique_ptr<TH2D> hGenInclusiveJetDefNinetyPtEtaFlipped;

    // Gen dijet histograms
    std::unique_ptr<TH1D> hGenDijetPtAve;
    std::unique_ptr<TH1D> hGenDijetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMForwardShifted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMBackwardShifted[nEtaShifts];

    std::unique_ptr<TH2D> hGenDijetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hGenDijetPtEtaLab;
    std::unique_ptr<TH2D> hGenDijetPtEtaCM;

    std::unique_ptr<TH2D> hGenDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaBackwardArr[nEtaCuts];

    // Gen dijet smeared histograms
    std::unique_ptr<TH2D> hGenDijetDefPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefPtEtaBackwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefExtraPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefExtraPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefExtraPtEtaBackwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefDoublePtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefDoublePtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefDoublePtEtaBackwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefNinetyPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefNinetyPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefNinetyPtEtaBackwardArr[nEtaCuts]; 

    // Unfolding histograms (response matrices and missing histograms)
    std::unique_ptr<THnSparseF> hGenDijetPtEtaCMVsRecoPtEtaCMArr[nEtaCuts];
    // std::unique_ptr<THnSparseF> hRefDijetPtEtaCMVsRecoPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMissArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMFakeArr[nEtaCuts];
    std::unique_ptr<THnSparseF> hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMissJerDefExtraArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMFakeJerDefExtraArr[nEtaCuts];
    std::unique_ptr<TH1D> hUnfoldingPairClassificationArr[nEtaCuts];
    std::unique_ptr<TH1D> hUnfoldingPairClassificationJerDefExtraArr[nEtaCuts];

    //
    // Ref-selected level histograms
    //
    std::unique_ptr<TH2D> hRefSelDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefSelDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefSelDijetPtEtaBackwardArr[nEtaCuts];

    //
    // Reco level histograms
    //

    // Reco jet histograms
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabUnflippedUnweighted;
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLab;
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaCM;
    std::unique_ptr<TH2D> hRecoInclusiveJetTrkMaxPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRecoInclusiveJetTrkMaxPtEtaLab;
    std::unique_ptr<TH2D> hRecoInclusiveJetTrkMaxPtEtaCM;
    std::unique_ptr<TH2D> hRecoInclusiveJetNoSelPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRecoInclusiveJetNoSelPtEtaLab;
    std::unique_ptr<TH2D> hRecoInclusiveJetNoSelPtEtaCM;
    // std::unique_ptr<TH1D> hRecoInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
    // std::unique_ptr<TH1D> hRecoInclusiveJetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaStdBins;
    std::unique_ptr<TH3D> hRecoInclusiveJetJESPtEta;
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDefNoSFPtEta;   // JER default variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDefDoublePtEta;   // JER default variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDefNinetyPtEta;   // JER default variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDefPtEta;   // JER default variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESUpPtEta;    // JER up variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDownPtEta;  // JER down variation for reco inclusive jet histograms
    std::unique_ptr<TH2D> hRecoInclusiveJetPtPtHat;

    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabMatched;
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabUnmatched;

    // Reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRecoDijetPtEtaLab;
    std::unique_ptr<TH2D> hRecoDijetPtEtaCM;

    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArr[nEtaCuts];
    
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefNoSF[nEtaCuts];             // JER no SF variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDefNoSF[nEtaCuts];        // JER no SF variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDefNoSF[nEtaCuts];       // JER no SF variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefDoubleSF[nEtaCuts];         // JER default variation with double SF for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDefDoubleSF[nEtaCuts];    // JER default variation with double SF for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDefDoubleSF[nEtaCuts];   // JER default variation with double SF for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefNinetySF[nEtaCuts];         // JER default variation with 90% SF for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDefNinetySF[nEtaCuts];    // JER default variation with 90% SF for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDefNinetySF[nEtaCuts];   // JER default variation with 90% SF for reco dijet histograms

    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerUp[nEtaCuts];         // JER up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerUp[nEtaCuts];    // JER up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerUp[nEtaCuts];   // JER up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDown[nEtaCuts];       // JER down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDown[nEtaCuts];  // JER down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDown[nEtaCuts]; // JER down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDef[nEtaCuts];        // JER default variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDef[nEtaCuts];   // JER default variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDef[nEtaCuts];  // JER default variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDownExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDownExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDownExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDefExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDefExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerUpExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerUpExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerUpExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefUnfold[nEtaCuts];        // JER default variation for reco dijet histograms for unfolding
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDefUnfold[nEtaCuts];   // JER default variation for reco dijet histograms for unfolding
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDefUnfold[nEtaCuts];  // JER default variation for reco dijet histograms for unfolding
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefExtraUnfold[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDefExtraUnfold[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDefExtraUnfold[nEtaCuts];

    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJeuUp[nEtaCuts];         // JEU up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJeuUp[nEtaCuts];    // JEU up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJeuUp[nEtaCuts];   // JEU up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJeuDown[nEtaCuts];       // JEU down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJeuDown[nEtaCuts];  // JEU down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJeuDown[nEtaCuts]; // JEU down variation for reco dijet histograms

    //
    // Ref level histograms
    //

    // Ref jet histograms
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaLab;
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaCM;
    std::unique_ptr<TH2D> hRefInclusiveJetTrkMaxPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefInclusiveJetTrkMaxPtEtaLab;
    std::unique_ptr<TH2D> hRefInclusiveJetTrkMaxPtEtaCM;
    std::unique_ptr<TH2D> hRefInclusiveJetNoSelPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefInclusiveJetNoSelPtEtaLab;
    std::unique_ptr<TH2D> hRefInclusiveJetNoSelPtEtaCM;
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaStdBins;

    // Ref dijet histograms

    std::unique_ptr<TH2D> hRefDijetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefDijetPtEtaLab;
    std::unique_ptr<TH2D> hRefDijetPtEtaCM;

    std::unique_ptr<TH2D> hRefDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArr[nEtaCuts];

    std::unique_ptr<TH2D> hRefDijetPtEtaCMArrJerUp[nEtaCuts];         // JER up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJerUp[nEtaCuts];    // JER up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJerUp[nEtaCuts];   // JER up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaCMArrJerDown[nEtaCuts];       // JER down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJerDown[nEtaCuts];  // JER down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJerDown[nEtaCuts]; // JER down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaCMArrJerDef[nEtaCuts];        // JER default variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJerDef[nEtaCuts];   // JER default variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJerDef[nEtaCuts];  // JER default variation for ref dijet histograms

};

//________________
// Event weight calculation for pPb8160 MC samples
/// @param ptHat: the ptHat
/// @param vz: the z position
/// @param fVzWeight: the z position weight function
/// @param eventsInSample: the number of events in the sample
/// @return the event weight
double eventWeight(const float &ptHat, const float &vz, TF1 &fVzWeight, Long64_t eventsInSample) {

    double genWeight = 1.0;
    double vzWeight = 1./ fVzWeight.Eval( vz );
    // vzWeight = 1.;
    // Magic numbers are (cross section x Nevents generated)
    if (ptHat > 15.0 && ptHat <= 30.)       { genWeight = 1.0404701e-06 * 961104 ; }
    else if (ptHat > 30. && ptHat <= 50.)   { genWeight = 7.7966624e-08 * 952110 ; }
    else if (ptHat > 50. && ptHat <= 80.)   { genWeight = 1.0016052e-08 * 952554 ; }
    else if (ptHat > 80. && ptHat <= 120.)  { genWeight = 1.3018269e-09 * 996844 ; }
    else if (ptHat > 120.&& ptHat <= 170.)  { genWeight = 2.2648493e-10 * 964681 ; }
    else if (ptHat > 170. && ptHat <= 220.) { genWeight = 4.0879112e-11 * 999260 ; }
    else if (ptHat > 220. && ptHat <= 280.) { genWeight = 1.1898939e-11 * 964336 ; }
    else if (ptHat > 280. && ptHat <= 370.) { genWeight = 3.3364433e-12 * 995036 ; }
    else if (ptHat > 370. && ptHat <= 460.) { genWeight = 7.6612402e-13 * 958160 ; }
    else if (ptHat > 460. && ptHat <= 540.) { genWeight = 2.1341026e-13 * 981427 ; }
    else if (ptHat > 540.)                  { genWeight = 7.9191586e-14 * 1000000; }
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
void createHistograms(Histograms &hs, const bool &isMc = false) {


    std::cout << "Creating histograms...";
    // Standard eta binning for L2L3 corrections
    double jetEtaL2L3StdVals[] = { -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, 
                                 -4.013, -3.839, -3.664, -3.489, -3.314, -3.139, 
                                 -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
                                 -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, 
                                 -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, 
                                 -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
                                 -0.435, -0.348, -0.261, -0.174, -0.087,  0.000, 
                                  0.087,  0.174,  0.261,  0.348,  0.435,  0.522, 
                                  0.609,  0.696,  0.783,  0.879,  0.957,  1.044, 
                                  1.131,  1.218,  1.305,  1.392,  1.479,  1.566, 
                                  1.653,  1.740,  1.830,  1.930,  2.043,  2.172, 
                                  2.322,  2.500,  2.650,  2.853,  2.964,  3.139, 
                                  3.314,  3.489,  3.664,  3.839,  4.013,  4.191, 
                                  4.363,  4.538,  4.716,  4.889,  5.191 };
    int jetEtaL2L3StdBins = sizeof(jetEtaL2L3StdVals)/sizeof(double) - 1;

    // Single jet pt binning
    const int nJetPtBins = 100;
    double jetPtBins[] = { 10., 1010. };
    const int nJetEtaBins = 160;
    double jetEtaBins[] = { -4.0, 4.0 };
    const int nJetJESBins = 50;
    double jetJESBins[] = { 0., 2. };

    // Dijet binning
    const int nDijetPtBins = 100;
    double dijetPtBins[] = { 0., 1000.};
    const int nDijetEtaBins = 72;
    double dijetEtaBins[] = { -3.6, 3.6 };
    const int nDijetEtaFBBins = 36;
    double dijetEtaFBBins[] = { 0., 3.6 };

    const int nDimResponse = 4;
    int responseBins[nDimResponse] = { nDijetPtBins, nDijetEtaBins, nDijetPtBins, nDijetEtaBins };
    double responseBinsMin[nDimResponse] = { dijetPtBins[0], dijetEtaBins[0], dijetPtBins[0], dijetEtaBins[0] };
    double responseBinsMax[nDimResponse] = { dijetPtBins[1], dijetEtaBins[1], dijetPtBins[1], dijetEtaBins[1] };

    //
    // Event level histograms
    //
    hs.hVzRaw = std::make_unique<TH1D>("hVzRaw", 
                                 "vz raw;vz (cm);dN/dvz", 
                                 60, -15., 15.);
    hs.hVzRaw->Sumw2();
    hs.hPtHatUnweighted = std::make_unique<TH1D>("hPtHatUnweighted", 
                                              "#hat{p}_{T} unweighted;#hat{p}_{T} (GeV);dN/d#hat{p}_{T}", 
                                              100, 0., 1000.);
    hs.hPtHatUnweighted->Sumw2();
    hs.hPtHat = std::make_unique<TH1D>("hPtHat", 
                                        "#hat{p}_{T};#hat{p}_{T} (GeV);dN/d#hat{p}_{T}", 
                                        100, 0., 1000.);
    hs.hPtHat->Sumw2();


    hs.hVzUnweighted = std::make_unique<TH1D>("hVzUnweighted", 
                                           "vz unweighted;vz (cm);dN/dvz", 
                                           60, -15., 15.);
    hs.hVzUnweighted->Sumw2();
    hs.hVz = std::make_unique<TH1D>("hVz", 
                                 "vz;vz (cm);dN/dvz", 
                                 60, -15., 15.);
    hs.hVz->Sumw2();

    if (isMc) {

        //
        // Event level histograms
        //
        hs.hGenDijetPtAveOverPtHatVsPtHat = std::make_unique<TH2D>("hGenDijetPtAveOverPtHatVsPtHat", 
                                            "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenDijetPtAveOverPtHatVsPtHat->Sumw2();
        hs.hGenDijetPtAveOverPtHatVsPtHatPass = std::make_unique<TH2D>("hGenDijetPtAveOverPtHatVsPtHatPass", 
                                            "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenDijetPtAveOverPtHatVsPtHatPass->Sumw2();
        hs.hGenLeadJetPtOverPtHatVsPtHat = std::make_unique<TH2D>("hGenLeadJetPtOverPtHatVsPtHat", 
                                            "Gen leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Gen p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenLeadJetPtOverPtHatVsPtHat->Sumw2();
        hs.hGenLeadJetPtOverPtHatVsPtHatPass = std::make_unique<TH2D>("hGenLeadJetPtOverPtHatVsPtHatPass", 
                                            "Gen leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Gen p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenLeadJetPtOverPtHatVsPtHatPass->Sumw2();
        hs.hRecoDijetPtAveOverPtHatVsPtHat = std::make_unique<TH2D>("hRecoDijetPtAveOverPtHatVsPtHat", 
                                            "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoDijetPtAveOverPtHatVsPtHat->Sumw2();
        hs.hRecoDijetPtAveOverPtHatVsPtHatPass = std::make_unique<TH2D>("hRecoDijetPtAveOverPtHatVsPtHatPass", 
                                            "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoDijetPtAveOverPtHatVsPtHatPass->Sumw2();
        hs.hRecoLeadJetPtOverPtHatVsPtHat = std::make_unique<TH2D>("hRecoLeadJetPtOverPtHatVsPtHat", 
                                            "Reco leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Reco p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoLeadJetPtOverPtHatVsPtHat->Sumw2();
        hs.hRecoLeadJetPtOverPtHatVsPtHatPass = std::make_unique<TH2D>("hRecoLeadJetPtOverPtHatVsPtHatPass", 
                                            "Reco leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T} (pass selection);#hat{p}_{T} (GeV);Reco p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoLeadJetPtOverPtHatVsPtHatPass->Sumw2();

        //
        // Gen level histograms
        //

        // Gen jet histograms
        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLabUnflippedUnweighted", 
                                                "Gen jet #eta (lab frame, unflipped, unweighted) vs p_{T};p_{T} (GeV);#eta",
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted->Sumw2();
        hs.hGenInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLabUnflipped", 
                                                "Gen jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta",
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hGenInclusiveJetPtEtaLab = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLab", 
                                                "Gen jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta", 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaLab->Sumw2();
        hs.hGenInclusiveJetPtEtaCM = std::make_unique<TH2D>("hGenInclusiveJetPtEtaCM", 
                                                "Gen jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM}", 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetPtEtaCM->Sumw2();
        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
                                                Form("Gen jet #eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
            hs.hGenInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%d", iShift), 
                                                Form("Gen jet #eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Sumw2();
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift) 


        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%.3f", etaShift[0]), 
                                                Form("Gen jet #eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0]->Sumw2();
        hs.hGenInclusiveJetEtaCMShifted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%.3f", etaShift[0]), 
                                                                    Form("Gen jet #eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShifted[0]->Sumw2();
        hs.hGenInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hGenInclusiveJetPtEtaStdBins", 
                                                                "Gen jet #eta (standard bins) vs p_{T};p_{T} (GeV);#eta",
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hs.hGenInclusiveJetPtEtaStdBins->Sumw2();
        hs.hGenInclusiveJetPtPtHat = std::make_unique<TH2D>("hGenInclusiveJetPtPtHat", 
                                                            "Gen jet #hat{p}_{T} vs p_{T};p_{T} (GeV);#hat{p}_{T} (GeV)",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtPtHat->Sumw2();

        hs.hGenInclusiveJetJESDefPtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefPtEta", 
                                                                "Gen jet JES default variation pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefPtEta->Sumw2();
        hs.hGenInclusiveJetJESDefExtraPtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefExtraPtEta", 
                                                                "Gen jet JES default variation with extra SF pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefExtraPtEta->Sumw2();
        hs.hGenInclusiveJetJESDefDoublePtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefDoublePtEta", 
                                                                "Gen jet JES default variation with double SF pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefDoublePtEta->Sumw2();
        hs.hGenInclusiveJetJESDefNinetyPtEta = std::make_unique<TH3D>("hGenInclusiveJetJESDefNinetyPtEta", 
                                                                "Gen jet JES default variation with 90%% SF pT vs #eta;p_{T}^{gen.smeared}/p_{T};p_{T}^{gen} (GeV);#eta",
                                                                nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetJESDefNinetyPtEta->Sumw2();

        hs.hGenInclusiveJetDefPtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefPtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefPtEtaFlipped->Sumw2();
        hs.hGenInclusiveJetDefExtraPtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefExtraPtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped, extra SF);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefExtraPtEtaFlipped->Sumw2();
        hs.hGenInclusiveJetDefDoublePtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefDoublePtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped, double SF);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefDoublePtEtaFlipped->Sumw2();
        hs.hGenInclusiveJetDefNinetyPtEtaFlipped = std::make_unique<TH2D>("hGenInclusiveJetDefNinetyPtEtaFlipped", 
                                                "Gen jet pT vs eta (flipped, 90%% SF);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetDefNinetyPtEtaFlipped->Sumw2();
        

        // Gen dijet histograms
        hs.hGenDijetPtAve = std::make_unique<TH1D>("hGenDijetPtAve", 
                                                    "Gen dijet pT average;p_{T}^{ave} (GeV);dN/dp_{T}^{ave}", 
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1]);
        hs.hGenDijetPtAve->Sumw2();
        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenDijetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMShiftedUnweighted_%d", iShift), 
                                                                                Form("Gen dijet eta (CM frame, shifted by %.3f, unweighted);#eta;dN/d#eta", etaShift[iShift]), 
                                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Sumw2();
            hs.hGenDijetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMShifted_%d", iShift), 
                                                                        Form("Gen dijet eta (CM frame, shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                                                        nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetEtaCMShifted[iShift]->Sumw2();
            hs.hGenDijetEtaCMForwardShifted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMForwardShifted_%d", iShift), 
                                                                            Form("Gen dijet eta (CM frame, forward shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetEtaCMForwardShifted[iShift]->Sumw2();
            hs.hGenDijetEtaCMBackwardShifted[iShift] = std::make_unique<TH1D>(Form("hGenDijetEtaCMBackwardShifted_%d", iShift), 
                                                                            Form("Gen dijet eta (CM frame, backward shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetEtaCMBackwardShifted[iShift]->Sumw2();
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        hs.hGenDijetPtEtaLabUnflipped = std::make_unique<TH2D>("hGenDijetPtEtaLabUnflipped", 
                                                               "Gen dijet #eta (lab frame, unflipped) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta",
                                                               nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                               nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hGenDijetPtEtaLabUnflipped->Sumw2();
        hs.hGenDijetPtEtaLab = std::make_unique<TH2D>("hGenDijetPtEtaLab", 
                                                      "Gen dijet #eta (lab frame) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta",
                                                      nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                      nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hGenDijetPtEtaLab->Sumw2();
        hs.hGenDijetPtEtaCM = std::make_unique<TH2D>("hGenDijetPtEtaCM", 
                                                     "Gen dijet #eta (CM frame) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", 
                                                     nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                     nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hGenDijetPtEtaCM->Sumw2();

        // Loop over eta cuts for dijet histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hGenDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaForward_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetPtEtaForwardArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaBackward_%d", iCut), 
                                                    Form("Gen dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetPtEtaBackwardArr[iCut]->Sumw2();

            // Smearing studies
            hs.hGenDijetDefPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (def smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetDefPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefPtEtaForward_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (def smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefPtEtaForwardArr[iCut]->Sumw2();
            hs.hGenDijetDefPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefPtEtaBackward_%d", iCut), 
                                                    Form("Gen dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (def smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefPtEtaBackwardArr[iCut]->Sumw2();

            hs.hGenDijetDefExtraPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefExtraPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (extra smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefExtraPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetDefExtraPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefExtraPtEtaForward_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (extra smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefExtraPtEtaForwardArr[iCut]->Sumw2();
            hs.hGenDijetDefExtraPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefExtraPtEtaBackward_%d", iCut), 
                                                    Form("Gen dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (extra smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefExtraPtEtaBackwardArr[iCut]->Sumw2();

            hs.hGenDijetDefDoublePtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefDoublePtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (def 2x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefDoublePtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetDefDoublePtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefDoublePtEtaForward_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (def 2x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefDoublePtEtaForwardArr[iCut]->Sumw2();
            hs.hGenDijetDefDoublePtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefDoublePtEtaBackward_%d", iCut), 
                                                    Form("Gen dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (def 2x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefDoublePtEtaBackwardArr[iCut]->Sumw2();
            hs.hGenDijetDefNinetyPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefNinetyPtEtaCM_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (def 0.9x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetDefNinetyPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetDefNinetyPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefNinetyPtEtaForward_%d", iCut), 
                                                Form("Gen dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (def 0.9x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefNinetyPtEtaForwardArr[iCut]->Sumw2();
            hs.hGenDijetDefNinetyPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetDefNinetyPtEtaBackward_%d", iCut), 
                                                    Form("Gen dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (def 0.9x smearing);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hGenDijetDefNinetyPtEtaBackwardArr[iCut]->Sumw2();

            // Response matrix and missing dijets histograms
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[iCut] = std::make_unique<THnSparseF>(Form("hGenDijetPtEtaCMVsRecoPtEtaCM_%d", iCut), 
                                                            Form("Gen dijet #eta_{CM} vs Reco dijet #eta_{CM} vs Gen dijet p_{T}^{ave} vs Reco dijet p_{T}^{ave};Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}"), 
                                                            nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[iCut]->Sumw2();
            // hs.hRefDijetPtEtaCMVsRecoPtEtaCMArr[iCut] = std::make_unique<THnSparseF>(Form("hRefDijetPtEtaCMVsRecoPtEtaCM_%d", iCut),
            //                                                 "Reference-matched dijet response;Ref p_{T}^{ave} (GeV);Ref #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
            //                                                 nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            // hs.hRefDijetPtEtaCMVsRecoPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMMissArr[iCut] = std::make_unique<TH2D>(Form("hGenDijetPtEtaCMMiss_%d", iCut), 
                                                                Form("Gen dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMissArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMFakeArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMFake_%d", iCut),
                                                                Form("Unmatched reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMFakeArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[iCut] = std::make_unique<THnSparseF>(
                Form("hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCM_%d", iCut),
                "Gen dijet vs reco dijet with eta-dependent default JER;Gen p_{T}^{ave} (GeV);Gen #eta_{CM};Reco p_{T}^{ave} (GeV);Reco #eta_{CM}",
                nDimResponse, responseBins, responseBinsMin, responseBinsMax);
            hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[iCut]->Sumw2();
            hs.hGenDijetPtEtaCMMissJerDefExtraArr[iCut] = std::make_unique<TH2D>(
                Form("hGenDijetPtEtaCMMissJerDefExtra_%d", iCut),
                Form("Gen dijet miss for eta-dependent default JER (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hGenDijetPtEtaCMMissJerDefExtraArr[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[iCut] = std::make_unique<TH2D>(
                Form("hRecoDijetPtEtaCMFakeJerDefExtra_%d", iCut),
                Form("Reco dijet fake for eta-dependent default JER (|#eta| < %.1f);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[iCut]->Sumw2();
            hs.hUnfoldingPairClassificationArr[iCut] = std::make_unique<TH1D>(
                Form("hUnfoldingPairClassification_%d", iCut),
                Form("Unfolding pair classification (|#eta| < %.1f);Pair category;Weighted events", etaCuts[iCut]),
                unfolding_diagnostics::nPairCategories, 0.5,
                unfolding_diagnostics::nPairCategories + 0.5);
            hs.hUnfoldingPairClassificationArr[iCut]->Sumw2();

            auto *classificationAxis = hs.hUnfoldingPairClassificationArr[iCut]->GetXaxis();
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMatchedDirect, "Matched direct");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMatchedSwapped, "Matched swapped");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kGenPassRecoFail, "Gen pass, reco fail");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kRecoPassGenFail, "Reco pass, gen fail");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kSelectedPairMismatch, "Selected-pair mismatch");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMissingOneRecoReference, "Missing one reco reference");
            classificationAxis->SetBinLabel(unfolding_diagnostics::kMissingBothRecoReferences, "Missing both reco references");
            hs.hUnfoldingPairClassificationJerDefExtraArr[iCut] = std::make_unique<TH1D>(
                Form("hUnfoldingPairClassificationJerDefExtra_%d", iCut),
                Form("Unfolding pair classification for eta-dependent default JER (|#eta| < %.1f);Pair category;Weighted events", etaCuts[iCut]),
                unfolding_diagnostics::nPairCategories, 0.5,
                unfolding_diagnostics::nPairCategories + 0.5);
            hs.hUnfoldingPairClassificationJerDefExtraArr[iCut]->Sumw2();
            auto *extraClassificationAxis = hs.hUnfoldingPairClassificationJerDefExtraArr[iCut]->GetXaxis();
            for (int category = 1; category <= unfolding_diagnostics::nPairCategories; ++category) {
                extraClassificationAxis->SetBinLabel(category, classificationAxis->GetBinLabel(category));
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)



        //
        // Ref level histograms
        //

        // Ref jet histograms
        hs.hRefInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetPtEtaLabUnflipped", 
                                                                      "Ref jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta;", 
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetPtEtaLab = std::make_unique<TH2D>("hRefInclusiveJetPtEtaLab", 
                                                              "Ref jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta;", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetPtEtaLab->Sumw2();
        hs.hRefInclusiveJetPtEtaCM = std::make_unique<TH2D>("hRefInclusiveJetPtEtaCM", 
                                                              "Ref jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM};", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetPtEtaCM->Sumw2();
        hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetTrkMaxPtEtaLabUnflipped", 
                                                                      "Ref jet track max p_{T} (lab frame, unflipped) vs jet p_{T} (TrkMax selection);p_{T}^{jet} (GeV);#eta;", 
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetTrkMaxPtEtaLab = std::make_unique<TH2D>("hRefInclusiveJetTrkMaxPtEtaLab", 
                                                                    "Ref jet track max p_{T} (lab frame) vs jet p_{T} (TrkMax selection);p_{T}^{jet} (GeV);#eta;", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetTrkMaxPtEtaLab->Sumw2();
        hs.hRefInclusiveJetTrkMaxPtEtaCM = std::make_unique<TH2D>("hRefInclusiveJetTrkMaxPtEtaCM", 
                                                                    "Ref jet track max p_{T} (CM frame) vs jet p_{T} (TrkMax selection);p_{T}^{jet} (GeV);#eta_{CM};", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetTrkMaxPtEtaCM->Sumw2();
        hs.hRefInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRefInclusiveJetPtEtaStdBins", 
                                                                "Ref jet #eta (standard bins) vs p_{T};p_{T} (GeV);#eta;",
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hs.hRefInclusiveJetPtEtaStdBins->Sumw2();
        hs.hRefInclusiveJetNoSelPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetNoSelPtEtaLabUnflipped", 
                                                                      "Ref jet #eta (lab frame, unflipped) vs p_{T} (no selection);p_{T} (GeV);#eta;", 
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetNoSelPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetNoSelPtEtaLab = std::make_unique<TH2D>("hRefInclusiveJetNoSelPtEtaLab", 
                                                              "Ref jet #eta (lab frame) vs p_{T} (no selection);p_{T} (GeV);#eta;", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetNoSelPtEtaLab->Sumw2();
        hs.hRefInclusiveJetNoSelPtEtaCM = std::make_unique<TH2D>("hRefInclusiveJetNoSelPtEtaCM", 
                                                              "Ref jet #eta (CM frame) vs p_{T} (no selection);p_{T} (GeV);#eta_{CM};", 
                                                              nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetNoSelPtEtaCM->Sumw2();

        // Ref dijet histograms

        hs.hRefDijetPtEtaLabUnflipped = std::make_unique<TH2D>("hRefDijetPtEtaLabUnflipped", 
                                                               "Ref dijet #eta (lab frame, unflipped) vs p_{T}^{ave} (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta", 
                                                               nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                               nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRefDijetPtEtaLabUnflipped->Sumw2();
        hs.hRefDijetPtEtaLab = std::make_unique<TH2D>("hRefDijetPtEtaLab", 
                                                      "Ref dijet #eta (lab frame) vs p_{T}^{ave} (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta", 
                                                      nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                      nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRefDijetPtEtaLab->Sumw2();
        hs.hRefDijetPtEtaCM = std::make_unique<TH2D>("hRefDijetPtEtaCM", 
                                                     "Ref dijet #eta (CM frame) vs p_{T}^{ave} (|#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta_{CM}", 
                                                     nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                     nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRefDijetPtEtaCM->Sumw2();

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCM_%d", iCut), 
                                                Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefDijetPtEtaCMArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaForward_%d", iCut), 
                                                Form("Ref dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRefDijetPtEtaForwardArr[iCut]->Sumw2();
            hs.hRefDijetPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaBackward_%d", iCut), 
                                                    Form("Ref dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRefDijetPtEtaBackwardArr[iCut]->Sumw2();

            if (isMc) {
                hs.hRefDijetPtEtaCMArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMJerUp_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                hs.hRefDijetPtEtaCMArrJerUp[iCut]->Sumw2();
                hs.hRefDijetPtEtaForwardArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaForwardJerUp_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaForwardArrJerUp[iCut]->Sumw2();
                hs.hRefDijetPtEtaBackwardArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaBackwardJerUp_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaBackwardArrJerUp[iCut]->Sumw2();
                hs.hRefDijetPtEtaCMArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMJerDown_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                hs.hRefDijetPtEtaCMArrJerDown[iCut]->Sumw2();
                hs.hRefDijetPtEtaForwardArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaForwardJerDown_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaForwardArrJerDown[iCut]->Sumw2();
                hs.hRefDijetPtEtaBackwardArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaBackwardJerDown_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaBackwardArrJerDown[iCut]->Sumw2();
                hs.hRefDijetPtEtaCMArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaCMJerDef_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                hs.hRefDijetPtEtaCMArrJerDef[iCut]->Sumw2();
                hs.hRefDijetPtEtaForwardArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaForwardJerDef_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaForwardArrJerDef[iCut]->Sumw2();
                hs.hRefDijetPtEtaBackwardArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaBackwardJerDef_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaBackwardArrJerDef[iCut]->Sumw2();
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    }

    //
    // Ref-selected level histograms
    //

    if (isMc) {
        // Loop over eta cuts for dijet histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefSelDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hRefSelDijetPtEtaCM_%d", iCut), 
                                                Form("Ref-selected dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRefSelDijetPtEtaCMArr[iCut]->Sumw2();
            hs.hRefSelDijetPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hRefSelDijetPtEtaForward_%d", iCut), 
                                                Form("Ref-selected dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRefSelDijetPtEtaForwardArr[iCut]->Sumw2();
            hs.hRefSelDijetPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hRefSelDijetPtEtaBackward_%d", iCut), 
                                                    Form("Ref-selected dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRefSelDijetPtEtaBackwardArr[iCut]->Sumw2();
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    }


    //
    // Reco level histograms
    //

    // Reco jet histograms
    hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnflippedUnweighted", 
                                                                    "Reco jet #eta (lab frame, unflipped, unweighted) vs p_{T};p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted->Sumw2();
    hs.hRecoInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnflipped", 
                                                                    "Reco jet #eta (lab frame, unflipped) vs p_{T};p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetPtEtaLab = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLab", 
                                                          "Reco jet #eta (lab frame) vs p_{T};p_{T} (GeV);#eta", 
                                                          nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                          nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaLab->Sumw2();
    hs.hRecoInclusiveJetPtEtaCM = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaCM", 
                                                          "Reco jet #eta (CM frame) vs p_{T};p_{T} (GeV);#eta_{CM}", 
                                                          nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                          nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetPtEtaCM->Sumw2();
    hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetTrkMaxPtEtaLabUnflipped", 
                                                                        "Reco jet #eta (lab frame, unflipped) vs p_{T} (TrkMax selection);p_{T} (GeV);#eta", 
                                                                        nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetTrkMaxPtEtaLab = std::make_unique<TH2D>("hRecoInclusiveJetTrkMaxPtEtaLab", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (TrkMax selection);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetTrkMaxPtEtaLab->Sumw2();
    hs.hRecoInclusiveJetTrkMaxPtEtaCM = std::make_unique<TH2D>("hRecoInclusiveJetTrkMaxPtEtaCM", 
                                                                    "Reco jet #eta (CM frame) vs p_{T} (TrkMax selection);p_{T} (GeV);#eta_{CM}", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetTrkMaxPtEtaCM->Sumw2();
    hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetNoSelPtEtaLabUnflipped", 
                                                                    "Reco jet #eta (lab frame, unflipped) vs p_{T} (no selection);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetNoSelPtEtaLab = std::make_unique<TH2D>("hRecoInclusiveJetNoSelPtEtaLab", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (no selection);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetNoSelPtEtaLab->Sumw2();
    hs.hRecoInclusiveJetNoSelPtEtaCM = std::make_unique<TH2D>("hRecoInclusiveJetNoSelPtEtaCM", 
                                                                    "Reco jet #eta (CM frame) vs p_{T} (no selection);p_{T} (GeV);#eta_{CM}", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetNoSelPtEtaCM->Sumw2();
    // for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
    //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
    //                                           Form("Reco jet #eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
    //                                           nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
    //     hs.hRecoInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShifted_%d", iShift), 
    //                                           Form("Reco jet #eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
    //                                           nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    //     hs.hRecoInclusiveJetEtaCMShifted[iShift]->Sumw2();
    // } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

    hs.hRecoInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaStdBins", 
                                                            "Reco jet #eta (standard bins) vs p_{T};p_{T} (GeV);#eta",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hs.hRecoInclusiveJetPtEtaStdBins->Sumw2();
    hs.hRecoInclusiveJetPtPtHat = std::make_unique<TH2D>("hRecoInclusiveJetPtPtHat", 
                                                            "Reco jet #hat{p}_{T} vs p_{T};p_{T} (GeV);#hat{p}_{T} (GeV)",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtPtHat->Sumw2();

    if (isMc) {
        hs.hRecoInclusiveJetPtEtaLabMatched = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabMatched", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (matched to ref);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetPtEtaLabMatched->Sumw2();
        hs.hRecoInclusiveJetPtEtaLabUnmatched = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnmatched", 
                                                                    "Reco jet #eta (lab frame) vs p_{T} (unmatched to ref);p_{T} (GeV);#eta", 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetPtEtaLabUnmatched->Sumw2();
        hs.hRecoInclusiveJetJESPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESPtEta", 
                                                                    "Reco jet JES (reco/ref) vs pT vs #eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESPtEta->Sumw2();

        hs.hRecoInclusiveJetJESDefNoSFPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefNoSFPtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins, no scaling);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefNoSFPtEta->Sumw2();
        hs.hRecoInclusiveJetJESDefDoublePtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefDoublePtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins, double scaling);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefDoublePtEta->Sumw2();
        hs.hRecoInclusiveJetJESDefNinetyPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefNinetyPtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins, 90% scaling);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefNinetyPtEta->Sumw2();
        hs.hRecoInclusiveJetJESDefPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDefPtEta", 
                                                                    "Reco jet JES default (reco/ref) vs pT vs #eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDefPtEta->Sumw2();
        hs.hRecoInclusiveJetJESUpPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESUpPtEta", 
                                                                    "Reco jet JES up (reco/ref) vs pT vs #eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESUpPtEta->Sumw2();
        hs.hRecoInclusiveJetJESDownPtEta = std::make_unique<TH3D>("hRecoInclusiveJetJESDownPtEta", 
                                                                        "Reco jet JES down (reco/ref) vs pT vs #eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                        nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                        nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESDownPtEta->Sumw2();
    }

    // Reco dijet histograms

    hs.hRecoDijetPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoDijetPtEtaLabUnflipped", 
                                                            "Reco dijet #eta (lab frame, unflipped) vs p_{T}^{ave} (with |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta",
                                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaLabUnflipped->Sumw2();
    hs.hRecoDijetPtEtaLab = std::make_unique<TH2D>("hRecoDijetPtEtaLab", 
                                                   "Reco dijet #eta (lab frame) vs p_{T}^{ave} (with |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta",
                                                   nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                   nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaLab->Sumw2();
    hs.hRecoDijetPtEtaCM = std::make_unique<TH2D>("hRecoDijetPtEtaCM", 
                                                  "Reco dijet #eta_{CM} (CM frame) vs p_{T}^{ave} (with |#eta_{CM}^{jet}|<1.9);p_{T}^{ave} (GeV);#eta_{CM}",
                                                  nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                  nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
    hs.hRecoDijetPtEtaCM->Sumw2();


    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        hs.hRecoDijetPtEtaCMArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCM_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
        hs.hRecoDijetPtEtaCMArr[iCut]->Sumw2();
        hs.hRecoDijetPtEtaForwardArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForward_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
        hs.hRecoDijetPtEtaForwardArr[iCut]->Sumw2();
        hs.hRecoDijetPtEtaBackwardArr[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackward_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave};p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
        hs.hRecoDijetPtEtaBackwardArr[iCut]->Sumw2();

        if (isMc) {
            hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefNoSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER no SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerDefNoSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerDefNoSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER no SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerDefNoSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerDefNoSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerDefNoSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER no SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerDefNoSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefDoubleSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default, double SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerDefDoubleSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerDefDoubleSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER default, double SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerDefDoubleSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerDefDoubleSF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerDefDoubleSF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER default, double SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerDefDoubleSF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefNinetySF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default, 90%% SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerDefNinetySF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerDefNinetySF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER default, 90%% SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerDefNinetySF[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerDefNinetySF[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerDefNinetySF_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER default, 90%% SF);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerDefNinetySF[iCut]->Sumw2();

            hs.hRecoDijetPtEtaCMArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDown[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerDown[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerDown[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDef_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDef[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerDef_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerDef[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerDef[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerDef_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER default);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerDef[iCut]->Sumw2();

            auto makeRecoJerExtraTriplet = [&](const char *variation, const char *title,
                                               std::unique_ptr<TH2D> *fullArr,
                                               std::unique_ptr<TH2D> *forwardArr,
                                               std::unique_ptr<TH2D> *backwardArr) {
                fullArr[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaCM%s_%d", variation, iCut),
                    Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (%s);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut], title),
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
                forwardArr[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaForward%s_%d", variation, iCut),
                    Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (%s);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut], title),
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                backwardArr[iCut] = std::make_unique<TH2D>(
                    Form("hRecoDijetPtEtaBackward%s_%d", variation, iCut),
                    Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (%s);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut], title),
                    nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                    nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                fullArr[iCut]->Sumw2();
                forwardArr[iCut]->Sumw2();
                backwardArr[iCut]->Sumw2();
            };
            makeRecoJerExtraTriplet("JerDownExtra", "JER down with eta-dependent scaling",
                hs.hRecoDijetPtEtaCMArrJerDownExtra, hs.hRecoDijetPtEtaForwardArrJerDownExtra,
                hs.hRecoDijetPtEtaBackwardArrJerDownExtra);
            makeRecoJerExtraTriplet("JerDefExtra", "JER default with eta-dependent scaling",
                hs.hRecoDijetPtEtaCMArrJerDefExtra, hs.hRecoDijetPtEtaForwardArrJerDefExtra,
                hs.hRecoDijetPtEtaBackwardArrJerDefExtra);
            makeRecoJerExtraTriplet("JerUpExtra", "JER up with eta-dependent scaling",
                hs.hRecoDijetPtEtaCMArrJerUpExtra, hs.hRecoDijetPtEtaForwardArrJerUpExtra,
                hs.hRecoDijetPtEtaBackwardArrJerUpExtra);
            makeRecoJerExtraTriplet("JerDefExtraUnfold", "JER default with eta-dependent scaling, unfolding",
                hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold, hs.hRecoDijetPtEtaForwardArrJerDefExtraUnfold,
                hs.hRecoDijetPtEtaBackwardArrJerDefExtraUnfold);

            hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJerDefUnfold_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JER default, unfolding);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJerDefUnfold[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJerDefUnfold_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JER default, unfolding);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJerDefUnfold[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJerDefUnfold[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJerDefUnfold_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JER default, unfolding);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJerDefUnfold[iCut]->Sumw2();
        } // if (isMc)

        if (!isMc) {
            hs.hRecoDijetPtEtaCMArrJeuUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJeuUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JEU up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJeuUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJeuUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJeuUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JEU up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJeuUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJeuUp[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJeuUp_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JEU up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJeuUp[iCut]->Sumw2();
            hs.hRecoDijetPtEtaCMArrJeuDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaCMJeuDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, |#eta| < %.1f) vs p_{T}^{ave} (JEU down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaBins, dijetEtaBins[0], dijetEtaBins[1]);
            hs.hRecoDijetPtEtaCMArrJeuDown[iCut]->Sumw2();
            hs.hRecoDijetPtEtaForwardArrJeuDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaForwardJeuDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JEU down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaForwardArrJeuDown[iCut]->Sumw2();
            hs.hRecoDijetPtEtaBackwardArrJeuDown[iCut] = std::make_unique<TH2D>(Form("hRecoDijetPtEtaBackwardJeuDown_%d", iCut), 
                                            Form("Reco dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JEU down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]), 
                                            nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                            nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
            hs.hRecoDijetPtEtaBackwardArrJeuDown[iCut]->Sumw2();
        } // if (!isMc)
    } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)


    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
} // for (int iCut{0}; iCut < nEtaCuts; ++iCut)

//________________
/// Set up the chains for the different trees
/// @param input: the input file or file list
/// @param hltTree: the HLT tree
/// @param eventTree: the event tree
/// @param skimTree: the skim tree
/// @param jetTree: the jet tree
void setupChains(const TString &input, TChain &hltTree, TChain &eventTree, TChain &skimTree, TChain &jetTree) {

    // Check input exists
    if ( input.Length() <= 0 ) {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        exit(0);
    }
    // Regular input
    else {
        // If input is a single ROOT file
        if ( input.Index(".root") > 0 ) {
            std::cout << Form( "Adding %s file to chains\n", input.Data() );
            hltTree.Add( input.Data() );
            eventTree.Add( input.Data() );
            skimTree.Add( input.Data() );
            jetTree.Add( input.Data() );
        }
        // Assuming that list of files is provided instead of a single file
        else {
            std::ifstream inputStream( input.Data() );

            if ( !inputStream ) std::cout << Form( "ERROR: Cannot open file list: %s\n", input.Data() );
            Int_t nFiles = 0;
            std::string file;
            size_t pos;
            while ( getline( inputStream, file ) ) {
                // NOTE: our external formatters may pass "file NumEvents"
                //       Take only the first part
                //cout << "DEBUG found " <<  file << endl;
                pos = file.find_first_of(" ");
                if ( pos != std::string::npos ) file.erase( pos, file.length() - pos );
                //cout << "DEBUG found [" <<  file << "]" << endl;

                // Check that file is of a correct name
                if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos &&
                     file.find("AOD") != std::string::npos ) {
                    
                    // Open file
                    TFile* ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if ( ftmp && !ftmp->IsZombie() && ftmp->GetNkeys() ) {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str() );
                        // Adding file to chains
                        hltTree.Add( file.c_str() );
                        eventTree.Add( file.c_str() );
                        skimTree.Add( file.c_str() );
                        jetTree.Add( file.c_str() );
                        ++nFiles;
                    } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp) {
                        ftmp->Close();
                    } //if (ftmp)
                } //if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos && file.find("AOD") != std::string::npos )
            } //while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
        } // else {   if file list
    } // else {   if normal input

    // Connect chains between each other
    hltTree.AddFriend(&eventTree);
    hltTree.AddFriend(&skimTree);
    hltTree.AddFriend(&jetTree);
}

//________________
/// Set up the branches for the main tree
/// @param mainTree: the main tree
/// @param isMc: true if the data is from Monte Carlo, false otherwise
void setupBranches(TChain &mainTree, const bool &isMc) {

    std::cout << "Setting up branches...";

    //
    // Disable all branches first
    //
    mainTree.SetBranchStatus("*", 0);

    //
    // Enable only the branches needed for this analysis
    //

    // Event level branches
    mainTree.SetBranchStatus("vz", 1);
    mainTree.SetBranchStatus("pthat", 1);

    // Enable skim/event filter branches
    mainTree.SetBranchStatus("pBeamScrapingFilter", 1);
    mainTree.SetBranchStatus("pPAprimaryVertexFilter", 1);
    mainTree.SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    mainTree.SetBranchStatus("phfCoincFilter", 1);
    mainTree.SetBranchStatus("pVertexFilterCutdz1p0", 1);
    mainTree.SetBranchStatus("pVertexFilterCutGplus", 1);
    mainTree.SetBranchStatus("pVertexFilterCutVtx1", 1);

    // Trigger branches
    mainTree.SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4", 1);
    mainTree.SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3", 1);
    mainTree.SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3", 1);

    // Jet level branches
    if (isMc) {
        mainTree.SetBranchStatus("ngen", 1);
        mainTree.SetBranchStatus("genpt", 1);
        mainTree.SetBranchStatus("geneta", 1);
        mainTree.SetBranchStatus("genphi", 1);
    }

    mainTree.SetBranchStatus("nref", 1);
    mainTree.SetBranchStatus("rawpt", 1);
    mainTree.SetBranchStatus("trackMax", 1);
    mainTree.SetBranchStatus("jteta", 1);
    mainTree.SetBranchStatus("jtphi", 1);
    mainTree.SetBranchStatus("jtPfNHF", 1);
    mainTree.SetBranchStatus("jtPfNEF", 1);
    mainTree.SetBranchStatus("jtPfCHF", 1);
    mainTree.SetBranchStatus("jtPfMUF", 1);
    mainTree.SetBranchStatus("jtPfCEF", 1);
    mainTree.SetBranchStatus("jtPfCHM", 1);
    mainTree.SetBranchStatus("jtPfCEM", 1);
    mainTree.SetBranchStatus("jtPfNHM", 1);
    mainTree.SetBranchStatus("jtPfNEM", 1);
    mainTree.SetBranchStatus("jtPfMUM", 1);

    if (isMc) {
        mainTree.SetBranchStatus("refpt", 1);
        mainTree.SetBranchStatus("refeta", 1);
        mainTree.SetBranchStatus("refphi", 1);
    }

    // Set branch addresses
    mainTree.SetBranchAddress("vz", &vz);
    mainTree.SetBranchAddress("pthat", &ptHat);

    mainTree.SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter);
    mainTree.SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter);
    mainTree.SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose);
    mainTree.SetBranchAddress("phfCoincFilter", &phfCoincFilter);
    mainTree.SetBranchAddress("pVertexFilterCutdz1p0", &pVertexFilterCutdz1p0);
    mainTree.SetBranchAddress("pVertexFilterCutGplus", &pVertexFilterCutGplus);
    mainTree.SetBranchAddress("pVertexFilterCutVtx1", &pVertexFilterCutVtx1);

    mainTree.SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4", &HLT_PAAK4PFJet60_Eta5p1_v4);
    mainTree.SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3", &HLT_PAAK4PFJet80_Eta5p1_v3);
    mainTree.SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3", &HLT_PAAK4PFJet100_Eta5p1_v3);

    if (isMc) {
        mainTree.SetBranchAddress("ngen",   &nGenJets);
        mainTree.SetBranchAddress("genpt",  &genJetPt);
        mainTree.SetBranchAddress("geneta", &genJetEta);
        mainTree.SetBranchAddress("genphi", &genJetPhi);
    }

    mainTree.SetBranchAddress("nref",     &nRecoJets);
    mainTree.SetBranchAddress("rawpt",    &recoJetPtRaw);
    mainTree.SetBranchAddress("trackMax", &recoJetTrackMaxPt);
    mainTree.SetBranchAddress("jteta",    &recoJetEta);
    mainTree.SetBranchAddress("jtphi",    &recoJetPhi);
    mainTree.SetBranchAddress("jtPfNHF",  &recoJetPfNHF);
    mainTree.SetBranchAddress("jtPfNEF",  &recoJetPfNEF);
    mainTree.SetBranchAddress("jtPfCHF",  &recoJetPfCHF);
    mainTree.SetBranchAddress("jtPfMUF",  &recoJetPfMUF);
    mainTree.SetBranchAddress("jtPfCEF",  &recoJetPfCEF);
    mainTree.SetBranchAddress("jtPfCHM",  &recoJetPfCHM);
    mainTree.SetBranchAddress("jtPfCEM",  &recoJetPfCEM);
    mainTree.SetBranchAddress("jtPfNHM",  &recoJetPfNHM);
    mainTree.SetBranchAddress("jtPfNEM",  &recoJetPfNEM);
    mainTree.SetBranchAddress("jtPfMUM",  &recoJetPfMUM);

    if (isMc) {
        mainTree.SetBranchAddress("refpt", &refJetPt);
        mainTree.SetBranchAddress("refeta", &refJetEta);
        mainTree.SetBranchAddress("refphi", &refJetPhi);
    }

    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
void setupInput(const TString &input, TChain &hltTree, TChain &eventTree, TChain &skimTree, TChain &jetTree, const bool &isMc) {
    setupChains(input, hltTree, eventTree, skimTree, jetTree);
    setupBranches(hltTree, isMc);
}

//________________
float etaCM(const float &etaLab, const float &etaShift, const bool &isPbGoing, const bool &isMc) {
    if (isMc) { // For embedding: Pb goes to negative, p goes to positive
        return isPbGoing ? -1.f * (etaLab + etaShift) : etaLab - etaShift;
    }
    return isPbGoing ? etaLab - etaShift : -1.f * (etaLab + etaShift);
}

//________________
float etaLab(const float &eta, const bool &isPbGoing, const bool &isMc) {
    if (isMc) { // For embedding: Pb goes to negative, p goes to positive
        return isPbGoing ? -eta : eta;
    }
    // For data: p goes to negative, Pb goes to positive in stored orientation.
    return isPbGoing ? eta : -eta;
}

//________________
/// Set the pT hat range for a given MC sample
/// @param ptHatSample: the ptHat sample
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
    std::cout << Form(": [%.1f, %.1f]", ptHatRange[0], ptHatRange[1]) << std::endl;
}

//________________
/// Check if event passes selection criteria
/// @param isPbGoing: true if Pb is going in the positive direction, false otherwise
/// @param isMc: true if the event is from Monte Carlo, false otherwise
/// @param triggerId: the trigger ID (0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100)
/// @return true if the event is good, false otherwise
bool isGoodEvent(const bool &isPbGoing, const bool &isMc, const int &triggerId) {

    bool goodTrigger{false};
    if (isMc) {
        if (ptHat < ptHatRange[0] || ptHat > ptHatRange[1]) {
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

    return ( ( std::abs(vz) <= 15. ) && 
             ( pBeamScrapingFilter == 1 ) && 
             ( pPAprimaryVertexFilter == 1 ) && 
             ( HBHENoiseFilterResultRun2Loose == 1 ) && 
             ( phfCoincFilter == 1 ) && 
             ( pVertexFilterCutdz1p0 == 1 )  && 
             goodTrigger );
}

//________________
void loadGenJets(std::vector<GenJet> &genJets) {
    // Clear the vector before filling it
    if ( !genJets.empty() ) { genJets.clear(); }

    for (int iGenJet{0}; iGenJet < nGenJets; ++iGenJet) {
        GenJet genJet{genJetPt[iGenJet], genJetEta[iGenJet], genJetPhi[iGenJet], -1.f};
        genJets.push_back(genJet);
    }

    // Sort the genJets vector in descending order of pT
    std::sort(genJets.begin(), genJets.end(), [](const GenJet &a, const GenJet &b) { return a.pt > b.pt; });
}

//________________
void loadRecoJets(std::vector<RecoJet> &recoJets, JetCorrector &jec) {
    // Clear the vector before filling it
    if ( !recoJets.empty() ) { recoJets.clear(); }
    float recoJetPtCorr{0.};

    // Loop over reco jets and add them to the vector, applying necessary corrections and selections
    for (int iRecoJet{0}; iRecoJet < nRecoJets; ++iRecoJet) {

        // Retrieve the raw pT and apply JEC to get the corrected pT
        jec.SetJetPT( recoJetPtRaw[iRecoJet] );
        jec.SetJetEta( recoJetEta[iRecoJet] );
        jec.SetJetPhi( recoJetPhi[iRecoJet] );
        recoJetPtCorr = jec.GetCorrectedPT();

        // Here we may precut reconstructed jets
        // if ( jec.GetCorrectedPT() < 10. ) { continue; } // Precut at 10 GeV to save computational time

        // Create reco jet object and add it to the vector
        RecoJet recoJet{recoJetPtRaw[iRecoJet], recoJetPtCorr, 
                        recoJetEta[iRecoJet], recoJetPhi[iRecoJet],
                        recoJetPfNHF[iRecoJet], recoJetPfNEF[iRecoJet], recoJetPfCHF[iRecoJet], 
                        recoJetPfMUF[iRecoJet], recoJetPfCEF[iRecoJet], recoJetPfCHM[iRecoJet], 
                        recoJetPfCEM[iRecoJet], recoJetPfNHM[iRecoJet], recoJetPfNEM[iRecoJet], 
                        recoJetPfMUM[iRecoJet], recoJetTrackMaxPt[iRecoJet], 
                        refJetPt[iRecoJet], refJetEta[iRecoJet], refJetPhi[iRecoJet],
                        -1., -1.,                  // JeuUp, JeuDown
                        -1.,                       // refPtSmeared
                        };

        recoJets.push_back(recoJet);
    } // for (int iRecoJet{0}; iRecoJet < nRecoJets; ++iRecoJet)

    // Sort the recoJets vector in descending order of pT
    std::sort(recoJets.begin(), recoJets.end(), [](const RecoJet &a, const RecoJet &b) { return a.recoPt > b.recoPt; });
}

//________________
void fillOverweightHistograms(std::vector<GenJet> &genJets, std::vector<RecoJet> &recoJets, const double &weight, 
                              Histograms &hs, const bool &isPythia, bool &isGenOverweight, bool &isRecoOverweight) {
    
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
void processGenJets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                    std::vector<GenJet> &genJets, Histograms &hs, const double &ptHat) {

    float genJetEtaLabFlipped{0.};
    float genJetEtaCM{0.};
    float pt{0.};
    
    // Loop over gen jets and fill histograms
    for (const auto &genJet : genJets) {
        pt = genJet.pt;
        genJetEtaLabFlipped = etaLab(genJet.eta, isPbGoing, isMc);
        genJetEtaCM = etaCM(genJet.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);

        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted->Fill(pt, genJet.eta, 1.);
        hs.hGenInclusiveJetPtEtaLabUnflipped->Fill(pt, genJet.eta, weight);
        hs.hGenInclusiveJetPtEtaLab->Fill(pt, genJetEtaLabFlipped, weight);
        hs.hGenInclusiveJetPtEtaCM->Fill(pt, genJetEtaCM, weight);

        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            genJetEtaCM = etaCM(genJet.eta, etaShift[iShift], isPbGoing, isMc);
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
void processGenDijets(const bool &isPbGoing, const bool &isMc, const double &weight, 
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
            float leadEtaLabFlipped = etaLab(leadingJet.eta, isPbGoing, isMc);
            float leadEtaCM = etaCM(leadingJet.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
            float subleadEtaLabUnflipped = subleadingJet.eta;
            float subleadEtaLabFlipped = etaLab(subleadingJet.eta, isPbGoing, isMc);
            float subleadEtaCM = etaCM(subleadingJet.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);

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
                float leadingEtaCMShifted = etaCM(leadingJet.eta, etaShift[iShift], isPbGoing, isMc);
                float subleadingEtaCMShifted = etaCM(subleadingJet.eta, etaShift[iShift], isPbGoing, isMc);
                float dijetEtaCMShifted = 0.5 * (leadingEtaCMShifted + subleadingEtaCMShifted);
                
                if (std::abs(leadingEtaCMShifted) > 1.9 || std::abs(subleadingEtaCMShifted) > 1.9) continue;

                // if (dijetPtAve < 60. || dijetPtAve > 80.) continue;
                hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Fill(dijetEtaCMShifted);
                hs.hGenDijetEtaCMShifted[iShift]->Fill(dijetEtaCMShifted, weight);
                (dijetEtaCMShifted >= 0.) ? hs.hGenDijetEtaCMForwardShifted[iShift]->Fill(dijetEtaCMShifted, weight) : 
                                            hs.hGenDijetEtaCMBackwardShifted[iShift]->Fill( std::abs(dijetEtaCMShifted), weight);
            } // for (int iShift = 0; iShift < nEtaShifts; ++iShift)
        } // if (std::abs(dphi) >= TMath::TwoPi() / 3.)
    } // if (leadingJet.pt >= 50.f && subleadingJet.pt >= 40.)

    // fillRecoDijetSystematic:
    // - `getPt`: callable taking `const RecoJet&` and returning the pT value
    //   used for ordering and dijet pT calculations (float-like).
    // - Finds the two highest-pT jets without sorting the full collection.
    // - Applies pT (lead >= 50, sublead >= 40) and |Δφ| > 2π/3 cuts, computes
    //   the dijet pT average using `getPt`, and fills the provided per-eta-cut
    //   2D histograms (`recoFull/Forward/Backward`).
    // - If optional ref arrays are provided and both jets have `refPt > 0`,
    //   the corresponding reference histograms are filled using the ref pT.
    // - Passing different accessors (e.g. `recoPt`, `recoPtJerUp`) lets this
    //   lambda fill nominal and systematic variations without duplicating code.
    std::vector<GenJet> genJetSortFallback;
    genJetSortFallback.reserve(genJets.size());
    auto fillGenDijetSystematic = [&](auto getPt,
                                      std::unique_ptr<TH2D> *genFullArr,
                                      std::unique_ptr<TH2D> *genForwardArr,
                                      std::unique_ptr<TH2D> *genBackwardArr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(genJets, getPt, genJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        // pT cuts for systematic variations (can be adjusted as needed)
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return; 
        float systDphi = deltaPhi(leadingJetSyst.phi, subleadingJetSyst.phi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaCM = etaCM(leadingJetSyst.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        float subleadEtaCM = etaCM(subleadingJetSyst.eta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;

            genFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
            if (systDijetEtaCM >= 0.) {
                genForwardArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
            }
            else {
                genBackwardArr[iCut]->Fill(systDijetPtAve, std::abs(systDijetEtaCM), weight);
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };
    
    fillGenDijetSystematic(
        [](const GenJet &jet) { return jet.pt; },
        hs.hGenDijetPtEtaCMArr,
        hs.hGenDijetPtEtaForwardArr,
        hs.hGenDijetPtEtaBackwardArr
    );

    // Smear gen jets and fill smeared-gen histograms for several extraScale factors.
    // Number of smearing iterations (runs) per scale.
    int nSmearRuns = 20; // adjust as needed

    // Scale index map: 0: 1.0f (Default), 1: 2.0f (Double), 2: 0.9f (Ninety), 3: 1.0f (Extra/eta-dependent)
    const std::array<float, 5> extraScales = {1.0f, 2.0f, 0.9f, 1.0f};

    for (int iScale = 0; iScale < static_cast<int>(extraScales.size()); ++iScale) {
        float extraScale = extraScales[iScale];
        // Choose the target histogram arrays based on the scale index
        std::unique_ptr<TH2D> *fullArr{nullptr};
        std::unique_ptr<TH2D> *forwardArr{nullptr};
        std::unique_ptr<TH2D> *backwardArr{nullptr};
        float genJetEtaLabFlipped = {-999.f};

        if (iScale == 0) {  // x1.0f (Default)
            fullArr = hs.hGenDijetDefPtEtaCMArr;
            forwardArr = hs.hGenDijetDefPtEtaForwardArr;
            backwardArr = hs.hGenDijetDefPtEtaBackwardArr;
        }
        else if (iScale == 1) {  // x2.0f (Double)
            fullArr = hs.hGenDijetDefDoublePtEtaCMArr;
            forwardArr = hs.hGenDijetDefDoublePtEtaForwardArr;
            backwardArr = hs.hGenDijetDefDoublePtEtaBackwardArr;
        }
        else if (iScale == 2) {  // x0.9f (Ninety)
            fullArr = hs.hGenDijetDefNinetyPtEtaCMArr;
            forwardArr = hs.hGenDijetDefNinetyPtEtaForwardArr;
            backwardArr = hs.hGenDijetDefNinetyPtEtaBackwardArr;
        }
        else if (iScale == 3) {  // x1.0f (Extra/eta-dependent)
            fullArr = hs.hGenDijetDefExtraPtEtaCMArr;
            forwardArr = hs.hGenDijetDefExtraPtEtaForwardArr;
            backwardArr = hs.hGenDijetDefExtraPtEtaBackwardArr;
        }

        // If any of the target histogram arrays weren't created, skip this scale
        if (!fullArr || !forwardArr || !backwardArr) continue;

        for (int iRun = 0; iRun < nSmearRuns; ++iRun) {
            // Apply smearing to each gen jet and store in smearedPt
            for (auto &jet : genJets) {
                jet.smearedPt = jet.pt;
                genJetEtaLabFlipped = etaLab(jet.eta, isPbGoing, isMc);

                // std::cout << "gen jet pt: " << jet.pt << ", eta: " << jet.eta << ", smeared pt: " << jet.smearedPt;

                if (iScale == 3) {  // x1.0f (Extra/eta-dependent) 
                    // Eta-dependent smearing: only smear jets with |eta| > 0.8, and use the fJERSmearVsEtaFunc to get 
                    // the scale factor for smearing based on eta
                    fJERSmearingHelper->smearMomentum(jet.pt, jet.eta, jet.smearedPt, 2);
                }
                else {
                    // Constant smearing: apply the same extraScale factor to all jets regardless of eta 
                    fJERSmearingHelper->smearMomentum(jet.pt, jet.eta, jet.smearedPt, 3, extraScale);
                }

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
                fullArr,
                forwardArr,
                backwardArr
            );
        } // for runs
    } // for iScale
}

//________________
void prepareUnfolding(const bool &isPbGoing, const bool &isMc, const float &weight, 
                      std::vector<GenJet> &genJets, std::vector<RecoJet> &recoJets, 
                      std::unique_ptr<THnSparseF> *responseArr,
                      std::unique_ptr<TH2D> *missArr,
                      std::unique_ptr<TH2D> *fakeArr,
                      std::unique_ptr<TH1D> *classificationArr = nullptr) {

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

    // Sort reco jets by the caller-prepared smeared pT so the selected measured
    // dijet matches the response, miss, fake, and classification outputs for
    // the requested JER variation.
    std::sort(recoJets.begin(), recoJets.end(), [](const RecoJet &a, const RecoJet &b) {
        return a.refPtSmeared > b.refPtSmeared;
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
                                   isPbGoing, isMc);
        genLeadingJetPhi = genLeadingJet.phi;

        genSubleadingJetPt = genSubleadingJet.pt;
        genSubleadingJetEtaCM = etaCM(genSubleadingJet.eta,
                                      analysis_contract::kNominalEtaShift,
                                      isPbGoing, isMc);
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

        recoLeadingJetPt = recoLeadingJet.refPtSmeared;
        recoLeadingJetEtaCM = etaCM(recoLeadingJet.recoEta,
                                    analysis_contract::kNominalEtaShift,
                                    isPbGoing, isMc);
        recoLeadingJetPhi = recoLeadingJet.recoPhi;

        recoSubleadingJetPt = recoSubleadingJet.refPtSmeared;
        recoSubleadingJetEtaCM = etaCM(recoSubleadingJet.recoEta,
                                       analysis_contract::kNominalEtaShift,
                                       isPbGoing, isMc);
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

        // Check if both gen and reco dijets pass selection
        const bool bothPairsPass = genPairPass && recoPairPass;
        // Check if the reco dijet matches the gen dijet (either direct or swapped)
        const bool matchedPair = bothPairsPass && goodMatching;
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

        // Fill the response matrix for matched pairs, and fill the appropriate 
        // histograms for gen-only or reco-only pairs
        if (matchedPair) {
            const double response[4] = {
                genDijetPtAve, genDijetEtaCM, recoDijetPtAve, recoDijetEtaCM
            };
            responseArr[iCut]->Fill(response, weight);
        }
        else if (genPairPass) {
            missArr[iCut]->Fill(
                genDijetPtAve, genDijetEtaCM, weight);
        }

        if (recoPairPass && !matchedPair) {
            fakeArr[iCut]->Fill(
                recoDijetPtAve, recoDijetEtaCM, weight);
        }

    } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
}


//________________
void processRecoJets(const bool &isPbGoing, const bool &isMc, const double &weight, 
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
        recoJetEtaLabFlipped = etaLab(recoJet.recoEta, isPbGoing, isMc);
        recoJetEtaCM = etaCM(recoJet.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        if (hasMatchingGenJet) {
            refJetEtaLabFlipped = etaLab(recoJet.refEta, isPbGoing, isMc);
            refJetEtaCM = etaCM(recoJet.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        }

        // Acceptance for reco jets without any selection applied
        if ( isGoodRecoJet(recoJet, 0) ) {
            hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.recoEta, weight);
            hs.hRecoInclusiveJetNoSelPtEtaLab->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            hs.hRecoInclusiveJetNoSelPtEtaCM->Fill(recoJet.recoPt, recoJetEtaCM, weight);
            if (hasMatchingGenJet) {
                hs.hRefInclusiveJetNoSelPtEtaLabUnflipped->Fill(recoJet.refPt, recoJet.refEta, weight);
                hs.hRefInclusiveJetNoSelPtEtaLab->Fill(recoJet.refPt, refJetEtaLabFlipped, weight);
                hs.hRefInclusiveJetNoSelPtEtaCM->Fill(recoJet.refPt, refJetEtaCM, weight);
            }
        }
        // Acceptance for reco jets with trkMax selection applied
        if ( isGoodRecoJet(recoJet, 1) ) {
            hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.recoEta, weight);
            hs.hRecoInclusiveJetTrkMaxPtEtaLab->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            hs.hRecoInclusiveJetTrkMaxPtEtaCM->Fill(recoJet.recoPt, recoJetEtaCM, weight);
            if (hasMatchingGenJet) {
                hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped->Fill(recoJet.refPt, recoJet.refEta, weight);
                hs.hRefInclusiveJetTrkMaxPtEtaLab->Fill(recoJet.refPt, refJetEtaLabFlipped, weight);
                hs.hRefInclusiveJetTrkMaxPtEtaCM->Fill(recoJet.refPt, refJetEtaCM, weight);
            }
        }
        // Acceptance for reco jets with jetId selection applied
        if ( isGoodRecoJet(recoJet, 2) ) {
            hs.hRecoInclusiveJetPtEtaLabUnflipped->Fill(recoJet.recoPt, recoJet.recoEta,  weight);
            hs.hRecoInclusiveJetPtEtaLab->Fill(recoJet.recoPt, recoJetEtaLabFlipped, weight);
            hs.hRecoInclusiveJetPtEtaCM->Fill(recoJet.recoPt, recoJetEtaCM, weight);
            if (hasMatchingGenJet) {
                hs.hRefInclusiveJetPtEtaLabUnflipped->Fill(recoJet.refPt, recoJet.refEta, weight);
                hs.hRefInclusiveJetPtEtaLab->Fill(recoJet.refPt, refJetEtaLabFlipped, weight);
                hs.hRefInclusiveJetPtEtaCM->Fill(recoJet.refPt, refJetEtaCM, weight);
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

        // std::cout << Form("Jet: eta = %.2f, phi = %.2f, rawPt = %.1f, recoPt = %.1f", 
        //                   recoJet.recoEta, recoJet.recoPhi, recoJet.recoRawPt, recoJet.recoPt) << std::endl;
        hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted->Fill(recoJet.recoPt, recoJet.recoEta, 1.);
        hs.hRecoInclusiveJetPtEtaStdBins->Fill(recoJet.recoPt, recoJet.recoEta, weight);

        // for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
        //     recoJetEtaCM = etaCM(recoJet.recoEta, etaShift[iShift], isPbGoing, isMc);
        //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Fill( recoJetEtaCM, 1.); // Unweighted
        //     hs.hRecoInclusiveJetEtaCMShifted[iShift]->Fill(recoJetEtaCM, weight);
        // } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        // Perform JES & JER smearing calculations
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
void processRecoDijets(const bool &isPbGoing, const bool &isMc, const float &weight, 
                       std::vector<RecoJet> &recoJets, std::vector<GenJet> &genJets,
                       Histograms &hs,
                       JERSmearingHelper *fJERSmearingHelper) {

    // std::cout << "processRecoDijets input parameters:" << std::endl;
    // std::cout << "  isPbGoing: " << (isPbGoing ? "true" : "false") << std::endl;
    // std::cout << "  isMc: " << (isMc ? "true" : "false") << std::endl;
    // std::cout << "  weight: " << weight << std::endl;
    // std::cout << "  recoJets.size(): " << recoJets.size() << std::endl;
    // std::cout << "  genJets.size(): " << genJets.size() << std::endl;
    // std::cout << "  hs: " << &hs << std::endl;
    // std::cout << "  fJERSmearingHelper: " << fJERSmearingHelper << std::endl;

    // Must be at least 2 jets to form a dijet system
    if (recoJets.size() < 2) return;

    std::sort(recoJets.begin(), recoJets.end(), [](const RecoJet &a, const RecoJet &b) { return a.recoPt > b.recoPt; });

    const auto& leadingJet = recoJets[0];
    const auto& subleadingJet = recoJets[1];
    if (leadingJet.recoPt >= 50.f && subleadingJet.recoPt >= 40.f) {
        float dphi = deltaPhi(leadingJet.recoPhi, subleadingJet.recoPhi);
        if (std::abs(dphi) >= TMath::TwoPi() / 3.) {
            float recoDijetPtAve = 0.5f * (leadingJet.recoPt + subleadingJet.recoPt);
            float recoLeadEtaLabUnflipped = leadingJet.recoEta;
            float recoLeadEtaLab = etaLab(leadingJet.recoEta, isPbGoing, isMc);
            float recoLeadEtaCM = etaCM(leadingJet.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
            float recoSubleadEtaLabUnflipped = subleadingJet.recoEta;
            float recoSubleadEtaLab = etaLab(subleadingJet.recoEta, isPbGoing, isMc);
            float recoSubleadEtaCM = etaCM(subleadingJet.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
            if (std::abs(recoLeadEtaCM) <= 1.9 && std::abs(recoSubleadEtaCM) <= 1.9) {
                float dijetEtaLabUnflipped = 0.5f * (recoLeadEtaLabUnflipped + recoSubleadEtaLabUnflipped);
                float dijetEtaLab = 0.5f * (recoLeadEtaLab + recoSubleadEtaLab);
                float dijetEtaCM = 0.5f * (recoLeadEtaCM + recoSubleadEtaCM);

                hs.hRecoDijetPtEtaLabUnflipped->Fill(recoDijetPtAve, dijetEtaLabUnflipped, weight);
                hs.hRecoDijetPtEtaLab->Fill(recoDijetPtAve, dijetEtaLab, weight);
                hs.hRecoDijetPtEtaCM->Fill(recoDijetPtAve, dijetEtaCM, weight);

                if (leadingJet.refPt > 0. && subleadingJet.refPt > 0.) {
                    float refLeadEtaLabUnflipped = leadingJet.refEta;
                    float refLeadEtaLab = etaLab(leadingJet.refEta, isPbGoing, isMc);
                    float refLeadEtaCM = etaCM(leadingJet.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
                    float refSubleadEtaLabUnflipped = subleadingJet.refEta;
                    float refSubleadEtaLab = etaLab(subleadingJet.refEta, isPbGoing, isMc);
                    float refSubleadEtaCM = etaCM(subleadingJet.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
                    
                    float refDijetPtAve = 0.5f * (leadingJet.refPt + subleadingJet.refPt);
                    float refDijetEtaLabUnflipped = 0.5f * (refLeadEtaLabUnflipped + refSubleadEtaLabUnflipped);
                    float refDijetEtaLab = 0.5f * (refLeadEtaLab + refSubleadEtaLab);
                    float refDijetEtaCM = 0.5f * (refLeadEtaCM + refSubleadEtaCM);

                    hs.hRefDijetPtEtaLabUnflipped->Fill(refDijetPtAve, refDijetEtaLabUnflipped, weight);
                    hs.hRefDijetPtEtaLab->Fill(refDijetPtAve, refDijetEtaLab, weight);
                    hs.hRefDijetPtEtaCM->Fill(refDijetPtAve, refDijetEtaCM, weight);
                } // if (leadingJet.refPt > 0. && subleadingJet.refPt > 0.)
            } // if (std::abs(recoLeadEtaCM) <= 1.9 && std::abs(recoSubleadEtaCM) <= 1.9)
        } // if (std::abs(dphi) >= TMath::TwoPi() / 3.)
    } // if (leadingJet.recoPt >= 50.f && subleadingJet.recoPt >= 40.f)

    // fillRecoDijetSystematic:
    // - `getPt`: callable taking `const RecoJet&` and returning the pT value
    //   used for ordering and dijet pT calculations (float-like).
    // - Finds the two highest-pT jets without sorting the full collection.
    // - Applies pT (lead >= 50, sublead >= 40) and |Δφ| > 2π/3 cuts, computes
    //   the dijet pT average using `getPt`, and fills the provided per-eta-cut
    //   2D histograms (`recoFull/Forward/Backward`).
    // - If optional ref arrays are provided and both jets have `refPt > 0`,
    //   the corresponding reference histograms are filled using the ref pT.
    // - Passing different accessors (e.g. `recoPt`, `recoPtJerUp`) lets this
    //   lambda fill nominal and systematic variations without duplicating code.
    std::vector<RecoJet> recoJetSortFallback;
    recoJetSortFallback.reserve(recoJets.size());
    auto fillRecoDijetSystematic = [&](auto getPt,
                                       std::unique_ptr<TH2D> *recoFullArr,
                                       std::unique_ptr<TH2D> *recoForwardArr,
                                       std::unique_ptr<TH2D> *recoBackwardArr,
                                       std::unique_ptr<TH2D> *refFullArr = nullptr,
                                       std::unique_ptr<TH2D> *refForwardArr = nullptr,
                                       std::unique_ptr<TH2D> *refBackwardArr = nullptr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(recoJets, getPt, recoJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        // pT cuts for systematic variations (can be adjusted as needed)
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return; 
        float systDphi = deltaPhi(leadingJetSyst.recoPhi, subleadingJetSyst.recoPhi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float systRefDijetPtAve = 0.5f * (leadingJetSyst.refPt + subleadingJetSyst.refPt);
        float leadEtaCM = etaCM(leadingJetSyst.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        float subleadEtaCM = etaCM(subleadingJetSyst.recoEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;
            recoFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
            if (refFullArr && leadingJetSyst.refPt > 0. && subleadingJetSyst.refPt > 0.) {
                refFullArr[iCut]->Fill(systRefDijetPtAve, systDijetEtaCM, weight);
            }

            if (systDijetEtaCM >= 0.) {
                recoForwardArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
                if (refForwardArr && leadingJetSyst.refPt > 0. && subleadingJetSyst.refPt > 0.) {
                    refForwardArr[iCut]->Fill(systRefDijetPtAve, systDijetEtaCM, weight);
                }
            }
            else {
                recoBackwardArr[iCut]->Fill(systDijetPtAve, std::abs(systDijetEtaCM), weight);
                if (refBackwardArr && leadingJetSyst.refPt > 0. && subleadingJetSyst.refPt > 0.) {
                    refBackwardArr[iCut]->Fill(systRefDijetPtAve, std::abs(systDijetEtaCM), weight);
                }
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };

    // Fill nominal dijet histograms
    if (isMc) {
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPt; },
                                hs.hRecoDijetPtEtaCMArr,
                                hs.hRecoDijetPtEtaForwardArr,
                                hs.hRecoDijetPtEtaBackwardArr,
                                hs.hRefDijetPtEtaCMArr,
                                hs.hRefDijetPtEtaForwardArr,
                                hs.hRefDijetPtEtaBackwardArr);
    }
    else {
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPt; },
                                hs.hRecoDijetPtEtaCMArr,
                                hs.hRecoDijetPtEtaForwardArr,
                                hs.hRecoDijetPtEtaBackwardArr);
    }

    // Fill JER systematic dijet histograms if requested
    if (isMc) {

        // Number of iterations for smearing
        const int nSmearRounds = 10;
        // 1.0, double, and 90% scale factors; nominal down/default/up; and
        // eta-dependent down/default/up variations.
        std::array<float, 9> jerScaleFactors = {1.f, 2.f, 0.9f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f};
        std::array<int, 9> jerVariations = {3, 3, 3, -1, 0, 1, 19, 20, 21};

        // Loop over smearing rounds
        for (int iIter = 0; iIter < nSmearRounds; ++iIter) {

            // Loop over the scale factors for JER smearing
            for (int iScale = 0; iScale < static_cast<int>(jerScaleFactors.size()); ++iScale) {
                
                // Loop over reconstructed jets and apply JER smearing
                for (auto &recoJet : recoJets) {
                    // Skip jets without a matching gen jet
                    if (recoJet.refPt < 0) continue;
                    // Initialize smeared pT with the reference pT
                    recoJet.refPtSmeared = recoJet.refPt;

                    // Apply the resolution smearing with the specified scale factor
                    if (iScale < 3) { // For the first three scales, apply smearing with the scale factor
                        fJERSmearingHelper->smearMomentum(recoJet.refPt, recoJet.recoEta, recoJet.refPtSmeared, 3, jerScaleFactors[iScale]);
                    } 
                    else {
                        fJERSmearingHelper->smearMomentum(recoJet.refPt, recoJet.recoEta,
                                                         recoJet.refPtSmeared, jerVariations[iScale]);
                    }
                    
                    // Skip filling histograms for the first 10 iterations to avoid overfilling
                    if (iIter > 9) continue; 

                    if (iScale == 0) { // JER def no SF (x1.0)
                        hs.hRecoInclusiveJetJESDefNoSFPtEta->Fill(recoJet.refPtSmeared / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 1) { // JER def double SF (x2.0)
                        hs.hRecoInclusiveJetJESDefDoublePtEta->Fill(recoJet.refPtSmeared / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 2) { // JER def ninety SF (x0.9)
                        hs.hRecoInclusiveJetJESDefNinetyPtEta->Fill(recoJet.refPtSmeared / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 3) { // JER down
                        hs.hRecoInclusiveJetJESDownPtEta->Fill(recoJet.refPtSmeared / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }   
                    if (iScale == 4) { // JER def
                        hs.hRecoInclusiveJetJESDefPtEta->Fill(recoJet.refPtSmeared / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                    if (iScale == 5) { // JER up
                        hs.hRecoInclusiveJetJESUpPtEta->Fill(recoJet.refPtSmeared / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
                    }
                } // for (const auto &recoJet : recoJets)

                if (iScale == 0) { // JER def no SF (x1.0)
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDefNoSF,
                                            hs.hRecoDijetPtEtaForwardArrJerDefNoSF,
                                            hs.hRecoDijetPtEtaBackwardArrJerDefNoSF);
                }
                else if (iScale == 1) { // JER def double SF (x2.0)
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDefDoubleSF,
                                            hs.hRecoDijetPtEtaForwardArrJerDefDoubleSF,
                                            hs.hRecoDijetPtEtaBackwardArrJerDefDoubleSF);
                }
                else if (iScale == 2) { // JER def ninety SF (x0.9)
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDefNinetySF,
                                            hs.hRecoDijetPtEtaForwardArrJerDefNinetySF,
                                            hs.hRecoDijetPtEtaBackwardArrJerDefNinetySF);
                }
                else if (iScale == 3) { // JER down
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDown,
                                            hs.hRecoDijetPtEtaForwardArrJerDown,
                                            hs.hRecoDijetPtEtaBackwardArrJerDown,
                                            hs.hRefDijetPtEtaCMArrJerDown,
                                            hs.hRefDijetPtEtaForwardArrJerDown,
                                            hs.hRefDijetPtEtaBackwardArrJerDown);
                }
                else if (iScale == 4) { // JER def
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDef,
                                            hs.hRecoDijetPtEtaForwardArrJerDef,
                                            hs.hRecoDijetPtEtaBackwardArrJerDef,
                                            hs.hRefDijetPtEtaCMArrJerDef,
                                            hs.hRefDijetPtEtaForwardArrJerDef,
                                            hs.hRefDijetPtEtaBackwardArrJerDef);

                    // For the first iteration, fill the nominal histograms with the default 
                    // JER smearing for unfolding purposes
                    if (iIter == 0) {
                        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                                hs.hRecoDijetPtEtaCMArrJerDefUnfold,
                                                hs.hRecoDijetPtEtaForwardArrJerDefUnfold,
                                                hs.hRecoDijetPtEtaBackwardArrJerDefUnfold);
                        prepareUnfolding(isPbGoing, isMc, weight, genJets, recoJets,
                                         hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr,
                                         hs.hGenDijetPtEtaCMMissArr,
                                         hs.hRecoDijetPtEtaCMFakeArr,
                                         hs.hUnfoldingPairClassificationArr);
                    }
                }
                else if (iScale == 5) { // JER up
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerUp,
                                            hs.hRecoDijetPtEtaForwardArrJerUp,
                                            hs.hRecoDijetPtEtaBackwardArrJerUp,
                                            hs.hRefDijetPtEtaCMArrJerUp,
                                            hs.hRefDijetPtEtaForwardArrJerUp,
                                            hs.hRefDijetPtEtaBackwardArrJerUp);
                }
                else if (iScale == 6) { // JER down with eta-dependent scaling
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDownExtra,
                                            hs.hRecoDijetPtEtaForwardArrJerDownExtra,
                                            hs.hRecoDijetPtEtaBackwardArrJerDownExtra);
                }
                else if (iScale == 7) { // JER default with eta-dependent scaling
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerDefExtra,
                                            hs.hRecoDijetPtEtaForwardArrJerDefExtra,
                                            hs.hRecoDijetPtEtaBackwardArrJerDefExtra);
                    if (iIter == 0) {
                        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                                hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold,
                                                hs.hRecoDijetPtEtaForwardArrJerDefExtraUnfold,
                                                hs.hRecoDijetPtEtaBackwardArrJerDefExtraUnfold);
                        prepareUnfolding(isPbGoing, isMc, weight, genJets, recoJets,
                                         hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr,
                                         hs.hGenDijetPtEtaCMMissJerDefExtraArr,
                                         hs.hRecoDijetPtEtaCMFakeJerDefExtraArr,
                                         hs.hUnfoldingPairClassificationJerDefExtraArr);
                    }
                }
                else if (iScale == 8) { // JER up with eta-dependent scaling
                    fillRecoDijetSystematic([](const RecoJet &jet) { return jet.refPtSmeared; },
                                            hs.hRecoDijetPtEtaCMArrJerUpExtra,
                                            hs.hRecoDijetPtEtaForwardArrJerUpExtra,
                                            hs.hRecoDijetPtEtaBackwardArrJerUpExtra);
                }
            } // for (int iScale = 0; iScale < static_cast<int>(jerScaleFactors.size()); ++iScale)
        } // for (int iIter = 0; iIter < nSmearRounds; ++iIter)

    } // if (isMc)
}

//________________
void processRefDijets(const bool &isPbGoing, const bool &isMc, const float &weight, 
                       std::vector<RecoJet> &recoJets, Histograms &hs) {

    // std::cout << "processRefDijets input parameters:" << std::endl;
    // std::cout << "  isPbGoing: " << (isPbGoing ? "true" : "false") << std::endl;
    // std::cout << "  isMc: " << (isMc ? "true" : "false") << std::endl;
    // std::cout << "  weight: " << weight << std::endl;
    // std::cout << "  recoJets.size(): " << recoJets.size() << std::endl;

    // Must be at least 2 jets to form a dijet system
    if (recoJets.size() < 2) return;

    // fillRecoDijetSystematic:
    // - `getPt`: callable taking `const RecoJet&` and returning the pT value
    //   used for ordering and dijet pT calculations (float-like).
    // - Finds the two highest-pT jets without sorting the full collection.
    // - Applies pT (lead >= 50, sublead >= 40) and |Δφ| > 2π/3 cuts, computes
    //   the dijet pT average using `getPt`, and fills the provided per-eta-cut
    //   2D histograms (`recoFull/Forward/Backward`).
    // - If optional ref arrays are provided and both jets have `refPt > 0`,
    //   the corresponding reference histograms are filled using the ref pT.
    // - Passing different accessors (e.g. `recoPt`, `recoPtJerUp`) lets this
    //   lambda fill nominal and systematic variations without duplicating code.
    std::vector<RecoJet> refJetSortFallback;
    refJetSortFallback.reserve(recoJets.size());
    auto fillRefDijetSystematic = [&](auto getPt,
                                       std::unique_ptr<TH2D> *refFullArr,
                                       std::unique_ptr<TH2D> *refForwardArr,
                                       std::unique_ptr<TH2D> *refBackwardArr) {
        const auto [leadingJetPtr, subleadingJetPtr] = selectLeadingTwo(recoJets, getPt, refJetSortFallback);
        const auto &leadingJetSyst = *leadingJetPtr;
        const auto &subleadingJetSyst = *subleadingJetPtr;
        // pT cuts for systematic variations (can be adjusted as needed)
        if (getPt(leadingJetSyst) < 50.f || getPt(subleadingJetSyst) < 40.f) return; 
        float systDphi = deltaPhi(leadingJetSyst.refPhi, subleadingJetSyst.refPhi);
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float leadEtaCM = etaCM(leadingJetSyst.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        float subleadEtaCM = etaCM(subleadingJetSyst.refEta, analysis_contract::kNominalEtaShift, isPbGoing, isMc);
        float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

        // Loop over eta cuts and fill the corresponding histograms
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;
            refFullArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
            if (systDijetEtaCM >= 0.) {
                refForwardArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
            }
            else {
                refBackwardArr[iCut]->Fill(systDijetPtAve, std::abs(systDijetEtaCM), weight);
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };

    // Fill nominal dijet histograms
    fillRefDijetSystematic([](const RecoJet &jet) { return jet.refPt; },
                           hs.hRefSelDijetPtEtaCMArr,
                           hs.hRefSelDijetPtEtaForwardArr,
                           hs.hRefSelDijetPtEtaBackwardArr);
}

//________________
void writeOutput(TString &oFileName, Histograms &hs, const bool &isMc) {

    std::cout << "Writing output to file: " << oFileName.Data();
    // Output file
    int compressionSetting = 208; // LZMA compression
    // TFile *fOut = TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting);
    auto fOut = std::unique_ptr<TFile>( TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting) );

    //
    // Event level histograms
    //
    hs.hVzRaw->Write();
    hs.hPtHatUnweighted->Write();
    hs.hPtHat->Write();
    hs.hVzUnweighted->Write();
    hs.hVz->Write();

    //
    // Gen level histograms
    //
    if (isMc) {
        // Event histograms
        hs.hGenDijetPtAveOverPtHatVsPtHat->Write();
        hs.hGenDijetPtAveOverPtHatVsPtHatPass->Write();
        hs.hGenLeadJetPtOverPtHatVsPtHat->Write();
        hs.hGenLeadJetPtOverPtHatVsPtHatPass->Write();
        hs.hRecoDijetPtAveOverPtHatVsPtHat->Write();
        hs.hRecoDijetPtAveOverPtHatVsPtHatPass->Write();
        hs.hRecoLeadJetPtOverPtHatVsPtHat->Write();
        hs.hRecoLeadJetPtOverPtHatVsPtHatPass->Write();

        // Gen jets
        hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted->Write();
        hs.hGenInclusiveJetPtEtaLabUnflipped->Write();
        hs.hGenInclusiveJetPtEtaLab->Write();
        hs.hGenInclusiveJetPtEtaCM->Write();
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Write();
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Write();
        }
        hs.hGenInclusiveJetPtEtaStdBins->Write();
        hs.hGenInclusiveJetPtPtHat->Write();

        // Gen dijets
        hs.hGenDijetPtAve->Write();
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Write();
            hs.hGenDijetEtaCMShifted[iShift]->Write();
            hs.hGenDijetEtaCMForwardShifted[iShift]->Write();
            hs.hGenDijetEtaCMBackwardShifted[iShift]->Write();
        }

        hs.hGenDijetPtEtaLabUnflipped->Write();
        hs.hGenDijetPtEtaLab->Write();
        hs.hGenDijetPtEtaCM->Write();

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hGenDijetPtEtaCMArr[iCut]->Write();
            hs.hGenDijetPtEtaForwardArr[iCut]->Write();
            hs.hGenDijetPtEtaBackwardArr[iCut]->Write();

            // Smearing studies
            hs.hGenDijetDefPtEtaCMArr[iCut]->Write();
            hs.hGenDijetDefPtEtaForwardArr[iCut]->Write();
            hs.hGenDijetDefPtEtaBackwardArr[iCut]->Write();
            hs.hGenDijetDefExtraPtEtaCMArr[iCut]->Write();
            hs.hGenDijetDefExtraPtEtaForwardArr[iCut]->Write();
            hs.hGenDijetDefExtraPtEtaBackwardArr[iCut]->Write();
            hs.hGenDijetDefDoublePtEtaCMArr[iCut]->Write();
            hs.hGenDijetDefDoublePtEtaForwardArr[iCut]->Write();
            hs.hGenDijetDefDoublePtEtaBackwardArr[iCut]->Write();
            hs.hGenDijetDefNinetyPtEtaCMArr[iCut]->Write();
            hs.hGenDijetDefNinetyPtEtaForwardArr[iCut]->Write();
            hs.hGenDijetDefNinetyPtEtaBackwardArr[iCut]->Write();

            // Response matrix 
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[iCut]->Write();
            // hs.hRefDijetPtEtaCMVsRecoPtEtaCMArr[iCut]->Write();
            hs.hGenDijetPtEtaCMMissArr[iCut]->Write();
            hs.hRecoDijetPtEtaCMFakeArr[iCut]->Write();
            hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[iCut]->Write();
            hs.hGenDijetPtEtaCMMissJerDefExtraArr[iCut]->Write();
            hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[iCut]->Write();
            hs.hUnfoldingPairClassificationArr[iCut]->Write();
            hs.hUnfoldingPairClassificationJerDefExtraArr[iCut]->Write();
        }

        hs.hGenInclusiveJetJESDefPtEta->Write();
        hs.hGenInclusiveJetJESDefExtraPtEta->Write();
        hs.hGenInclusiveJetJESDefDoublePtEta->Write();
        hs.hGenInclusiveJetJESDefNinetyPtEta->Write();

        hs.hGenInclusiveJetDefPtEtaFlipped->Write();
        hs.hGenInclusiveJetDefExtraPtEtaFlipped->Write();
        hs.hGenInclusiveJetDefDoublePtEtaFlipped->Write();
        hs.hGenInclusiveJetDefNinetyPtEtaFlipped->Write();

        //
        // Ref level histograms
        //
        
        // Ref jets
        hs.hRefInclusiveJetNoSelPtEtaLabUnflipped->Write();
        hs.hRefInclusiveJetNoSelPtEtaLab->Write();
        hs.hRefInclusiveJetNoSelPtEtaCM->Write();
        hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped->Write();
        hs.hRefInclusiveJetTrkMaxPtEtaLab->Write();
        hs.hRefInclusiveJetTrkMaxPtEtaCM->Write();
        hs.hRefInclusiveJetPtEtaLabUnflipped->Write();
        hs.hRefInclusiveJetPtEtaLab->Write();
        hs.hRefInclusiveJetPtEtaCM->Write();

        hs.hRefInclusiveJetPtEtaStdBins->Write();

        // Ref dijets

        hs.hRefDijetPtEtaLabUnflipped->Write();
        hs.hRefDijetPtEtaLab->Write();
        hs.hRefDijetPtEtaCM->Write();

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefDijetPtEtaCMArr[iCut]->Write();
            hs.hRefDijetPtEtaForwardArr[iCut]->Write();
            hs.hRefDijetPtEtaBackwardArr[iCut]->Write();
        }

        //
        // Ref-selected level histograms
        //
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefSelDijetPtEtaCMArr[iCut]->Write();
            hs.hRefSelDijetPtEtaForwardArr[iCut]->Write();
            hs.hRefSelDijetPtEtaBackwardArr[iCut]->Write();
        }

    } // if (isMc)

    //
    // Reco level histograms
    //
    
    // Reco jets
    hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetNoSelPtEtaLab->Write();
    hs.hRecoInclusiveJetNoSelPtEtaCM->Write();
    hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetTrkMaxPtEtaLab->Write();
    hs.hRecoInclusiveJetTrkMaxPtEtaCM->Write();
    hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted->Write();
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetPtEtaLab->Write();
    hs.hRecoInclusiveJetPtEtaCM->Write();

    // for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
    //     hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Write();
    //     hs.hRecoInclusiveJetEtaCMShifted[iShift]->Write();
    // }
    hs.hRecoInclusiveJetPtEtaStdBins->Write();
    if (isMc) {
        hs.hRecoInclusiveJetJESPtEta->Write();
        hs.hRecoInclusiveJetPtEtaLabMatched->Write();
        hs.hRecoInclusiveJetPtEtaLabUnmatched->Write();
        hs.hRecoInclusiveJetJESDefNoSFPtEta->Write();
        hs.hRecoInclusiveJetJESDefDoublePtEta->Write();
        hs.hRecoInclusiveJetJESDefNinetyPtEta->Write();
        hs.hRecoInclusiveJetJESDefPtEta->Write();
        hs.hRecoInclusiveJetJESUpPtEta->Write();
        hs.hRecoInclusiveJetJESDownPtEta->Write();
        hs.hRecoInclusiveJetPtPtHat->Write();
    }
    
    // Reco dijets
    hs.hRecoDijetPtEtaLabUnflipped->Write();
    hs.hRecoDijetPtEtaLab->Write();
    hs.hRecoDijetPtEtaCM->Write();

    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        hs.hRecoDijetPtEtaCMArr[iCut]->Write();
        hs.hRecoDijetPtEtaForwardArr[iCut]->Write();
        hs.hRecoDijetPtEtaBackwardArr[iCut]->Write();
        if (isMc) {
            // Reco
            if (hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut]) hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDefNoSF[iCut]) hs.hRecoDijetPtEtaForwardArrJerDefNoSF[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDefNoSF[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDefNoSF[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut]) hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDefDoubleSF[iCut]) hs.hRecoDijetPtEtaForwardArrJerDefDoubleSF[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDefDoubleSF[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDefDoubleSF[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut]) hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDefNinetySF[iCut]) hs.hRecoDijetPtEtaForwardArrJerDefNinetySF[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDefNinetySF[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDefNinetySF[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerUp[iCut]) hs.hRecoDijetPtEtaCMArrJerUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerUp[iCut]) hs.hRecoDijetPtEtaForwardArrJerUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerUp[iCut]) hs.hRecoDijetPtEtaBackwardArrJerUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDown[iCut]) hs.hRecoDijetPtEtaCMArrJerDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDown[iCut]) hs.hRecoDijetPtEtaForwardArrJerDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDown[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDef[iCut]) hs.hRecoDijetPtEtaCMArrJerDef[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDef[iCut]) hs.hRecoDijetPtEtaForwardArrJerDef[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDef[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDef[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDownExtra[iCut]) hs.hRecoDijetPtEtaCMArrJerDownExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDownExtra[iCut]) hs.hRecoDijetPtEtaForwardArrJerDownExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDownExtra[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDownExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDefExtra[iCut]) hs.hRecoDijetPtEtaCMArrJerDefExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDefExtra[iCut]) hs.hRecoDijetPtEtaForwardArrJerDefExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDefExtra[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDefExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerUpExtra[iCut]) hs.hRecoDijetPtEtaCMArrJerUpExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerUpExtra[iCut]) hs.hRecoDijetPtEtaForwardArrJerUpExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerUpExtra[iCut]) hs.hRecoDijetPtEtaBackwardArrJerUpExtra[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut]) hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDefUnfold[iCut]) hs.hRecoDijetPtEtaForwardArrJerDefUnfold[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDefUnfold[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDefUnfold[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold[iCut]) hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDefExtraUnfold[iCut]) hs.hRecoDijetPtEtaForwardArrJerDefExtraUnfold[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDefExtraUnfold[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDefExtraUnfold[iCut]->Write();
            // Ref
            if (hs.hRefDijetPtEtaCMArrJerUp[iCut]) hs.hRefDijetPtEtaCMArrJerUp[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJerUp[iCut]) hs.hRefDijetPtEtaForwardArrJerUp[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJerUp[iCut]) hs.hRefDijetPtEtaBackwardArrJerUp[iCut]->Write();
            if (hs.hRefDijetPtEtaCMArrJerDown[iCut]) hs.hRefDijetPtEtaCMArrJerDown[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJerDown[iCut]) hs.hRefDijetPtEtaForwardArrJerDown[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJerDown[iCut]) hs.hRefDijetPtEtaBackwardArrJerDown[iCut]->Write();
            if (hs.hRefDijetPtEtaCMArrJerDef[iCut]) hs.hRefDijetPtEtaCMArrJerDef[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJerDef[iCut]) hs.hRefDijetPtEtaForwardArrJerDef[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJerDef[iCut]) hs.hRefDijetPtEtaBackwardArrJerDef[iCut]->Write();
        }

        if (!isMc) {
            if (hs.hRecoDijetPtEtaCMArrJeuUp[iCut]) hs.hRecoDijetPtEtaCMArrJeuUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJeuUp[iCut]) hs.hRecoDijetPtEtaForwardArrJeuUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJeuUp[iCut]) hs.hRecoDijetPtEtaBackwardArrJeuUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaCMArrJeuDown[iCut]) hs.hRecoDijetPtEtaCMArrJeuDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJeuDown[iCut]) hs.hRecoDijetPtEtaForwardArrJeuDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJeuDown[iCut]) hs.hRecoDijetPtEtaBackwardArrJeuDown[iCut]->Write();
        }
    }

    fOut->Close();
    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
void processEvents(const bool &isPbGoing, const bool &isMc, const bool &isPythia, TChain &mainTree, 
                   JetCorrector &jec, JetUncertainty &jeu,
                   Histograms &hs,
                   const int &triggerId, const int &jetSelectionMethod) {

    std::cout << "Start event processing..." << std::endl;

    // Vz weight function for MC, derived from data/MC ratio of vz distribution in minimum bias events
    auto fVzWeight = std::make_unique<TF1>("fVzWeight", "pol8", -15.1, 15.1);
    fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);

    // JER smearing helper for calculating the smearing factor based on eta and pT
    auto fJerSmearingHelper = std::make_unique<JERSmearingHelper>();

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

    int nEventsProcessed{0};
    int nGoodEvents{0};
    bool isGenOverweight{false};
    bool isRecoOverweight{false};
    double weight{1.};

    // Vectors to hold gen and reco jets collections in the current event
    std::vector<GenJet> genJets{};
    std::vector<RecoJet> recoJets{};

    int printEvery{100000};
    Long64_t nextPrintAt{printEvery};
    auto startTime = std::chrono::steady_clock::now();

    Long64_t nEntries = mainTree.GetEntries();
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {

        mainTree.GetEntry(iEntry);
        nEventsProcessed++;
        Long64_t entriesPassed = iEntry + 1;
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
        // std::cout << Form("Event: %d ptHat = %.1f GeV, vz = %.2f cm", iEntry, ptHat, vz) << std::endl;
        // std::cout << Form("%.1f <= ptHat <= %.1f ? %s", ptHatRange[0], ptHatRange[1], 
        //                  (ptHat >= ptHatRange[0] && ptHat <= ptHatRange[1]) ? "true" : "false") << std::endl;
        // std::cout << Form("|vz| <= 15 cm ? %s", (std::abs(vz) <= 15.) ? "true" : "false") << std::endl;

        // Check the event is satisfies basic event selection
        if (!isGoodEvent(isPbGoing, isMc, triggerId)) continue;

        if ( isMc ) {
            weight = eventWeight(ptHat, vz, *fVzWeight, nEntries);
        }

        // Load reconstructed jets
        loadRecoJets( recoJets, jec );

        if (isMc) {
            // Load generated jets
            loadGenJets( genJets );
            // Fill histograms for overweight events and skip them
            fillOverweightHistograms(genJets, recoJets, weight, hs, isPythia, isGenOverweight, isRecoOverweight);
            if (isGenOverweight) continue;
            if (isRecoOverweight) continue;
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
            processGenJets(isPbGoing, isMc, weight, genJets, hs, ptHat);
            processGenDijets(isPbGoing, isMc, weight, genJets, hs, fJerSmearingHelper.get());
            processRefDijets(isPbGoing, isMc, weight, recoJets, hs);
        }

        //
        // Reco level processing
        //
        processRecoJets(isPbGoing, isMc, weight, recoJets, hs, jeu, jetSelectionMethod, ptHat);
        processRecoDijets(isPbGoing, isMc, weight, recoJets, genJets, hs, 
                          fJerSmearingHelper.get());
    } // for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry)

    std::cout << Form("Processed %d events, of which %d passed selection.", nEventsProcessed, nGoodEvents) << std::endl;
}

//________________
/**
    List of input parameters for the analysis
    
    input: input file name (or file list)
    output: output file name
    mcType (default 2): 0 - data, 1 - embedding, 2 - pythia
    isPbGoingDir (default 1): 0 - p-going, 1 - Pb-going
    ptHatSample (default 30, for MC only): 15, 30, 50, 80, 120, 170, 220, 280, 370, 460, 540
    triggerId (default 0): 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100
    recoJetSelMethod (default 2): 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetIdTightLeptVeto, 3 - jetIdLoose
*/
#if defined(__CINT__) || defined(__CLING__)
void processForestSimple( const char* input = "", const char* output = "", 
                          int mcType = 2, int isPbGoingDir = 1, int ptHatSample = 50,
                          int triggerId = 0, int recoJetSelMethod = 2 ) {

    const char *path2auxFiles = "./aux_files/pPb_8160";
#else
int main(int argc, char* argv[]) {

    const char* input = "";
    const char* output = "";
    int mcType = {2};
    int isPbGoingDir = {1};
    int ptHatSample = {50};
    int triggerId = {0};
    int recoJetSelMethod = {2};

    if (argc != 8) {
        std::cerr << RED << "Incorrect number of input parameters. Terminating." << RESET << std::endl;
        usage();
    }

    input = argv[1];
    output = argv[2];
    mcType = std::atoi(argv[3]); 
    isPbGoingDir = std::atoi(argv[4]);
    ptHatSample = std::atoi(argv[5]);
    triggerId = std::atoi(argv[6]);
    recoJetSelMethod = std::atoi(argv[7]);

    const char *path2auxFiles = "./aux_files/pPb_8160";
#endif

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

    bool isMc = (mcType != 0);
    bool isPythia = ((mcType != 0) && (mcType == 2));
    bool isEmbedding = (isMc == 1);
    bool isPbGoing = (isPbGoingDir == 1);

    // Define the standard eta bins for L2L3 corrections
    TString generator = isPythia ? "pythia" : "embedding";
    TString tag = isPythia ? "unembedded" : "embedded";
    TString directionStr = isPbGoing ? "Pbgoing" : "pgoing";

    std::vector< std::string > fJECFileNames;
    std::string fJEUFileName;

    if (isMc) {
        if (isPbGoingDir) {
            if (isEmbedding) {
                // PYTHIA+EPOS
                fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_Pbgoing_Embedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
            }
            else {
                // PYTHIA
                fJECFileNames.push_back( Form("%s/JEC/Autumn16_HI_pPb_Pbgoing_Unembedded_MC_L2Relative_AK4PF.txt", path2auxFiles) );
            }
        }
        else {
            if (isEmbedding) {
                // PYTHIA+EPOS
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
        
        fJEUFileName = Form("%s/JEC/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF.txt", path2auxFiles);
    } // else


    // Set pT hat range for MC samples
    if (isMc) setPtHatRange(ptHatSample);

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
    std::cout << "  JEC files:" << std::endl;
    for (const auto &fn : fJECFileNames) std::cout << "    " << GREEN << fn.c_str() << RESET << std::endl;

    // TFile *f = TFile::Open( fname );
    auto mainTree = std::make_unique<TChain>("hltanalysis/HltTree");
    auto eventTree = std::make_unique<TChain>("hiEvtAnalyzer/HiTree");
    auto skimTree = std::make_unique<TChain>("skimanalysis/HltTree");
    auto jetTree = std::make_unique<TChain>("ak4PFJetAnalyzer/t");
    setupInput(fname, *mainTree, *eventTree, *skimTree, *jetTree, isMc);

    ///////////////////////
    // Create histograms //
    ///////////////////////
    Histograms hs;
    createHistograms(hs, isMc);

    //////////////////////////////////
    // Setup Jet Energy Corrections //
    //////////////////////////////////
    auto jec = std::make_unique<JetCorrector>(fJECFileNames);
    auto jeu = std::make_unique<JetUncertainty>(fJEUFileName);

    //////////////////////////////
    //        Processing        //
    //////////////////////////////

    std::cout << Form("%s%s%s %s%s%s direction, ptHat%d sample, number of entries:%s %lld %s", 
                       GREEN, (mcType ? ( (isPythia) ? "Pythia" : "Embedding") : "Data" ), RESET,
                       GREEN, (isPbGoing ? "Pb-going" : "p-going"), RESET,
                       ptHatSample, 
                       GREEN, mainTree->GetEntries(), RESET) << std::endl;

    // Process the events
    processEvents(isPbGoing, isMc, isPythia, *mainTree, *jec, *jeu, hs,
                  triggerId, recoJetSelMethod);

    // Write output
    writeOutput(oFileName, hs, isMc);
}
