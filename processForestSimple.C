#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"

#include "JetCorrector.h"
#include "JetUncertainty.h"

#include <iostream>
#include <array>
#include <limits>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
#include <random>

const char* RED    = "\033[1;31m";
const char* GREEN  = "\033[1;32m";
const char* RESET  = "\033[0m";

// Eta shifts for pPb 8.16 TeV collisions, used for etaCM calculation
const int nEtaShifts = 13;
static constexpr std::array<float, nEtaShifts> etaShift{0.463, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.470, 0.475, 0.480, 0.485, 0.490, 0.495 };
const int nEtaCuts = 7;
static constexpr std::array<float, nEtaCuts> etaCuts{1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.5};

// JER correction factors
std::vector<float> fJerEtaLow{-5.191, -3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139};
std::vector<float> fJerEtaHi{-3.139, -2.964, -2.853, -2.500, -2.322, -2.043, -1.930, -1.740, -1.305, -1.131, -0.783, -0.522, 0.000, 0.522, 0.783, 1.131, 1.305, 1.740, 1.930, 2.043, 2.322, 2.500, 2.853, 2.964, 3.139, 5.191};
std::vector<float> fJerDef{1.1922, 1.1869, 1.7788, 1.3418, 1.2963, 1.1512, 1.1426, 1.1000, 1.1278, 1.1609, 1.1464, 1.1948, 1.15958, 1.15958, 1.1948, 1.1464, 1.1609, 1.1278, 1.1000, 1.1426, 1.1512, 1.2963, 1.3418, 1.7788, 1.1869, 1.1922};
std::vector<float> fJerLow{1.0434, 1.0626, 1.578, 1.1327, 1.0592, 1.0372, 1.0212, 0.9921, 1.0292, 1.0584, 1.0832, 1.1296, 1.095, 1.095, 1.1296, 1.0832, 1.0584, 1.0292, 0.9921, 1.0212, 1.0372, 1.0592, 1.1327, 1.578, 1.0626, 1.0434};
std::vector<float> fJerHi{1.341, 1.3112, 1.9796, 1.5509, 1.5334, 1.2652, 1.264, 1.2079, 1.2264, 1.2634, 1.2096, 1.26, 1.224, 1.224, 1.26, 1.2096, 1.2634, 1.2264, 1.2079, 1.264, 1.2652, 1.5334, 1.5509, 1.9796, 1.3112, 1.341};

std::vector<float> etaBinEdges = {-3.60f, -3.46f, -3.31f, -3.17f, -3.02f, -2.88f, -2.74f, -2.59f, -2.45f, -2.30f, -2.16f, -2.02f, -1.87f, -1.73f, -1.58f, -1.44f, -1.30f, -1.15f, -1.01f, -0.86f, -0.72f, -0.58f, -0.43f, -0.29f, -0.14f, 0.00f, 0.14f, 0.29f, 0.43f, 0.58f, 0.72f, 0.86f, 1.01f, 1.15f, 1.30f, 1.44f, 1.58f, 1.73f, 1.87f, 2.02f, 2.16f, 2.30f, 2.45f, 2.59f, 2.74f, 2.88f, 3.02f, 3.17f, 3.31f, 3.46f, 3.60f};
std::vector<float> jerValues_pt_30_50;
std::vector<float> jerValues_pt_50_80;
std::vector<float> jerValues_pt_80_120;
std::vector<float> jerValues_pt_120_180;
std::vector<float> jerValues_pt_180_250;
std::vector<float> jerValues_pt_250_540;
std::vector<float> jerValues_pt_540_1000;

//________________
void setValuesForExtraJerSmearing(const bool &isPythia, const bool &isPbGoing) {
    if (isPythia) {
        if (isPbGoing) {   
            jerValues_pt_30_50   = {0.0000f, 0.0000f, 13.7345f, 1.3745f, 1.2328f, 1.1388f, 1.1334f, 1.1234f, 1.0900f, 1.0901f, 1.1103f, 1.1412f, 1.1773f, 1.2527f, 1.3295f, 1.2763f, 1.1989f, 1.1162f, 1.0561f, 1.0208f, 1.0130f, 0.9996f, 0.9934f, 1.0115f, 1.0008f, 0.9990f, 0.9960f, 1.0007f, 1.0033f, 0.9828f, 1.0066f, 1.0633f, 1.1404f, 1.2356f, 1.3077f, 1.3434f, 1.2351f, 1.2070f, 1.0926f, 1.0523f, 1.0200f, 1.0320f, 1.0592f, 1.0659f, 1.1426f, 1.2017f, 1.2244f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_50_80   = {0.0000f, 0.0000f, 0.0000f, 1.3350f, 1.2188f, 1.0559f, 0.9962f, 0.9962f, 0.9993f, 1.0202f, 1.0286f, 1.0638f, 1.0926f, 1.1791f, 1.2838f, 1.2706f, 1.1911f, 1.1229f, 1.0478f, 1.0224f, 1.0134f, 1.0070f, 0.9916f, 0.9932f, 1.0059f, 1.0171f, 0.9920f, 1.0057f, 0.9918f, 0.9824f, 1.0173f, 1.0547f, 1.1322f, 1.2136f, 1.2657f, 1.2993f, 1.1727f, 1.1202f, 1.0125f, 0.9830f, 0.9358f, 0.9371f, 0.9533f, 0.9597f, 1.0290f, 1.2117f, 1.6779f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_80_120  = {0.0000f, 0.0000f, 0.0000f, 1.1694f, 1.1969f, 1.0171f, 0.8922f, 0.9015f, 0.9247f, 0.9322f, 0.9726f, 0.9797f, 1.0343f, 1.1071f, 1.2688f, 1.2548f, 1.1693f, 1.1005f, 1.0609f, 1.0153f, 0.9981f, 1.0037f, 1.0055f, 1.0030f, 1.0084f, 1.0032f, 0.9974f, 0.9918f, 0.9901f, 0.9987f, 1.0044f, 1.0628f, 1.1322f, 1.1729f, 1.2624f, 1.2621f, 1.1252f, 1.0541f, 0.9627f, 0.9253f, 0.8844f, 0.8813f, 0.8627f, 0.8680f, 1.0044f, 1.2168f, 1.5782f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_120_180 = {0.0000f, 0.0000f, 0.0000f, 0.8883f, 1.3277f, 1.1296f, 0.9198f, 0.9107f, 0.9106f, 0.9154f, 0.9467f, 0.9408f, 0.9666f, 1.0231f, 1.2140f, 1.2396f, 1.1744f, 1.0822f, 1.0538f, 1.0345f, 0.9974f, 1.0006f, 1.0064f, 1.0111f, 1.0187f, 0.9923f, 0.9940f, 0.9944f, 0.9956f, 0.9894f, 1.0156f, 1.0568f, 1.1039f, 1.1542f, 1.2445f, 1.1984f, 1.0383f, 1.0049f, 0.9137f, 0.8678f, 0.8493f, 0.8473f, 0.8883f, 0.8660f, 0.9470f, 1.3471f, 0.0000f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_180_250 = {0.0000f, 0.0000f, 0.0000f, 1.4487f, 1.3431f, 1.1334f, 0.9285f, 0.8936f, 0.8987f, 0.8931f, 0.9289f, 0.9272f, 0.9495f, 1.0064f, 1.1732f, 1.2297f, 1.1701f, 1.0813f, 1.0572f, 1.0288f, 0.9966f, 1.0068f, 0.9998f, 1.0117f, 1.0077f, 1.0052f, 0.9837f, 0.9928f, 0.9972f, 0.9984f, 1.0066f, 1.0446f, 1.1084f, 1.1651f, 1.2555f, 1.1642f, 1.0002f, 0.9872f, 0.9008f, 0.8854f, 0.8356f, 0.8424f, 0.8525f, 0.9187f, 1.1300f, 1.3784f, 0.0000f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_250_540 = {0.0000f, 0.0000f, 0.0000f, 4.7637f, 1.4480f, 1.1482f, 0.9554f, 0.8845f, 0.8741f, 0.8996f, 0.9184f, 0.9192f, 0.9376f, 0.9782f, 1.1415f, 1.2620f, 1.1770f, 1.0821f, 1.0673f, 1.0421f, 1.0063f, 1.0028f, 0.9996f, 1.0121f, 1.0135f, 1.0033f, 0.9850f, 0.9879f, 0.9883f, 1.0012f, 1.0098f, 1.0675f, 1.1063f, 1.2033f, 1.2706f, 1.1560f, 0.9754f, 0.9611f, 0.9038f, 0.8594f, 0.8350f, 0.8536f, 0.8813f, 7.7637f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_540_1000 = {0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.8128f, 0.8985f, 0.8034f, 0.8353f, 0.8745f, 0.8878f, 0.8986f, 0.9261f, 0.9603f, 1.1329f, 1.3425f, 1.2176f, 1.1056f, 1.0734f, 1.0475f, 1.0021f, 1.0085f, 0.9956f, 1.0118f, 1.0130f, 0.9951f, 0.9791f, 0.9882f, 1.0011f, 1.0054f, 1.0338f, 1.0822f, 1.1335f, 1.2241f, 1.3479f, 1.1544f, 0.9906f, 1.0144f, 0.9211f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f};
        }
        else {
            jerValues_pt_30_50   = {0.0000f, 0.0000f, 0.0000f, 1.3057f, 1.2777f, 1.1098f, 1.0934f, 1.0520f, 1.0481f, 1.0650f, 1.1003f, 1.1364f, 1.1688f, 1.2295f, 1.3074f, 1.2788f, 1.1880f, 1.1266f, 1.0505f, 1.0304f, 0.9879f, 1.0064f, 0.9979f, 1.0079f, 0.9996f, 1.0000f, 0.9952f, 1.0017f, 1.0021f, 1.0015f, 1.0125f, 1.0679f, 1.1485f, 1.2262f, 1.3013f, 1.3304f, 1.2411f, 1.1804f, 1.1131f, 1.0779f, 1.0332f, 1.0329f, 1.1035f, 1.1078f, 1.1595f, 1.2561f, 1.4631f, 14.4011f, 0.0000f, 0.0000f};
            jerValues_pt_50_80   = {0.0000f, 0.0000f, 0.0000f, 1.2374f, 1.1827f, 1.0609f, 0.9603f, 0.9814f, 0.9795f, 0.9783f, 1.0200f, 1.0567f, 1.0951f, 1.1871f, 1.2997f, 1.2551f, 1.1685f, 1.1091f, 1.0676f, 1.0312f, 0.9934f, 1.0039f, 1.0033f, 0.9944f, 1.0174f, 1.0007f, 0.9853f, 0.9979f, 0.9889f, 1.0149f, 1.0288f, 1.0638f, 1.1256f, 1.2032f, 1.3000f, 1.2863f, 1.1855f, 1.1207f, 1.0487f, 0.9896f, 0.9766f, 0.9528f, 0.9896f, 0.9902f, 1.0771f, 1.2325f, 1.5838f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_80_120  = {0.0000f, 0.0000f, 0.0000f, 0.8003f, 1.1630f, 0.9473f, 0.8376f, 0.8730f, 0.8849f, 0.9438f, 0.9400f, 0.9723f, 1.0111f, 1.0796f, 1.2702f, 1.2451f, 1.1609f, 1.1025f, 1.0347f, 1.0257f, 0.9944f, 0.9890f, 1.0156f, 1.0097f, 1.0198f, 1.0107f, 0.9892f, 0.9847f, 0.9968f, 0.9900f, 1.0040f, 1.0557f, 1.1282f, 1.1964f, 1.2557f, 1.2756f, 1.1091f, 1.0371f, 0.9716f, 0.9176f, 0.8816f, 0.8826f, 0.8843f, 0.8860f, 1.0542f, 1.2131f, 1.3742f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_120_180 = {0.0000f, 0.0000f, 0.0000f, 0.0000f, 1.2451f, 1.0228f, 0.8464f, 0.8722f, 0.8682f, 0.8736f, 0.9120f, 0.9372f, 0.9575f, 1.0206f, 1.2075f, 1.2270f, 1.1465f, 1.0646f, 1.0586f, 1.0271f, 0.9979f, 1.0121f, 0.9950f, 1.0183f, 1.0095f, 1.0027f, 0.9872f, 0.9872f, 0.9888f, 1.0014f, 1.0079f, 1.0583f, 1.0964f, 1.1657f, 1.2337f, 1.2069f, 1.0516f, 1.0168f, 0.9275f, 0.8864f, 0.8736f, 0.8659f, 0.8813f, 0.9030f, 1.1733f, 1.4208f, 1.6644f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_180_250 = {0.0000f, 0.0000f, 0.0000f, 0.0000f, 1.8657f, 1.0699f, 0.9055f, 0.9014f, 0.8441f, 0.8694f, 0.8811f, 0.9018f, 0.9372f, 0.9800f, 1.1602f, 1.2364f, 1.1530f, 1.0751f, 1.0485f, 1.0292f, 0.9985f, 1.0036f, 0.9974f, 1.0125f, 1.0131f, 0.9983f, 1.0021f, 0.9905f, 0.9902f, 0.9937f, 1.0264f, 1.0515f, 1.0969f, 1.1647f, 1.2597f, 1.1992f, 1.0123f, 0.9921f, 0.9069f, 0.8885f, 0.8528f, 0.8535f, 0.8691f, 0.9238f, 1.1601f, 1.3991f, 2.4023f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_250_540 = {0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 1.4558f, 0.9426f, 0.8956f, 0.8674f, 0.8764f, 0.9144f, 0.9149f, 0.9739f, 1.1459f, 1.2684f, 1.1859f, 1.0912f, 1.0700f, 1.0306f, 1.0042f, 1.0080f, 1.0024f, 1.0146f, 1.0138f, 0.9977f, 0.9866f, 0.9809f, 0.9941f, 0.9976f, 1.0135f, 1.0768f, 1.1152f, 1.1865f, 1.2699f, 1.1495f, 0.9993f, 0.9786f, 0.9055f, 0.8844f, 0.8665f, 0.8623f, 0.8732f, 0.9463f, 1.1497f, 1.4829f, 1.6876f, 0.0000f, 0.0000f, 0.0000f};
            jerValues_pt_540_1000 = {0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 1.0878f, 0.9699f, 0.9457f, 1.0988f, 1.3651f, 1.2052f, 1.1082f, 1.0793f, 1.0550f, 1.0039f, 1.0100f, 1.0017f, 1.0029f, 1.0054f, 1.0000f, 0.9883f, 0.9879f, 0.9945f, 1.0055f, 1.0317f, 1.0830f, 1.1328f, 1.2299f, 1.3554f, 1.1411f, 0.9855f, 0.9586f, 0.8950f, 0.8668f, 0.8396f, 0.8115f, 0.8223f, 0.8614f, 1.0673f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f};
        }
    }
    else {
        if (isPbGoing) {

        }
        else {

        }
    }
}

//________________
int findEtaBinInVectorOfEtaBins(const float &eta) {
    for (size_t i = 0; i + 1 < etaBinEdges.size(); ++i) {
        if (etaBinEdges[i] < eta && eta < etaBinEdges[i + 1]) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

//________________
float jerExtraScalingForEta(const float &eta, const float &pt, const bool &isMc) {
    float extraScaling = 1.f;
    if ( !isMc) {
        return extraScaling;
    }

    // if ( std::abs(eta) > 0.8) {
    //     int etaBin = findEtaBinInVectorOfEtaBins(eta);
    //     if (etaBin >= 0 ) {
    //         if (30.f <= pt && pt < 50.f) {
    //             extraScaling = jerValues_pt_30_50[etaBin];
    //         }
    //         else if (50.f <= pt && pt < 80.f) {
    //             extraScaling = jerValues_pt_50_80[etaBin];
    //         }
    //         else if (80.f <= pt && pt < 120.f) {
    //             extraScaling = jerValues_pt_80_120[etaBin];
    //         }
    //         else if (120.f <= pt && pt < 180.f) {
    //             extraScaling = jerValues_pt_120_180[etaBin];
    //         }
    //         else if (180.f <= pt && pt < 250.f) {
    //             extraScaling = jerValues_pt_180_250[etaBin];
    //         }
    //         else if (250.f <= pt && pt < 540.f) {
    //             extraScaling = jerValues_pt_250_540[etaBin];
    //         }
    //         else if (540.f <= pt && pt < 1000.f) {
    //             extraScaling = jerValues_pt_540_1000[etaBin];
    //         }
    //     } // if (etaBin >= 0 )
    // } // if ( std::abs(eta) > 0.8)
    return extraScaling;
}

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
void usage() {
    std::cout << GREEN << "Usage: ./processForestSimple  [options]\n"
                << "Options:\n"
                << "  <input_file_or_list>  Input ROOT file or file list (default: empty)\n"
                << "  <output_file>         Output ROOT file name (default: empty)\n"
                << "  <0|1|2>               MC type (default 2): 0 - data, 1 - embedding, 2 - pythia\n"
                << "  <0|1>                 Direction (default 1): 0 - p-going, 1 - Pb-going\n"
                << "  <value>               pT hat sample for MC (default: 30)\n"
                << "  <-1|0|1>              JEU systematic variation (default 0): -1 - JEU-, 0 - not applied, 1 - JEU+\n"
                << "  <-99|0|1|-1>          JER systematic variation (default -99): -99 - only JEC applied, -1 - JER-, 0 - std.JER, 1 - JER+ \n"
                << "  <value>               Trigger ID (default 0): 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100\n"
                << "  <value>               Reco jet selection method (default 2): 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetIdTightLeptVeto, 3 - jetIdLoose\n"
                << RESET;
    exit (0);
}

//________________
struct GenJet {
    float pt;
    float eta;
    float phi;
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

    float recoPtJeuUp;
    float recoPtJeuDown;
    float recoPtJerDef;
    float recoPtJerUp;
    float recoPtJerDown;
};

//________________
/// Define the histograms for the analysis
struct Histograms {
    //
    // Event level histograms
    //
    std::unique_ptr<TH1D> hPtHatUnweighted;
    std::unique_ptr<TH1D> hPtHat;
    std::unique_ptr<TH1D> hVzUnweighted;
    std::unique_ptr<TH1D> hVz;

    std::unique_ptr<TH2D> hGenDijetPtAveOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hGenLeadJetPtOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hRecoDijetPtAveOverPtHatVsPtHat;
    std::unique_ptr<TH2D> hRecoLeadJetPtOverPtHatVsPtHat;

    //
    // Gen level histograms
    //

    // Gen jet histograms
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaLabUnflipped;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaLab;
    std::unique_ptr<TH1D> hGenInclusiveJetPt;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hGenInclusiveJetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH2D> hGenInclusiveJetPtEtaStdBins;
    std::unique_ptr<TH2D> hGenInclusiveJetPtPtHat;
    std::unique_ptr<TH1D> hGenInclusiveJetEtaPt80To120;

    // Gen dijet histograms
    std::unique_ptr<TH1D> hGenDijetPtAve;
    std::unique_ptr<TH1D> hGenDijetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMForwardShifted[nEtaShifts];
    std::unique_ptr<TH1D> hGenDijetEtaCMBackwardShifted[nEtaShifts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaBackwardArr[nEtaCuts];

    //
    // Reco level histograms
    //

    // Reco jet histograms
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaLabUnflipped;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaLab;
    std::unique_ptr<TH1D> hRecoInclusiveJetPt;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaCMShifted[nEtaShifts];
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaStdBins;
    std::unique_ptr<TH3D> hRecoInclusiveJetJESPtEtaStdBins;
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDefNoSFPtEtaStdBins;   // JER default variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDefPtEtaStdBins;   // JER default variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESUpPtEtaStdBins;    // JER up variation for reco inclusive jet histograms
    std::unique_ptr<TH3D> hRecoInclusiveJetJESDownPtEtaStdBins;  // JER down variation for reco inclusive jet histograms
    std::unique_ptr<TH2D> hRecoInclusiveJetPtPtHat;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaPt80To120;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaPtMatched80To120;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaPtDefNoSF80To120;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaPtDef80To120;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaPtDefExtraNoSF80To120;
    std::unique_ptr<TH1D> hRecoInclusiveJetEtaPtDefExtra80To120;

    // Reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArr[nEtaCuts];
    
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerUp[nEtaCuts];    // JER up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerUp[nEtaCuts];   // JER up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDown[nEtaCuts];  // JER down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDown[nEtaCuts]; // JER down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJerDef[nEtaCuts];   // JER default variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJerDef[nEtaCuts];  // JER default variation for reco dijet histograms

    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJeuUp[nEtaCuts];    // JEU up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJeuUp[nEtaCuts];   // JEU up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaForwardArrJeuDown[nEtaCuts];  // JEU down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaBackwardArrJeuDown[nEtaCuts]; // JEU down variation for reco dijet histograms

    //
    // Ref level histograms
    //

    // Ref jet histograms
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH1D> hRefInclusiveJetEtaLabUnflipped;
    std::unique_ptr<TH1D> hRefInclusiveJetEtaLab;
    std::unique_ptr<TH1D> hRefInclusiveJetPt;
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaStdBins;
    std::unique_ptr<TH1D> hRefInclusiveJetEtaPt80To120;

    // Ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArr[nEtaCuts];

    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJerUp[nEtaCuts];    // JER up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJerUp[nEtaCuts];   // JER up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJerDown[nEtaCuts];  // JER down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJerDown[nEtaCuts]; // JER down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJerDef[nEtaCuts];   // JER default variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJerDef[nEtaCuts];  // JER default variation for ref dijet histograms

    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJeuUp[nEtaCuts];    // JEU up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJeuUp[nEtaCuts];   // JEU up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaForwardArrJeuDown[nEtaCuts];   // JEU down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaBackwardArrJeuDown[nEtaCuts];  // JEU down variation for ref dijet histograms
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
// Create the histograms for the analysis
/// @param hs: the histograms object
/// @param isMc: true if the data is from Monte Carlo, false otherwise
void createHistograms(Histograms &hs, const bool &isMc = false, 
                      const bool &jerSyst = false, const bool &jeuSyst = false) {


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
    double jetPtBins[] = { 0., 1000. };
    const int nJetEtaBins = 144;
    double jetEtaBins[] = { -3.6, 3.6 };
    const int nJetJESBins = 100;
    double jetJESBins[] = { 0., 2. };

    // Dijet binning
    const int nDijetPtBins = 47;
    double dijetPtBins[] = { 30., 500.};
    const int nDijetEtaBins = 60;
    double dijetEtaBins[] = { -3.0, 3.0 };
    const int nDijetEtaFBBins = 30;
    double dijetEtaFBBins[] = { 0., 3. };

    //
    // Event level histograms
    //
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
        hs.hGenLeadJetPtOverPtHatVsPtHat = std::make_unique<TH2D>("hGenLeadJetPtOverPtHatVsPtHat", 
                                            "Gen leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Gen p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hGenLeadJetPtOverPtHatVsPtHat->Sumw2();
        hs.hRecoDijetPtAveOverPtHatVsPtHat = std::make_unique<TH2D>("hRecoDijetPtAveOverPtHatVsPtHat", 
                                            "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoDijetPtAveOverPtHatVsPtHat->Sumw2();
        hs.hRecoLeadJetPtOverPtHatVsPtHat = std::make_unique<TH2D>("hRecoLeadJetPtOverPtHatVsPtHat", 
                                            "Reco leading jet p_{T}/#hat{p}_{T} vs #hat{p}_{T};#hat{p}_{T} (GeV);Reco p_{T}^{lead}/#hat{p}_{T}", 
                                            100, 15., 1015.,
                                            200, 0., 4.);
        hs.hRecoLeadJetPtOverPtHatVsPtHat->Sumw2();

        //
        // Gen level histograms
        //

        // Gen jet histograms
        hs.hGenInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hGenInclusiveJetPtEtaLabUnflipped", 
                                                "Gen jet pT vs eta (lab frame, unflipped);#eta;p_{T} (GeV)", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hGenInclusiveJetEtaLabUnflipped = std::make_unique<TH1D>("hGenInclusiveJetEtaLabUnflipped", 
                                                "Gen jet eta (lab frame, unflipped);#eta;dN/d#eta", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaLabUnflipped->Sumw2();
        hs.hGenInclusiveJetEtaLab = std::make_unique<TH1D>("hGenInclusiveJetEtaLab", 
                                                "Gen jet eta (lab frame);#eta;dN/d#eta", 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaLab->Sumw2();
        hs.hGenInclusiveJetPt = std::make_unique<TH1D>("hGenInclusiveJetPt", 
                                                "Gen jet p_{T};p_{T} (GeV);dN/dp_{T}", 
                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPt->Sumw2();
        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
            hs.hGenInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%d", iShift), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[iShift]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Sumw2();
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift) 


        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShiftedUnweighted_%.3f", etaShift[0]), 
                                                Form("Gen jet eta (CM frame, shifted by %.3f, unweighted);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShiftedUnweighted[0]->Sumw2();
        hs.hGenInclusiveJetEtaCMShifted[0] = std::make_unique<TH1D>(Form("hGenInclusiveJetEtaCMShifted_%.3f", etaShift[0]), 
                                                                    Form("Gen jet eta (CM frame, shifted by %.3f);#eta_{CM};dN/d#eta_{CM}", etaShift[0]), 
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaCMShifted[0]->Sumw2();
        hs.hGenInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hGenInclusiveJetPtEtaStdBins", 
                                                                "Gen jet pT vs eta (standard bins);#eta;p_{T} (GeV)",
                                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals,
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtEtaStdBins->Sumw2();
        hs.hGenInclusiveJetPtPtHat = std::make_unique<TH2D>("hGenInclusiveJetPtPtHat", 
                                                            "Gen jet #hat{p}_{T} vs p_{T};p_{T} (GeV);#hat{p}_{T} (GeV)",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hGenInclusiveJetPtPtHat->Sumw2();
        hs.hGenInclusiveJetEtaPt80To120 = std::make_unique<TH1D>("hGenInclusiveJetEtaPt80To120", 
                                                                "Gen jet eta (p_{T} 80-120 GeV);#eta;dN/d#eta",
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hGenInclusiveJetEtaPt80To120->Sumw2();

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
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)

        //
        // Ref level histograms
        //

        // Ref jet histograms
        hs.hRefInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRefInclusiveJetPtEtaLabUnflipped", 
                                                                    "Ref jet pT vs eta (lab frame, unflipped);#eta;p_{T} (GeV)", 
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hRefInclusiveJetPtEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetEtaLabUnflipped = std::make_unique<TH1D>("hRefInclusiveJetEtaLabUnflipped", 
                                                                    "Ref jet eta (lab frame, unflipped);#eta;dN/d#eta", 
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetEtaLabUnflipped->Sumw2();
        hs.hRefInclusiveJetEtaLab = std::make_unique<TH1D>("hRefInclusiveJetEtaLab", 
                                              "Ref jet eta (lab frame);#eta;dN/d#eta", 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetEtaLab->Sumw2();
        hs.hRefInclusiveJetPt = std::make_unique<TH1D>("hRefInclusiveJetPt", 
                                                        "Ref jet p_{T};p_{T} (GeV);dN/dp_{T}", 
                                                        nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hRefInclusiveJetPt->Sumw2();
        hs.hRefInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRefInclusiveJetPtEtaStdBins", 
                                                                "Ref jet pT vs eta (standard bins);#eta;p_{T} (GeV)",
                                                                jetEtaL2L3StdBins, jetEtaL2L3StdVals,
                                                                nJetPtBins, jetPtBins[0], jetPtBins[1]);
        hs.hRefInclusiveJetPtEtaStdBins->Sumw2();
        hs.hRefInclusiveJetEtaPt80To120 = std::make_unique<TH1D>("hRefInclusiveJetEtaPt80To120", 
                                                                "Ref jet eta (p_{T} 80-120 GeV);#eta;dN/d#eta",
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRefInclusiveJetEtaPt80To120->Sumw2();

        // Ref dijet histograms
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

            if (jerSyst && isMc) {
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

            if (jeuSyst && !isMc) {
                hs.hRefDijetPtEtaForwardArrJeuUp[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaForwardJeuUp_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JEU up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaForwardArrJeuUp[iCut]->Sumw2();
                hs.hRefDijetPtEtaBackwardArrJeuUp[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaBackwardJeuUp_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JEU up);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaBackwardArrJeuUp[iCut]->Sumw2();
                hs.hRefDijetPtEtaForwardArrJeuDown[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaForwardJeuDown_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, forward, |#eta| < %.1f) vs p_{T}^{ave} (JEU down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaForwardArrJeuDown[iCut]->Sumw2();
                hs.hRefDijetPtEtaBackwardArrJeuDown[iCut] = std::make_unique<TH2D>(Form("hRefDijetPtEtaBackwardJeuDown_%d", iCut),
                                                Form("Ref dijet #eta_{CM} (CM frame, backward, |#eta| < %.1f) vs p_{T}^{ave} (JEU down);p_{T}^{ave} (GeV);#eta_{CM}", etaCuts[iCut]),
                                                nDijetPtBins, dijetPtBins[0], dijetPtBins[1],
                                                nDijetEtaFBBins, dijetEtaFBBins[0], dijetEtaFBBins[1]);
                hs.hRefDijetPtEtaBackwardArrJeuDown[iCut]->Sumw2();
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    }

    //
    // Reco level histograms
    //

    // Reco jet histograms
    hs.hRecoInclusiveJetPtEtaLabUnflipped = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaLabUnflipped", 
                                                                    "Reco jet pT vs eta (lab frame, unflipped);#eta;p_{T} (GeV)", 
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1], 
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetEtaLabUnflipped = std::make_unique<TH1D>("hRecoInclusiveJetEtaLabUnflipped", 
                                                                "Reco jet eta (lab frame, unflipped);#eta;dN/d#eta", 
                                                                nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaLabUnflipped->Sumw2();
    hs.hRecoInclusiveJetEtaLab = std::make_unique<TH1D>("hRecoInclusiveJetEtaLab", 
                                                        "Reco jet eta (lab frame);#eta;dN/d#eta", 
                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaLab->Sumw2();
    hs.hRecoInclusiveJetPt = std::make_unique<TH1D>("hRecoInclusiveJetPt", 
                                                    "Reco jet p_{T};p_{T} (GeV);dN/dp_{T}", 
                                                    nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPt->Sumw2();
    for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
        hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShiftedUnweighted_%d", iShift), 
                                              Form("Reco jet eta (CM frame, shifted by %.3f, unweighted);#eta;dN/d#eta", etaShift[iShift]), 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Sumw2();
        hs.hRecoInclusiveJetEtaCMShifted[iShift] = std::make_unique<TH1D>(Form("hRecoInclusiveJetEtaCMShifted_%d", iShift), 
                                              Form("Reco jet eta (CM frame, shifted by %.3f);#eta;dN/d#eta", etaShift[iShift]), 
                                              nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetEtaCMShifted[iShift]->Sumw2();
    } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

    hs.hRecoInclusiveJetPtEtaStdBins = std::make_unique<TH2D>("hRecoInclusiveJetPtEtaStdBins", 
                                                            "Reco jet pT vs eta (standard bins);#eta;p_{T} (GeV)",
                                                            jetEtaL2L3StdBins, jetEtaL2L3StdVals,
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtEtaStdBins->Sumw2();
    hs.hRecoInclusiveJetPtPtHat = std::make_unique<TH2D>("hRecoInclusiveJetPtPtHat", 
                                                            "Reco jet #hat{p}_{T} vs p_{T};p_{T} (GeV);#hat{p}_{T} (GeV)",
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                            nJetPtBins, jetPtBins[0], jetPtBins[1]);
    hs.hRecoInclusiveJetPtPtHat->Sumw2();

    if (isMc) {
        hs.hRecoInclusiveJetJESPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESPtEtaStdBins", 
                                                                    "Reco jet JES (reco/ref) vs pT vs eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                    nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                    nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                    nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
        hs.hRecoInclusiveJetJESPtEtaStdBins->Sumw2();
        if (jerSyst) {
            hs.hRecoInclusiveJetJESDefNoSFPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESDefNoSFPtEtaStdBins", 
                                                                        "Reco jet JES default (reco/ref) vs pT vs eta (standard bins, no scaling);JES;p_{T} (GeV);#eta",
                                                                        nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                        nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hRecoInclusiveJetJESDefNoSFPtEtaStdBins->Sumw2();
            hs.hRecoInclusiveJetJESDefPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESDefPtEtaStdBins", 
                                                                        "Reco jet JES default (reco/ref) vs pT vs eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                        nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                        nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                        nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hRecoInclusiveJetJESDefPtEtaStdBins->Sumw2();
            hs.hRecoInclusiveJetJESUpPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESUpPtEtaStdBins", 
                                                                      "Reco jet JES up (reco/ref) vs pT vs eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                      nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                      nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                      nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hRecoInclusiveJetJESUpPtEtaStdBins->Sumw2();
            hs.hRecoInclusiveJetJESDownPtEtaStdBins = std::make_unique<TH3D>("hRecoInclusiveJetJESDownPtEtaStdBins", 
                                                                         "Reco jet JES down (reco/ref) vs pT vs eta (standard bins);JES;p_{T} (GeV);#eta",
                                                                         nJetJESBins, jetJESBins[0], jetJESBins[1],
                                                                         nJetPtBins, jetPtBins[0], jetPtBins[1],
                                                                         nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
            hs.hRecoInclusiveJetJESDownPtEtaStdBins->Sumw2();
        }
    }
    hs.hRecoInclusiveJetEtaPt80To120 = std::make_unique<TH1D>("hRecoInclusiveJetEtaPt80To120", 
                                                            "Reco jet eta (p_{T} 80-120 GeV);#eta;dN/d#eta",
                                                            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaPt80To120->Sumw2();
    hs.hRecoInclusiveJetEtaPtMatched80To120 = std::make_unique<TH1D>("hRecoInclusiveJetEtaPtMatched80To120", 
                                                            "Reco jet eta (p_{T} 80-120 GeV, gen-reco matched);#eta;dN/d#eta",
                                                            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaPtMatched80To120->Sumw2();
    hs.hRecoInclusiveJetEtaPtDefNoSF80To120 = std::make_unique<TH1D>("hRecoInclusiveJetEtaPtDefNoSF80To120", 
                                                            "Reco jet eta (p_{T} 80-120 GeV, JER default, no scaling);#eta;dN/d#eta",
                                                            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaPtDefNoSF80To120->Sumw2();
    hs.hRecoInclusiveJetEtaPtDef80To120 = std::make_unique<TH1D>("hRecoInclusiveJetEtaPtDef80To120", 
                                                            "Reco jet eta (p_{T} 80-120 GeV, JER default);#eta;dN/d#eta",
                                                            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaPtDef80To120->Sumw2();
    hs.hRecoInclusiveJetEtaPtDefExtraNoSF80To120 = std::make_unique<TH1D>("hRecoInclusiveJetEtaPtDefExtraNoSF80To120", 
                                                            "Reco jet eta (p_{T} 80-120 GeV, JER default, no scaling);#eta;dN/d#eta",
                                                            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaPtDefExtraNoSF80To120->Sumw2();
    hs.hRecoInclusiveJetEtaPtDefExtra80To120 = std::make_unique<TH1D>("hRecoInclusiveJetEtaPtDefExtra80To120", 
                                                            "Reco jet eta (p_{T} 80-120 GeV, JER default);#eta;dN/d#eta",
                                                            nJetEtaBins, jetEtaBins[0], jetEtaBins[1]);
    hs.hRecoInclusiveJetEtaPtDefExtra80To120->Sumw2();

    // Reco dijet histograms
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

        if (jerSyst && isMc) {
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
        } // if (jerSyst && isMc)

        if (jeuSyst && !isMc) {
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
        } // if (jeuSyst && !isMc)
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
/// @param isMc: true if the data is from Monte Carlo, false otherwise
void setupChains(const TString &input, TChain &hltTree, TChain &eventTree, TChain &skimTree, 
                 TChain &jetTree, const bool &isMc) {

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

    mainTree.SetBranchStatus("refpt", 1);
    mainTree.SetBranchStatus("refeta", 1);

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

    mainTree.SetBranchAddress("refpt", &refJetPt);
    mainTree.SetBranchAddress("refeta", &refJetEta);
    mainTree.SetBranchAddress("refphi", &refJetPhi);

    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
void setupInput(const TString &input, TChain &hltTree, TChain &eventTree, TChain &skimTree, 
                TChain &jetTree, const bool &isMc) {
    setupChains(input, hltTree, eventTree, skimTree, jetTree, isMc);
    setupBranches(hltTree, isMc);
}

//________________
float etaCM(const float &etaLab, const float &etaShift, const bool &isPbGoing, const bool &isMc) {
    float etaCM;
    if (isPbGoing) {
        etaCM = -1.0 * (etaLab + etaShift);
    }
    else {
        etaCM = etaLab - etaShift;
    }
    return etaCM;
}

//________________
float etaLab(const float &eta, const bool &isPbGoing, const bool &isMc) {
    float etaL = eta;
    if ( isMc ) { // For embedding: Pb goes to negative, p goes to positive
        if (isPbGoing) {
            etaL = -etaL;
        }
    }
    else { // For data: p goes to negative, Pb goes to positive
        if (isPbGoing) {
        }
        else {
            etaL = -etaL;
        }
    }
    return etaL;
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
    else if (ptHatSample == 540) { ptHatRange[0] = 540.; ptHatRange[1] = std::numeric_limits<float>::max(); }
    std::cout << Form(": [%.1f, %.1f]", ptHatRange[0], ptHatRange[1]) << std::endl;
}

//________________
/// Check if event passes selection criteria
/// @param isPbGoing: true if Pb is going in the positive direction, false otherwise
/// @param isMc: true if the event is from Monte Carlo, false otherwise
/// @param ptHatSample: the ptHat sample
/// @param triggerId: the trigger ID (0 - MB, 1 - Jet60, 2 - Jet80, 3 - Jet100)
/// @return true if the event is good, false otherwise
bool isGoodEvent(const bool &isPbGoing, const bool &isMc, const int &ptHatSample, const int &triggerId) {

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
        GenJet genJet{genJetPt[iGenJet], genJetEta[iGenJet], genJetPhi[iGenJet]};
        genJets.push_back(genJet);
    }

    if ( !genJets.empty() && genJets.size() >=2 ) {
        // Sort gen jets by pT in descending order
        std::sort(genJets.begin(), genJets.end(),
                [](const GenJet &a, const GenJet &b) { return a.pt > b.pt; });
    }
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
                        refJetPt[iRecoJet], refJetEta[iRecoJet],
                        -1., -1., 
                        -1., -1., -1.};
        recoJets.push_back(recoJet);
    } // for (int iRecoJet{0}; iRecoJet < nRecoJets; ++iRecoJet)

    if ( !recoJets.empty() && recoJets.size() >= 2 ) {
        // Sort reco jets by recoPt in descending order
        std::sort(recoJets.begin(), recoJets.end(),
                [](const RecoJet &a, const RecoJet &b) { return a.recoPt > b.recoPt; });
    }
}

//________________
void fillOverweightHistograms(std::vector<GenJet> &genJets, std::vector<RecoJet> &recoJets, 
                              const double &weight, Histograms &hs,
                              bool &isGenOverweight, bool &isRecoOverweight) {
    
    isGenOverweight = {true};
    isRecoOverweight = {true};

    float dijetPtAve{0.};
    if ( genJets.size() >= 2 ) {
        dijetPtAve = 0.5 * (genJets[0].pt + genJets[1].pt);
        hs.hGenDijetPtAveOverPtHatVsPtHat->Fill(ptHat, dijetPtAve / ptHat, weight);
        hs.hGenLeadJetPtOverPtHatVsPtHat->Fill(ptHat, genJets[0].pt / ptHat, weight);
    }

    if ( recoJets.size() >= 2 ) {
        dijetPtAve = 0.5 * (recoJets[0].recoPt + recoJets[1].recoPt);
        hs.hRecoDijetPtAveOverPtHatVsPtHat->Fill(ptHat, dijetPtAve / ptHat, weight);
        hs.hRecoLeadJetPtOverPtHatVsPtHat->Fill(ptHat, recoJets[0].recoPt / ptHat, weight);
    }
}

//________________
void processGenJets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                    std::vector<GenJet> &genJets, Histograms &hs, const double &ptHat) {

    float genJetEtaLabFlipped{0.};
    float genJetEtaCM{0.};
    
    // Loop over gen jets and fill histograms
    for (const auto &genJet : genJets) {
        genJetEtaLabFlipped = etaLab(genJet.eta, isPbGoing, isMc);

        hs.hGenInclusiveJetPtEtaLabUnflipped->Fill(genJet.eta, genJet.pt, weight);
        hs.hGenInclusiveJetEtaLabUnflipped->Fill(genJet.eta, weight);
        hs.hGenInclusiveJetEtaLab->Fill(genJetEtaLabFlipped, weight);
        hs.hGenInclusiveJetPt->Fill(genJet.pt, weight);
        if (genJet.pt > 80. && genJet.pt < 120.) {
            hs.hGenInclusiveJetEtaPt80To120->Fill(genJet.eta, weight);
        }

        for (int iShift{0}; iShift < nEtaShifts; ++iShift) {
            genJetEtaCM = etaCM(genJet.eta, etaShift[iShift], isPbGoing, isMc);
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Fill( genJetEtaCM, 1.); // Unweighted
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Fill(genJetEtaCM, weight);
        } // for (int iShift{0}; iShift < nEtaShifts; ++iShift)

        hs.hGenInclusiveJetPtEtaStdBins->Fill(genJet.eta, genJet.pt, weight);
        hs.hGenInclusiveJetPtPtHat->Fill(genJet.pt, ptHat, weight);
    }
}

//________________
void processGenDijets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                      std::vector<GenJet> &genJets, Histograms &hs) {

    // Must be at least 2 jets to form a dijet system
    if (genJets.size() < 2) return;

    // Sort jets by pT to identify leading and subleading jets
    std::sort(genJets.begin(), genJets.end(), [](const GenJet &a, const GenJet &b) { return a.pt > b.pt; });

    GenJet leadingJet = genJets[0];
    GenJet subleadingJet = genJets[1];

    // Make a dijet cut
    if (leadingJet.pt < 50. || subleadingJet.pt < 40.) return;
    float dphi = leadingJet.phi - subleadingJet.phi;
    if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
    if (std::abs(dphi) < TMath::TwoPi() / 3.) return;
    float dijetPtAve = 0.5 * (leadingJet.pt + subleadingJet.pt);

    // Check shifts
    for (int iShift = 0; iShift < nEtaShifts; ++iShift) {

        float leadingEtaCMShifted = etaCM(leadingJet.eta, etaShift[iShift], isPbGoing, isMc);
        float subleadingEtaCMShifted = etaCM(subleadingJet.eta, etaShift[iShift], isPbGoing, isMc);
        float dijetEtaCMShifted = 0.5 * (leadingEtaCMShifted + subleadingEtaCMShifted);
        

        if (std::abs(leadingEtaCMShifted) > 1.9 || std::abs(subleadingEtaCMShifted) > 1.9) continue;
        
        hs.hGenDijetPtAve->Fill(dijetPtAve, weight);

        // if (dijetPtAve < 60. || dijetPtAve > 80.) continue;
        hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Fill(dijetEtaCMShifted);
        hs.hGenDijetEtaCMShifted[iShift]->Fill(dijetEtaCMShifted, weight);
        (dijetEtaCMShifted >= 0.) ? hs.hGenDijetEtaCMForwardShifted[iShift]->Fill(dijetEtaCMShifted, weight) : 
                                    hs.hGenDijetEtaCMBackwardShifted[iShift]->Fill( std::abs(dijetEtaCMShifted), weight);
    } // for (int iShift = 0; iShift < nEtaShifts; ++iShift)

    // For each eta cut, fill corresponding histograms
    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        float leadEtaCM = etaCM(leadingJet.eta, 0.465f, isPbGoing, isMc);
        float subleadEtaCM = etaCM(subleadingJet.eta, 0.465f, isPbGoing, isMc);
        if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;
        float dijetEtaCM = 0.5 * (leadEtaCM + subleadEtaCM);
        hs.hGenDijetPtEtaCMArr[iCut]->Fill(dijetPtAve, dijetEtaCM, weight);
        (dijetEtaCM >= 0.) ? hs.hGenDijetPtEtaForwardArr[iCut]->Fill(dijetPtAve, dijetEtaCM, weight) :
                             hs.hGenDijetPtEtaBackwardArr[iCut]->Fill(dijetPtAve, std::abs(dijetEtaCM), weight);
    } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
}


//________________
/// Check if a reconstructed jet passes the selection criteria
/// @param recoJet: the reconstructed jet object
/// @param jetSelectionMethod: the jet selection method (0 - no selection, 1 - trackMax, 2 - jetIdTightLept, 3 - jetIdLoose)
/// @return true if the jet is good, false otherwise
bool isGoodRecoJet(const RecoJet &recoJet, const int &jetSelectionMethod) {
    // jetSelectionMethod: 0 - no selection, 1 - trackMax, 2 - jetIdTightLept, 3 - jetIdLoose

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

        // std::cout << Form("\n---------\nJet eta = %.2f, NHF = %.3f, NEF = %.3f, CHF = %.3f, MUF = %.3f, CEF = %.3f, CHM = %d, CEM = %d, NHM = %d, NEM = %d, MUM = %d \t [chargedMult: %d, neutralMult: %d]\n",
        //                   recoJet.recoEta, recoJet.recoPfNHF, recoJet.recoPfNEF, 
        //                   recoJet.recoPfCHF, recoJet.recoPfMUF, recoJet.recoPfCEF, 
        //                   recoJet.recoPfCHM, recoJet.recoPfCEM, recoJet.recoPfNHM, 
        //                   recoJet.recoPfNEM, recoJet.recoPfMUM, chargedMult, neutralMult);


                // Check cuts depending on jet pseudorapdity
        if ( std::abs( recoJet.recoEta ) <= 2.7 ) {
            
            passNHF = ( recoJet.recoPfNHF < neutFracCut ) ? true : false;
            passNEF = ( recoJet.recoPfNEF < neutFracCut ) ? true : false;
            passNumOfConstituents = ( numberOfConstituents > 1 ) ? true : false;

            if ( jetSelectionMethod == 2 ) { 
                passMuonFrac = ( recoJet.recoPfMUF < 0.8 ) ? true : false; 
            } // if ( jetSelectionMethod == 2 )

            if( std::abs( recoJet.recoEta ) <= 2.4 ) {
                passChargedFrac = ( recoJet.recoPfCHF > 0 ) ? true : false;
                passChargedMult = ( chargedMult > 0 ) ? true : false;
                passChargedEmFrac = ( recoJet.recoPfCEF < chargedEmFracCut ) ? true : false;
            } // if( std::abs( recoJetEta[jetIndex] ) <= 2.4 )

        } // if ( std::abs( recoJetEta[jetIndex] ) <= 2.7 )
        else if ( std::abs( recoJet.recoEta ) <= 3.0) {

            passNEF = ( recoJet.recoPfNEF > 0.01 ) ? true : false;
            passNHF = ( recoJet.recoPfNHF < 0.98 ) ? true : false;
            passNeutralMult = ( neutralMult > 2 ) ? true : false;

        } // else if ( std::abs( recoJetEta[jetIndex] ) <= 3.0)
        else  {
            passNEF = ( recoJet.recoPfNEF < 0.9 ) ? true : false;
            passNeutralMult = (neutralMult > 10 ) ? true : false; // CAUTION: JET MET says it should be >10
        } // else 

        isGood = passNHF && passNEF && passNumOfConstituents && passMuonFrac && passChargedFrac && 
                 passChargedMult && passChargedEmFrac && passNeutralMult;

        // std::cout << Form("Jet %d: passNHF: %d, passNEF: %d, passNumOfConstituents: %d, passMuonFrac: %d, passChargedFrac: %d, passChargedMult: %d, passChargedEmFrac: %d, passNeutralMult: %d \t [isGood: %d]\n",
        //                   jetIndex, passNHF, passNEF, passNumOfConstituents, passMuonFrac, passChargedFrac, 
        //                   passChargedMult, passChargedEmFrac, passNeutralMult, isGood);
    } // else if (jetSelectionMethod == 2 || jetSelectionMethod == 3)

    return isGood;
}

//________________
/// Calculate smeared momentum addition due to JER (JER: defatul, up and down)
/// @param ptDef: default pt (initial value is refPt), returned value is smeared pt component (not total smeared pt)
/// @param ptUp: up pt (initial value is refPt), returned value is smeared pt component (not total smeared pt)
/// @param ptDown: down pt (initial value is refPt), returned value is smeared pt component (not total smeared pt)
/// @param ptDefNoScale: default pt without scaling (initial value is refPt), returned value is smeared pt component (not total smeared pt)
/// @param eta: reco jet pseudorapidity
/// @param fJERSmearFunc: JER smearing function vs pT
/// @param fRndm: random number generator (can be nullptr to use std::random)
void calculateResolutionSmearing(float &ptDef, float &ptUp, float &ptDown, float &ptDefNoScale, 
                                 const float &eta, TF1 *fJERSmearFunc, TRandom3 *fRndm) {
    // Find eta bin for JER constants: fJerEtaLow <= recoEta < fJerEtaHi
    int iEtaBin = -1;
    float ptRef = ptDef; // Use default pt as reference for smearing

    for (size_t i = 0; i < fJerEtaLow.size(); ++i) {
        if ( fJerEtaLow[i] <= eta && eta < fJerEtaHi[i] ) {
            iEtaBin = static_cast<int>(i);
            break;
        }
    }

    if ( iEtaBin < 0 ) {
        std::cerr << RED << Form("No JER bin found for eta = %f", eta) << RESET<< std::endl;
        return;
    }

    double scaleFactorDef = std::sqrt( std::max(fJerDef[iEtaBin] * fJerDef[iEtaBin] - 1., 0.) );
    double scaleFactorUp = std::sqrt( std::max(fJerHi[iEtaBin] * fJerHi[iEtaBin] - 1., 0.) );
    double scaleFactorDown = std::sqrt( std::max(fJerLow[iEtaBin] * fJerLow[iEtaBin] - 1., 0.) );

    // std::cout << Form("JER scale factors for eta = %f: etaBin: %d, sfDef = %f, sfUp = %f, sfDown = %f, fJerDef = %f, fJerHi = %f, fJerLow = %f", 
    //                   eta, iEtaBin,
    //                   scaleFactorDef, scaleFactorUp, scaleFactorDown, 
    //                   fJerDef[iEtaBin], fJerHi[iEtaBin], fJerLow[iEtaBin]) << std::endl;

    double sigmaSmearDef{0.};
    double sigmaSmearUp{0.};
    double sigmaSmearDown{0.};
    double sigmaSmearDefNoScale{0.};
    double evalValue{0.};
    if ( ptDef <= 30.) {
        evalValue = fJERSmearFunc->Eval( 31. );
    }
    else if ( ptDef >= 800 ) {
        evalValue = fJERSmearFunc->Eval( 799. );
    }
    else {
        evalValue = fJERSmearFunc->Eval( ptDef );
    }
    sigmaSmearDef = scaleFactorDef * evalValue;
    sigmaSmearDefNoScale = evalValue;
    sigmaSmearUp = scaleFactorUp * evalValue;
    sigmaSmearDown = scaleFactorDown * evalValue;

    double extraCorrDef = fRndm->Gaus( 1., sigmaSmearDef );
    double extraCorrDefNoScale = fRndm->Gaus( 1., sigmaSmearDefNoScale );
    double extraCorrUp = fRndm->Gaus( 1., sigmaSmearUp );
    double extraCorrDown = fRndm->Gaus( 1., sigmaSmearDown );

    if (fRndm) {
        ptDef = ptDef * (extraCorrDef - 1.);
        ptDefNoScale = ptDefNoScale * (extraCorrDefNoScale - 1.);
        ptUp = ptUp * (extraCorrUp - 1.);
        ptDown = ptDown * (extraCorrDown - 1.);
    }
    else {
        // Use C++ std random as an alternative to ROOT TRandom3
        static thread_local std::mt19937_64 rng(
            static_cast<unsigned long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count())
        );
        std::normal_distribution<double> gaus(1.0, sigmaSmearDef);
        ptDef = ptDef * (gaus(rng) - 1.);
        gaus.param(std::normal_distribution<double>::param_type(1.0, sigmaSmearDefNoScale));
        ptDefNoScale = ptDefNoScale * (gaus(rng) - 1.);
        gaus.param(std::normal_distribution<double>::param_type(1.0, sigmaSmearUp));
        ptUp = ptUp * (gaus(rng) - 1.);
        gaus.param(std::normal_distribution<double>::param_type(1.0, sigmaSmearDown));
        std::normal_distribution<double> gDown(1.0, sigmaSmearDown);
        ptDown = ptDown * (gDown(rng) - 1.);
    }
    // std::cout << Form("refPt: %f, eta: %f jerBin: %d eval: %f", ptRef, eta, iEtaBin, evalValue) << std::endl;
    // std::cout << Form("Default: %f sf: %f, sigma: %f, extraCorr: %f, ptSmeared: %f", fJerDef[iEtaBin], scaleFactorDef, sigmaSmearDef, extraCorrDef, ptDef) << std::endl;
    // std::cout << Form("Default no scale: %f sf: %f, sigma: %f, extraCorr: %f, ptSmeared: %f", fJerDef[iEtaBin], 1.f, sigmaSmearDefNoScale, extraCorrDefNoScale, ptDefNoScale) << std::endl;
    // std::cout << Form("Up: %f sf: %f, sigma: %f, extraCorr: %f, ptSmeared: %f", fJerHi[iEtaBin], scaleFactorUp, sigmaSmearUp, extraCorrUp, ptUp) << std::endl;
    // std::cout << Form("Down: %f sf: %f, sigma: %f, extraCorr: %f, ptSmeared: %f", fJerLow[iEtaBin], scaleFactorDown, sigmaSmearDown, extraCorrDown, ptDown) << std::endl;
    // std::cout << Form("refPt: %f, smeared ptDef: %f, ptDefNoScale: %f, ptUp: %f, ptDown: %f", ptRef, ptRef+ptDef, ptRef+ptDefNoScale, ptRef+ptUp, ptRef+ptDown) << std::endl;
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
/// Calculate momentum smearing (JER: defatul, up and down)
/// @param ptDef: default pt (initial value is refPt)
/// @param ptUp: up pt (initial value is refPt)
/// @param ptDown: down pt (initial value is refPt)
/// @param ptDefNoScale: default pt without scaling (initial value is refPt)
/// @param eta: reco jet pseudorapidity
/// @param fJERSmearFunc: JER smearing function vs pT
/// @param fRndm: random number generator (can be nullptr to use std::random)
void calculateSmearedMomentum(float &ptDef, float &ptUp, float &ptDown, float &ptDefNoScale, 
                              float &eta, TF1 *fJERSmearFunc, TRandom3 *fRndm) {

    // The next variables will store only the smearing component (i.e., the additional momentum due to JER), 
    // which will be added to the scaled pt later
    float ptSmearDef = ptDef;
    float ptSmearUp = ptUp;
    float ptSmearDown = ptDown;
    float ptSmearDefNoScale = ptDefNoScale;

    calculateResolutionSmearing(ptSmearDef, ptSmearUp, ptSmearDown, ptSmearDefNoScale, eta, fJERSmearFunc, fRndm);
    calculateScaleSmearing(ptDef, eta);
    calculateScaleSmearing(ptUp, eta);
    calculateScaleSmearing(ptDown, eta);
    calculateScaleSmearing(ptDefNoScale, eta);

    ptDef += ptSmearDef;
    ptUp += ptSmearUp;
    ptDown += ptSmearDown;
    ptDefNoScale += ptSmearDefNoScale;
}

//________________
void processRecoJets(const bool &isPbGoing, const bool &isMc, const double &weight, 
                     std::vector<RecoJet> &recoJets, Histograms &hs, 
                     JetCorrector &jec, JetUncertainty &jeu,
                     const int &jetSelectionMethod, const double &ptHat,
                     TF1 *fJERSmearFunc, TRandom3 *fRndm) {
    
    float recoJetEtaLabFlipped{0.};
    float recoJetEtaCM{0.};
    float recoJetJES{0.};
    float recoJetPtJeuUp{0.};
    float recoJetPtJeuDown{0.};
    float recoJetPtJerDef{0.};
    float recoJetPtJerUp{0.};
    float recoJetPtJerDown{0.};
    float recoJerPtJerDefNoSF{0.};
    bool  hasMatchingGenJet{false};
    float extraCorrFactor{1.};

    // Loop over reconstructed jets
    for (auto &recoJet : recoJets) {
        // Here we can apply a preselection on the reconstructed jets if needed
        // if (recoJet.recoPt() < 10.f) continue; // Precut at 10 GeV to save computational time

        if ( !isGoodRecoJet(recoJet, jetSelectionMethod) ) continue;

        hasMatchingGenJet = recoJet.refPt > 0;

        // std::cout << Form("Jet: eta = %.2f, phi = %.2f, rawPt = %.1f, recoPt = %.1f", 
        //                   recoJet.recoEta, recoJet.recoPhi, recoJet.recoRawPt, recoJet.recoPt) << std::endl;

        // Original correction
        if (recoJet.recoPt > 80.f && recoJet.recoPt < 120.f) {
            hs.hRecoInclusiveJetEtaPt80To120->Fill(recoJet.recoEta, weight);
            if (recoJet.refPt > 0) {
                hs.hRecoInclusiveJetEtaPtMatched80To120->Fill(recoJet.recoEta, weight);
                hs.hRefInclusiveJetEtaPt80To120->Fill(recoJet.refEta, weight);
            }
        }

        // Perform JES & JER smearing calculations
        if (isMc && hasMatchingGenJet) {
            recoJetPtJerDef = recoJet.refPt;
            recoJetPtJerUp = recoJet.refPt;
            recoJetPtJerDown = recoJet.refPt;
            recoJerPtJerDefNoSF = recoJet.refPt;
            calculateSmearedMomentum(recoJetPtJerDef, recoJetPtJerUp, 
                                     recoJetPtJerDown, recoJerPtJerDefNoSF,
                                     recoJet.recoEta, 
                                     fJERSmearFunc, fRndm);

            recoJet.recoPtJerDef = recoJetPtJerDef;
            recoJet.recoPtJerUp = recoJetPtJerUp;
            recoJet.recoPtJerDown = recoJetPtJerDown;

            hs.hRecoInclusiveJetJESPtEtaStdBins->Fill(recoJet.recoPt / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
            hs.hRecoInclusiveJetJESDefNoSFPtEtaStdBins->Fill(recoJerPtJerDefNoSF / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
            hs.hRecoInclusiveJetJESDefPtEtaStdBins->Fill(recoJetPtJerDef / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
            hs.hRecoInclusiveJetJESUpPtEtaStdBins->Fill(recoJetPtJerUp / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);
            hs.hRecoInclusiveJetJESDownPtEtaStdBins->Fill(recoJetPtJerDown / recoJet.refPt, recoJet.refPt, recoJet.refEta, weight);

            // Smeared, but without JER scaling factor (SF)
            if (recoJerPtJerDefNoSF > 80.f && recoJerPtJerDefNoSF < 120.f) {
                hs.hRecoInclusiveJetEtaPtDefNoSF80To120->Fill(recoJet.recoEta, weight);
            }
            // Smeared with JER scaling factor (SF), a.k.a. standard JER 
            if (recoJetPtJerDef > 80.f && recoJetPtJerDef < 120.f) {
                hs.hRecoInclusiveJetEtaPtDef80To120->Fill(recoJet.recoEta, weight);
            }

            // if ( std::abs( recoJet.recoEta) > 0.8f && recoJetPtJerDef > 15.f) {
            //     extraCorrFactor = jerExtraScalingForEta( recoJet.recoEta, recoJet.refPt, isMc);
            //     // std::cout << Form("eta = %.2f, refPt = %.1f, recoJetPtJerDef = %.1f, extraCorrFactor: %.3f", 
            //                          recoJet.recoEta, recoJet.refPt, recoJetPtJerDef, extraCorrFactor) << std::endl;
            //     // recoJerPtJerDefNoSF *= extraCorrFactor;
            //     // recoJetPtJerDef *= extraCorrFactor;
            //     // recoJetPtJerUp *= extraCorrFactor;
            //     // recoJetPtJerDown *= extraCorrFactor;
            // }
            // // Smeared + eta dependence, but without JER scaling factor (SF)
            // if (recoJerPtJerDefNoSF * extraCorrFactor > 80.f && recoJerPtJerDefNoSF * extraCorrFactor < 120.f) {
            //     hs.hRecoInclusiveJetEtaPtDefExtraNoSF80To120->Fill(recoJet.recoEta, weight);
            // }
            // // Smeared + eta dependence with JER scaling factor (SF)
            // if (recoJetPtJerDef * extraCorrFactor > 80.f && recoJetPtJerDef * extraCorrFactor < 120.f) {
            //     hs.hRecoInclusiveJetEtaPtDefExtra80To120->Fill(recoJet.recoEta, weight);
            // }

            // Fill reference jet histograms
            hs.hRefInclusiveJetPtEtaLabUnflipped->Fill(recoJet.refEta, recoJet.refPt, weight);
            hs.hRefInclusiveJetEtaLabUnflipped->Fill(recoJet.refEta, weight);
            hs.hRefInclusiveJetEtaLab->Fill(etaLab(recoJet.refEta, isPbGoing, isMc), weight);
            hs.hRefInclusiveJetPt->Fill(recoJet.refPt, weight);
            hs.hRefInclusiveJetPtEtaStdBins->Fill(recoJet.refEta, recoJet.refPt, weight);
        } // if (isMc && hasMatchingGenJet)

        // Perform JEU calculations
        if ( !isMc ) {
            jeu.SetJetPT( recoJet.recoPt );
            jeu.SetJetEta( recoJet.recoEta );
            jeu.SetJetPhi( recoJet.recoPhi );
            recoJet.recoPtJeuUp = recoJet.recoPt * (1. + jeu.GetUncertainty().first);
            recoJet.recoPtJeuDown = recoJet.recoPt * (1. - jeu.GetUncertainty().second);
        } // if ( !isMc )

        recoJetEtaLabFlipped = etaLab(recoJet.recoEta, isPbGoing, isMc);
        // recoJetEtaCM = etaCM(recoJetEta[iRecoJet], 0.465f, isPbGoing, isMc);

        // std::cout << Form("Reco jet rawpT = %.1f GeV, corrected pT = %.1f GeV, ref pT: %.1f, eta = %.2f", 
        //                   recoJet.recoPtRaw, recoJet.recoPt, recoJet.refPt, recoJet.recoEta) << std::endl;

        hs.hRecoInclusiveJetPtEtaLabUnflipped->Fill(recoJet.recoEta, recoJet.recoPt, weight);
        hs.hRecoInclusiveJetEtaLabUnflipped->Fill(recoJet.recoEta, weight);
        hs.hRecoInclusiveJetEtaLab->Fill(recoJetEtaLabFlipped, weight);
        hs.hRecoInclusiveJetPt->Fill(recoJet.recoPt, weight);
        hs.hRecoInclusiveJetPtPtHat->Fill(recoJet.recoPt, ptHat, weight);
        // hs.hRecoInclusiveJetEtaCMShiftedUnweighted[nEtaShifts];
        // hs.hRecoInclusiveJetEtaCMShifted[nEtaShifts];

        hs.hRecoInclusiveJetPtEtaStdBins->Fill(recoJet.recoEta, recoJet.recoPt, weight);
    } // for (const auto &recoJet : recoJets)
}

//________________
void processRecoDijets(const bool &isPbGoing, const bool &isMc, const float &weight, 
                       std::vector<RecoJet> &recoJets, Histograms &hs,
                       const bool &doJerSyst = false, const bool &doJeuSyst = false) {

    // std::cout << "processRecoDijets input parameters:" << std::endl;
    // std::cout << "  isPbGoing: " << (isPbGoing ? "true" : "false") << std::endl;
    // std::cout << "  isMc: " << (isMc ? "true" : "false") << std::endl;
    // std::cout << "  weight: " << weight << std::endl;
    // std::cout << "  recoJets.size(): " << recoJets.size() << std::endl;
    // std::cout << "  doJerSyst: " << (doJerSyst ? "true" : "false") << std::endl;
    // std::cout << "  doJeuSyst: " << (doJeuSyst ? "true" : "false") << std::endl;

    // Must be at least 2 jets to form a dijet system
    if (recoJets.size() < 2) return;

    // Sort jets by pT to identify leading and subleading jets
    std::sort(recoJets.begin(), recoJets.end(), [](const RecoJet &a, const RecoJet &b) { return a.recoPt > b.recoPt; });
    const auto &leadingJet = recoJets[0];
    const auto &subleadingJet = recoJets[1];
    float dphi = leadingJet.recoPhi - subleadingJet.recoPhi;
    if (dphi > TMath::Pi()) dphi -= TMath::TwoPi();
    if (dphi < -TMath::Pi()) dphi += TMath::TwoPi();

    // std::cout << Form("Leading jet: pT = %.1f GeV, eta = %.2f, phi = %.2f; Subleading jet: pT = %.1f GeV, eta = %.2f, phi = %.2f; dphi = %.2f", 
    //                   leadingJet.recoPt, leadingJet.recoEta, leadingJet.recoPhi, 
    //                   subleadingJet.recoPt, subleadingJet.recoEta, subleadingJet.recoPhi, dphi) << std::endl;

    if (std::abs(dphi) < TMath::TwoPi() / 3.) return; // Dijet azimuthal angle cut
    float dijetPtAve = 0.5 * (leadingJet.recoPt + subleadingJet.recoPt);

    bool leadHasMatchingRef{false}; 
    bool subleadHasMatchingRef{false}; 
    if (isMc) {
        leadHasMatchingRef = (leadingJet.refPt > 0.);
        subleadHasMatchingRef = (subleadingJet.refPt > 0.);
    }

    // std::cout << Form("Dijet system: lead pT = %.1f GeV, sublead pT = %.1f GeV, dijet pT ave = %.1f GeV, lead eta = %.2f, sublead eta = %.2f", 
    //                   leadingJet.recoPt, subleadingJet.recoPt, dijetPtAve, leadingJet.recoEta, subleadingJet.recoEta) << std::endl;

    // Fill for the given eta cuts
    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        float leadEtaCM = etaCM(leadingJet.recoEta, 0.465f, isPbGoing, isMc);
        float subleadEtaCM = etaCM(subleadingJet.recoEta, 0.465f, isPbGoing, isMc);
        if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;
        float dijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);

        // std::cout << Form("Reco dijet: lead pT = %.1f GeV, sublead pT = %.1f GeV, dijet pT ave = %.1f GeV, lead eta = %.2f, sublead eta = %.2f, dijet etaCM = %.2f", 
        //                   leadingJet.recoPt, subleadingJet.recoPt, dijetPtAve, leadingJet.recoEta, subleadingJet.recoEta, dijetEtaCM) << std::endl;

        hs.hRecoDijetPtEtaCMArr[iCut]->Fill(dijetPtAve, dijetEtaCM, weight);
        (dijetEtaCM >= 0.) ? hs.hRecoDijetPtEtaForwardArr[iCut]->Fill(dijetPtAve, dijetEtaCM, weight) :
                             hs.hRecoDijetPtEtaBackwardArr[iCut]->Fill(dijetPtAve, std::abs(dijetEtaCM), weight);

        if (isMc && leadHasMatchingRef && subleadHasMatchingRef) {
            float leadRefEtaCM = etaCM(leadingJet.refEta, 0.465f, isPbGoing, isMc);
            float subleadRefEtaCM = etaCM(subleadingJet.refEta, 0.465f, isPbGoing, isMc);
            float dijetRefEtaCM = 0.5f * (leadRefEtaCM + subleadRefEtaCM);
            float dijetRefPtAve = 0.5f * (leadingJet.refPt + subleadingJet.refPt);
            hs.hRefDijetPtEtaCMArr[iCut]->Fill(dijetPtAve, dijetRefEtaCM, weight);
            (dijetRefEtaCM >= 0.) ? hs.hRefDijetPtEtaForwardArr[iCut]->Fill(dijetRefPtAve, dijetRefEtaCM, weight) :
                                    hs.hRefDijetPtEtaBackwardArr[iCut]->Fill(dijetRefPtAve, std::abs(dijetRefEtaCM), weight);
        } // if (isMc && leadHasMatchingRef && subleadHasMatchingRef)
    } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)

    auto fillRecoDijetSystematic = [&](auto getPt,
                                       std::unique_ptr<TH2D> *forwardArr,
                                       std::unique_ptr<TH2D> *backwardArr,
                                       std::unique_ptr<TH2D> *refForwardArr = nullptr,
                                       std::unique_ptr<TH2D> *refBackwardArr = nullptr) {
        std::vector<RecoJet> sortedJets = recoJets;
        std::sort(sortedJets.begin(), sortedJets.end(), [&](const RecoJet &a, const RecoJet &b) {
            return getPt(a) > getPt(b);
        });

        const auto &leadingJetSyst = sortedJets[0];
        const auto &subleadingJetSyst = sortedJets[1];
        float systDphi = leadingJetSyst.recoPhi - subleadingJetSyst.recoPhi;
        if (systDphi > TMath::Pi()) systDphi -= TMath::TwoPi();
        if (systDphi < -TMath::Pi()) systDphi += TMath::TwoPi();
        if (std::abs(systDphi) < TMath::TwoPi() / 3.) return;

        float systDijetPtAve = 0.5f * (getPt(leadingJetSyst) + getPt(subleadingJetSyst));
        float systRefDijetPtAve = 0.5f * (leadingJetSyst.refPt + subleadingJetSyst.refPt);

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            float leadEtaCM = etaCM(leadingJetSyst.recoEta, 0.465f, isPbGoing, isMc);
            float subleadEtaCM = etaCM(subleadingJetSyst.recoEta, 0.465f, isPbGoing, isMc);
            if (std::abs(leadEtaCM) > etaCuts[iCut] || std::abs(subleadEtaCM) > etaCuts[iCut]) continue;

            float systDijetEtaCM = 0.5f * (leadEtaCM + subleadEtaCM);
            if (systDijetEtaCM >= 0.) {
                forwardArr[iCut]->Fill(systDijetPtAve, systDijetEtaCM, weight);
                if (refForwardArr && leadingJetSyst.refPt > 0. && subleadingJetSyst.refPt > 0.) {
                    refForwardArr[iCut]->Fill(systRefDijetPtAve, systDijetEtaCM, weight);
                }
            }
            else {
                backwardArr[iCut]->Fill(systDijetPtAve, std::abs(systDijetEtaCM), weight);
                if (refBackwardArr && leadingJetSyst.refPt > 0. && subleadingJetSyst.refPt > 0.) {
                    refBackwardArr[iCut]->Fill(systRefDijetPtAve, std::abs(systDijetEtaCM), weight);
                }
            }
        } // for (int iCut{0}; iCut < nEtaCuts; ++iCut)
    };

    if (isMc && doJerSyst) {
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJerDef; },
                                hs.hRecoDijetPtEtaForwardArrJerDef,
                                hs.hRecoDijetPtEtaBackwardArrJerDef,
                                hs.hRefDijetPtEtaForwardArrJerDef,
                                hs.hRefDijetPtEtaBackwardArrJerDef);
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJerUp; },
                                hs.hRecoDijetPtEtaForwardArrJerUp,
                                hs.hRecoDijetPtEtaBackwardArrJerUp,
                                hs.hRefDijetPtEtaForwardArrJerUp,
                                hs.hRefDijetPtEtaBackwardArrJerUp);
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJerDown; },
                                hs.hRecoDijetPtEtaForwardArrJerDown,
                                hs.hRecoDijetPtEtaBackwardArrJerDown,
                                hs.hRefDijetPtEtaForwardArrJerDown,
                                hs.hRefDijetPtEtaBackwardArrJerDown);
    }

    if (!isMc && doJeuSyst) {
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJeuUp; },
                                hs.hRecoDijetPtEtaForwardArrJeuUp,
                                hs.hRecoDijetPtEtaBackwardArrJeuUp,
                                hs.hRefDijetPtEtaForwardArrJeuUp,
                                hs.hRefDijetPtEtaBackwardArrJeuUp);
        fillRecoDijetSystematic([](const RecoJet &jet) { return jet.recoPtJeuDown; },
                                hs.hRecoDijetPtEtaForwardArrJeuDown,
                                hs.hRecoDijetPtEtaBackwardArrJeuDown,
                                hs.hRefDijetPtEtaForwardArrJeuDown,
                                hs.hRefDijetPtEtaBackwardArrJeuDown);
    }
}

//________________
void writeOutput(TString &oFileName, Histograms &hs, const bool &isMc, const bool &doJerSyst, const bool &doJeuSyst) {

    std::cout << "Writing output to file: " << oFileName.Data();
    // Output file
    int compressionSetting = 208; // LZMA compression
    // TFile *fOut = TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting);
    auto fOut = std::unique_ptr<TFile>( TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting) );

    //
    // Event level histograms
    //
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
        hs.hGenLeadJetPtOverPtHatVsPtHat->Write();
        hs.hRecoDijetPtAveOverPtHatVsPtHat->Write();
        hs.hRecoLeadJetPtOverPtHatVsPtHat->Write();

        // Gen jets
        hs.hGenInclusiveJetPtEtaLabUnflipped->Write();
        hs.hGenInclusiveJetEtaLabUnflipped->Write();
        hs.hGenInclusiveJetEtaLab->Write();
        hs.hGenInclusiveJetPt->Write();
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift]->Write();
            hs.hGenInclusiveJetEtaCMShifted[iShift]->Write();
        }
        hs.hGenInclusiveJetPtEtaStdBins->Write();
        hs.hGenInclusiveJetPtPtHat->Write();
        hs.hGenInclusiveJetEtaPt80To120->Write();

        // Gen dijets
        hs.hGenDijetPtAve->Write();
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            hs.hGenDijetEtaCMShiftedUnweighted[iShift]->Write();
            hs.hGenDijetEtaCMShifted[iShift]->Write();
            hs.hGenDijetEtaCMForwardShifted[iShift]->Write();
            hs.hGenDijetEtaCMBackwardShifted[iShift]->Write();
        }

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hGenDijetPtEtaCMArr[iCut]->Write();
            hs.hGenDijetPtEtaForwardArr[iCut]->Write();
            hs.hGenDijetPtEtaBackwardArr[iCut]->Write();
        }

        //
        // Ref level histograms
        //
        
        // Ref jets
        hs.hRefInclusiveJetPtEtaLabUnflipped->Write();
        hs.hRefInclusiveJetEtaLabUnflipped->Write();
        hs.hRefInclusiveJetEtaLab->Write();
        hs.hRefInclusiveJetPt->Write();
        hs.hRefInclusiveJetPtEtaStdBins->Write();
        hs.hRefInclusiveJetEtaPt80To120->Write();

        // Ref dijets
        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            hs.hRefDijetPtEtaCMArr[iCut]->Write();
            hs.hRefDijetPtEtaForwardArr[iCut]->Write();
            hs.hRefDijetPtEtaBackwardArr[iCut]->Write();
        }

    } // if (isMc)

    //
    // Reco level histograms
    //
    
    // Reco jets
    hs.hRecoInclusiveJetPtEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetEtaLabUnflipped->Write();
    hs.hRecoInclusiveJetEtaLab->Write();
    hs.hRecoInclusiveJetPt->Write();
    for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
        hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift]->Write();
        hs.hRecoInclusiveJetEtaCMShifted[iShift]->Write();
    }
    hs.hRecoInclusiveJetPtEtaStdBins->Write();
    if (isMc) {
        hs.hRecoInclusiveJetJESPtEtaStdBins->Write();
        if (doJerSyst) {
            hs.hRecoInclusiveJetJESDefNoSFPtEtaStdBins->Write();
            hs.hRecoInclusiveJetJESDefPtEtaStdBins->Write();
            hs.hRecoInclusiveJetJESUpPtEtaStdBins->Write();
            hs.hRecoInclusiveJetJESDownPtEtaStdBins->Write();
        }
    }
    hs.hRecoInclusiveJetPtPtHat->Write();
    hs.hRecoInclusiveJetEtaPt80To120->Write();
    hs.hRecoInclusiveJetEtaPtMatched80To120->Write();
    hs.hRecoInclusiveJetEtaPtDefNoSF80To120->Write();
    hs.hRecoInclusiveJetEtaPtDef80To120->Write();
    hs.hRecoInclusiveJetEtaPtDefExtraNoSF80To120->Write();
    hs.hRecoInclusiveJetEtaPtDefExtra80To120->Write();

    // Reco dijets
    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        hs.hRecoDijetPtEtaCMArr[iCut]->Write();
        hs.hRecoDijetPtEtaForwardArr[iCut]->Write();
        hs.hRecoDijetPtEtaBackwardArr[iCut]->Write();
        if (isMc && doJerSyst) {
            if (hs.hRefDijetPtEtaForwardArrJerDef[iCut]) hs.hRefDijetPtEtaForwardArrJerDef[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJerDef[iCut]) hs.hRefDijetPtEtaBackwardArrJerDef[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerUp[iCut]) hs.hRecoDijetPtEtaForwardArrJerUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerUp[iCut]) hs.hRecoDijetPtEtaBackwardArrJerUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJerDown[iCut]) hs.hRecoDijetPtEtaForwardArrJerDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJerDown[iCut]) hs.hRecoDijetPtEtaBackwardArrJerDown[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJerUp[iCut]) hs.hRefDijetPtEtaForwardArrJerUp[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJerUp[iCut]) hs.hRefDijetPtEtaBackwardArrJerUp[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJerDown[iCut]) hs.hRefDijetPtEtaForwardArrJerDown[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJerDown[iCut]) hs.hRefDijetPtEtaBackwardArrJerDown[iCut]->Write();
        }
        if (doJeuSyst && !isMc) {
            if (hs.hRecoDijetPtEtaForwardArrJeuUp[iCut]) hs.hRecoDijetPtEtaForwardArrJeuUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJeuUp[iCut]) hs.hRecoDijetPtEtaBackwardArrJeuUp[iCut]->Write();
            if (hs.hRecoDijetPtEtaForwardArrJeuDown[iCut]) hs.hRecoDijetPtEtaForwardArrJeuDown[iCut]->Write();
            if (hs.hRecoDijetPtEtaBackwardArrJeuDown[iCut]) hs.hRecoDijetPtEtaBackwardArrJeuDown[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJeuUp[iCut]) hs.hRefDijetPtEtaForwardArrJeuUp[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJeuUp[iCut]) hs.hRefDijetPtEtaBackwardArrJeuUp[iCut]->Write();
            if (hs.hRefDijetPtEtaForwardArrJeuDown[iCut]) hs.hRefDijetPtEtaForwardArrJeuDown[iCut]->Write();
            if (hs.hRefDijetPtEtaBackwardArrJeuDown[iCut]) hs.hRefDijetPtEtaBackwardArrJeuDown[iCut]->Write();
        }
    }

    fOut->Close();
    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}

//________________
void processEvents(const bool &isPbGoing, const bool &isMc, TChain &mainTree, 
                   JetCorrector &jec, JetUncertainty &jeu,
                   Histograms &hs, const int &ptHatSample,
                   const bool &doJeuSyst, const bool &doJerSyst,
                   const int &triggerId, const int &jetSelectionMethod) {

    std::cout << "Start event processing..." << std::endl;

    // Vz weight function for MC, derived from data/MC ratio of vz distribution in minimum bias events
    auto fVzWeight = std::make_unique<TF1>("fVzWeight", "pol8", -15.1, 15.1);
    fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);

    // JER smearing function (parametrization of sigma_smear vs pT)
    auto fJERSmearFunc = std::make_unique<TF1>("fJERSmearFunc", "sqrt(max(0., [0] + [1]/x))", 30., 800.);
    if (isPbGoing) {
        fJERSmearFunc->SetParameters(0.0018, 0.9352); // -1.6 < eta < 1.6
        // fJERSmearFunc->SetParameters(0.00183, 0.75500); // -0.8 < eta < 0.8
    }
    else {
        fJERSmearFunc->SetParameters(0.00176, 0.76438); // -0.8 < eta < 0.8
    }

    auto fRndm = std::make_unique<TRandom3>(0);

    int nEventsProcessed{0};
    int nGoodEvents{0};
    bool isGenOverweight{false};
    bool isRecoOverweight{false};
    double weight{1.};

    // Vectors to hold gen and reco jets collections in the current event
    std::vector<GenJet> genJets;
    std::vector<RecoJet> recoJets;

    int printEvery{50000};
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

        // Load Monte Carlo (MC) truth information if processing MC
        if (isMc) {
            loadGenJets( genJets );
        }
        // Load reconstructed jets
        loadRecoJets( recoJets, jec );

        // std::cout << "\n========================================\n";
        // std::cout << Form("Event: %d ptHat = %.1f GeV, vz = %.2f cm", iEntry, ptHat, vz) << std::endl;
        // std::cout << Form("%.1f <= ptHat <= %.1f ? %s", ptHatRange[0], ptHatRange[1], 
        //                  (ptHat >= ptHatRange[0] && ptHat <= ptHatRange[1]) ? "true" : "false") << std::endl;
        // std::cout << Form("|vz| <= 15 cm ? %s", (std::abs(vz) <= 15.) ? "true" : "false") << std::endl;

        // Check the event is satisfies basic event selection

        if (!isGoodEvent(isPbGoing, isMc, ptHatSample, triggerId)) continue;

        if ( isMc ) {
            if (isPbGoing) vz = -vz; 
            weight = eventWeight(ptHat, vz, *fVzWeight, nEntries);
            fillOverweightHistograms(genJets, recoJets, weight, hs, isGenOverweight, isRecoOverweight);
            if (isGenOverweight) continue;
            if (isRecoOverweight) continue;
        } // if ( isMc )


        nGoodEvents++;

        hs.hPtHatUnweighted->Fill(ptHat);
        hs.hPtHat->Fill(ptHat, weight);
        hs.hVzUnweighted->Fill(vz);
        hs.hVz->Fill(vz, weight);

        //
        // Gen level processing
        //
        if (isMc) {
            // processGenJets(isPbGoing, isMc, weight, genJets, hs, ptHat);
            // processGenDijets(isPbGoing, isMc, weight, genJets, hs);
        }

        //
        // Reco level processing
        //
        processRecoJets(isPbGoing, isMc, weight, recoJets, hs, jec, jeu, jetSelectionMethod, 
                        ptHat, fJERSmearFunc.get(), fRndm.get());
        // processRecoDijets(isPbGoing, isMc, weight, recoJets, hs, doJerSyst, doJeuSyst);
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
    jeuSyst (default 0): 0 - default/no systematics, 1 - use systematics (up and down)
    jerSyst (default 0): 0 - do not use systematics, 1 - use systematcis (default, up and down)
    triggerId (default 0): 0 - no trigger (or MB), 1 - jet60, 2 - jet80, 3 - jet100
    recoJetSelMethod (default 2): 0 - no selection, 1 - trkMaxPt/RawPt, 2 - jetIdTightLeptVeto, 3 - jetIdLoose
*/
#if defined(__CINT__) || defined(__CLING__)
void processForestSimple( const char* input = "", const char* output = "", 
                          int mcType = 2, int isPbGoingDir = 1, int ptHatSample = 50, bool jeuSyst = false,
                          bool jerSyst = true, int triggerId = 0, int recoJetSelMethod = 2 ) {

    const char *path2auxFiles = "./aux_files/pPb_8160";
#else
int main(int argc, char* argv[]) {

    const char* input = "";
    const char* output = "";
    int mcType = {2};
    int isPbGoingDir = {1};
    int ptHatSample = {50};
    bool jeuSyst = {false};
    bool jerSyst = {false};
    int triggerId = {0};
    int recoJetSelMethod = {2};

    if (argc <= 1 || argc != 10) {
        std::cerr << RED << "No input parameters provided. Terminating." << RESET << std::endl;
        usage();
    }

    input = argv[1];
    output = argv[2];
    mcType = std::atoi(argv[3]); 
    isPbGoingDir = std::atoi(argv[4]);
    ptHatSample = std::atoi(argv[5]);
    int tmpJeuSyst = std::atoi(argv[6]);
    if (tmpJeuSyst == 0) {
        jeuSyst = false;
    }
    else if (tmpJeuSyst == 1) {
        jeuSyst = true;
    }
    else {
        std::cerr << RED << "Invalid jeuSyst parameter. Must be 0 (default/no systematics) or 1 (use systematics). Terminating." << RESET << std::endl;
        usage();
    }
    int tmpJerSyst = std::atoi(argv[7]);
    if (tmpJerSyst == 0) {
        jerSyst = false;
    }
    else if (tmpJerSyst == 1) {
        jerSyst = true;
    }
    else {
        std::cerr << RED << "Invalid jerSyst parameter. Must be 0 (do not use systematics) or 1 (use systematics). Terminating." << RESET << std::endl;
        usage();
    }
    triggerId = std::atoi(argv[8]);
    recoJetSelMethod = std::atoi(argv[9]);

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
    std::cout << "  JEU sys: " << GREEN << (jeuSyst?"true":"false") << RESET << std::endl;
    std::cout << "  JER sys: " << GREEN << (jerSyst?"true":"false") << RESET << std::endl;
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
    createHistograms(hs, isMc, jerSyst, jeuSyst);
    if (isMc) {
        setValuesForExtraJerSmearing(isPythia, isPbGoing);
    }

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
    processEvents(isPbGoing, isMc, *mainTree, *jec, *jeu, hs, ptHatSample, 
                  jeuSyst, jerSyst, triggerId, recoJetSelMethod);

    // Write output
    writeOutput(oFileName, hs, isMc, jerSyst, jeuSyst);

}
