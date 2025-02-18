/**
 * @file ForestAODReader.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief CMS ForestAOD reader
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "ForestAODReader.h"

// ROOT headers
#include "TFile.h"

// C++ headers
#include <cstddef>
#include <cstring>
#include <fstream>

//_________________
ForestAODReader::ForestAODReader() : fEvent{nullptr}, fInFileName{nullptr}, fEvents2Read{0}, fEventsProcessed{0},
    fIsMc{false}, fCorrectCentMC{false}, fUseHltBranch{kTRUE}, fUseSkimmingBranch{kTRUE}, 
    fUseRecoJetBranch{kTRUE}, 
    fUseTrackBranch{false}, fUseGenTrackBranch{false},
    fHltTree{nullptr}, fSkimTree{nullptr}, fEventTree{nullptr},
    fRecoJetTree{nullptr}, fTrkTree{nullptr}, fGenTrkTree{nullptr},
    fRecoJetTreeName{"akCs4PFJetAnalyzer"},
    fJEC{nullptr}, fJECFiles{}, fJECPath{}, fJEU{nullptr}, fJEUInputFileName{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{false}, 
    fFixJetArrays{false}, fEventCut{nullptr}, fJetCut{nullptr},
    fRecoJet2GenJetId{}, fGenJet2RecoJet{}, 
    fUseExtraJECforAk4Cs{false}, fJECScaleCorr{nullptr}, fUseJEU{0},
    fUseJERSystematics{0}, fAlphaJER{0.0415552}, fBetaJER{0.960013},
    fJERSmearFunc{nullptr}, fRndm{nullptr},
    fVerbose{false} {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::ForestAODReader()" << std::endl;
    }
    // Initialize many variables
    fRndm = new TRandom3(0);
    clearVariables();
    setJERSystParams();
}

//_________________
ForestAODReader::ForestAODReader(const char* inputStream,  
                                 const bool& useHltBranch, const bool& useSkimmingBranch, 
                                 const bool& useRecoJetBranch, 
                                 const bool& useTrackBranch, const bool& useGenTrackBranch, 
                                 const bool& isMc) : 
    fEvent{nullptr}, fInFileName{inputStream}, fEvents2Read{0}, 
    fEventsProcessed{0}, fIsMc{isMc}, fCorrectCentMC{false}, 
    fUseManualJEC{false}, fIsPbGoing{true},
    fUseHltBranch{useHltBranch}, fUseSkimmingBranch{useSkimmingBranch}, 
    fUseRecoJetBranch{useRecoJetBranch}, 
    fUseTrackBranch{useTrackBranch}, fUseGenTrackBranch{useGenTrackBranch},
    fRecoJetTreeName{"akCs4PFJetAnalyzer"},
    fJEC{nullptr}, fJECFiles{}, fJECPath{}, fJEU{nullptr}, fJEUInputFileName{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{false}, 
    fFixJetArrays{false}, fEventCut{nullptr}, fJetCut{nullptr},
    fJECScaleCorr{nullptr}, fUseJEU{0}, fUseJERSystematics{0}, 
    fAlphaJER{0.0415552}, fBetaJER{0.960013}, fJERSmearFunc{nullptr}, 
    fVerbose{false} {
    // Initialize many variables
    fRndm = new TRandom3(0);
    clearVariables();

    if ( fVerbose ) {
        std::cout << "ForestAODReader::ForestAODReader()" << std::endl;
    }
    setJERSystParams();
}

//_________________
ForestAODReader::~ForestAODReader() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::~ForestAODReader()";
    }
    if (fEvent) delete fEvent;
    if (fHltTree) delete fHltTree;
    if (fSkimTree) delete fSkimTree;
    if (fEventTree) delete fEventTree;
    if (fRecoJetTree) delete fRecoJetTree;
    if (fTrkTree) delete fTrkTree;
    if (fGenTrkTree) delete fGenTrkTree;
    if (fJEC) delete fJEC;
    if (fJEU) delete fJEU;
    if (fEventCut) delete fEventCut;
    if (fJetCut) delete fJetCut;
    if (fJECScaleCorr) delete fJECScaleCorr;
    if (fJERSmearFunc) delete fJERSmearFunc;
}

//________________
void ForestAODReader::setJERSystParams() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::setJERSystParams";
    }
    fJerEtaLow.clear(); fJerEtaHi.clear(); fJerDef.clear(); fJerLow.clear(); fJerHi.clear();
    fJerEtaLow.push_back(-5.191); fJerEtaHi.push_back(-3.139); fJerDef.push_back(1.1922); fJerLow.push_back(1.0434); fJerHi.push_back(1.341);
    fJerEtaLow.push_back(-3.139); fJerEtaHi.push_back(-2.964); fJerDef.push_back(1.1869); fJerLow.push_back(1.0626); fJerHi.push_back(1.3112);
    fJerEtaLow.push_back(-2.964); fJerEtaHi.push_back(-2.853); fJerDef.push_back(1.7788); fJerLow.push_back(1.578); fJerHi.push_back(1.9796);
    fJerEtaLow.push_back(-2.853); fJerEtaHi.push_back(-2.500); fJerDef.push_back(1.3418); fJerLow.push_back(1.1327); fJerHi.push_back(1.5509);
    fJerEtaLow.push_back(-2.500); fJerEtaHi.push_back(-2.322); fJerDef.push_back(1.2963); fJerLow.push_back(1.0592); fJerHi.push_back(1.5334);
    fJerEtaLow.push_back(-2.322); fJerEtaHi.push_back(-2.043); fJerDef.push_back(1.1512); fJerLow.push_back(1.0372); fJerHi.push_back(1.2652);
    fJerEtaLow.push_back(-2.043); fJerEtaHi.push_back(-1.930); fJerDef.push_back(1.1426); fJerLow.push_back(1.0212); fJerHi.push_back(1.264);
    fJerEtaLow.push_back(-1.930); fJerEtaHi.push_back(-1.740); fJerDef.push_back(1.1000); fJerLow.push_back(0.9921); fJerHi.push_back(1.2079);
    fJerEtaLow.push_back(-1.740); fJerEtaHi.push_back(-1.305); fJerDef.push_back(1.1278); fJerLow.push_back(1.0292); fJerHi.push_back(1.2264);
    fJerEtaLow.push_back(-1.305); fJerEtaHi.push_back(-1.131); fJerDef.push_back(1.1609); fJerLow.push_back(1.0584); fJerHi.push_back(1.2634);
    fJerEtaLow.push_back(-1.131); fJerEtaHi.push_back(-0.783); fJerDef.push_back(1.1464); fJerLow.push_back(1.0832); fJerHi.push_back(1.2096);
    fJerEtaLow.push_back(-0.783); fJerEtaHi.push_back(-0.522); fJerDef.push_back(1.1948); fJerLow.push_back(1.1296); fJerHi.push_back(1.26);
    fJerEtaLow.push_back(-0.522); fJerEtaHi.push_back(0.000); fJerDef.push_back(1.15958); fJerLow.push_back(1.095); fJerHi.push_back(1.224);
    fJerEtaLow.push_back(0.000); fJerEtaHi.push_back(0.522); fJerDef.push_back(1.15958); fJerLow.push_back(1.095); fJerHi.push_back(1.224);
    fJerEtaLow.push_back(0.522); fJerEtaHi.push_back(0.783); fJerDef.push_back(1.1948); fJerLow.push_back(1.1296); fJerHi.push_back(1.26);
    fJerEtaLow.push_back(0.783); fJerEtaHi.push_back(1.131); fJerDef.push_back(1.1464); fJerLow.push_back(1.0832); fJerHi.push_back(1.2096);
    fJerEtaLow.push_back(1.131); fJerEtaHi.push_back(1.305); fJerDef.push_back(1.1609); fJerLow.push_back(1.0584); fJerHi.push_back(1.2634);
    fJerEtaLow.push_back(1.305); fJerEtaHi.push_back(1.740); fJerDef.push_back(1.1278); fJerLow.push_back(1.0292); fJerHi.push_back(1.2264);
    fJerEtaLow.push_back(1.740); fJerEtaHi.push_back(1.930); fJerDef.push_back(1.1000); fJerLow.push_back(0.9921); fJerHi.push_back(1.2079);
    fJerEtaLow.push_back(1.930); fJerEtaHi.push_back(2.043); fJerDef.push_back(1.1426); fJerLow.push_back(1.0212); fJerHi.push_back(1.264);
    fJerEtaLow.push_back(2.043); fJerEtaHi.push_back(2.322); fJerDef.push_back(1.1512); fJerLow.push_back(1.0372); fJerHi.push_back(1.2652);
    fJerEtaLow.push_back(2.322); fJerEtaHi.push_back(2.500); fJerDef.push_back(1.2963); fJerLow.push_back(1.0592); fJerHi.push_back(1.5334);
    fJerEtaLow.push_back(2.500); fJerEtaHi.push_back(2.853); fJerDef.push_back(1.3418); fJerLow.push_back(1.1327); fJerHi.push_back(1.5509);
    fJerEtaLow.push_back(2.853); fJerEtaHi.push_back(2.964); fJerDef.push_back(1.7788); fJerLow.push_back(1.578); fJerHi.push_back(1.9796);
    fJerEtaLow.push_back(2.964); fJerEtaHi.push_back(3.139); fJerDef.push_back(1.1869); fJerLow.push_back(1.0626); fJerHi.push_back(1.3112);
    fJerEtaLow.push_back(3.139); fJerEtaHi.push_back(5.191); fJerDef.push_back(1.1922); fJerLow.push_back(1.0434); fJerHi.push_back(1.341);

    fJERSmearFunc = new TF1("fJERSmearFunc","sqrt( [0] * [0] + [1] * [1] / x )", 30., 800.);
    fJERSmearFunc->SetParameter(0, fAlphaJER);
    fJERSmearFunc->SetParameter(1, fBetaJER);

    if ( fVerbose ) {
        std::cout << "\t[DONE]\n";
    }
}

//________________
double ForestAODReader::retrieveResolutionFactor(const double& eta) {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::retrieveResolutionFactor\n";
    }

    double val{1.};
    double res{0.};

    // Search for the bin index
    for (int i{0}; i<fJerEtaLow.size(); i++) {
        if ( eta>=fJerEtaLow.at(i) && eta<fJerEtaHi.at(i) ) {
            if ( fUseJERSystematics == -1 ) {
                val = fJerLow.at(i);
            }
            else if ( fUseJERSystematics == 0 ) {
                val = fJerDef.at(i);
            }
            else if ( fUseJERSystematics == 1 ) {
                val = fJerHi.at(i);
            }
            else {
                val = {1.};
            }
            break;
        }
    } // for (int i{0}; i<fJerEtaLow.size(); i++)

    res = TMath::Sqrt( TMath::Max(val * val - 1., 0.) );

    if ( fVerbose ) {
        std::cout << "eta: " << eta << "JER val: " << val << " Resolution factor: " << res << std::endl;
        std::cout << "\t[DONE]\n";
    }
    return res;
}

//________________
double ForestAODReader::extraJERCorr(const double &ptCorr, const double &eta) {
    // This factor applied to the reco jet pT after JEC to match data.
    // By default JEC is not fully cover data/MC JEC difference.
    // Extra correction should be applied to MC only.
    if ( fVerbose ) {
        std::cout << "ForestAODReader::extraJERCorr\n";
    }

    double res = retrieveResolutionFactor(eta);
    double sigmaSmear{0.};
    if ( ptCorr <= 30.) {
        sigmaSmear = res * fJERSmearFunc->Eval( 31. );
    }
    else if ( ptCorr >= 800 ) {
        sigmaSmear = res * fJERSmearFunc->Eval( 799. );
    }
    else {
        sigmaSmear = res * fJERSmearFunc->Eval( ptCorr );
    }

    double extraCorr = fRndm->Gaus( 1., sigmaSmear );

    if ( fVerbose ) {
        std::cout << "Resolution factor: " << res << " sigma: " << sigmaSmear 
                  << " correction factor: " << extraCorr << std::endl;
        std::cout << "\t[DONE]\n";
    }
    return extraCorr;
}

//________________
int ForestAODReader::findBinIndex(const double &val, double *array, int nBins) {
    int binIndex = -1;
    for (int i{0}; i<nBins; i++) {
        if ( val >= array[i] && val < array[i+1] ) {
            binIndex = i;
            break;
        }
    }
    return binIndex;
}

//________________
double ForestAODReader::jecManualCorrection(const double &pt, const double &eta) {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::jecManualCorrection - begin\n";
    }

    // Array of eta bin edges (left edge included, right edge excluded)
    double jetEtaL2L3StdVals[] = { 
        -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, 
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
         4.363,  4.538,  4.716,  4.889,  5.191 
    };
    int jetNEtaL2L3StdBins = sizeof(jetEtaL2L3StdVals)/sizeof(double)-1;

    // Find eta bin
    int iEta = findBinIndex(eta, jetEtaL2L3StdVals, jetNEtaL2L3StdBins);

    // Array of pt bin edges (left edge included, right edge excluded)
    double jetPtVals[] = {
        10, 15, 20, 25, 35, 45, 55, 65, 75, 85,
        95, 110, 130, 150, 180, 210, 250, 300, 360, 420,
        500, 600, 720, 860, 1100, 6500
    };
    int jetNPtBins = sizeof(jetPtVals) / sizeof(double) - 1;

    // Find pt bin
    int iPt = findBinIndex(pt, jetPtVals, jetNPtBins);

    double retVal = {0.};
    if ( (fIsMc && fIsPbGoing) || (!fIsMc && !fIsPbGoing) ) {

        // [eta][pt] binning
        double corrFactor[82][25] = {
            {0.57482, 0.582115, 0.59952, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242, 0.594242},
            {0.703931, 0.70757, 0.73298, 0.759685, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971, 0.89971},
            {0.811999, 0.811747, 0.81234, 0.833431, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214, 0.909214},
            {0.871698, 0.868552, 0.865951, 0.891889, 0.896865, 0.869295, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849, 0.85849},
            {0.906874, 0.904491, 0.901818, 0.918155, 0.908792, 0.913063, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384, 0.925384},
            {0.914368, 0.905878, 0.903684, 0.918102, 0.920212, 0.938205, 0.933365, 0.922843, 0.930435, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95},
            {0.904319, 0.898593, 0.903269, 0.918633, 0.939525, 0.939721, 0.946444, 0.959347, 0.958318, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297, 0.943297},
            {0.886983, 0.887344, 0.891294, 0.90894, 0.930444, 0.934165, 0.934344, 0.935785, 0.95786, 0.972574, 0.982899, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257, 1.04257},
            {0.858551, 0.866679, 0.876018, 0.888182, 0.901724, 0.914908, 0.92877, 0.938289, 0.953861, 0.965993, 0.977378, 1.02978, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336, 0.982336},
            {0.820293, 0.833248, 0.849359, 0.86332, 0.883825, 0.902238, 0.924544, 0.940937, 0.954956, 0.966211, 0.973933, 0.996924, 0.994166, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497, 0.946497},
            {0.759024, 0.772135, 0.792746, 0.810997, 0.838511, 0.851707, 0.860012, 0.86899, 0.867794, 0.882711, 0.888496, 0.882199, 0.877673, 0.890845, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115, 0.806115},
            {0.734224, 0.752294, 0.777953, 0.795172, 0.812652, 0.816657, 0.824728, 0.827757, 0.829353, 0.843665, 0.853913, 0.865692, 0.873347, 0.903245, 0.91371, 0.894636, 0.894636, 0.894636, 0.894636, 0.894636, 0.894636, 0.894636, 0.894636, 0.894636, 0.894636},
            {0.748996, 0.789777, 0.814819, 0.838483, 0.856947, 0.864053, 0.870441, 0.871457, 0.879189, 0.884046, 0.896514, 0.908253, 0.917991, 0.937306, 0.947071, 0.97646, 0.999385, 0.999385, 0.999385, 0.999385, 0.999385, 0.999385, 0.999385, 0.999385, 0.999385},
            {0.744229, 0.787811, 0.812309, 0.836799, 0.868536, 0.880613, 0.883306, 0.884194, 0.888097, 0.890519, 0.897243, 0.907659, 0.917343, 0.929254, 0.941875, 0.94456, 0.952827, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95},
            {0.767884, 0.794637, 0.815549, 0.834948, 0.863727, 0.874473, 0.87766, 0.878835, 0.883159, 0.884608, 0.887539, 0.900561, 0.90868, 0.918621, 0.931252, 0.940644, 0.946885, 0.943503, 0.975997, 0.975997, 0.975997, 0.975997, 0.975997, 0.975997, 0.975997},
            {0.813609, 0.825887, 0.842149, 0.857991, 0.873684, 0.879108, 0.882586, 0.883468, 0.885675, 0.88802, 0.890639, 0.90019, 0.905306, 0.914549, 0.923712, 0.932465, 0.939781, 0.946591, 0.946303, 0.939709, 0.939709, 0.939709, 0.939709, 0.939709, 0.939709},
            {0.844793, 0.852509, 0.866578, 0.87486, 0.881496, 0.889668, 0.89063, 0.891436, 0.89118, 0.892192, 0.893757, 0.897367, 0.90359, 0.908406, 0.919135, 0.922543, 0.92965, 0.938575, 0.94522, 0.943902, 0.95, 0.95, 0.95, 0.95, 0.95},
            {0.853491, 0.862398, 0.874441, 0.876187, 0.87933, 0.886506, 0.890389, 0.89253, 0.891753, 0.890737, 0.891319, 0.894437, 0.899539, 0.901012, 0.907424, 0.917592, 0.921025, 0.929369, 0.937366, 0.944268, 0.94073, 0.95, 0.95, 0.95, 0.95},
            {0.856132, 0.860919, 0.869023, 0.876172, 0.876997, 0.883097, 0.886667, 0.884827, 0.888161, 0.885008, 0.884785, 0.889836, 0.89114, 0.891741, 0.896316, 0.90818, 0.91137, 0.918411, 0.924729, 0.932909, 0.941706, 0.93501, 0.93501, 0.93501, 0.93501},
            {0.858556, 0.856003, 0.863769, 0.876096, 0.869381, 0.88155, 0.885117, 0.882606, 0.885985, 0.884136, 0.882465, 0.885306, 0.885782, 0.887055, 0.891621, 0.897159, 0.902556, 0.910226, 0.911705, 0.921682, 0.936152, 0.906843, 0.95, 0.95, 0.95},
            {0.860712, 0.860642, 0.869634, 0.8756, 0.871006, 0.887041, 0.892709, 0.889888, 0.891664, 0.892173, 0.89161, 0.890436, 0.891836, 0.890279, 0.8927, 0.897091, 0.903371, 0.909968, 0.909997, 0.924544, 0.931199, 0.928998, 0.905673, 0.905673, 0.905673},
            {0.865682, 0.872019, 0.879477, 0.882553, 0.888354, 0.900963, 0.909071, 0.900627, 0.904993, 0.906072, 0.906937, 0.903275, 0.904221, 0.90191, 0.90276, 0.906886, 0.910983, 0.91669, 0.917798, 0.928444, 0.932879, 0.932273, 0.923942, 0.923942, 0.923942},
            {0.869049, 0.876698, 0.884423, 0.886074, 0.893357, 0.892461, 0.906503, 0.900938, 0.903626, 0.902597, 0.900906, 0.899551, 0.89941, 0.896231, 0.895977, 0.896572, 0.899126, 0.904265, 0.908788, 0.9134, 0.918533, 0.919211, 0.925982, 0.888711, 0.888711},
            {0.860163, 0.870272, 0.885658, 0.8785, 0.885649, 0.885831, 0.891503, 0.89315, 0.892911, 0.888299, 0.886594, 0.886638, 0.883719, 0.878853, 0.879257, 0.880705, 0.877577, 0.882467, 0.885171, 0.888587, 0.890682, 0.896046, 0.902811, 0.88801, 0.88801},
            {0.856347, 0.864583, 0.879364, 0.874246, 0.879624, 0.891516, 0.890058, 0.893464, 0.892401, 0.8896, 0.888561, 0.885526, 0.88115, 0.880614, 0.879652, 0.884933, 0.879916, 0.884567, 0.887319, 0.887951, 0.890288, 0.893453, 0.907409, 0.884799, 0.884799},
            {0.864263, 0.867942, 0.871548, 0.875808, 0.880086, 0.901296, 0.893638, 0.895728, 0.894771, 0.896165, 0.893718, 0.890841, 0.886555, 0.888713, 0.888168, 0.889292, 0.886367, 0.890458, 0.893518, 0.894739, 0.895119, 0.898422, 0.908131, 0.901982, 0.901982},
            {0.87185, 0.873684, 0.879729, 0.885391, 0.887567, 0.909136, 0.902472, 0.899966, 0.901059, 0.900462, 0.899847, 0.89859, 0.895116, 0.896433, 0.896511, 0.893167, 0.889912, 0.894997, 0.896781, 0.898924, 0.899546, 0.905104, 0.910844, 0.910732, 0.910732},
            {0.884165, 0.887649, 0.897122, 0.902822, 0.907577, 0.917078, 0.915001, 0.913183, 0.91498, 0.911904, 0.911208, 0.908993, 0.907069, 0.90629, 0.905734, 0.903386, 0.901116, 0.905529, 0.906481, 0.907069, 0.908161, 0.912039, 0.914755, 0.915972, 0.95},
            {0.899685, 0.904422, 0.908476, 0.91884, 0.927859, 0.930524, 0.926789, 0.927416, 0.929553, 0.926641, 0.921703, 0.921895, 0.919651, 0.917357, 0.917048, 0.915268, 0.914024, 0.914531, 0.915744, 0.91526, 0.916558, 0.919475, 0.923286, 0.923154, 0.932135},
            {0.906829, 0.91344, 0.919385, 0.931278, 0.935478, 0.938793, 0.936216, 0.936044, 0.937381, 0.937587, 0.933076, 0.933457, 0.928083, 0.927827, 0.924515, 0.921596, 0.921794, 0.919312, 0.922051, 0.920202, 0.922647, 0.922305, 0.926632, 0.930074, 0.94197},
            {0.917465, 0.926967, 0.928176, 0.932138, 0.941622, 0.950893, 0.947266, 0.944575, 0.942143, 0.944411, 0.943228, 0.942044, 0.936596, 0.937411, 0.93112, 0.927466, 0.92866, 0.925284, 0.926489, 0.92616, 0.927391, 0.927326, 0.930818, 0.935448, 0.958197},
            {0.929892, 0.93385, 0.930345, 0.933143, 0.943582, 0.959411, 0.955181, 0.952755, 0.95032, 0.950357, 0.946829, 0.94636, 0.942348, 0.941378, 0.938615, 0.93673, 0.93466, 0.932391, 0.932416, 0.931965, 0.933158, 0.935215, 0.936483, 0.939099, 0.949301},
            {0.93275, 0.935751, 0.937092, 0.945199, 0.95048, 0.962747, 0.957389, 0.95928, 0.957049, 0.956174, 0.952769, 0.950625, 0.948691, 0.9431, 0.942425, 0.941271, 0.938097, 0.936923, 0.937065, 0.93654, 0.937992, 0.938773, 0.941076, 0.942366, 0.948131},
            {0.931234, 0.938765, 0.939203, 0.94749, 0.951945, 0.957084, 0.958062, 0.959414, 0.955565, 0.957024, 0.955492, 0.952362, 0.95084, 0.945241, 0.944363, 0.940305, 0.940834, 0.938331, 0.938811, 0.939013, 0.938787, 0.939947, 0.943059, 0.943231, 0.949202},
            {0.929919, 0.93362, 0.934625, 0.940154, 0.945615, 0.951103, 0.954179, 0.954548, 0.952409, 0.957051, 0.952736, 0.949979, 0.948758, 0.946523, 0.944937, 0.941115, 0.940778, 0.938768, 0.939372, 0.937839, 0.938367, 0.939353, 0.941673, 0.942235, 0.948839},
            {0.926425, 0.930257, 0.933318, 0.939296, 0.949651, 0.958204, 0.951933, 0.950001, 0.950602, 0.95351, 0.951685, 0.948507, 0.947626, 0.946805, 0.943121, 0.942736, 0.940358, 0.938827, 0.937841, 0.937245, 0.938746, 0.938624, 0.940523, 0.941694, 0.948813},
            {0.923086, 0.931436, 0.94064, 0.946594, 0.957596, 0.961701, 0.949957, 0.951761, 0.952167, 0.949609, 0.950942, 0.948178, 0.948331, 0.948692, 0.943726, 0.94357, 0.940913, 0.939389, 0.938881, 0.939339, 0.938875, 0.939915, 0.941791, 0.943617, 0.948834},
            {0.925373, 0.933405, 0.946619, 0.94957, 0.956479, 0.958508, 0.951082, 0.954553, 0.95279, 0.951167, 0.951971, 0.950229, 0.950766, 0.949465, 0.944941, 0.944136, 0.941089, 0.941049, 0.940735, 0.939532, 0.938692, 0.940554, 0.942358, 0.944399, 0.947306},
            {0.925175, 0.932754, 0.943141, 0.946366, 0.947584, 0.958889, 0.955454, 0.955048, 0.950661, 0.951793, 0.9523, 0.949666, 0.949536, 0.946691, 0.945127, 0.943941, 0.942175, 0.940709, 0.938856, 0.938097, 0.937857, 0.939385, 0.941322, 0.942675, 0.946575},
            {0.921766, 0.933286, 0.941517, 0.947031, 0.946804, 0.956855, 0.953394, 0.951857, 0.950639, 0.950768, 0.951708, 0.947906, 0.947265, 0.946654, 0.945407, 0.942874, 0.941262, 0.93903, 0.937375, 0.937957, 0.937953, 0.93895, 0.941454, 0.942782, 0.945528},
            {0.921546, 0.936538, 0.943111, 0.948462, 0.950126, 0.953943, 0.952004, 0.94972, 0.951843, 0.947759, 0.94949, 0.947816, 0.946263, 0.946425, 0.943087, 0.94204, 0.939867, 0.937608, 0.93423, 0.936919, 0.937493, 0.937757, 0.940234, 0.941594, 0.945929},
            {0.923778, 0.933629, 0.940721, 0.945476, 0.94788, 0.949386, 0.949861, 0.952269, 0.951428, 0.946301, 0.947422, 0.944936, 0.94492, 0.944014, 0.940224, 0.940405, 0.938251, 0.935395, 0.932192, 0.93423, 0.934254, 0.934611, 0.935801, 0.937327, 0.944298},
            {0.926523, 0.932347, 0.936109, 0.939456, 0.94283, 0.944124, 0.947854, 0.951199, 0.950396, 0.949339, 0.948712, 0.944816, 0.944051, 0.943789, 0.940574, 0.937463, 0.937458, 0.934054, 0.933349, 0.932618, 0.93294, 0.933341, 0.934808, 0.937396, 0.940945},
            {0.926951, 0.93206, 0.934191, 0.938051, 0.939648, 0.945827, 0.949307, 0.945298, 0.951754, 0.950695, 0.948106, 0.947216, 0.945441, 0.943893, 0.942036, 0.937197, 0.937266, 0.935037, 0.934377, 0.932974, 0.932632, 0.934045, 0.935875, 0.938246, 0.94086},
            {0.925066, 0.929221, 0.93542, 0.941229, 0.94261, 0.950786, 0.949256, 0.942318, 0.950121, 0.947967, 0.946722, 0.947833, 0.947317, 0.94335, 0.94165, 0.938559, 0.936779, 0.936515, 0.933656, 0.93325, 0.932779, 0.934099, 0.934923, 0.936671, 0.939808},
            {0.922393, 0.926257, 0.937437, 0.943322, 0.948451, 0.95039, 0.948676, 0.947195, 0.947902, 0.945519, 0.944897, 0.945852, 0.945038, 0.942307, 0.93983, 0.938266, 0.93627, 0.934939, 0.931901, 0.931261, 0.932148, 0.932822, 0.93365, 0.935883, 0.938299},
            {0.923258, 0.925352, 0.933459, 0.938394, 0.94997, 0.949103, 0.952609, 0.952939, 0.948635, 0.950072, 0.948409, 0.948037, 0.94514, 0.941317, 0.94076, 0.9383, 0.936114, 0.933141, 0.931885, 0.931237, 0.93177, 0.933176, 0.934218, 0.936354, 0.939314},
            {0.927843, 0.931885, 0.932284, 0.940535, 0.952525, 0.953069, 0.954464, 0.954417, 0.951632, 0.952936, 0.951699, 0.948143, 0.946925, 0.942917, 0.942168, 0.938168, 0.936091, 0.933725, 0.933773, 0.932565, 0.932854, 0.934316, 0.935473, 0.937398, 0.941134},
            {0.927302, 0.933088, 0.93455, 0.941783, 0.953305, 0.94904, 0.952652, 0.952167, 0.951449, 0.94828, 0.94853, 0.944404, 0.943307, 0.940829, 0.937684, 0.935189, 0.932043, 0.931498, 0.930185, 0.928238, 0.928811, 0.930968, 0.930638, 0.932523, 0.937181},
            {0.921299, 0.924443, 0.929628, 0.936245, 0.946091, 0.943727, 0.946723, 0.946646, 0.947267, 0.944318, 0.945173, 0.941398, 0.938312, 0.934903, 0.931313, 0.930246, 0.926783, 0.925833, 0.923019, 0.922469, 0.922644, 0.925197, 0.924406, 0.925405, 0.929573},
            {0.912881, 0.91598, 0.924893, 0.933549, 0.939278, 0.941682, 0.938588, 0.940401, 0.944138, 0.939772, 0.939, 0.936369, 0.934478, 0.930717, 0.928166, 0.925177, 0.921385, 0.920521, 0.919088, 0.918392, 0.919292, 0.920838, 0.921172, 0.922085, 0.926359},
            {0.901594, 0.907253, 0.917307, 0.921955, 0.928131, 0.934009, 0.933296, 0.933952, 0.937425, 0.931818, 0.929389, 0.930459, 0.929998, 0.925485, 0.922731, 0.919989, 0.91714, 0.916277, 0.915188, 0.915287, 0.915682, 0.917603, 0.919772, 0.921005, 0.929527},
            {0.892243, 0.898839, 0.903858, 0.907253, 0.913877, 0.921306, 0.925899, 0.926673, 0.928693, 0.924366, 0.92147, 0.921669, 0.920618, 0.917347, 0.91322, 0.912233, 0.911774, 0.911352, 0.909972, 0.910736, 0.911134, 0.914663, 0.915622, 0.91876, 0.930501},
            {0.88168, 0.884951, 0.891669, 0.89713, 0.901483, 0.908227, 0.913951, 0.914426, 0.917569, 0.913631, 0.912022, 0.909239, 0.906609, 0.906086, 0.902123, 0.900374, 0.900749, 0.901473, 0.900868, 0.900987, 0.903212, 0.906269, 0.906566, 0.911626, 0.920865},
            {0.868464, 0.8699, 0.878551, 0.887057, 0.890815, 0.899204, 0.899746, 0.900059, 0.902395, 0.90109, 0.899652, 0.897681, 0.894468, 0.893239, 0.889246, 0.88736, 0.888835, 0.88981, 0.891327, 0.890887, 0.893158, 0.89575, 0.898119, 0.902009, 0.910647},
            {0.856672, 0.860289, 0.866061, 0.871267, 0.88084, 0.889826, 0.889982, 0.890098, 0.889854, 0.889512, 0.891467, 0.889505, 0.885347, 0.884032, 0.882851, 0.881392, 0.882398, 0.884536, 0.8862, 0.885962, 0.888317, 0.891108, 0.894014, 0.895895, 0.902491},
            {0.853296, 0.858885, 0.861521, 0.867, 0.877836, 0.885973, 0.889453, 0.88745, 0.884896, 0.883736, 0.885363, 0.883799, 0.879991, 0.878132, 0.880564, 0.876877, 0.878019, 0.877546, 0.878897, 0.880113, 0.883407, 0.885038, 0.888965, 0.890574, 0.8943},
            {0.856108, 0.86167, 0.862785, 0.873942, 0.878943, 0.883511, 0.888853, 0.889615, 0.886509, 0.886173, 0.884748, 0.883875, 0.88074, 0.876493, 0.876835, 0.874072, 0.876029, 0.877069, 0.879594, 0.88475, 0.88966, 0.892534, 0.89749, 0.8995, 0.910292},
            {0.855886, 0.862059, 0.866568, 0.878002, 0.883085, 0.884878, 0.891916, 0.896168, 0.894038, 0.895925, 0.896025, 0.894953, 0.894455, 0.89107, 0.889442, 0.889468, 0.891923, 0.897139, 0.900877, 0.906741, 0.91298, 0.918988, 0.924001, 0.929363, 0.938814},
            {0.852622, 0.858203, 0.867369, 0.873156, 0.884993, 0.898054, 0.897615, 0.898586, 0.900412, 0.903576, 0.900234, 0.901431, 0.904816, 0.900892, 0.901691, 0.902242, 0.906725, 0.912154, 0.916182, 0.920153, 0.927494, 0.93392, 0.938716, 0.947209, 0.950588},
            {0.84754, 0.849652, 0.857777, 0.861461, 0.875276, 0.889211, 0.889374, 0.890347, 0.893456, 0.89556, 0.892559, 0.893533, 0.895922, 0.891636, 0.895512, 0.899478, 0.904895, 0.912682, 0.9166, 0.922158, 0.929837, 0.935917, 0.941444, 0.947896, 0.949134},
            {0.840396, 0.841654, 0.845696, 0.853359, 0.860066, 0.871988, 0.877934, 0.881299, 0.881537, 0.880246, 0.881164, 0.880694, 0.882288, 0.881915, 0.887997, 0.89318, 0.90094, 0.909381, 0.914779, 0.921169, 0.929337, 0.935437, 0.942085, 0.945655, 0.936801},
            {0.837563, 0.839664, 0.842331, 0.852869, 0.861145, 0.871957, 0.873048, 0.874716, 0.875938, 0.873635, 0.87398, 0.87441, 0.876716, 0.879169, 0.883285, 0.888837, 0.897046, 0.907683, 0.91335, 0.92114, 0.927745, 0.933649, 0.941522, 0.943241, 0.946034},
            {0.834312, 0.840804, 0.847212, 0.856137, 0.869422, 0.869825, 0.87192, 0.871744, 0.875146, 0.870701, 0.871906, 0.876017, 0.876444, 0.878476, 0.885093, 0.89148, 0.899047, 0.909824, 0.915424, 0.922687, 0.928605, 0.933131, 0.93937, 0.940964, 0.955387},
            {0.82749, 0.836114, 0.846143, 0.854376, 0.864868, 0.866446, 0.871287, 0.872776, 0.872728, 0.86959, 0.872074, 0.87709, 0.879215, 0.881955, 0.891497, 0.898165, 0.906703, 0.913823, 0.921397, 0.926482, 0.931729, 0.936762, 0.941942, 0.940316, 0.977434},
            {0.798727, 0.811329, 0.825937, 0.840088, 0.85249, 0.861464, 0.865201, 0.867967, 0.866826, 0.870611, 0.87275, 0.879619, 0.887348, 0.891821, 0.901216, 0.910637, 0.919343, 0.926146, 0.933957, 0.93915, 0.940655, 0.945646, 0.950549, 0.958888, 0.958888},
            {0.750664, 0.776088, 0.799341, 0.824472, 0.847457, 0.858774, 0.861689, 0.863012, 0.866252, 0.870988, 0.873657, 0.883651, 0.89403, 0.900845, 0.912083, 0.923064, 0.934918, 0.942083, 0.950042, 0.956554, 0.957819, 0.960896, 0.955484, 0.993795, 0.993795},
            {0.727904, 0.767305, 0.796241, 0.826597, 0.853432, 0.862331, 0.868468, 0.871718, 0.874697, 0.877645, 0.884476, 0.894723, 0.905319, 0.916957, 0.931131, 0.942501, 0.957351, 0.964718, 0.976277, 0.980206, 0.983721, 0.985208, 0.944389, 0.944389, 0.944389},
            {0.735326, 0.770118, 0.794546, 0.821073, 0.844857, 0.84912, 0.856365, 0.859633, 0.862814, 0.86764, 0.878031, 0.888032, 0.901543, 0.916302, 0.935703, 0.950917, 0.967611, 0.978797, 0.994767, 0.99744, 1.00651, 1.03522, 1.03522, 1.03522, 1.03522},
            {0.721929, 0.742997, 0.761763, 0.778214, 0.798507, 0.808196, 0.813732, 0.820565, 0.829211, 0.837779, 0.846043, 0.85687, 0.866658, 0.880013, 0.89958, 0.918697, 0.933206, 0.948137, 0.961954, 0.967642, 0.99875, 0.99875, 0.99875, 0.99875, 0.99875},
            {0.743616, 0.763511, 0.783621, 0.795149, 0.817341, 0.84384, 0.853726, 0.867559, 0.876775, 0.887924, 0.894724, 0.906817, 0.915243, 0.925656, 0.931936, 0.944475, 0.949425, 0.954666, 0.941933, 0.920053, 1.00027, 1.00027, 1.00027, 1.00027, 1.00027},
            {0.799568, 0.814258, 0.828413, 0.843243, 0.866928, 0.89526, 0.910002, 0.925434, 0.934961, 0.947253, 0.953964, 0.968984, 0.982373, 0.997563, 1.00652, 1.01669, 1.03223, 1.05799, 1.04972, 1.06541, 1.06541, 1.06541, 1.06541, 1.06541, 1.06541},
            {0.837726, 0.845097, 0.854288, 0.863395, 0.88765, 0.901856, 0.917987, 0.929902, 0.935907, 0.944513, 0.953801, 0.966785, 0.97818, 0.989839, 1.00173, 1.01927, 1.0226, 1.06066, 1.08178, 1.08178, 1.08178, 1.08178, 1.08178, 1.08178, 1.08178},
            {0.868547, 0.869971, 0.878015, 0.885353, 0.906472, 0.912, 0.927955, 0.938178, 0.94506, 0.950694, 0.959732, 0.970447, 0.984816, 0.990019, 0.989299, 0.999745, 0.969905, 1.09008, 1.09008, 1.09008, 1.09008, 1.09008, 1.09008, 1.09008, 1.09008},
            {0.886294, 0.884283, 0.891502, 0.902626, 0.914903, 0.919929, 0.936225, 0.945307, 0.95668, 0.962859, 0.968369, 0.982199, 0.993436, 0.99793, 0.994226, 0.998818, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85},
            {0.885855, 0.8827, 0.884913, 0.896107, 0.910476, 0.916955, 0.931315, 0.937832, 0.944953, 0.948434, 0.958823, 0.969698, 0.984913, 0.981042, 1.00248, 0.996619, 0.996619, 0.996619, 0.996619, 0.996619, 0.996619, 0.996619, 0.996619, 0.996619, 0.996619},
            {0.874466, 0.86992, 0.870944, 0.884772, 0.899985, 0.909074, 0.919161, 0.928075, 0.931581, 0.931238, 0.946911, 0.940742, 0.957358, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904, 1.02904},
            {0.844132, 0.840483, 0.844242, 0.854454, 0.87354, 0.885239, 0.898579, 0.905473, 0.914208, 0.923553, 0.928632, 0.925785, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234, 0.884234},
            {0.781699, 0.781321, 0.784587, 0.795758, 0.816979, 0.822684, 0.847535, 0.863272, 0.88921, 0.904231, 0.889566, 0.913563, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85},
            {0.673517, 0.670273, 0.677134, 0.698893, 0.701921, 0.734922, 0.766415, 0.761047, 0.794727, 0.842985, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747, 0.855747},
            {0.552067, 0.555755, 0.564438, 0.581354, 0.596681, 0.631931, 0.659867, 0.597055, 0.608149, 0.843495, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85, 0.85},
            {0.449336, 0.456211, 0.473186, 0.487885, 0.530754, 0.578856, 0.578505, 0.578544, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55}
        };

        // L2L3Residual
        if ( iEta >= 0 && iEta < jetNEtaL2L3StdBins && iPt >= 0 && iPt < jetNPtBins ) {
            retVal = 1./corrFactor[iEta][iPt];
        }
    }

    if ( fVerbose ) {
        std::cout << Form("eta range: %f - %f, pT range: %f - %f \t", jetEtaL2L3StdVals[iEta], jetEtaL2L3StdVals[iEta+1], jetPtVals[iPt], jetPtVals[iPt+1]);
        std::cout << Form("iEta = %d, iPt = %d, retVal = %f\n", iEta, iPt, retVal);
        std::cout << "ForestAODReader::jecManualCorrection - end" << std::endl;
    }
    
    return retVal;
}

//_________________
void ForestAODReader::clearVariables() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::clearVariables " << std::endl;
    }
    fRunId = {0};
    fEventId = {0};
    fLumi = {0};
    fVertexZ = {-999.f};
    fHiBin = {-1};
    fPtHatWeight = {-1.f};
    fPtHat = {-1.f};

    // bad jets and multiplicity to be added

    fNRecoJets = {0};
    fNGenJets = {0};
    fNTracks = {0};

    fHLT_HIAK4CaloJet60_v1 = {0};
    fHLT_HIAK4CaloJet80_v1 = {0};
    fHLT_PAAK4CaloJet60_Eta5p1_v3 = {0};
    fHLT_PAAK4CaloJet80_Eta5p1_v3 = {0};
    fHLT_PAAK4CaloJet100_Eta5p1_v3 = {0};
    fHLT_PAAK4PFJet60_Eta5p1_v4 = {0};
    fHLT_PAAK4PFJet80_Eta5p1_v3 = {0};
    fHLT_PAAK4PFJet100_Eta5p1_v3 = {0};
    fHLT_PAAK4PFJet120_Eta5p1_v2 = {0};

    fHLT_HIAK4PFJet15_v1 = {0};
    fHLT_HIAK4PFJet15_v1_Prescl = {0};
    fHLT_HIAK4PFJet30_v1 = {0};
    fHLT_HIAK4PFJet30_v1_Prescl = {0};
    fHLT_HIAK4PFJet40_v1 = {0};
    fHLT_HIAK4PFJet40_v1_Prescl = {0};
    fHLT_HIAK4PFJet60_v1 = {0};
    fHLT_HIAK4PFJet60_v1_Prescl = {0};
    fHLT_HIAK4PFJet80_v1 = {0};
    fHLT_HIAK4PFJet80_v1_Prescl = {0};
    fHLT_HIAK4PFJet120_v1 = {0};
    fHLT_HIAK4PFJet120_v1_Prescl = {0};

    fHLT_HIAK8PFJet15_v1 = {0};
    fHLT_HIAK8PFJet15_v1_Prescl = {0};
    fHLT_HIAK8PFJet25_v1 = {0};
    fHLT_HIAK8PFJet25_v1_Prescl = {0};
    fHLT_HIAK8PFJet40_v1 = {0};
    fHLT_HIAK8PFJet40_v1_Prescl = {0};
    fHLT_HIAK8PFJet60_v1 = {0};
    fHLT_HIAK8PFJet60_v1_Prescl = {0};
    fHLT_HIAK8PFJet80_v1 = {0};
    fHLT_HIAK8PFJet80_v1_Prescl = {0};
    fHLT_HIAK8PFJet140_v1 = {0};
    fHLT_HIAK8PFJet140_v1_Prescl = {0};

    fHLT_HIPFJet25_v1 = {0};
    fHLT_HIPFJet25_v1_Prescl = {0};
    fHLT_HIPFJet140_v1 = {0};
    fHLT_HIPFJet140_v1_Prescl = {0};

    fHLT_HIPuAK4CaloJet80Eta5p1_v1 = {0};
    fHLT_HIPuAK4CaloJet100Eta5p1_v1 = {0};

    fHBHENoiseFilterResultRun2Loose = {0};
    fHBHENoiseFilterResultRun2Tight = {0};
    fHBHEIsoNoiseFilterResult = {0};
    fCollisionEventSelectionAODv2 = {0};
    fPhfCoincFilter2Th4 = {0};
    fPPAprimaryVertexFilter = {0};
    fPBeamScrapingFilter = {0};
    fPprimaryVertexFilter = {0};
    fPVertexFilterCutG = {0};
    fPVertexFilterCutGloose = {0};
    fPVertexFilterCutGtight = {0};
    fPVertexFilterCutE = {0};
    fPVertexFilterCutEandG = {0};
    fPClusterCompatibilityFilter = {0};

    fPhfCoincFilter = {0};
    fPVertexFilterCutdz1p0 = {0};
    fPVertexFilterCutGplus = {0};
    fPVertexFilterCutVtx1 = {0};

    // Loop over jets and tracks
    for (short i{0}; i<20000; i++) {

        // Jet variables
        if (i<10000) {
            fRecoJetPt[i] = {0.f};
            fRecoJetEta[i] = {0.f};
            fRecoJetPhi[i] = {0.f};
            fRecoJetWTAEta[i] = {0.f};
            fRecoJetWTAPhi[i] = {0.f};
            fRecoJetTrackMax[i] = {0.f};
            fRefJetPt[i] = {0.f};
            fRefJetEta[i] = {0.f};
            fRefJetPhi[i] = {0.f};
            fRefJetWTAEta[i] = {0.f};
            fRefJetWTAPhi[i] = {0.f};
            fRefJetPartonFlavor[i] = {-999};
            fRefJetPartonFlavorForB[i] = {-99};
            fGenJetPt[i] = {0.f};
            fGenJetEta[i] = {0.f};
            fGenJetPhi[i] = {0.f};
            fGenJetWTAEta[i] = {0.f};
            fGenJetWTAPhi[i] = {0.f};
        } // if (i<100)

        // Track variables
        fTrackPt[i] = {0.f};
        fTrackEta[i] = {0.f};
        fTrackPhi[i] = {0.f};
        fTrackPtErr[i] = {0.f};
        fTrackDcaXY[i] = {0.f};
        fTrackDcaZ[i] = {0.f};
        fTrackDcaXYErr[i] = {0.f};
        fTrackDcaZErr[i] = {0.f};
        fTrackChi2[i] = {0.f};
        fTrackNDOF[i] = {0};
        fTrackPartFlowEcal[i] = {0.f};
        fTrackPartFlowHcal[i] = {0.f};
        fTrackMVA[i] = {0.f};
        fTrackAlgo[i] = {0};
        fTrackCharge[i] = {0};
        fTrackNHits[i] = {0};
        fTrackNLayers[i] = {0};
        fTrackHighPurity[i] = {false};
    } // for (short i{0}; i<9999; i++)

    fGenTrackPt.clear();
    fGenTrackEta.clear();
    fGenTrackPhi.clear();
    fGenTrackCharge.clear();
    fGenTrackPid.clear();
    fGenTrackSube.clear();

    if (fIsMc) {
        fRecoJet2GenJetId.clear();
        fGenJet2RecoJet.clear();
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//_________________
int ForestAODReader::init() {
    if ( fVerbose ) {
        std::cout << "ForestAODReader::init()" << std::endl;
    }
    int status = 0;
    // Setup chains to read
    status = setupChains();
    // Setup branches to read
    setupBranches();
    // Setup jet energy correction files and pointer
    setupJEC();
    // Setup jet energy uncertainty files and pointer
    setupJEU();
    if ( fIsMc ) {
        if ( TMath::Abs(fUseJERSystematics)<=1 ) {
            setJERSystParams();
        }
    }
    if ( fVerbose ) {
        std::cout << "ForestAODReader::init() is finished " << std::endl;
    }
    return status;
}

//________________
void ForestAODReader::setupJEC() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupJEC()" << std::endl;
    }

    // If no path to the aux_file
    if ( fJECPath.Length() <= 0 ) {
        // Set default values
        std::cout << "[WARNING] Default path to JEC files will be used" << std::endl;
        setPath2JetAnalysis();
    }

    if ( fJECFiles.empty() ) {
        std::cout << "[WARNING] Default JEC file with parameters will be used" << std::endl;
        addJECFile();
    }

    // Need to add path in front of the file name
    std::vector< std::string > tmp;
    for (unsigned int i{0}; i<fJECFiles.size(); i++) {
        tmp.push_back( Form( "%s/aux_files/%s_%i/JEC/%s", 
                             fJECPath.Data(), fCollidingSystem.Data(),
                             fCollidingEnergyGeV, fJECFiles.at(i).c_str() ) );
    }
        
    fJECFiles.clear();
    fJECFiles = tmp;

    std::cout << "JEC files added: " << std::endl;
    for (unsigned int i{0}; i<fJECFiles.size(); i++) {
        std::cout << i << fJECFiles.at(i) << std::endl;
    }
	
	fJEC = new JetCorrector( fJECFiles );

    if ( fUseExtraJECforAk4Cs ) {
        createExtraJECScaleCorrFunction();
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//________________
void ForestAODReader::setupJEU() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupJEU()" << std::endl;
    }

    // Next part is needed only if JEU correction is applied
    if ( fUseJEU == 0) return;

    // If no path to the aux_file
    if ( fJECPath.Length() <= 0 ) {
        // Set default values
        std::cout << "[WARNING] Default path to JEU files will be used" << std::endl;
        setPath2JetAnalysis();
    }

    // If no correction file is specified
    if ( fJEUInputFileName.Length() <= 0 ) {
        std::cout << "[WARNING] Default JEU file with parameters will be used" << std::endl;
        setJEUFileName();
    }

    TString tmp = Form( "%s/aux_files/%s_%i/JEC/%s", 
                        fJECPath.Data(), fCollidingSystem.Data(),
                        fCollidingEnergyGeV, fJEUInputFileName.Data() );
    fJEUInputFileName = tmp;

    fJEU = new JetUncertainty( fJEUInputFileName.Data() );
    std::cout << "JEU file: " << fJEUInputFileName.Data() << std::endl;

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }   
}

//________________
void ForestAODReader::createExtraJECScaleCorrFunction() {
    fJECScaleCorr = new TF1("JetScaleCorrection","[3] + ([0]-[3]) / ( 1.0 + pow( x/[2],[1] ) )", 30, 800);
    fJECScaleCorr->SetParameters(2.27008e+00, 9.18625e-01, 1.43067e+00, 1.00002e+00);
}

//________________
double ForestAODReader::evalCentralityWeight(const double& x) {
    double weight{1.};
    double p0{4.363352};
    double p1{-8.957467e-02};
    double p2{7.301890e-04};
    double p3{-2.885492e-06};
    double p4{4.741175e-09};
    double p5{0.};
    //double p5{-1.407975e-09};

    weight = p0 + p1 * x + p2 * TMath::Power(x, 2) +
             p3 * TMath::Power(x, 3) + p4 * TMath::Power(x, 4) +
             p5 * TMath::Power(x, 5);

    return weight;
}

//_________________
void ForestAODReader::finish() {
    // Nothing to do here
}

//_________________
int ForestAODReader::setupChains() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupChains()";
    }

    // Setup chains (0-good, 1-bad)
    int returnStatus = 1;

    // Setup chains to read

    std::cout << "Setting chains... ";

    // Use event branch
    fEventTree = new TChain("hiEvtAnalyzer/HiTree");

    // Use HLT branch
    if ( fUseHltBranch ) {
        fHltTree = new TChain("hltanalysis/HltTree");
    }
    // Use skimming branch
    if ( fUseSkimmingBranch ) {
        fSkimTree = new TChain("skimanalysis/HltTree");
    }
    // Use particle flow jet branch
    if ( fUseRecoJetBranch ) {
        fRecoJetTree = new TChain( Form( "%s/t", fRecoJetTreeName.Data() ) );
    }
    // Use reconstructed track branch
    if ( fUseTrackBranch ) {
        fTrkTree = new TChain("ppTrack/trackTree");
    }
    // Use generated track branch
    if ( fIsMc && fUseGenTrackBranch ) {
        fGenTrkTree = new TChain("HiGenParticleAna/hi");
    }
    std::cout << "\t[DONE]\n";

    // Initialize input file name (should switch to const char* processing later)
    TString input(fInFileName);

    // Check input exists
    if (  input.Length()<= 0 ) {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        returnStatus = 1;
        exit(0);
    }
    // Normail input
    else {
        // If input is a single ROOT file
        if ( input.Index(".root") > 0 ) {
            std::cout << Form( "Adding %s file to chains\n", input.Data() );
            fEventTree->Add( input.Data() );
            if ( fUseHltBranch ) fHltTree->Add( input.Data() );
            if ( fUseSkimmingBranch ) fSkimTree->Add( input.Data() );
            if ( fUseRecoJetBranch ) fRecoJetTree->Add( input.Data() );
            if ( fUseTrackBranch ) fTrkTree->Add( input.Data() );
            if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( input.Data() );

            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
            Long64_t fEvents2Read2 = fRecoJetTree->GetEntries();
            std::cout << Form("Total number of events to read2: %lld\n", fEvents2Read2 );
        }
        // Assuming that list of files is provided instead of a single file
        else {
            std::ifstream inputStream( input.Data() );

            if ( !inputStream ) std::cout << Form( "ERROR: Cannot open file list: %s\n", input.Data() );
            int nFiles = 0;
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
                if ( file.find(".root") != std::string::npos /* && file.find("Forest") != std::string::npos &&
                     file.find("AOD") != std::string::npos */ ) {
                    
                    // Open file
                    TFile* ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if ( ftmp && !ftmp->IsZombie() && ftmp->GetNkeys() ) {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str() );
                        // Adding file to chains
                        fEventTree->Add( file.c_str() );
                        if ( fUseHltBranch ) fHltTree->Add( file.c_str() );
                        if ( fUseSkimmingBranch ) fSkimTree->Add( file.c_str() );
                        if ( fUseRecoJetBranch ) fRecoJetTree->Add( file.c_str() );
                        if ( fUseTrackBranch ) fTrkTree->Add( file.c_str() );
                        if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( file.c_str() );
                        ++nFiles;
                    } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp) {
                        ftmp->Close();
                    } //if (ftmp)
                } //if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos && file.find("AOD") != std::string::npos )
            } //while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
        } // else {   if file list
        returnStatus = 0;
    } // else {   if normal input

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
    return returnStatus;
}

//_________________
void ForestAODReader::setupBranches() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::setupChains()";
    }

    // Disable all branches - this is important while reading big files
    fEventTree->SetBranchStatus("*", 0);
    if ( fUseHltBranch ) fHltTree->SetBranchStatus("*", 0);
    if ( fUseSkimmingBranch ) fSkimTree->SetBranchStatus("*", 0);
    if ( fUseRecoJetBranch ) fRecoJetTree->SetBranchStatus("*", 0);
    if ( fUseTrackBranch ) fTrkTree->SetBranchStatus("*", 0);
    if ( fUseGenTrackBranch && fIsMc ) fGenTrkTree->SetBranchStatus("*", 0);


    // enable branches of interest -> see definition of each variables above

    // Event quantities
    fEventTree->SetBranchStatus("run", 1);
    fEventTree->SetBranchStatus("evt", 1);
    fEventTree->SetBranchStatus("lumi", 1);
    fEventTree->SetBranchStatus("vz", 1);
    fEventTree->SetBranchStatus("hiBin", 1); //centrality only for PbPb and XeXe
    fEventTree->SetBranchAddress("run", &fRunId);
    fEventTree->SetBranchAddress("evt", &fEventId);
    fEventTree->SetBranchAddress("lumi", &fLumi);
    fEventTree->SetBranchAddress("vz", &fVertexZ);
    fEventTree->SetBranchAddress("hiBin", &fHiBin); 

    if( fIsMc ) {
        fEventTree->SetBranchStatus("weight", 1);
        fEventTree->SetBranchStatus("pthat", 1); 
        fEventTree->SetBranchAddress("weight", &fPtHatWeight);
        fEventTree->SetBranchAddress("pthat", &fPtHat);
    }

    if ( fUseHltBranch ) {

        // Status
        fHltTree->SetBranchStatus("HLT_HIAK4CaloJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4CaloJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4CaloJet60_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4CaloJet80_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4CaloJet100_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet60_Eta5p1_v4", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet80_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet100_Eta5p1_v3", 1);
        fHltTree->SetBranchStatus("HLT_PAAK4PFJet120_Eta5p1_v2", 1);

        fHltTree->SetBranchStatus("HLT_HIAK4PFJet15_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet15_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet30_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet30_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet40_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet40_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet60_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet80_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet120_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet120_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIAK8PFJet15_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet15_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet25_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet25_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet40_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet40_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet60_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet80_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet140_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet140_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIPFJet25_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet25_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet140_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet140_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1", 1);

        // Address
        fHltTree->SetBranchAddress("HLT_HIAK4CaloJet60_v1", &fHLT_HIAK4CaloJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4CaloJet80_v1", &fHLT_HIAK4CaloJet80_v1);
        fHltTree->SetBranchAddress("HLT_PAAK4CaloJet60_Eta5p1_v3", &fHLT_PAAK4CaloJet60_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4CaloJet80_Eta5p1_v3", &fHLT_PAAK4CaloJet80_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4CaloJet100_Eta5p1_v3", &fHLT_PAAK4CaloJet100_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet60_Eta5p1_v4", &fHLT_PAAK4PFJet60_Eta5p1_v4);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet80_Eta5p1_v3", &fHLT_PAAK4PFJet80_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet100_Eta5p1_v3", &fHLT_PAAK4PFJet100_Eta5p1_v3);
        fHltTree->SetBranchAddress("HLT_PAAK4PFJet120_Eta5p1_v2", &fHLT_PAAK4PFJet120_Eta5p1_v2);

        fHltTree->SetBranchAddress("HLT_HIAK4PFJet15_v1", &fHLT_HIAK4PFJet15_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet15_v1_Prescl", &fHLT_HIAK4PFJet15_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet30_v1", &fHLT_HIAK4PFJet30_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet30_v1_Prescl", &fHLT_HIAK4PFJet30_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet40_v1", &fHLT_HIAK4PFJet40_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet40_v1_Prescl", &fHLT_HIAK4PFJet40_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet60_v1", &fHLT_HIAK4PFJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet60_v1_Prescl", &fHLT_HIAK4PFJet60_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet80_v1", &fHLT_HIAK4PFJet80_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet80_v1_Prescl", &fHLT_HIAK4PFJet80_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet120_v1", &fHLT_HIAK4PFJet120_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet120_v1_Prescl", &fHLT_HIAK4PFJet120_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIAK8PFJet15_v1", &fHLT_HIAK8PFJet15_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet15_v1_Prescl", &fHLT_HIAK8PFJet15_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet25_v1", &fHLT_HIAK8PFJet25_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet25_v1_Prescl", &fHLT_HIAK8PFJet25_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet40_v1", &fHLT_HIAK8PFJet40_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet40_v1_Prescl", &fHLT_HIAK8PFJet40_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet60_v1", &fHLT_HIAK8PFJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet60_v1_Prescl", &fHLT_HIAK8PFJet60_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet80_v1", &fHLT_HIAK8PFJet80_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet80_v1_Prescl", &fHLT_HIAK8PFJet80_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet140_v1", &fHLT_HIAK8PFJet140_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet140_v1_Prescl", &fHLT_HIAK8PFJet140_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIPFJet25_v1", &fHLT_HIPFJet25_v1);
        fHltTree->SetBranchAddress("HLT_HIPFJet25_v1_Prescl", &fHLT_HIPFJet25_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIPFJet140_v1", &fHLT_HIPFJet140_v1);
        fHltTree->SetBranchAddress("HLT_HIPFJet140_v1_Prescl", &fHLT_HIPFJet140_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1", &fHLT_HIPuAK4CaloJet80Eta5p1_v1);
        fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1", &fHLT_HIPuAK4CaloJet100Eta5p1_v1);
    } // if ( fUseHltBranch )

    // Skimming quantities
    if ( fUseSkimmingBranch ) {
        // Status
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Tight", 1);
        fSkimTree->SetBranchStatus("HBHEIsoNoiseFilterResult", 1);
        fSkimTree->SetBranchStatus("collisionEventSelectionAODv2", 1);
        fSkimTree->SetBranchStatus("phfCoincFilter2Th4", 1);
        fSkimTree->SetBranchStatus("pPAprimaryVertexFilter", 1);
        fSkimTree->SetBranchStatus("pBeamScrapingFilter", 1);
        fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutG", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGloose", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGtight", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutE", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutEandG", 1);
        fSkimTree->SetBranchStatus("pclusterCompatibilityFilter", 1);

        fSkimTree->SetBranchStatus("phfCoincFilter", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutdz1p0", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGplus", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutVtx1", 1);


        // Address
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &fHBHENoiseFilterResultRun2Loose);
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &fHBHENoiseFilterResultRun2Tight);
        fSkimTree->SetBranchAddress("HBHEIsoNoiseFilterResult", &fHBHEIsoNoiseFilterResult);
        fSkimTree->SetBranchAddress("collisionEventSelectionAODv2", &fCollisionEventSelectionAODv2);
        fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fPhfCoincFilter2Th4);
        fSkimTree->SetBranchAddress("pPAprimaryVertexFilter", &fPPAprimaryVertexFilter);
        fSkimTree->SetBranchAddress("pBeamScrapingFilter", &fPBeamScrapingFilter);
        fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPprimaryVertexFilter);
        fSkimTree->SetBranchAddress("pVertexFilterCutG", &fPVertexFilterCutG);
        fSkimTree->SetBranchAddress("pVertexFilterCutGloose", &fPVertexFilterCutGloose);
        fSkimTree->SetBranchAddress("pVertexFilterCutGtight", &fPVertexFilterCutGtight);
        fSkimTree->SetBranchAddress("pVertexFilterCutE", &fPVertexFilterCutE);
        fSkimTree->SetBranchAddress("pVertexFilterCutEandG", &fPVertexFilterCutEandG);
        fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fPClusterCompatibilityFilter);

        fSkimTree->SetBranchAddress("phfCoincFilter", &fPhfCoincFilter);
        fSkimTree->SetBranchAddress("pVertexFilterCutdz1p0", &fPVertexFilterCutdz1p0);
        fSkimTree->SetBranchAddress("pVertexFilterCutGplus", &fPVertexFilterCutGplus);
        fSkimTree->SetBranchAddress("pVertexFilterCutVtx1", &fPVertexFilterCutVtx1);
    } // if ( fUseSkimmingBranch )

    // Jet quantities
    if ( fUseRecoJetBranch ) {

        fRecoJetTree->SetBranchStatus("nref", 1);
        fRecoJetTree->SetBranchStatus("rawpt", 1);
        fRecoJetTree->SetBranchStatus("trackMax", 1);
        fRecoJetTree->SetBranchStatus("jteta", 1);
        fRecoJetTree->SetBranchStatus("jtphi", 1);
        fRecoJetTree->SetBranchStatus("WTAeta", 1);
        fRecoJetTree->SetBranchStatus("WTAphi", 1);

        fRecoJetTree->SetBranchAddress("nref", &fNRecoJets);
        fRecoJetTree->SetBranchAddress("rawpt", &fRecoJetPt);
        fRecoJetTree->SetBranchAddress("trackMax", &fRecoJetTrackMax);
        fRecoJetTree->SetBranchAddress("jteta", &fRecoJetEta);
        fRecoJetTree->SetBranchAddress("jtphi", &fRecoJetPhi);
        fRecoJetTree->SetBranchAddress("WTAeta", &fRecoJetWTAEta);
        fRecoJetTree->SetBranchAddress("WTAphi", &fRecoJetWTAPhi);

        fRecoJetTree->SetBranchStatus("jtPfNHF", 1);
        fRecoJetTree->SetBranchStatus("jtPfNEF", 1);
        fRecoJetTree->SetBranchStatus("jtPfCHF", 1);
        fRecoJetTree->SetBranchStatus("jtPfMUF", 1);
        fRecoJetTree->SetBranchStatus("jtPfCEF", 1);
        fRecoJetTree->SetBranchStatus("jtPfCHM", 1);
        fRecoJetTree->SetBranchStatus("jtPfCEM", 1);
        fRecoJetTree->SetBranchStatus("jtPfNHM", 1);
        fRecoJetTree->SetBranchStatus("jtPfNEM", 1);
        fRecoJetTree->SetBranchStatus("jtPfMUM", 1);

      	fRecoJetTree->SetBranchAddress("jtPfNHF", &fRecoJtPfNHF);    
 	    fRecoJetTree->SetBranchAddress("jtPfNEF", &fRecoJtPfNEF);
      	fRecoJetTree->SetBranchAddress("jtPfCHF", &fRecoJtPfCHF);
      	fRecoJetTree->SetBranchAddress("jtPfMUF", &fRecoJtPfMUF);
      	fRecoJetTree->SetBranchAddress("jtPfCEF", &fRecoJtPfCEF);
      	fRecoJetTree->SetBranchAddress("jtPfCHM", &fRecoJtPfCHM);
      	fRecoJetTree->SetBranchAddress("jtPfCEM", &fRecoJtPfCEM);
      	fRecoJetTree->SetBranchAddress("jtPfNHM", &fRecoJtPfNHM);
      	fRecoJetTree->SetBranchAddress("jtPfNEM", &fRecoJtPfNEM);
  	    fRecoJetTree->SetBranchAddress("jtPfMUM", &fRecoJtPfMUM);

        // Gen jet quantities
        if ( fIsMc ) {
            fRecoJetTree->SetBranchStatus("ngen", 1);
            fRecoJetTree->SetBranchStatus("genpt", 1);
            fRecoJetTree->SetBranchStatus("geneta", 1);
            fRecoJetTree->SetBranchStatus("genphi", 1);
            fRecoJetTree->SetBranchStatus("WTAgeneta", 1);
            fRecoJetTree->SetBranchStatus("WTAgenphi", 1);
            fRecoJetTree->SetBranchAddress("ngen", &fNGenJets);
            fRecoJetTree->SetBranchAddress("genpt", &fGenJetPt);
            fRecoJetTree->SetBranchAddress("geneta", &fGenJetEta);
            fRecoJetTree->SetBranchAddress("genphi", &fGenJetPhi);
            fRecoJetTree->SetBranchAddress("WTAgeneta", &fGenJetWTAEta);
            fRecoJetTree->SetBranchAddress("WTAgenphi", &fGenJetWTAPhi);
        }

        // Jet-matching quantities
        if ( fIsMc ) {
            fRecoJetTree->SetBranchStatus("refpt", 1);
            fRecoJetTree->SetBranchStatus("refeta", 1);
            fRecoJetTree->SetBranchStatus("refphi", 1);
            fRecoJetTree->SetBranchStatus("refWTAeta", 1);
            fRecoJetTree->SetBranchStatus("refWTAphi", 1);
            fRecoJetTree->SetBranchStatus("refparton_flavor", 1);
            fRecoJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fRecoJetTree->SetBranchAddress("refpt", &fRefJetPt);
            fRecoJetTree->SetBranchAddress("refeta", &fRefJetEta);
            fRecoJetTree->SetBranchAddress("refphi", &fRefJetPhi);
            fRecoJetTree->SetBranchAddress("refWTAeta", &fRefJetWTAEta);
            fRecoJetTree->SetBranchAddress("refWTAphi", &fRefJetWTAPhi);
            fRecoJetTree->SetBranchAddress("refparton_flavor", &fRefJetPartonFlavor);
            fRecoJetTree->SetBranchAddress("refparton_flavorForB", &fRefJetPartonFlavorForB);
        }

    } // if ( fUseRecoJetBranch )

    // Track quantities
    if ( fUseTrackBranch ) {
        fTrkTree->SetBranchStatus("nTrk", 1);
        fTrkTree->SetBranchStatus("trkPt", 1);
        fTrkTree->SetBranchStatus("trkEta", 1);
        fTrkTree->SetBranchStatus("trkPhi", 1);
        fTrkTree->SetBranchStatus("trkPtError", 1);
        fTrkTree->SetBranchStatus("trkDxy1", 1);
        fTrkTree->SetBranchStatus("trkDxyError1", 1);
        fTrkTree->SetBranchStatus("trkDz1", 1);
        fTrkTree->SetBranchStatus("trkDzError1", 1);
        fTrkTree->SetBranchStatus("trkChi2", 1);
        fTrkTree->SetBranchStatus("trkNdof", 1);
        fTrkTree->SetBranchStatus("trkCharge", 1);
        fTrkTree->SetBranchStatus("trkNHit", 1);
        fTrkTree->SetBranchStatus("trkNlayer", 1);
        fTrkTree->SetBranchStatus("highPurity", 1);
        fTrkTree->SetBranchStatus("pfEcal", 1);
        fTrkTree->SetBranchStatus("pfHcal", 1);
        fTrkTree->SetBranchStatus("trkMVA", 1);
        fTrkTree->SetBranchStatus("trkAlgo", 1);

        fTrkTree->SetBranchAddress("nTrk", &fNTracks);
        fTrkTree->SetBranchAddress("trkPt", &fTrackPt);
        fTrkTree->SetBranchAddress("trkEta", &fTrackEta);
        fTrkTree->SetBranchAddress("trkPhi", &fTrackPhi);
        fTrkTree->SetBranchAddress("trkPtError", &fTrackPtErr);
        fTrkTree->SetBranchAddress("trkDxy1", &fTrackDcaXY);
        fTrkTree->SetBranchAddress("trkDxyError1", &fTrackDcaXYErr);
        fTrkTree->SetBranchAddress("trkDz1", &fTrackDcaZ);
        fTrkTree->SetBranchAddress("trkDzError1", &fTrackDcaZErr);
        fTrkTree->SetBranchAddress("trkChi2", &fTrackChi2);
        fTrkTree->SetBranchAddress("trkNdof", &fTrackNDOF);
        fTrkTree->SetBranchAddress("trkCharge", &fTrackCharge);
        fTrkTree->SetBranchAddress("trkNHit", &fTrackNHits);
        fTrkTree->SetBranchAddress("trkNlayer", &fTrackNLayers);
        fTrkTree->SetBranchAddress("highPurity", &fTrackHighPurity);
        fTrkTree->SetBranchAddress("pfEcal", &fTrackPartFlowEcal);
        fTrkTree->SetBranchAddress("pfHcal", &fTrackPartFlowHcal);
        fTrkTree->SetBranchAddress("trkMVA", &fTrackMVA);
        fTrkTree->SetBranchAddress("trkAlgo", &fTrackAlgo);
    } // if ( fUseTrackBranch ) 

    // Gen particle quantities
    if( fIsMc && fUseGenTrackBranch ) {
        fGenTrkTree->SetBranchStatus("pt", 1);
        fGenTrkTree->SetBranchStatus("eta", 1);
        fGenTrkTree->SetBranchStatus("phi", 1);
        fGenTrkTree->SetBranchStatus("chg", 1);
        fGenTrkTree->SetBranchStatus("pdg", 1);
        fGenTrkTree->SetBranchStatus("sube", 1);

        fGenTrkTree->SetBranchAddress("pt", &fGenTrackPt);
        fGenTrkTree->SetBranchAddress("eta", &fGenTrackEta);
        fGenTrkTree->SetBranchAddress("phi", &fGenTrackPhi);
        fGenTrkTree->SetBranchAddress("chg", &fGenTrackCharge);
        fGenTrkTree->SetBranchAddress("pdg", &fGenTrackPid);
        fGenTrkTree->SetBranchAddress("sube", &fGenTrackSube);
    }

    if ( fVerbose ) {
        std::cout << "\t[DONE]" << std::endl;
    }
}

//_________________
void ForestAODReader::report() {
    std::cout << "ForestAODReader::reporting" << std::endl;
}

//_________________
void ForestAODReader::readEvent() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::readEvent()\n";
    }

    if (fIsMc) {
        fRecoJet2GenJetId.clear();
        fGenJet2RecoJet.clear();
    }

    // Or one can call the clearVariables() function (will take more time)
    if (fUseGenTrackBranch && fIsMc) {
        fGenTrackPt.clear();
        fGenTrackEta.clear();
        fGenTrackPhi.clear();
        fGenTrackCharge.clear();
        fGenTrackPid.clear();
        fGenTrackSube.clear();
    }

    if ( fEventsProcessed >= fEvents2Read ) { 
        std::cerr << "ForestAODReader::readEvent() out of entry numbers\n"; 
        fReaderStatus = 2; // End of input stream
    }
    fEventTree->GetEntry( fEventsProcessed );
    if (fUseHltBranch) fHltTree->GetEntry(fEventsProcessed);
    if (fUseSkimmingBranch) fSkimTree->GetEntry(fEventsProcessed);
    if (fUseRecoJetBranch) fRecoJetTree->GetEntry(fEventsProcessed);
    if (fUseTrackBranch) fTrkTree->GetEntry(fEventsProcessed);
    if (fUseGenTrackBranch) fGenTrkTree->GetEntry(fEventsProcessed);
    fEventsProcessed++;

    if ( fVerbose ) {
        std::cout << "Events processed: " << fEventsProcessed << std::endl;
        std::cout << "ForestAODReader::readEvent() \t[DONE]" << std::endl;
    }
}

//________________
void ForestAODReader::fixIndices() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::fixIndices()";
    }

    if (fUseRecoJetBranch) {

        // Loop over reconstructed jets
        for (int iRecoJet{0}; iRecoJet<fNRecoJets; iRecoJet++) {

            if ( fNGenJets <= 0 ) {
                fRecoJet2GenJetId.push_back(-1);
                continue;
            }

            // Must have a gen-matched jet
            if ( fRefJetPt[iRecoJet] < 0) {
                fRecoJet2GenJetId.push_back(-1);
                continue;
            }

            for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++) {
                // Skip Ref and Gen jets that do not match on pT within computational precision
                //std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fGenJetPt[iGenJet] - fRefJetPt[iRecoJet]) << std::endl;

                //std::cout << "iGen: " << iGenJet << " genPt: " << fGenJetPt[iGenJet] << std::endl;
                if ( TMath::Abs(fGenJetPt[iGenJet] - fRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON ) {
                    // If it is the last one
                    if ( iGenJet == (fNGenJets - 1) ) fRecoJet2GenJetId.push_back(-1);
                    continue;
                }
                else {
                    fRefJetEta[iRecoJet] = fGenJetEta[iGenJet];
                    fRefJetPhi[iRecoJet] = fGenJetPhi[iGenJet];
                    fRefJetWTAEta[iRecoJet] = fGenJetWTAEta[iGenJet];
                    fRefJetWTAPhi[iRecoJet] = fGenJetWTAPhi[iGenJet];
                    fRecoJet2GenJetId.push_back(iGenJet);
                    break;
                }

            }
        } //for (int iRecoJet=0; iRecoJet<fNRecoJets; iRecoJet++)

        // Fill the corresponding index in the reco vector and fill the gen
        for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++) {
            std::vector<int>::iterator it=std::find(fRecoJet2GenJetId.begin(), fRecoJet2GenJetId.end(), iGenJet);
            if (it != fRecoJet2GenJetId.end()) {
                fGenJet2RecoJet.push_back( std::distance(fRecoJet2GenJetId.begin(), it) );
            }
            else {
                fGenJet2RecoJet.push_back(-1);
            }
        }
    } // if (fUseRecoJetBranch)

    if ( fVerbose ) {
        std::cout << "\t[DONE]\n";
    }
}

//_________________
Event* ForestAODReader::returnEvent() {

    if ( fVerbose ) {
        std::cout << "ForestAODReader::returnEvent() \n";
    }

    //std::cout << "ForestAODReader::returnEvent" << std::endl;
    readEvent();

    int nBadRecoJets{0};

    if (fFixJetArrays && fIsMc) {
        fixIndices();
    }

    fEvent = new Event();

    // Remove UPC bins
    if ( fIsMc && fCorrectCentMC && fHiBin<10) {
        delete fEvent;
        fEvent = nullptr;
        return fEvent;
    }

    fEvent->setRunId( fRunId );
    fEvent->setEventId( fEventId );
    fEvent->setLumi( fLumi );
    fEvent->setVz( fVertexZ );
    float centW{1.f};
    if ( fIsMc && fCorrectCentMC) {
        // To handle centrality weight
        fEvent->setHiBin( fHiBin - 10 );
        centW = evalCentralityWeight( fHiBin - 10 );
    }
    else if ( fIsMc ) {
        fEvent->setHiBin( fHiBin );
    }
    else {
        fEvent->setHiBin( fHiBin );
    }
    fEvent->setCentralityWeight( centW );
    
    if ( fIsMc ) {
        fEvent->setPtHat( fPtHat );
        fEvent->setPtHatWeight( fPtHatWeight );
    }
    else {
        fEvent->setPtHat( 1. );
        fEvent->setPtHatWeight( 1. );        
    }

    // Fill HLT branch
    if ( fUseHltBranch ) {

        fEvent->trigAndSkim()->setHLT_HIAK4CaloJet60_v1(fHLT_HIAK4CaloJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4CaloJet80_v1(fHLT_HIAK4CaloJet80_v1);
        fEvent->trigAndSkim()->setHLT_PAAK4CaloJet60_Eta5p1_v3(fHLT_PAAK4CaloJet60_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4CaloJet80_Eta5p1_v3(fHLT_PAAK4CaloJet80_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4CaloJet100_Eta5p1_v3(fHLT_PAAK4CaloJet100_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet60_Eta5p1_v4(fHLT_PAAK4PFJet60_Eta5p1_v4);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet80_Eta5p1_v3(fHLT_PAAK4PFJet80_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet100_Eta5p1_v3(fHLT_PAAK4PFJet100_Eta5p1_v3);
        fEvent->trigAndSkim()->setHLT_PAAK4PFJet120_Eta5p1_v2(fHLT_PAAK4PFJet120_Eta5p1_v2);

        fEvent->trigAndSkim()->setHLT_HIAK4PFJet15_v1(fHLT_HIAK4PFJet15_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet15_v1_Prescl(fHLT_HIAK4PFJet15_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet30_v1(fHLT_HIAK4PFJet30_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet30_v1_Prescl(fHLT_HIAK4PFJet30_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet40_v1(fHLT_HIAK4PFJet40_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet40_v1_Prescl(fHLT_HIAK4PFJet40_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet60_v1(fHLT_HIAK4PFJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet60_v1_Prescl(fHLT_HIAK4PFJet60_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet80_v1(fHLT_HIAK4PFJet80_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet80_v1_Prescl(fHLT_HIAK4PFJet80_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet120_v1(fHLT_HIAK4PFJet120_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet120_v1_Prescl(fHLT_HIAK4PFJet120_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIAK8PFJet15_v1(fHLT_HIAK8PFJet15_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet15_v1_Prescl(fHLT_HIAK8PFJet15_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet25_v1(fHLT_HIAK8PFJet25_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet25_v1_Prescl(fHLT_HIAK8PFJet25_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet40_v1(fHLT_HIAK8PFJet40_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet40_v1_Prescl(fHLT_HIAK8PFJet40_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet60_v1(fHLT_HIAK8PFJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet60_v1_Prescl(fHLT_HIAK8PFJet60_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet80_v1(fHLT_HIAK8PFJet80_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet80_v1_Prescl(fHLT_HIAK8PFJet80_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet140_v1(fHLT_HIAK8PFJet140_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet140_v1_Prescl(fHLT_HIAK8PFJet140_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIPFJet25_v1(fHLT_HIPFJet25_v1);
        fEvent->trigAndSkim()->setHLT_HIPFJet25_v1_Prescl(fHLT_HIPFJet25_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIPFJet140_v1(fHLT_HIPFJet140_v1);
        fEvent->trigAndSkim()->setHLT_HIPFJet140_v1_Prescl(fHLT_HIPFJet140_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIPuAK4CaloJet80Eta5p1_v1(fHLT_HIPuAK4CaloJet80Eta5p1_v1);
        fEvent->trigAndSkim()->setHLT_HIPuAK4CaloJet100Eta5p1_v1(fHLT_HIPuAK4CaloJet100Eta5p1_v1);
    }

    // Fill skimming branch
    if ( fUseSkimmingBranch ) {
        fEvent->trigAndSkim()->setHBHENoiseFilterResultRun2Loose(fHBHENoiseFilterResultRun2Loose);
        fEvent->trigAndSkim()->setHBHENoiseFilterResultRun2Tight(fHBHENoiseFilterResultRun2Tight);
        fEvent->trigAndSkim()->setHBHEIsoNoiseFilterResult(fHBHEIsoNoiseFilterResult);
        fEvent->trigAndSkim()->setCollisionEventSelectionAODv2(fCollisionEventSelectionAODv2);
        fEvent->trigAndSkim()->setPhfCoincFilter2Th4(fPhfCoincFilter2Th4);
        fEvent->trigAndSkim()->setPPAprimaryVertexFilter(fPPAprimaryVertexFilter);
        fEvent->trigAndSkim()->setPBeamScrapingFilter(fPBeamScrapingFilter);
        fEvent->trigAndSkim()->setPprimaryVertexFilter(fPprimaryVertexFilter);
        fEvent->trigAndSkim()->setPVertexFilterCutG(fPVertexFilterCutG);
        fEvent->trigAndSkim()->setPVertexFilterCutGloose(fPVertexFilterCutGloose);
        fEvent->trigAndSkim()->setPVertexFilterCutGtight(fPVertexFilterCutGtight);
        fEvent->trigAndSkim()->setPVertexFilterCutE(fPVertexFilterCutE);
        fEvent->trigAndSkim()->setPVertexFilterCutEandG(fPVertexFilterCutEandG);
        fEvent->trigAndSkim()->setPClusterCompatibilityFilter(fPClusterCompatibilityFilter);

        fEvent->trigAndSkim()->setPhfCoincFilter(fPhfCoincFilter);
        fEvent->trigAndSkim()->setPVertexFilterCutdz1p0(fPVertexFilterCutdz1p0);
        fEvent->trigAndSkim()->setPVertexFilterCutGplus(fPVertexFilterCutGplus);
        fEvent->trigAndSkim()->setPVertexFilterCutVtx1(fPVertexFilterCutVtx1);
    }

    //fEvent->print();
    
    // Create particle flow jet instances
    if ( fUseRecoJetBranch ) {

        if ( fVerbose ) {
            std::cout << "Use PF branch \n";
        }

        // Loop over generated jets
        if ( fIsMc && !fEvent->isGenJetCollectionFilled() ) {

            if ( fVerbose ) {
                std::cout << "nGenRecoJets: " << fNGenJets << std::endl;
            }

            if ( fVerbose ) {
                std::cout << "Filling GenJets: " << fNGenJets << std::endl;
            }
            for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++) {
                GenJet *jet = new GenJet{};
                jet->setPt( fGenJetPt[iGenJet] );
                jet->setEta( fGenJetEta[iGenJet] );
                jet->setPhi( fGenJetPhi[iGenJet] );
                jet->setWTAEta( fGenJetWTAEta[iGenJet] );
                jet->setWTAPhi( fGenJetWTAPhi[iGenJet] );
                jet->setFlavor( fRefJetPartonFlavor[fGenJet2RecoJet.at(iGenJet)] );
                jet->setFlavorForB( fRefJetPartonFlavorForB[fGenJet2RecoJet.at(iGenJet)] );
                jet->setPtWeight( 1. );
                if ( fVerbose ) {
                    jet->print();
                }
                
                fEvent->genJetCollection()->push_back( jet );
            } // for (int iGenJet{0}; iGenJet<fNGenJets; iGenJet++)
            if ( fVerbose) {
                std::cout << "GenJetCollection filled" << std::endl;
            }

            // Projection from filling the collection several times
            fEvent->setGenJetCollectionIsFilled();
        } // if ( fIsMc )
        

        // Loop over reconstructed jets
        
        if ( fVerbose ) {
            std::cout << "nRecoJets: " << fNRecoJets << std::endl;
        }

        for (int iJet{0}; iJet<fNRecoJets; iJet++) {

            // Create a new jet instance
            RecoJet *jet = new RecoJet{};

            if ( fIsMc ) {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if ( fRecoJetPt[iJet] > fPtHat ) {
                    nBadRecoJets++;
                }

                // Add index of the matched GenJet
                jet->setGenJetId( fRecoJet2GenJetId.at(iJet) );
            } // if ( fIsMc )

            // Reco
            jet->setPt( fRecoJetPt[iJet] );
            jet->setEta( fRecoJetEta[iJet] );
            jet->setPhi( fRecoJetPhi[iJet] );
            jet->setWTAEta( fRecoJetWTAEta[iJet] );
            jet->setWTAPhi( fRecoJetWTAPhi[iJet] );
            jet->setRawPt( fRecoJetPt[iJet] );
            jet->setTrackMaxPt( fRecoJetTrackMax[iJet] );
            if ( fJEC ) {
                fJEC->SetJetPT( fRecoJetPt[iJet] );
                fJEC->SetJetEta( fRecoJetEta[iJet] );
                fJEC->SetJetPhi( fRecoJetPhi[iJet] );
                double pTcorr = fJEC->GetCorrectedPT();
                if ( fUseManualJEC ) {
                    pTcorr = fRecoJetPt[iJet] * jecManualCorrection( fRecoJetEta[iJet], fRecoJetEta[iJet] );
                }
                if ( fVerbose ) {
                    std::cout << "pTCorr: " << pTcorr << std::endl; 
                }
                if ( fUseExtraJECforAk4Cs ) {
                    if ( fJECScaleCorr ){
                        pTcorr *= fJECScaleCorr->Eval( pTcorr );
                    }
                    else {
                        std::cerr << "No extra correction for ak4cs exists!" << std::endl;
                    }
                }

                if ( fVerbose && fUseExtraJECforAk4Cs ) {
                    std::cout << "pTCorr after extra ak4cs correction for the jetType: " << pTcorr << std::endl; 
                }

                // To check JEC fUseJECSystematics should be outside [-1, 1] range
                if ( fIsMc && ( TMath::Abs( fUseJERSystematics ) <= 1 ) ) {
                    // pTcorr *= extraJERCorr( pTcorr, fRecoJetEta[iJet]);
                    if ( jet->hasMatching() ) {
                        pTcorr *= extraJERCorr( fEvent->genJetCollection()->at( fRecoJet2GenJetId.at(iJet) )->pt(), 
                                                fRecoJetEta[iJet]);
                    }
                }

                // JEU correction for the real data for systematic uncertainty calculation
                if ( fUseJEU !=0 && fJEU && !fIsMc ) {
                    fJEU->SetJetPT( pTcorr );
			        fJEU->SetJetEta( fRecoJetEta[iJet] );
			        fJEU->SetJetPhi( fRecoJetPhi[iJet] );
                    if ( fUseJEU > 0 ) {
                        pTcorr *= (1. + fJEU->GetUncertainty().first);
                    }
                    else {
                        pTcorr *= (1. - fJEU->GetUncertainty().second);
                    }

                    if ( fVerbose ) {
                        std::cout << "pTCorr after JEU: " << pTcorr << std::endl; 
                    }
                }
                jet->setPtJECCorr( pTcorr );
            }
            else { // If no JEC available
                if ( fVerbose ) {
                    std::cout << "No JEC available" << std::endl;
                }
                jet->setPtJECCorr( -999.f );
            }
            jet->setJtPfNHF( fRecoJtPfNHF[iJet] );
            jet->setJtPfNEF( fRecoJtPfNEF[iJet] );
            jet->setJtPfCHF( fRecoJtPfCHF[iJet] );
            jet->setJtPfMUF( fRecoJtPfMUF[iJet] );
            jet->setJtPfCEF( fRecoJtPfCEF[iJet] );
            jet->setJtPfCHM( fRecoJtPfCHM[iJet] );
            jet->setJtPfCEM( fRecoJtPfCEM[iJet] );
            jet->setJtPfNHM( fRecoJtPfNHM[iJet] );
            jet->setJtPfNEM( fRecoJtPfNEM[iJet] );
            jet->setJtPfMUM( fRecoJtPfMUM[iJet] );

            if ( fVerbose ) {
                jet->print();
            }

            // Check fronе-loaded cut
            if ( fJetCut && !fJetCut->pass(jet) ) {
                delete jet;
                continue;
            }

            fEvent->recoJetCollection()->push_back( jet );
        } // for (int iJet{0}; iJet<fNRecoJets; iJet++)
    } // if ( fUseRecoJetBranch )

    if ( fEventCut && !fEventCut->pass(fEvent) ) {
        delete fEvent;
        fEvent = nullptr;
    }

    return fEvent;
}