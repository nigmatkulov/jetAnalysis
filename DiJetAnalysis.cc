/**
 * @file DiJetAnalysis.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Dijet analysis
 * @version 0.1
 * @date 2024-01-09
 * 
 * @copyright Copyright (c) 2024
 * 
 */

// ROOT headers
#include "TF1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TMath.h"

// C++ headers
#include <iostream>

// Jet analysis headers
#include "DiJetAnalysis.h"

ClassImp(DiJetAnalysis)

//________________
DiJetAnalysis::DiJetAnalysis() : BaseAnalysis(), 
    fVzWeight{nullptr}, fDijetPtAveWeight{nullptr},
    fUseCentralityWeight{}, fHM{nullptr},
    fEtaShift{0}, fIsMc{kFALSE}, fIsPPb{kTRUE},
    fLeadJetPtLow{50.}, fSubleadJetPtLow{40.},
    fDijetPhiCut{ 5. * TMath::Pi() / 6},
    fIsPbGoingDir{kFALSE}, fVerbose{kFALSE},
    fNEventsInSample{1000000},
    fIsDijetFound{kFALSE}, fIsDijetJetIdFound{kFALSE},
    fUseMcReweighting{0}, fJetPtBins{75}, fJetPtLow{20},
    fJetPtHi{1520}, fJetPtStep{20}, fSelectJetsInCMFrame{kFALSE},
    fMcReweight{1},
    fEventCounter{0}, fCycleCounter{0},
    fTotalCounter{0} {

    fJetEtaLab[0] = -3.; fJetEtaLab[1] = 3.;
    fJetEtaCM[0] = -2.5; fJetEtaCM[1] = 2.5;
    fPtHatRange[0] = {15.};
    fPtHatRange[1] = {30.};
    for (Int_t i=0; i<fJetPtBins; i++) {
        for (Int_t j=0; j<fJetPtBins; j++) {
            fJetPtLeadPtSubleadReweightMatrix[i][j] = 1;
        }
    } // for (Int_t i=0; i<fJetPtBins; i++)  
}

//________________
DiJetAnalysis::~DiJetAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void DiJetAnalysis::init() {
    // Initialize analysis
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::init" << std::endl;
        print();
    }

    // pT leading, pT subleading weighting matrix
    if ( fUseMcReweighting != 0 ) {
        // Minimum bias
        if ( fUseMcReweighting == 1 ) {

            Double_t nCorr[75][75] = {
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.995499, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.995079, 1.08477, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.89656, 1.06268, 1.11711, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.808471, 0.971749, 1.08511, 1.13102, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.703757, 0.868933, 0.987011, 1.10313, 1.14135, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.729863, 0.775657, 0.894124, 1.01259, 1.09236, 1.10959, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.74028, 0.713858, 0.846946, 0.914762, 1.04291, 1.12351, 1.15514, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.749974, 0.697858, 0.749715, 0.872775, 0.965266, 1.05935, 1.12325, 1.22079, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.853884, 0.674756, 0.718812, 0.839934, 0.860951, 0.974492, 0.984198, 1.12301, 1.11968, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.586357, 0.81878, 0.716049, 0.71069, 0.810955, 0.957883, 0.997449, 1.09737, 1.13382, 1.21784, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1.06215, 0.615578, 0.629749, 0.73391, 0.911308, 0.860416, 0.98769, 0.976659, 1.12483, 1.21906, 1.35041, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 0.77268, 1.07192, 0.853045, 0.673274, 0.856571, 0.681838, 0.803027, 1.13872, 1.03461, 1.06382, 1.17336, 1.13749, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 2.35589, 0.450158, 0.835726, 0.66213, 0.62903, 0.677643, 0.999199, 0.819487, 0.838809, 0.792467, 1.23365, 1.2534, 0.899042, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 3.94298, 1.58952, 0.875866, 0.36284, 0.712022, 0.829763, 0.754013, 0.725412, 0.895415, 0.631429, 1.01474, 1.06823, 1.32973, 1.1453, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 3.35083, 1, 0.6962, 0.892187, 0.897281, 0.689577, 0.527882, 0.688625, 0.969758, 0.682462, 0.920315, 1.09502, 1.3731, 0.641771, 1.19591, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 4.52111, 1.63045, 0.78503, 0.636677, 0.646492, 0.55895, 0.631746, 0.793179, 0.843699, 0.984156, 1.17543, 0.870027, 0.926256, 1.33214, 1.06207, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 6.40666, 8.60684, 1, 1, 0.970381, 0.540695, 0.872672, 0.808735, 0.69977, 1.13454, 0.653863, 0.790035, 1.31893, 1.06561, 1.22553, 1.18426, 1.44756, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 10.2633, 1, 1, 2.01195, 1.63043, 0.503434, 0.833373, 0.609546, 0.630575, 0.577548, 0.703876, 0.906971, 0.746172, 0.959578, 1.32437, 1.195, 1.11067, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 27.0457, 1, 1, 1, 0.533968, 0.602534, 1.16266, 1.48179, 1.11069, 0.911364, 0.496766, 0.779063, 1.21689, 1.25155, 1.04052, 1.2863, 1.18335, 1.06295, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 12.6829, 1, 1, 0.797632, 0.510305, 0.647367, 1.10465, 0.444706, 0.781271, 1.1073, 1.53029, 0.756157, 0.826162, 0.969245, 0.78877, 1.06466, 0.870943, 0.981236, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1.64232, 0.918636, 1.25108, 0.852246, 0.965515, 0.597329, 1.30839, 0.50989, 0.484329, 1.04866, 0.765476, 1.50645, 0.456065, 0.726962, 1.34799, 0.333439, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 6.46282, 1, 1.64593, 1.92037, 1.35669, 1.06628, 1, 1, 2.50355, 0.68806, 1.6023, 1.75169, 1.04021, 0.684093, 1.23812, 0.99222, 1.22494, 0.447758, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1.46784, 1.10842, 1.53332, 1, 2.32399, 1.64494, 1, 0.957691, 0.430072, 0.395075, 0.691688, 0.930347, 0.283549, 0.268405, 1.66517, 2.4681, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 4.3991, 2.53117, 1, 1, 0.951636, 1, 1, 0.742418, 1, 0.642109, 0.581719, 1, 1, 2.43, 1.51133, 1.0808, 0.368298, 0.827042, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 66.1378, 1, 1, 1, 1, 1, 1, 4.24306, 1, 1.24425, 1, 1, 1, 1.81441, 1, 1.53319, 2.06676, 1, 1.10515, 0.979178, 0.480501, 0.493127, 5.43513, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 5.22266, 1, 1, 1.76797, 5.02972, 1, 1.35264, 1.24874, 1.21702, 1, 1.04315, 0.900611, 0.791002, 2.82633, 1.32585, 0.624901, 1.29338, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 87.4172, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1.75749, 1, 1, 1.57032, 4.4692, 2.62525, 1, 1, 1.84138, 0.847514, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5.57903, 1, 1, 1, 1, 2.43812, 1, 2.23438, 2.07921, 1.94394, 1, 1, 4.06176, 1.21309, 1.11192, 1, 1.09271, 2.4338, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 252.874, 1, 1, 349.366, 1, 1, 1, 1, 1, 1, 1, 4.72186, 4.19863, 3.2906, 1, 1, 1, 2.93128, 1, 1, 8.85792, 3.92696, 1.74554, 1, 1, 2.75486, 1.43012, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 70.9991, 1, 1, 1, 1, 1, 1, 1, 4.45301, 1, 4.24824, 1, 1, 1, 1, 1, 1, 1, 1.96315, 1.77972, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 31.8754, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 40.3095, 1, 21.6394, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3.54046, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 34.0973, 1, 1, 1, 1, 1, 1, 1, 1, 10.7918, 1, 1, 1, 1, 1, 1, 1, 1, 8.96026, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5.9751, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 14.8499, 1, 1, 1, 1, 1, 1, 1, 6.80772, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 26.6325, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 388.715, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 226.403, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 664.624, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
            };

            // Copy matrix above to the one that will be used for correction
            for (Int_t i=0; i<fJetPtBins; i++) {
                for (Int_t j=0; j<fJetPtBins; j++) {
                    fJetPtLeadPtSubleadReweightMatrix[i][j] = nCorr[i][j];
                }
            } // for (Int_t i=0; i<fJetPtBins; i++)
        } // else if ( fUseMcReweighting == 1 )
        else {
            // Copy matrix above to the one that will be used for correction
            for (Int_t i=0; i<fJetPtBins; i++) {
                for (Int_t j=0; j<fJetPtBins; j++) {
                    fJetPtLeadPtSubleadReweightMatrix[i][j] = 1;
                }
            } // for (Int_t i=0; i<fJetPtBins; i++)       
        }
    }

    // For MC
    if ( fIsMc ) {
        // If weighting function does not exist
        if ( !fVzWeight ) {
            if ( fIsPPb ) { // Assumed to be pPb8160
                fVzWeight = new TF1("fVzWeight", "pol8", -15.1, 15.1);
                fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);
            }
            else { // Assummed to be pp5020
                fVzWeight = new TF1("fVzWeight", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x", -20., 20.);
                fVzWeight->SetParameters(0.973941, 0.00310622, 0.000711664, -1.83098e-06, 6.9346e-07, 0., 0.);
            }
        }
        if ( fVerbose ) {
            std::cout << "Vz weight function: ";
            fVzWeight->Print();
        }
    }
}

//________________
void DiJetAnalysis::print() {
    std::cout << "----------------------------------------\n";
    std::cout << "DiJetAnalysis parameters:\n";
    std::cout << "Use centrality weight       : " << fUseCentralityWeight << std::endl
              << "Histogram manager           : " << fHM << std::endl
              << "Is MC                       : " << fIsMc << std::endl
              << "Is pPb                      : " << fIsPPb << std::endl
              << "Is Pb-going direction       : " << fIsPbGoingDir << std::endl
              << "eta shift                   : " << fEtaShift << std::endl
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl
              << "Leading jet pT              : " << fLeadJetPtLow << std::endl
              << "SubLeading jet pT           : " << fSubleadJetPtLow << std::endl
              << "Dijet phi cut               : " << fDijetPhiCut << std::endl
              << "Select jets in CM frame     : " << fSelectJetsInCMFrame << std::endl;
    std::cout << "----------------------------------------\n";
}

//________________
Double_t DiJetAnalysis::eventWeight(const Bool_t& isMc, const Bool_t& isPPb, 
                                    const Double_t& ptHat, const Double_t& vz) {
    Double_t weight{1.};
    Double_t genWeight{1.};
    Double_t vzWeight{1.};

    // For Monte Carlo samples
    if ( isMc ) {
    
        if ( fVzWeight ) {
            vzWeight = fVzWeight->Eval( vz );
        }

        // In case of pPb (assumed to be pPb8160)
        if ( isPPb ) {

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
            genWeight /= fNEventsInSample;
            
            // Vz weighting
            vzWeight = 1. / vzWeight;
        } // if ( isPPb )
        else { // Assuming pp5020
           // Only vz weighting will be used. pT hat weighting is
        }

        weight = genWeight * vzWeight;

    } // if ( isMc )


    if ( fVerbose) {
        std::cout << "fNEventsInSample: " << fNEventsInSample << " genWeight: " 
                  << genWeight << " vzWeight: " << vzWeight 
                  << " weight: " << weight << std::endl;
    }

    return weight;
}

//________________
Double_t DiJetAnalysis::deltaPhi(const Double_t& phi1, const Double_t &phi2) {
    Double_t dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ) dphi -= TMath::TwoPi();
    if ( dphi < -TMath::Pi() ) dphi += TMath::TwoPi();
    return dphi;
}

//________________
Bool_t DiJetAnalysis::isGoodGenJet(const GenJet* jet) {
    Bool_t goodJet{kFALSE};
    Double_t etaCut[2] {fJetEtaLab[0], fJetEtaLab[1]}; 
    Double_t eta = jet->eta();

    if ( fSelectJetsInCMFrame ) {

        eta = boostEta2CM( eta );

        etaCut[0] = fJetEtaCM[0];
        etaCut[1] = fJetEtaCM[1];
    }
    else {
        eta = etaLab( eta );
    }

    if ( jet->pt() > 20. && etaCut[0] < eta && eta < etaCut[1] ) {
        goodJet = {kTRUE};
    }
    
    if ( fVerbose ) {
        std::cout << Form("Gen jet cut %s\n", goodJet ? "\t[passed]" : "\t[failed]" );
    }

    return goodJet;
}

//________________
Double_t DiJetAnalysis::dijetEtaInFrame(const Double_t& eta1, const Double_t& eta2, Bool_t isCM) {
    Double_t etaDijet = 0.5 * (eta1 + eta2);
    if ( isCM ) {
        etaDijet = boostEta2CM( etaDijet );
    }
    else {
        etaDijet = etaLab( etaDijet );
    }
    return etaDijet;
}

//________________
Bool_t DiJetAnalysis::isGoodTrkMax(const RecoJet* jet) {
    Bool_t goodTrackMax = {kTRUE};
    Double_t rawPt = jet->rawPt();
    Double_t trackMaxPt = jet->trackMaxPt();
    if ( TMath::Abs( jet->eta() ) < 2.4 && 
         ( trackMaxPt/rawPt < 0.01 ||
           trackMaxPt/rawPt > 0.98) ) {
        goodTrackMax = {kFALSE};
    }

    if ( fVerbose ) {
        std::cout << "TrackMaxPt/rawPt: " << trackMaxPt/rawPt << ( (goodTrackMax) ? " [passed]" : " [failed]" ) 
                  << ( (trackMaxPt/rawPt < 0.01) ? " too low value " : "" ) << ( (trackMaxPt/rawPt > 0.98) ? " too large value " : "" )
                  << std::endl;
    }

    return goodTrackMax;
}

//_________________
Bool_t DiJetAnalysis::isGoodJetId(const RecoJet* jet) {

	Bool_t passJetId = {kFALSE};

    int chm = jet->jtPfCHM();
    int cem = jet->jtPfCEM();
    int mum = jet->jtPfMUM();
    int nhm = jet->jtPfNHM();
    int nem = jet->jtPfNEM();

    float chf = jet->jtPfCHF();
    float cef = jet->jtPfCEF();
    float nhf = jet->jtPfNHF();
    float nef = jet->jtPfNEF();
    float muf = jet->jtPfMUF();

    int chargedMult = chm + cem + mum;
    int neutralMult = nhm + nem;
    int numberOfConstituents = chargedMult + neutralMult;

    float eta = jet->eta();
	
	float chargedEmFracCut{1.}, neutFracCut{1.};
    if ( !fIsLooseJetIdCut ) {
        chargedEmFracCut = {0.9};
        neutFracCut = {0.9};
    }
    else {
        chargedEmFracCut = {0.99};
        neutFracCut = {0.99};
    }

    Bool_t passNHF{kFALSE};
    Bool_t passNEF{kFALSE};
    Bool_t passNumOfConstituents{kTRUE};
    Bool_t passMuonFrac{kTRUE};
    Bool_t passChargedFrac{kTRUE};
    Bool_t passChargedMult{kTRUE};
    Bool_t passChargedEmFrac{kTRUE};
    Bool_t passNeutralMult{kTRUE};
	
    // Check cuts depending on jet pseudorapdity
    if ( TMath::Abs( eta ) <= 2.7 ) {
        
        passNHF = ( nhf < neutFracCut ) ? kTRUE : kFALSE;
        passNEF = ( nhf < neutFracCut ) ? kTRUE : kFALSE;
        passNumOfConstituents = ( numberOfConstituents > 1 ) ? kTRUE : kFALSE;

        if ( !fIsLooseJetIdCut ) { 
            passMuonFrac = ( muf < 0.8 ) ? kTRUE : kFALSE; 
        } // if ( !fIsLooseJetIdCut )

        if( TMath::Abs( eta ) <= 2.4 ) {
            passChargedFrac = ( chf > 0 ) ? kTRUE : kFALSE;
            passChargedMult = ( chargedMult > 0 ) ? kTRUE : kFALSE;
            passChargedEmFrac = ( cef < chargedEmFracCut ) ? kTRUE : kFALSE;
        } // if( TMath::Abs( eta ) <= 2.4 )

    } // if ( TMath::Abs( eta ) <= 2.7 )
    else if ( TMath::Abs( eta ) <= 3.0) {

        passNEF = ( nef > 0.01 ) ? kTRUE : kFALSE;
        passNHF = ( nhf < 0.98 ) ? kTRUE : kFALSE;
        passNeutralMult = ( neutralMult > 2 ) ? kTRUE : kFALSE;

    } // else if ( TMath::Abs( eta ) <= 3.0)
    else  {
        passNEF = ( nef < 0.9 ) ? kTRUE : kFALSE;
        passNeutralMult = (neutralMult > 10 ) ? kTRUE : kFALSE; // CAUTION: JET MET says it should be >10
    } // else 

    passJetId = passNHF && passNEF && passNumOfConstituents && passMuonFrac && 
                passChargedFrac && passChargedMult && passChargedEmFrac && passNeutralMult;

    if ( fVerbose ) {
        std::cout << "JetId selection results: " << ( (passJetId) ? "[passed]" : "[failed]" ) << " Reasons ";
        std::cout << Form("passNHF: %d \tpassNEF: %d \tpassNumConst: %d \tpassMuonFrac: %d \tpassChFrac: %d \tpassChMult: %d \tpassChEmFrac: %d \tpassNeutMult: %d\n", 
                          passNHF, passNEF, passNumOfConstituents, passMuonFrac, passChargedFrac, 
                          passChargedMult , passChargedEmFrac , passNeutralMult);
    }
		
	return passJetId;
}

//________________
Double_t DiJetAnalysis::boostEta2CM(const Double_t &eta) {
    Double_t etaCM = eta;
    // Apply lab frame boost to CM for the pPb 
    if ( fIsPPb ) {
        if ( fIsMc ) { // For embedding: Pb goes to negative, p goes to positive
            if ( fIsPbGoingDir ) {
                etaCM += fEtaShift;
                etaCM = -etaCM;
            }
            else {
                etaCM -= fEtaShift;
            }
        }
        else { // For data: p goes to negative, Pb goes to positive
            if ( fIsPbGoingDir ) {
                etaCM -= fEtaShift;
            }
            else {
                etaCM += fEtaShift;
                etaCM = -etaCM;
            }
        }
    } // if ( fIsPPb )
    else { // if pp
        
    }
    return etaCM;
}

//________________
Double_t DiJetAnalysis::etaLab(const Double_t &eta) {
    Double_t etaL = eta;
    if ( fIsPPb) {
        if ( fIsMc ) { // For embedding: Pb goes to negative, p goes to positive
            if (fIsPbGoingDir) {
                etaL = -etaL;
            }
        }
        else { // For data: p goes to negative, Pb goes to positive
            if (fIsPbGoingDir) {
            }
            else {
                etaL = -etaL;
            }
        }
    }
    else { // For pp apply eta shift
        etaL += fEtaShift;
    }
    return etaL;
}
    
//________________
Bool_t DiJetAnalysis::isGoodRecoJet(const RecoJet* jet) {
    Bool_t goodJet{kFALSE};
    Bool_t goodKine{kFALSE};
    Bool_t hasMatching{kFALSE};

    Double_t etaCut[2] {fJetEtaLab[0], fJetEtaLab[1]};

    Double_t eta = jet->eta();

    if ( fSelectJetsInCMFrame ) {

        eta = boostEta2CM( eta );

        etaCut[0] = fJetEtaCM[0]; 
        etaCut[1] = fJetEtaCM[1];
    }
    else {
        eta = etaLab( eta );
    }

    if ( jet->ptJECCorr() > 20 && etaCut[0] < eta && eta < etaCut[1] ) {
        goodKine = {kTRUE};
    }

    if ( fIsMc ) {
        if ( jet->hasMatching() ) {
            hasMatching = {kTRUE};
        }
    }
    else {
        hasMatching = {kTRUE};
    }

    goodJet = goodKine && hasMatching;

    if ( fVerbose ) {
        std::cout << Form("Reco jet cut %s", goodJet ? "\t[passed]" : "\t[failed]"); 
        if ( goodJet ) {
            std::cout << std::endl;
        }
        else {
            std::cout << Form("\t goodKine: %d hasMatching: %d\n", goodKine, hasMatching);
        }
    } // if ( fVerbose )

    return goodJet;
}

//________________
void DiJetAnalysis::processGenJets(const Event* event, Double_t ptHatW) {

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processGenJets" << std::endl;
    }

    fMcReweight = {1};

    Double_t ptLead{-1.}, ptSubLead{-1.}, etaLead{0.}, etaSubLead{0.},
             phiLead{0.},  phiSubLead{0.};
    Int_t idLead{-1}, idSubLead{-1};
    Bool_t isDijetFound{ kFALSE };

    GenJetIterator genJetIter;
    Int_t counter{0};

    // Loop over generated jets
    for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

        Double_t pt = (*genJetIter)->pt();
        Double_t eta = (*genJetIter)->eta();
        Double_t phi = (*genJetIter)->phi();

        if ( fVerbose ) {
            std::cout << "Gen jet #" << counter << " ";
            (*genJetIter)->print();
        }

        // Apply single-jet selection to gen jets
        //if ( !isGoodGenJet( *genJetIter ) ) continue;
        
        if ( pt > ptLead ) {
            ptSubLead = ptLead;
            etaSubLead = etaLead;
            phiSubLead = phiLead;
            idSubLead = idLead;
            ptLead = pt;
            etaLead = eta;
            phiLead = phi;
            idLead = counter;
        }
        else if ( pt > ptSubLead ) {
            ptSubLead = pt;
            etaSubLead = eta;
            phiSubLead = phi;
            idSubLead = counter;
        }
        // Fill inclusive jet pt
        fHM->hGenInclusiveJetPt->Fill(pt, ptHatW);
        fHM->hGenInclusiveJetPtEta->Fill(eta, pt, ptHatW);

        if ( fVerbose ) {
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idLead: %d idSubLead: %d\n", 
                              ptLead, ptSubLead, idLead, idSubLead);
        }

        if ( pt > 30. ) {
            fHM->hGenGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), ptHatW);
            fHM->hGenGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), ptHatW );
        }

        counter++;
    } // for ( genJetIter = event->genJetCollection()->begin();

    //
    // Check for gen dijet
    //
    if ( idLead>=0 && idSubLead>=0 ) {

        Bool_t goodLeadJet = isGoodGenJet( event->genJetCollection()->at( idLead ) );
        Bool_t goodSubLeadJet = isGoodGenJet( event->genJetCollection()->at( idSubLead ) );
        Bool_t goodDijet = isGoodDijet( ptLead, ptSubLead, TMath::Abs( deltaPhi(phiLead, phiSubLead) ) );
        isDijetFound = goodLeadJet && goodSubLeadJet && goodDijet;

        // Analyze gen dijets
        if ( isDijetFound ) {

            fHM->hGenPtLeadPtSublead->Fill( ptLead, ptSubLead, ptHatW );
            fHM->hGenEtaLeadEtaSublead->Fill( etaLab(etaLead), etaLab(etaSubLead), ptHatW );
            fHM->hGenEtaCMLeadEtaCMSublead->Fill( boostEta2CM(etaLead), boostEta2CM(etaSubLead), ptHatW );
            fHM->hGenPtLeadPtSubleadMcReweight->Fill( ptLead, ptSubLead, ptHatW * fMcReweight );
            fHM->hGenEtaLeadEtaSubleadMcReweight->Fill( etaLab(etaLead), etaLab(etaSubLead), ptHatW * fMcReweight );

            Double_t dijetPt = 0.5 * (ptLead + ptSubLead);
            Double_t dijetEta = dijetEtaInFrame(etaLead, etaSubLead, kFALSE);
            Double_t dijetDphi = deltaPhi(phiLead, phiSubLead);
            Double_t dijetEtaCM = dijetEtaInFrame(etaLead, etaSubLead, kTRUE);

            Double_t genDijetLeadSublead[9] {dijetPt, dijetEta, dijetDphi, 
                                             ptLead, etaLab(etaLead), phiLead, 
                                             ptSubLead, etaLab(etaSubLead), phiSubLead };
            fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(genDijetLeadSublead);
            fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(genDijetLeadSublead, ptHatW * fMcReweight );
            fHM->hGenDijetEta->Fill(dijetEta, ptHatW * fMcReweight );
            fHM->hGenDijetPtEtaDphi->Fill(dijetPt, dijetEta, dijetDphi, 1.);
            fHM->hGenDijetPtEtaDphiWeighted->Fill(dijetPt, dijetEta, dijetDphi, ptHatW * fMcReweight );
            fHM->hGenDijetEtaCM->Fill(dijetEtaCM, ptHatW * fMcReweight );
            fHM->hGenDijetPtEtaDphiCM->Fill(dijetPt, dijetEtaCM, dijetDphi, 1.);
            fHM->hGenDijetPtEtaDphiCMWeighted->Fill(dijetPt, dijetEtaCM, dijetDphi, ptHatW * fMcReweight );

            (dijetEta >= 0) ? fHM->hGenDijetPtEtaForward->Fill(dijetPt, dijetEta) : fHM->hGenDijetPtEtaBackward->Fill(dijetPt, TMath::Abs(dijetEta));
            (dijetEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForward->Fill(dijetPt, dijetEtaCM) : fHM->hGenDijetPtEtaCMBackward->Fill(dijetPt, TMath::Abs(dijetEtaCM));
            (dijetEta >= 0) ? fHM->hGenDijetPtEtaForwardWeighted->Fill(dijetPt, dijetEta, ptHatW * fMcReweight) : fHM->hGenDijetPtEtaBackwardWeighted->Fill(dijetPt, TMath::Abs(dijetEta), ptHatW * fMcReweight);
            (dijetEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForwardWeighted->Fill(dijetPt, dijetEtaCM, ptHatW * fMcReweight) : fHM->hGenDijetPtEtaCMBackwardWeighted->Fill(dijetPt, TMath::Abs(dijetEtaCM), ptHatW * fMcReweight);
        } // if ( isDijetFound )
    } // if ( idLead>=0 && idSubLead>=0 )

    if ( fVerbose ) {
        std::cout << "Gen dijet found: " << ( (isDijetFound) ? "[true]" : "[false]" ) << std::endl;
        std::cout << "Reporting from DiJetAnalysis::processGenJets - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::findMcWeight(const Double_t& ptLead, const Double_t& ptSublead) {

    if ( fUseMcReweighting !=0 ) {
        if ( fVerbose ) {
            std::cout << Form("DiJetAnalysis::findMcWeight - ptLead: %5.1f ptSublead: %5.1f\n", ptLead, ptSublead);
        }

        Int_t ptLeadBin{-1}; 
        if ( ptLead >= fJetPtLow && ptLead <= fJetPtHi ) {
            ptLeadBin = ( ptLead - fJetPtLow ) / fJetPtStep;
        }
        Int_t ptSubleadBin{-1};
        if ( ptSublead >= fJetPtLow && ptSublead <= fJetPtHi ) {
            ptSubleadBin = ( ptSublead - fJetPtLow ) / fJetPtStep;
        }
        Double_t val = ( ptLeadBin >=0 && ptSubleadBin >= 0 ) ? 
                    fJetPtLeadPtSubleadReweightMatrix[ptLeadBin][ptSubleadBin] : 1.;

        if ( fVerbose ) {
            std::cout << Form("\t ptLeadBin: %d ptSubleadBin: %d weight: %6.3f\n", ptLeadBin, ptSubleadBin, val);
        }    

        fMcReweight = val;
    }
    else {
        fMcReweight = {1};
    }
}

//________________
void DiJetAnalysis::processRecoJets(const Event* event, Double_t ptHatW) {

    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processRecoJets" << std::endl;
    }

    fMcReweight = {1.};

    // Define variables
    Double_t ptRecoLead{-1.}, ptRecoSubLead{-1.},
             ptRawRecoLead{-1.}, ptRawRecoSubLead{-1.},
             etaRecoLead{0.}, etaRecoSubLead{0.},
             phiRecoLead{0.},  phiRecoSubLead{0.}, 
             ptRefLead{-1.}, ptRefSubLead{-1.},
             etaRefLead{0.}, etaRefSubLead{0.},
             phiRefLead{0.}, phiRefSubLead{0.},
             ptRecoLeadJetId{-1.}, ptRecoSubLeadJetId{-1.},
             etaRecoLeadJetId{0.}, etaRecoSubLeadJetId{0.},
             phiRecoLeadJetId{0.},  phiRecoSubLeadJetId{0.},
             ptRefLeadJetId{-1.}, ptRefSubLeadJetId{-1.},
             etaRefLeadJetId{0.}, etaRefSubLeadJetId{0.},
             phiRefLeadJetId{0.}, phiRefSubLeadJetId{0.};
    Int_t  idRecoLead{-1}, idRecoSubLead{-1}, idRecoLeadJetId{-1}, idRecoSubLeadJetId{-1};

    // Loop over reconstructed jets
    RecoJetIterator recoJetIter;
    Int_t counter{0};
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        Double_t pt = (*recoJetIter)->ptJECCorr();
        Double_t eta = (*recoJetIter)->eta();
        Double_t phi = (*recoJetIter)->phi();
        Double_t ptRaw = (*recoJetIter)->pt();

        if ( fVerbose ) {
            std::cout << "Reco jet #" << counter << " ";
            (*recoJetIter)->print();
        }

        // JetId parameters
        int chargedMult = (*recoJetIter)->jtPfCHM() + (*recoJetIter)->jtPfCEM() + (*recoJetIter)->jtPfMUM();
        int neutralMult = (*recoJetIter)->jtPfNHM() + (*recoJetIter)->jtPfNEM();
        int numberOfConstituents = chargedMult + neutralMult;

        Int_t dummyIter{0};
        if ( TMath::Abs( eta ) <= 2.4 ) { dummyIter = {0}; }
        else if ( TMath::Abs( eta ) <= 2.7 ) { dummyIter = {1}; }
        else if ( TMath::Abs( eta ) <= 3.0 ) { dummyIter = {2}; }
        else { dummyIter = {3}; }

        //std::cout << "pT: " << pt << " eta: " << eta << " ptHatW: " << ptHatW << std::endl;
        fHM->hRecoInclusiveJetPt->Fill(pt, ptHatW);
        fHM->hRecoInclusiveAllJetPtVsEta->Fill(eta, pt, ptHatW);
        fHM->hRecoInclusiveJetPtVsEtaKineCut->Fill(eta, pt, ptHatW);

        // JetId histograms
        fHM->hNHF[dummyIter]->Fill( (*recoJetIter)->jtPfNHF(), ptHatW );
        fHM->hNEmF[dummyIter]->Fill( (*recoJetIter)->jtPfNEF(), ptHatW );
        fHM->hNumOfConst[dummyIter]->Fill( numberOfConstituents, ptHatW );
        fHM->hMUF[dummyIter]->Fill( (*recoJetIter)->jtPfMUF(), ptHatW );
        fHM->hCHF[dummyIter]->Fill( (*recoJetIter)->jtPfCHF(), ptHatW );
        fHM->hChargedMult[dummyIter]->Fill( chargedMult, ptHatW );
        fHM->hCEmF[dummyIter]->Fill( (*recoJetIter)->jtPfCEF(), ptHatW );
        fHM->hNumOfNeutPart[dummyIter]->Fill( neutralMult, ptHatW );

        // Local variables for the current jet analysis
        Bool_t passTrkMax{kFALSE};
        Bool_t passJetId{kFALSE};
        Bool_t hasMatching{kFALSE};
        GenJet *matchedJet{nullptr};
        Double_t genPt{-999};
        Double_t genEta{-999};
        Double_t genPhi{-999};
        Double_t JES {-999};
        Double_t dEta{-999};
        Double_t dPhi{-999};
        Double_t res[4] { JES, genPt, genEta, genPhi };

        // Check selection criteria
        passTrkMax = isGoodTrkMax( (*recoJetIter) );
        passJetId = isGoodJetId( (*recoJetIter) );

        // On MC check reco jet matching to gen
        if ( fIsMc ) {
            hasMatching = (*recoJetIter)->hasMatching();
            if ( hasMatching ) {
                matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
                genPt = matchedJet->pt();
                genEta = matchedJet->eta();
                genPhi = matchedJet->phi();

                JES = pt/genPt;
                dEta = eta - genEta;
                dPhi = phi - genPhi;
                res[0] = JES;
                res[1] = genPt; 
                res[2] = genEta;
                res[3] = genPhi;

                fHM->hRefInclusiveJetPt->Fill(genPt, ptHatW);
                fHM->hRefInclusiveJetPtEta->Fill(genEta, genPt, ptHatW);
                fHM->hRecoMatchedPtEta->Fill(eta, pt, ptHatW);

                Double_t correl[5] { pt, ptRaw, genPt, eta, genEta };
                fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, ptHatW);
                fHM->hJESInclusiveJetPtEtaPhi->Fill(res);
                fHM->hJESInclusiveJetPtEtaPhiWeighted->Fill(res, ptHatW);

                fHM->hRecoInclusiveMatchedJetPtVsEta->Fill(eta, pt, ptHatW);
                fHM->hRecoInclusiveMatchedJetPtVsEtaKineCut->Fill(eta, pt, ptHatW);
                fHM->hRecoInclusiveJetRefPtVsEtaKineCut->Fill(genEta, genPt, ptHatW);
                fHM->hRecoInclusiveJetJESPtEtaPhiKineCut->Fill(res, ptHatW);
                fHM->hRecoInclusiveJetDEtaPtEtaKineCut->Fill(dEta, genPt, genEta, ptHatW);
            } // if ( hasMatching )
            else {
                // Fill unmatched jets
                fHM->hRecoInclusiveUnmatchedJetPtVsEta->Fill(eta, pt, ptHatW);
                fHM->hRecoInclusiveUnmatchedJetPtVsEtaKineCut->Fill(eta, pt, ptHatW);
            } // else
        } // if ( fIsMc )


        //
        // For track max selection
        //
        if ( passTrkMax ) {

            //
            // Inclusive jet part
            //

            fHM->hRecoInclusiveJetPtVsEtaTrkMaxCut->Fill(eta, pt, ptHatW);

            // For MC only
            if ( fIsMc ) {
                if ( !hasMatching ) {
                    fHM->hRecoInclusiveUnmatchedJetPtVsEtaTrkMaxCut->Fill(eta, pt, ptHatW);
                }
                else {
                    fHM->hRecoInclusiveMatchedJetPtVsEtaTrkMaxCut->Fill(eta, pt, ptHatW);
                    fHM->hRecoInclusiveJetRefPtVsEtaTrkMaxCut->Fill(genEta, genPt, ptHatW);
                    fHM->hRecoInclusiveJetJESPtEtaPhiTrkMaxCut->Fill(res, ptHatW);
                    fHM->hRecoInclusiveJetDEtaPtEtaTrkMaxCut->Fill(dEta, genPt, genEta, ptHatW);
                } // 
            } // if ( fIsMc )

            //
            // Dijet part
            //

            // Check for leading and subleading jets
            if ( pt > ptRecoLead ) {
                ptRecoSubLead = ptRecoLead;
                etaRecoSubLead = etaRecoLead;
                phiRecoSubLead = phiRecoLead;
                ptRawRecoSubLead = ptRawRecoLead;
                idRecoSubLead = idRecoLead;
                ptRecoLead = pt;
                ptRawRecoLead = ptRaw;
                etaRecoLead = eta;
                phiRecoLead = phi;
                idRecoLead = counter;

                if ( fIsMc ) {
                    if ( hasMatching ) {
                        ptRefSubLead = ptRefLead;
                        etaRefSubLead = etaRefLead;
                        phiRefSubLead = phiRefLead;
                        ptRefLead = genPt;
                        etaRefLead = genEta;
                        phiRefLead = genPhi;
                    }
                    else {
                        ptRefSubLead = ptRefLead;
                        etaRefSubLead = etaRefLead;
                        phiRefSubLead = phiRefLead;
                        ptRefLead = {-999.};
                        etaRefLead = {-999.};
                        phiRefLead = {-999.};
                    }
                } // if ( fIsMc )
            } // if ( pt > ptRecoLead )
            else if ( pt > ptRecoSubLead ) {
                ptRecoSubLead = pt;
                ptRawRecoSubLead = ptRaw;
                etaRecoSubLead = eta;
                phiRecoSubLead = phi;
                idRecoSubLead = counter;

                if ( fIsMc ) {
                    if ( hasMatching ) {
                        ptRefSubLead = genPt;
                        etaRefSubLead = genEta;
                        phiRefSubLead = genPhi;
                    }
                    else {
                        ptRefSubLead = {-999.};
                        etaRefSubLead = {-999.};
                        phiRefSubLead = {-999.};
                    }
                } // if ( fIsMc )
            } // else if ( pt > ptRecoSubLead )

            if ( pt > 30. ) {
                fHM->hRecoGoodInclusiveJetEtaLabFrame->Fill(eta, ptHatW);
                fHM->hRecoGoodInclusiveJetEtaCMFrame->Fill(boostEta2CM(eta), ptHatW);
            }
        } // if ( passTrkMax )

        //
        // For jetId selection
        //
        if ( passJetId ) {

            //
            // Inclusive jet part
            //

            fHM->hRecoInclusiveJetPtVsEtaJetIdCut->Fill(eta, pt, ptHatW);

            // For MC only
            if ( fIsMc ) { 
                // If not matched
                if ( !hasMatching ) {
                    fHM->hRecoInclusiveUnmatchedJetPtVsEtaJetIdCut->Fill(eta, pt, ptHatW);
                }
                else { // If matched
                    fHM->hRecoInclusiveMatchedJetPtVsEtaJetIdCut->Fill(eta, pt, ptHatW);
                    fHM->hRecoInclusiveJetRefPtVsEtaJetIdCut->Fill(genEta, genPt, ptHatW);
                    fHM->hRecoInclusiveJetJESPtEtaPhiJetIdCut->Fill(res, ptHatW);
                    fHM->hRecoInclusiveJetDEtaPtEtaJetIdCut->Fill(dEta, genPt, genEta, ptHatW);
                }
            } // if ( fIsMc )

            //
            // Dijet part
            //

            if ( pt > ptRecoLeadJetId ) {
                ptRecoSubLeadJetId = ptRecoLeadJetId;
                etaRecoSubLeadJetId = etaRecoLeadJetId;
                phiRecoSubLeadJetId = phiRecoLeadJetId;
                idRecoSubLeadJetId = idRecoLeadJetId;
                ptRecoLeadJetId = pt;
                etaRecoLeadJetId = eta;
                phiRecoLeadJetId = phi;
                idRecoLeadJetId = counter;

                if ( fIsMc ) { 
                    if (hasMatching ) {
                        ptRefSubLeadJetId = ptRefLeadJetId;
                        etaRefSubLeadJetId = etaRefLeadJetId;
                        phiRefSubLeadJetId = phiRefLeadJetId;
                        ptRefLeadJetId = genPt;
                        etaRefLeadJetId = genEta;
                        phiRefLeadJetId = genPhi;
                    }
                    else {
                        ptRefSubLeadJetId = ptRefLeadJetId;
                        etaRefSubLeadJetId = etaRefLeadJetId;
                        phiRefSubLeadJetId = phiRefLeadJetId;
                        ptRefLeadJetId = {-999.};
                        etaRefLeadJetId = {-999.};
                        phiRefLeadJetId = {-999.};                       
                    }
                } // if ( fIsMc )
            } // if ( pt > ptRecoLeadJetId )
            else if ( pt > ptRecoSubLeadJetId ) {
                ptRecoSubLeadJetId = pt;
                etaRecoSubLeadJetId = eta;
                phiRecoSubLeadJetId = phi;
                idRecoSubLeadJetId = counter;

                if ( fIsMc ) {
                    if ( hasMatching ) {
                        ptRefSubLeadJetId = genPt;
                        etaRefSubLeadJetId = genEta;
                        phiRefSubLeadJetId = genPhi;
                    }
                    else {
                        ptRefLeadJetId = {-999.};
                        etaRefLeadJetId = {-999.};
                        phiRefLeadJetId = {-999.}; 
                    }
                } // if ( fIsMc )
            } // else if ( pt > ptRecoSubLeadJetId )
        } // if ( passJetId )

        if ( fVerbose ) {
            std::cout << Form("TrkMax selection --> Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d\n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead);
            std::cout << Form("JetId selection  --> Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d\n", 
                              ptRecoLeadJetId, ptRecoSubLeadJetId, idRecoLeadJetId, idRecoSubLeadJetId);
        }

        // Increment counter
        counter++;
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )


    //
    // TrkMax dijets
    //

    if ( fVerbose ) {
        std::cout << "Start checking TrkMax dijets\n";
    }

    fIsDijetFound = {kFALSE};

    if ( idRecoLead>=0 && idRecoSubLead>=0 ) {

        fHM->hRecoLeadJetAllPtVsEta->Fill(etaRecoLead, ptRecoLead, ptHatW);
        fHM->hRecoSubLeadJetAllPtVsEta->Fill(etaRecoSubLead, ptRecoSubLead, ptHatW);   

        if ( fIsMc ) {
            // Check leading jet matching to gen
            if ( event->recoJetCollection()->at( idRecoLead )->hasMatching() ) {
                fHM->hRecoLeadJetMatchedPtVsEta->Fill( etaLab(etaRecoLead), ptRecoLead, ptHatW);
            }
            else {
                fHM->hRecoLeadJetUnmatchedPtVsEta->Fill( etaLab(etaRecoLead), ptRecoLead, ptHatW);
            }

            // Check subleading jet matching to gen
            if ( event->recoJetCollection()->at( idRecoSubLead )->hasMatching() ) {
                fHM->hRecoSubLeadJetMatchedPtVsEta->Fill( etaLab(etaRecoSubLead), ptRecoSubLead, ptHatW);
            }
            else {
                fHM->hRecoSubLeadJetUnmatchedPtVsEta->Fill( etaLab(etaRecoSubLead), ptRecoSubLead, ptHatW);
            }
        } // if ( fIsMc )

        Bool_t goodLeadJet = isGoodRecoJet( event->recoJetCollection()->at( idRecoLead ) );
        if ( fVerbose ) {
            std::cout << Form("Leading jet is %s\n", ((goodLeadJet) ? "good" : "bad") );  
        }
        Bool_t goodSubLeadJet = isGoodRecoJet( event->recoJetCollection()->at( idRecoSubLead ) );
        if ( fVerbose ) {
            std::cout << Form("SubLeading jet is %s\n", ((goodSubLeadJet) ? "good" : "bad") );  
        }
        Bool_t goodDijet = isGoodDijet( ptRecoLead, ptRecoSubLead, TMath::Abs( deltaPhi(phiRecoLead, phiRecoSubLead) ) );
        if ( fVerbose ) {
            std::cout << Form("Dijet is %s\n", ((goodDijet) ? "good" : "bad") );  
        }
        fIsDijetFound = goodLeadJet && goodSubLeadJet && goodDijet;

        if ( fVerbose ) {
            std::cout << Form("Dijet status %s\n", ((fIsDijetFound) ? "[GOOD]" : "[BAD]") );  
        }

        // Analyze trkMax dijets
        if ( fIsDijetFound ) {

            if ( fVerbose ) {
                std::cout << "Start filling dijet histograms";
            }

            // Dijet analysis
            Double_t dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
            Double_t dijetRecoEta = dijetEtaInFrame(etaRecoLead, etaRecoSubLead, kFALSE);
            Double_t dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);
            Double_t dijetRecoEtaCM = dijetEtaInFrame(etaRecoLead, etaRecoSubLead, kTRUE);


            // Correlation between leading and subleading
            fHM->hRecoPtLeadPtSublead->Fill( ptRecoLead, ptRecoSubLead, ptHatW );
            fHM->hRecoEtaLeadEtaSublead->Fill( etaLab(etaRecoLead), etaLab(etaRecoSubLead), ptHatW );
            fHM->hRecoEtaCMLeadEtaCMSublead->Fill( boostEta2CM(etaRecoLead), boostEta2CM(etaRecoSubLead), ptHatW );
            fHM->hRecoPtLeadPtSubleadMcReweight->Fill( ptRecoLead, ptRecoSubLead, ptHatW * fMcReweight );
            fHM->hRecoEtaLeadEtaSubleadMcReweight->Fill( etaLab(etaRecoLead), etaLab(etaRecoSubLead), ptHatW * fMcReweight );

            Double_t dijetRecoInfo[9] { dijetRecoPt, dijetRecoEta, dijetRecoDphi,
                                        ptRecoLead, etaLab(etaRecoLead), phiRecoLead,
                                        ptRecoSubLead, etaLab(etaRecoSubLead), phiRecoSubLead };
            fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(dijetRecoInfo);
            fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(dijetRecoInfo, ptHatW * fMcReweight);
            fHM->hRecoDijetEta->Fill( dijetRecoEta, ptHatW * fMcReweight);
            fHM->hRecoDijetPtEta->Fill( dijetRecoPt, dijetRecoEta, ptHatW * fMcReweight);
            fHM->hRecoDijetPtEtaDphi->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, 1. );
            fHM->hRecoDijetPtEtaDphiWeighted->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, ptHatW * fMcReweight);
            fHM->hRecoDijetEtaCM->Fill( dijetRecoEtaCM, ptHatW * fMcReweight);
            fHM->hRecoDijetPtEtaDphiCM->Fill( dijetRecoPt, dijetRecoEtaCM, dijetRecoDphi, 1. );
            fHM->hRecoDijetPtEtaDphiCMWeighted->Fill( dijetRecoPt, dijetRecoEtaCM, dijetRecoDphi, ptHatW * fMcReweight);

            if (dijetRecoEta >= 0) {
                fHM->hRecoDijetPtEtaForward->Fill(dijetRecoPt, dijetRecoEta);
                fHM->hRecoDijetPtEtaForwardWeighted->Fill(dijetRecoPt, dijetRecoEta, ptHatW * fMcReweight);
            }
            else {
                fHM->hRecoDijetPtEtaBackward->Fill(dijetRecoPt, TMath::Abs(dijetRecoEta));
                fHM->hRecoDijetPtEtaBackwardWeighted->Fill(dijetRecoPt, TMath::Abs(dijetRecoEta), ptHatW * fMcReweight);
            }
                
            if (dijetRecoEtaCM >= 0) {
                fHM->hRecoDijetPtEtaCMForward->Fill(dijetRecoPt, dijetRecoEtaCM);
                fHM->hRecoDijetPtEtaCMForwardWeighted->Fill(dijetRecoPt, dijetRecoEtaCM, ptHatW * fMcReweight);
            }
            else {
                fHM->hRecoDijetPtEtaCMBackward->Fill(dijetRecoPt, TMath::Abs(dijetRecoEtaCM));
                fHM->hRecoDijetPtEtaCMBackwardWeighted->Fill(dijetRecoPt, TMath::Abs(dijetRecoEtaCM), ptHatW * fMcReweight);
            }

            if ( fIsMc ) {

                Double_t dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
                Double_t dijetRefEta = dijetEtaInFrame( etaRefLead, etaRefSubLead, kFALSE);
                Double_t dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);
                Double_t dijetRefEtaCM = dijetEtaInFrame( etaRefLead, etaRefSubLead, kTRUE);

                // Leading jet information
                Double_t correl[5] { ptRecoLead, ptRawRecoLead, ptRefLead, etaLab(etaRecoLead), etaLab(etaRefLead) };
                fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, ptHatW * fMcReweight );

                // Subleading jet information
                correl[0] = ptRecoSubLead;
                correl[1] = ptRawRecoSubLead;
                correl[2] = ptRefSubLead;
                correl[3] = etaLab(etaRecoSubLead); 
                correl[4] = etaLab(etaRefSubLead);
                fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, ptHatW * fMcReweight);

                fHM->hRefPtLeadPtSublead->Fill( ptRefLead, ptRefSubLead, ptHatW );
                fHM->hRefEtaLeadEtaSublead->Fill( ptRefLead, ptRefSubLead, ptHatW );
                fHM->hRefEtaCMLeadEtaCMSublead->Fill( boostEta2CM(etaRefLead), boostEta2CM(etaRefSubLead), ptHatW );
                fHM->hRefPtLeadPtSubleadMcReweight->Fill( ptRefLead, ptRefSubLead, ptHatW * fMcReweight );
                fHM->hRefEtaLeadEtaSubleadMcReweight->Fill( ptRefLead, ptRefSubLead, ptHatW * fMcReweight );

                Double_t dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                                 ptRecoLead, etaLab(etaRecoLead),
                                                 ptRecoSubLead, etaLab(etaRecoSubLead),
                                                 dijetRefPt, dijetRefEta,
                                                 ptRefLead, etaLab(etaRefLead),
                                                 ptRefSubLead, etaLab(etaRefSubLead) };

                Double_t dijetUnfold[4] = { dijetRecoPt, dijetRecoEta, dijetRefPt, dijetRefEta };

                fHM->hRecoDijetPtEtaRefDijetPtEta->Fill(dijetUnfold, 1.);
                fHM->hRecoDijetPtEtaRefDijetPtEtaWeighted->Fill(dijetUnfold, ptHatW * fMcReweight);

                fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Fill(dijetRecoUnfold);
                fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, ptHatW * fMcReweight );
                fHM->hRefDijetEta->Fill( dijetRefEta, ptHatW * fMcReweight );
                fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEta, dijetRefEta, ptHatW * fMcReweight );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, 1.);
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, ptHatW * fMcReweight );
                fHM->hRefDijetPtEtaDphi->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, 1. );
                fHM->hRefDijetPtEtaDphiWeighted->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW * fMcReweight );

                fHM->hRefDijetEtaCM->Fill( dijetRefEtaCM, ptHatW );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPt, 1.);
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPt, ptHatW * fMcReweight );
                fHM->hRefDijetPtEtaDphiCM->Fill( dijetRefPt, dijetRefEtaCM, dijetRefDphi, 1. );
                fHM->hRefDijetPtEtaDphiCMWeighted->Fill( dijetRefPt, dijetRefEtaCM, dijetRefDphi, ptHatW * fMcReweight );

                (dijetRefEta >= 0) ? fHM->hRefDijetPtEtaForward->Fill(dijetRefPt, dijetRefEta) : fHM->hRefDijetPtEtaBackward->Fill(dijetRefPt, TMath::Abs(dijetRefEta));
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForward->Fill(dijetRefPt, dijetRefEtaCM) : fHM->hRefDijetPtEtaCMBackward->Fill(dijetRefPt, TMath::Abs(dijetRefEtaCM));
                (dijetRefEta >= 0) ? fHM->hRefDijetPtEtaForwardWeighted->Fill(dijetRefPt, dijetRefEta, ptHatW * fMcReweight) : fHM->hRefDijetPtEtaBackwardWeighted->Fill(dijetRefPt, TMath::Abs(dijetRefEta), ptHatW * fMcReweight);
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForwardWeighted->Fill(dijetRefPt, dijetRefEtaCM, ptHatW * fMcReweight) : fHM->hRefDijetPtEtaCMBackwardWeighted->Fill(dijetRefPt, TMath::Abs(dijetRefEtaCM), ptHatW * fMcReweight);
                        
            } // if ( fIsMc )

            if ( fVerbose ) {
                std::cout << "\t[DONE]\n";
            }
        } // if ( fIsDijetFound )
    } // if ( idRecoLead>=0 && idRecoSubLead>=0 )


    //
    // JetId dijets
    //

    if ( fVerbose ) {
        std::cout << "Start checking JetId dijets\n";
    }

    fIsDijetJetIdFound = {kFALSE};

    if ( idRecoLeadJetId>=0 && idRecoSubLeadJetId>=0 ) {

        fHM->hRecoLeadJetAllPtVsEtaJetIdCut->Fill( etaLab(etaRecoLeadJetId), ptRecoLeadJetId, ptHatW);
        fHM->hRecoSubLeadJetAllPtVsEtaJetIdCut->Fill( etaLab(etaRecoSubLeadJetId), ptRecoSubLeadJetId, ptHatW);   

        if ( fIsMc ) {
            // Check leading jet matching to gen
            if ( event->recoJetCollection()->at( idRecoLeadJetId )->hasMatching() ) {
                fHM->hRecoLeadJetMatchedPtVsEtaJetIdCut->Fill( etaLab(etaRecoLeadJetId), ptRecoLeadJetId, ptHatW);
            }
            else {
                fHM->hRecoLeadJetUnmatchedPtVsEtaJetIdCut->Fill( etaLab(etaRecoLeadJetId), ptRecoLeadJetId, ptHatW);
            }

            // Check subleading jet matching to gen
            if ( event->recoJetCollection()->at( idRecoSubLeadJetId )->hasMatching() ) {
                fHM->hRecoSubLeadJetMatchedPtVsEtaJetIdCut->Fill( etaLab(etaRecoSubLeadJetId), ptRecoSubLeadJetId, ptHatW);
            }
            else {
                fHM->hRecoSubLeadJetUnmatchedPtVsEtaJetIdCut->Fill( etaLab(etaRecoSubLeadJetId), ptRecoSubLeadJetId, ptHatW);
            }
        } // if ( fIsMc )
        
        Bool_t goodLeadJet = isGoodRecoJet( event->recoJetCollection()->at( idRecoLeadJetId ) );
        if ( fVerbose ) {
            std::cout << Form("Leading jet is %s\n", ((goodLeadJet) ? "good" : "bad") );  
        }
        Bool_t goodSubLeadJet = isGoodRecoJet( event->recoJetCollection()->at( idRecoSubLeadJetId ) );
        if ( fVerbose ) {
            std::cout << Form("SubLeading jet is %s\n", ((goodSubLeadJet) ? "good" : "bad") );  
        }
        Bool_t goodDijet = isGoodDijet( ptRecoLeadJetId, ptRecoSubLeadJetId, TMath::Abs( deltaPhi(phiRecoLeadJetId, phiRecoSubLeadJetId) ) );
        if ( fVerbose ) {
            std::cout << Form("Dijet is %s\n", ((goodDijet) ? "good" : "bad") );  
        }
        fIsDijetJetIdFound = goodLeadJet && goodSubLeadJet && goodDijet;

        if ( fVerbose ) {
            std::cout << Form("Dijet status %s\n", ((fIsDijetJetIdFound) ? "[GOOD]" : "[BAD]") );  
        }

        // Analyze jetId dijets
        if ( fIsDijetJetIdFound ) {

            if ( fVerbose ) {
                std::cout << "Start filling dijet histograms";
            }

            // Dijet analysis
            Double_t dijetRecoPt = 0.5 * (ptRecoLeadJetId + ptRecoSubLeadJetId);
            Double_t dijetRecoEta = dijetEtaInFrame( etaRecoLeadJetId, etaRecoSubLeadJetId, kFALSE );
            Double_t dijetRecoDphi = deltaPhi(phiRecoLeadJetId, phiRecoSubLeadJetId);
            Double_t dijetRecoEtaCM = dijetEtaInFrame( etaRecoLeadJetId, etaRecoSubLeadJetId, kTRUE );

            // Apply lab frame boost to CM for the pPb 
            dijetRecoEta = etaLab( dijetRecoEta );
            dijetRecoEtaCM = boostEta2CM( dijetRecoEtaCM );

            fHM->hRecoDijetPtEtaDphiJetId->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, ptHatW );

            if ( fIsMc ) {

                Double_t dijetRefPt = 0.5 * (ptRefLeadJetId + ptRefSubLeadJetId);
                Double_t dijetRefEta = dijetEtaInFrame( etaRefLeadJetId, etaRefSubLeadJetId, kFALSE);
                Double_t dijetRefDphi = deltaPhi(phiRefLeadJetId, phiRefSubLeadJetId);
                Double_t dijetRefEtaCM = dijetEtaInFrame( etaRefLeadJetId, etaRefSubLeadJetId, kTRUE);

                dijetRefEta = etaLab( dijetRefEta );
                dijetRefEtaCM = boostEta2CM( dijetRefEtaCM );


                fHM->hRefDijetPtEtaDphiJetId->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtJetId->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, ptHatW);
            } // if ( fIsMc )

            if ( fVerbose ) {
                std::cout << "\t[DONE]\n";
            }

        } // if ( fIsDijetJetIdFound )
    } // if ( idRecoLeadJetId>=0 && idRecoSubLeadJetId>=0 )

    // Fill matching between the histograms
    if ( fIsDijetFound && fIsDijetJetIdFound ) {
        if ( (idRecoLead == idRecoLeadJetId) && (idRecoSubLead == idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(1., ptHatW);
        }
        else if ( (idRecoLead == idRecoLeadJetId) && (idRecoSubLead != idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(2., ptHatW);
        }
        else if ( (idRecoLead != idRecoLeadJetId) && (idRecoSubLead == idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(3., ptHatW);
        }
        else if ( (idRecoLead != idRecoLeadJetId) && (idRecoSubLead != idRecoSubLeadJetId) ) {
            fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(4., ptHatW);
        }
    }
    else if ( fIsDijetFound && !fIsDijetJetIdFound ) {
        fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(5., ptHatW);
    }
    else if ( !fIsDijetFound && fIsDijetJetIdFound ) {
        fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(6., ptHatW);   
    }
    else if ( !fIsDijetFound && !fIsDijetJetIdFound ) {
        fHM->hRecoTrkMaxToJetIdDijetMatching->Fill(7., ptHatW);
    }

    if ( fVerbose ) {
        std::cout << "TrkMax dijet found: " << ( (fIsDijetFound) ? "[true]" : "[false]" ) << std::endl;
        // std::cout << "JetId  dijet found: " << ( (fIsDijetJetIdFound) ? "[true]" : "[false]" ) << std::endl;
        std::cout << "Reporting from DiJetAnalysis::processRecoJets - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRefJets(const Event* event, Double_t ptHatW) {
    if ( fVerbose ) {
        std::cout << "Reporting from DiJetAnalysis::processRefJets" << std::endl;
    }

    fMcReweight = {1.};

    Double_t ptRecoLead{-1.}, ptRecoSubLead{-1.},
             ptRawRecoLead{-1.}, ptRawRecoSubLead{-1.},
             etaRecoLead{0.}, etaRecoSubLead{0.},
             phiRecoLead{0.},  phiRecoSubLead{0.}, 
             ptRefLead{-1.}, ptRefSubLead{-1.},
             etaRefLead{0.}, etaRefSubLead{0.},
             phiRefLead{0.}, phiRefSubLead{0.};
    Bool_t isDijetFound{kFALSE};
    Int_t  idRecoLead{-1}, idRecoSubLead{-1};

    // Loop over reconstructed jets
    RecoJetIterator recoJetIter;
    Int_t counter{0};
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        if ( !(*recoJetIter)->hasMatching() ) continue;

        GenJet *matchedJet{nullptr};
        Double_t genPt{999.};
        Double_t genEta{-999.};
        Double_t genPhi{-999.};

        Double_t pt = (*recoJetIter)->ptJECCorr();
        Double_t eta = (*recoJetIter)->eta();
        Double_t phi = (*recoJetIter)->phi();
        Double_t ptRaw = (*recoJetIter)->pt();

        matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
        genPt = matchedJet->pt();
        genEta = matchedJet->eta();
        genPhi = matchedJet->phi();

        if ( fVerbose ) {
            std::cout << "Ref jet info for reco jet #" << counter;
            matchedJet->print();
            std::cout << "Reco jet #" << counter << " ";
            (*recoJetIter)->print();
        }
        
        if ( genPt > ptRefLead ) {
            ptRecoSubLead = ptRecoLead;
            etaRecoSubLead = etaRecoLead;
            phiRecoSubLead = phiRecoLead;
            ptRawRecoSubLead = ptRawRecoLead;
            idRecoSubLead = idRecoLead;
            ptRecoLead = pt;
            ptRawRecoLead = ptRaw;
            etaRecoLead = eta;
            phiRecoLead = phi;
            idRecoLead = counter;

            ptRefSubLead = ptRefLead;
            etaRefSubLead = etaRefLead;
            phiRefSubLead = phiRefLead;
            ptRefLead = genPt;
            etaRefLead = genEta;
            phiRefLead = genPhi;
        }
        else if ( genPt > ptRefSubLead ) {
            ptRecoSubLead = pt;
            ptRawRecoSubLead = ptRaw;
            etaRecoSubLead = eta;
            phiRecoSubLead = phi;
            idRecoSubLead = counter;

            ptRefSubLead = genPt;
            etaRefSubLead = genEta;
            phiRefSubLead = genPhi;
        }

        if ( fVerbose ) {
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d genId: %d \n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead, (*recoJetIter)->genJetId());
        }

        // Increment counter
        counter++;
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    //
    // Check if leading and subleading jets were found
    //
    if (idRecoLead>=0 && idRecoSubLead>=0) {
        if ( fVerbose ) {
            std::cout << Form("Checking dijet Lead pT: %5.2f SubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d Lead genId: %d SubLead genId: %d\n", 
                              ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead, 
                              event->recoJetCollection()->at( idRecoLead )->genJetId(), 
                              event->recoJetCollection()->at( idRecoSubLead )->genJetId() );
        }
        Bool_t goodLeadJet{kFALSE};
        Bool_t goodSubLeadJet{kFALSE};
        Bool_t goodDijet{kFALSE};
        if ( event->recoJetCollection()->at( idRecoLead )->hasMatching() && 
             event->recoJetCollection()->at( idRecoSubLead )->hasMatching() ) {
            goodLeadJet = isGoodGenJet( event->genJetCollection()->at( event->recoJetCollection()->at( idRecoLead )->genJetId() ) );
            goodSubLeadJet = isGoodGenJet( event->genJetCollection()->at( event->recoJetCollection()->at( idRecoSubLead )->genJetId() ) );
            goodDijet = isGoodDijet( ptRefLead, ptRefSubLead, TMath::Abs( deltaPhi(phiRefLead, phiRefSubLead) ) );
        }
        isDijetFound = goodLeadJet && goodSubLeadJet && goodDijet;

        // Analyze trkMax dijets
        if ( isDijetFound ) {

            // Dijet analysis
            Double_t dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
            Double_t dijetRecoEta = dijetEtaInFrame(etaRecoLead, etaRecoSubLead, kFALSE);
            Double_t dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);
            Double_t dijetRecoEtaCM = dijetEtaInFrame(etaRecoLead, etaRecoSubLead, kTRUE);

            Double_t dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
            Double_t dijetRefEta = dijetEtaInFrame(etaRefLead, etaRefSubLead, kFALSE);
            Double_t dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);
            Double_t dijetRefEtaCM = dijetEtaInFrame(etaRefLead, etaRefSubLead, kTRUE);

            // Dijet reco vs ref for unfolding
            Double_t dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                             ptRecoLead, etaLab(etaRecoLead),
                                             ptRecoSubLead, etaLab(etaRecoSubLead),
                                             dijetRefPt, dijetRefEta,
                                             ptRefLead, etaLab(etaRefLead),
                                             ptRefSubLead, etaLab(etaRefSubLead) };    

            fHM->hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, ptHatW * fMcReweight );
            fHM->hRefSelDijetEta->Fill(dijetRefEta, ptHatW * fMcReweight );
            fHM->hRefSelDijetPtEtaDphi->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, 1.);
            fHM->hRefSelDijetPtEtaDphiWeighted->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, ptHatW * fMcReweight );
            fHM->hRefSelDijetEtaCM->Fill(dijetRefEtaCM, ptHatW * fMcReweight );
            fHM->hRefSelDijetPtEtaDphiCM->Fill(dijetRefPt, dijetRefEtaCM, dijetRefDphi, 1.);
            fHM->hRefSelDijetPtEtaDphiCMWeighted->Fill(dijetRefPt, dijetRefEtaCM, dijetRefDphi, ptHatW * fMcReweight );
        } // if ( isDijetFound )
    } // if (idRecoLead>=0 && idRecoSubLead>=0)

    if ( fVerbose ) {
        std::cout << "Dijet found: " << isDijetFound << std::endl;
        std::cout << "Reporting from DiJetAnalysis::processRefJets - [DONE]" << std::endl;
    }
}

//________________
Bool_t DiJetAnalysis::isGoodDijet(const Double_t& ptLead, const Double_t& ptSublead, const Double_t& dphi) {
    Bool_t isGood = ( ptLead > fLeadJetPtLow &&
                      ptSublead > fSubleadJetPtLow && 
                      dphi > fDijetPhiCut );
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::isGoodDijet " << isGood << " ";
        std::cout << Form("pTlead: %5.2f pTsub: %5.2f dphi: %4.2f\n", ptLead, ptSublead, dphi);
    }

    // Check reweight
    if ( fIsMc && fUseMcReweighting != 0 ) {
        findMcWeight(ptLead, ptSublead);
    }
    else {
        fMcReweight = {1.};
    }

    if ( fIsMc ) {
        if ( (0.5 * (ptLead + ptSublead) ) > (2 * fPtHatRange[1]) ) {
            isGood = {kFALSE};
        }
    }

    return isGood;
}

//________________
void DiJetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    if ( fVerbose ) {
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "DiJetAnalysis::processEvent" << std::endl;
    }

    fTotalCounter++;
    fEventCounter++;
    if ( fEventCounter >= 50000 ) {
        fCycleCounter++;
        std::cout << Form("DiJetAnalysis::processEvent [INFO] Events processed: %d Sample fraction: %3.2f%%", 
                          fCycleCounter * 50000, (Double_t)(fCycleCounter * 50000) / fNEventsInSample * 100. )
                  << std::endl;
        fEventCounter = {0};
    }

    if ( !fHM ) {
        std::cout << "[Warning] No histogram manager connected to the DiJetAnalysis\n";
    }

    // Must be flushed for each event !!!!
    fIsDijetFound = {kFALSE};
    fIsDijetJetIdFound = {kFALSE};

    //
    // Event quantities
    //

    // ptHat
    Double_t ptHat = event->ptHat();
    Double_t vz = event->vz();
    // ptHat weight (a.k.a. event weight that includes the ptHat and vz)
    Double_t ptHatW{1.};
    // Check correct MC sample
    if ( fIsPPb ) { // Assume pPb8160
        if ( fIsMc ) {
             
            // Skip events with ptHat that is outside the ranged embedded
            if ( ptHat <= fPtHatRange[0] || ptHat > fPtHatRange[1] ) {
                if ( fVerbose ) {
                    std::cout << Form("[WARNING] Bad ptHat value: %4.1f < %4.1f <= %4.1f\n", fPtHatRange[0], ptHat, fPtHatRange[1]);
                }
                return;
            }

            // For MC we need to flip the direction of Pb-going to properly reweight distributions
            if ( fIsPbGoingDir ) {
                vz = -vz;
            }
            ptHatW = eventWeight(fIsMc, fIsPPb, ptHat, vz);
        }
        else {
            if ( !fIsPbGoingDir ) {
                vz = -vz;
            }
        }
    } // if ( fIsPPb )
    else { // Assume pp5020
        ptHatW = event->ptHatWeight() * eventWeight(fIsMc, fIsPPb, ptHat, vz);
    }

    // Process and analyze reco jets
    processRecoJets(event, ptHatW);

    if ( fIsMc ) {
        // Process and analyze gen jets
        processGenJets(event, ptHatW);
        processRefJets(event, ptHatW);
    }

    Double_t centW = event->centralityWeight();
    centW = {1.}; // Do not apply weight for pPb or pp

    //std::cout << "centrality weight: " << centW << std::endl;

    // For dijet analysis and reweighting purposes it is important 
    // to fill event histograms only when dijet is found
    if ( fIsDijetFound ) {
        fHM->hHiBin->Fill( event->hiBin() );

        fHM->hVz->Fill( vz,  centW );
        fHM->hVzWeighted->Fill( vz, ptHatW * centW );

        fHM->hHiBinWeighted->Fill( event->hiBin(), ptHatW * centW );

        fHM->hPtHatWeighted->Fill( ptHat, ptHatW * centW );
        fHM->hPtHatWeight->Fill( ptHatW, centW );

        Double_t vzPtHat[2] = { vz, ptHat };
        fHM->hVzPtHat->Fill( vzPtHat, centW );
        fHM->hVzPtHatWeighted->Fill( vzPtHat, ptHatW * centW );
    } // if ( fIsDijetFound ) 

    // Fill this outside to be sure that MC is calculated properly,
    // i.e. different ptHat samples are taken with the proper weight
    fHM->hPtHat->Fill( ptHat, centW );

    if ( fVerbose ) {
        std::cout << "Event quantities were read properly" << std::endl;
        //event->print();
        std::cout << "DiJetAnalysis::processEvent - [DONE]" << std::endl;
    }
}

//________________
void DiJetAnalysis::finish() {
    // Save data and close files
    fCycleCounter++;
    std::cout << Form("DiJetAnalysis::processEvent [INFO] Total events processed: %d Sample fraction: %3.2f%%", 
                      (fCycleCounter * 50000) + fEventCounter, 
                      (Double_t)(fCycleCounter * 50000 + fEventCounter) / fNEventsInSample )
              << std::endl;
    std::cout << Form("DiJetAnalysis::processEvent [INFO]: Total number of events processed: %d", fTotalCounter) << std::endl;
    std::cout << "DiJetAnalysis::finish" << std::endl;
}

//________________
void DiJetAnalysis::report() {
    // Force to report everyone
}

//________________
TList* DiJetAnalysis::getOutputList() {
    TList *outputList = new TList();

    // Add list of settings for cuts
    return outputList;
}
