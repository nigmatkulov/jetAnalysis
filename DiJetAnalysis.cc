/**
 * @file DiJetAnalysis.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Dijet analysis
 * @version 1.1
 * @date 2025-01-09
 * 
 * @copyright Copyright (c) 2025
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

//________________
DiJetAnalysis::DiJetAnalysis() : BaseAnalysis(), 
    fVzWeight{nullptr}, fDijetPtAveWeight{nullptr},
    fUseCentralityWeight{}, fHM{nullptr},
    fEtaShift{0.465}, fIsMc{false}, fCollisionSystem{1}, fCollisionEnergy{8160},
    fLeadJetPtLow{50.}, fSubleadJetPtLow{40.},
    fDijetPhiCut{ 2. * TMath::Pi() / 3},
    fIsPbGoingDir{false}, fVerbose{false},
    fNEventsInSample{1000000},
    fUseJetIdSelection{false}, fIsLooseJetIdCut{false}, 
    fIsGenDijetLabFound{false}, fIsGenDijetCMFound{false},
    fIsRecoDijetLabFound{false}, fIsRecoDijetCMFound{false},
    fIsRefSelDijetLabFound{false}, fIsRefSelDijetCMFound{false},
    fUseMcReweighting{0}, fJetPtBins{75}, fJetPtLow{20},
    fJetPtHi{1520}, fJetPtStep{20}, fSelectJetsInCMFrame{false},
    fMcReweight{1}, fPtAveBins{}, fPtAveOldBins{} {

    fJetEtaLab[0] = -3.; fJetEtaLab[1] = 3.;
    fJetEtaCM[0] = -2.5; fJetEtaCM[1] = 2.5;
    fPtHatRange[0] = {0};
    fPtHatRange[1] = {100000000};
    for (int i=0; i<fJetPtBins; i++) {
        for (int j=0; j<fJetPtBins; j++) {
            fJetPtLeadPtSubleadReweightMatrix[i][j] = 1;
        }
    } // for (int i=0; i<fJetPtBins; i++)

    double dijetPtVals[17] {  50.,  60.,   70.,  80.,  90.,
                              100., 110.,  120., 130., 140.,
                              150., 160.,  180., 200., 250., 
                              300., 500.};
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    fPtAveBins.assign(dijetPtVals, dijetPtVals + sizeOfPtVals);

    double dijetPtOldVals[7] {25., 55., 75., 95., 115., 150., 400.};
    int sizeOfPtOldVals = sizeof(dijetPtOldVals)/sizeof(dijetPtOldVals[0]);
    fPtAveOldBins.assign(dijetPtOldVals, dijetPtOldVals + sizeOfPtOldVals);
}

//________________
DiJetAnalysis::~DiJetAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
}

//________________
void DiJetAnalysis::init() {
    // Initialize analysis
    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::init -- begin" << std::endl;
        print();
    }

    // pT leading, pT subleading weighting matrix
    if ( fUseMcReweighting != 0 ) {
        // Minimum bias
        if ( fUseMcReweighting == 1 ) {

            double nCorr[75][75] = {
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
            for (int i=0; i<fJetPtBins; i++) {
                for (int j=0; j<fJetPtBins; j++) {
                    fJetPtLeadPtSubleadReweightMatrix[i][j] = nCorr[i][j];
                }
            } // for (int i=0; i<fJetPtBins; i++)
        } // else if ( fUseMcReweighting == 1 )
        else {
            // Copy matrix above to the one that will be used for correction
            for (int i=0; i<fJetPtBins; i++) {
                for (int j=0; j<fJetPtBins; j++) {
                    fJetPtLeadPtSubleadReweightMatrix[i][j] = 1;
                }
            } // for (int i=0; i<fJetPtBins; i++)       
        }
    }

    // For MC
    if ( fIsMc ) {
        // Initialize vz weight function
        initVzWeightFunction();
    }

    // Initialize vz weight function
    initVzWeightFunction();

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::init -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::initVzWeightFunction() {

    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::initVzWeightFunction -- begin" << std::endl;
    }

    // Check if Vz weight function exists
    if ( !fVzWeight ) {
        if ( fCollisionSystem == 0 ) { // Assume pp 5020
            fVzWeight = new TF1("fVzWeight", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x", -20., 20.);
            fVzWeight->SetParameters(0.973941, 0.00310622, 0.000711664, -1.83098e-06, 6.9346e-07, 0., 0.);
        }
        else if ( fCollisionSystem == 1 ) { // Assume pPb 8160
            fVzWeight = new TF1("fVzWeight", "pol8", -15.1, 15.1);
            fVzWeight->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);
        }
        else if ( fCollisionSystem == 2 ) { // Assume PbPb 5020
            fVzWeight = new TF1("fVzWeight", "pol0", -15.1, 15.1);
            fVzWeight->SetParameter(0, 1.);
        }
        else { // Unknown collision system
            fVzWeight = new TF1("fVzWeight", "pol0", -200.1, 200.1);
            fVzWeight->SetParameter(0, 1.); 
        }
    }

    if ( fVerbose ) {
        std::cout << "Vz weight function: ";
        fVzWeight->Print();

        std::cout << "DiJetAnalysis::initVzWeightFunction -- end" << std::endl;
    }
}

//________________
TString DiJetAnalysis::collisionSystem() const {
    TString collSys = "PbPb";
    if ( fCollisionSystem == 0 ) {
        collSys = "pp";
    }
    else if ( fCollisionSystem == 1 ) {
        collSys = "pPb";
    }
    else if ( fCollisionSystem == 2 ) {
        collSys = "PbPb";
    }
    else {
        collSys = "Unknown";
    }

    return collSys;
}

//________________
void DiJetAnalysis::print() {
    std::cout << "----------------------------------------\n";
    std::cout << "DiJetAnalysis parameters:\n";
    std::cout << "Use centrality weight       : " << fUseCentralityWeight << std::endl
              << "Histogram manager           : " << fHM << std::endl
              << "Is MC                       : " << fIsMc << std::endl
              << "Collision system            : " << collisionSystem().Data() << std::endl
              << "Collision energy (GeV)      : " << fCollisionEnergy << std::endl
              << "Is Pb-going direction       : " << fIsPbGoingDir << std::endl
              << "eta shift                   : " << fEtaShift << std::endl
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl
              << "Use jetId selection         : " << fUseJetIdSelection << std::endl
              << "Use loose jetId cut         : " << fIsLooseJetIdCut << std::endl
              << "Leading jet pT              : " << fLeadJetPtLow << std::endl
              << "SubLeading jet pT           : " << fSubleadJetPtLow << std::endl
              << "Dijet phi cut               : " << fDijetPhiCut << std::endl
              << "Select jets in CM frame     : " << fSelectJetsInCMFrame << std::endl;
    std::cout << "----------------------------------------\n";
}

//________________
int DiJetAnalysis::findDijetPtAveBin(const double &ptAve) {
    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::findDijetPtAveBin -- begin" << std::endl;
    }
    int bin{-1};
    if ( fPtAveBins[0] < ptAve && ptAve < fPtAveBins.at( fPtAveBins.size()-1 ) ) {
        for (unsigned int i=0; i<fPtAveBins.size()-1; i++) {
            if ( fPtAveBins[i] <= ptAve && ptAve < fPtAveBins[i+1] ) {
                bin = i;
                break;
            }
        }
    }

    if ( fVerbose ) {
        std::cout << Form("ptAve: %5.2f bin: %d\n", ptAve, bin);
        std::cout << "DiJetAnalysis::findDijetPtAveBin -- end" << std::endl;
    }
    return bin;
}

//________________
int DiJetAnalysis::findDijetPtAveOldBin(const double &ptAve) {
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::findDijetPtAveOldBin -- begin" << std::endl;
    }
    int bin{-1};
    if ( fPtAveOldBins[0] < ptAve && ptAve < fPtAveOldBins.at( fPtAveOldBins.size()-1 ) ) {
        for (unsigned int i=0; i<fPtAveOldBins.size()-1; i++) {
            if ( fPtAveOldBins[i] <= ptAve && ptAve < fPtAveOldBins[i+1] ) {
                bin = i;
                break;
            }
        }
    }

    if ( fVerbose ) {
        std::cout << Form("ptAve: %5.2f bin: %d\n", ptAve, bin);
        std::cout << "DiJetAnalysis::findDijetPtAveOldBin -- end" << std::endl;
    }
    return bin;
}

//________________
double DiJetAnalysis::eventWeight(const double& ptHat, const double& vz, 
                                  const double& centWeight, const double& ptHatW) {

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::eventWeight -- begin" << std::endl;
    }

    // Calculate event weight
    double weight{1.};
    double genWeight{1.};
    double vzWeight{1.};

    if ( fIsMc ) {
        // In case of pPb (assumed to be pPb8160)
        if ( fCollisionSystem == 0 ) {     // Assuming pp 5020
            weight = ptHatW;
            if ( fVzWeight ) {
                vzWeight = fVzWeight->Eval( vz );
            }
            weight *= vzWeight;
        }
        else if ( fCollisionSystem == 1 ) { // Assuming pPb 8160

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
            if ( fVzWeight ) {
                vzWeight = fVzWeight->Eval( vz );
            }
            vzWeight = 1. / vzWeight;

            weight = genWeight * vzWeight;
        } 
        else if ( fCollisionSystem == 2 ) { // Assuming PbPb 5020
            weight = ptHatW;
            if ( fUseCentralityWeight ) {
                weight *= centWeight;
            }
        }
        else {
            weight = 1.;
        }
    } // if ( fIsMc )

    if ( fVerbose ) {
        std::cout << Form("Input parameters:\nptHat: %5.2f vz: %5.2f centWeight: %5.2f ptHatW: %5.2f\n", 
                          ptHat, vz, centWeight, ptHatW);
        std::cout << Form("Calculated parameters:\nweight: %5.2f genWeight: %5.2f vzWeight: %5.2f\n", weight, genWeight, vzWeight);
        if ( fCollisionSystem == 1 ) {
            std::cout << Form("fNEventsInSample: %d\n", fNEventsInSample);
        }
        std::cout << "DiJetAnalysis::eventWeight -- end" << std::endl;
    }

    return weight;
}

//________________
void DiJetAnalysis::findLeadSubleadJets(const double &pt, const int &counter,
                                        double &ptLead, double &ptSublead,
                                        int &idLead, int &idSubLead) {
    // Find leading and subleading jets
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::findLeadSubleadJets -- begin" << std::endl;
    }

    if ( pt > ptLead ) {
        ptSublead = ptLead;
        idSubLead = idLead;
        ptLead = pt;
        idLead = counter;
    }
    else if ( pt > ptSublead ) {
        ptSublead = pt;
        idSubLead = counter;
    }

    if ( fVerbose ) {
        std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f Lead id: %d SubLead id: %d\n", 
                          ptLead, ptSublead, idLead, idSubLead);
        std::cout << "DiJetAnalysis::findLeadSubleadJets - end" << std::endl;
    }
}

//________________
double DiJetAnalysis::deltaPhi(const double& phi1, const double &phi2) {
    double dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ) dphi -= TMath::TwoPi();
    if ( dphi < -TMath::Pi() ) dphi += TMath::TwoPi();
    return dphi;
}

//________________
bool DiJetAnalysis::isGoodGenJet(const GenJet* jet) {
    bool goodJet{false};
    double etaCut[2] {fJetEtaLab[0], fJetEtaLab[1]}; 
    double eta = jet->eta();

    if ( fSelectJetsInCMFrame ) {

        eta = boostEta2CM( eta );

        etaCut[0] = fJetEtaCM[0];
        etaCut[1] = fJetEtaCM[1];
    }
    else {
        eta = etaLab( eta );
    }

    if ( jet->pt() > 20. && etaCut[0] < eta && eta < etaCut[1] ) {
        goodJet = {true};
    }
    
    if ( fVerbose ) {
        std::cout << Form("Gen jet cut %s\n", goodJet ? "\t[passed]" : "\t[failed]" );
    }

    return goodJet;
}

//________________
double DiJetAnalysis::dijetEtaInFrame(const double& eta1, const double& eta2, bool isCM) {
    double etaDijet = 0.5 * (eta1 + eta2);
    if ( isCM ) {
        etaDijet = boostEta2CM( etaDijet );
    }
    else {
        etaDijet = etaLab( etaDijet );
    }
    return etaDijet;
}

//________________
bool DiJetAnalysis::isGoodTrkMax(const RecoJet* jet) {
    bool goodTrackMax = {true};
    double rawPt = jet->rawPt();
    double trackMaxPt = jet->trackMaxPt();
    if ( TMath::Abs( jet->eta() ) < 2.4 &&
         ( trackMaxPt/rawPt < 0.01 ||
           trackMaxPt/rawPt > 0.98) ) {
        goodTrackMax = {false};
    }

    if ( fVerbose ) {
        std::cout << "TrackMaxPt/rawPt: " << trackMaxPt/rawPt << ( (goodTrackMax) ? " [passed]" : " [failed]" ) 
                  << ( (trackMaxPt/rawPt < 0.01) ? " too low value " : "" ) << ( (trackMaxPt/rawPt > 0.98) ? " too large value " : "" )
                  << std::endl;
    }

    return goodTrackMax;
}

//_________________
bool DiJetAnalysis::isGoodJetId(const RecoJet* jet) {

	bool passJetId = {false};

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

    bool passNHF{false};
    bool passNEF{false};
    bool passNumOfConstituents{true};
    bool passMuonFrac{true};
    bool passChargedFrac{true};
    bool passChargedMult{true};
    bool passChargedEmFrac{true};
    bool passNeutralMult{true};
	
    // Check cuts depending on jet pseudorapdity
    if ( TMath::Abs( eta ) <= 2.7 ) {
        
        passNHF = ( nhf < neutFracCut ) ? true : false;
        passNEF = ( nhf < neutFracCut ) ? true : false;
        passNumOfConstituents = ( numberOfConstituents > 1 ) ? true : false;

        if ( !fIsLooseJetIdCut ) { 
            passMuonFrac = ( muf < 0.8 ) ? true : false; 
        } // if ( !fIsLooseJetIdCut )

        if( TMath::Abs( eta ) <= 2.4 ) {
            passChargedFrac = ( chf > 0 ) ? true : false;
            passChargedMult = ( chargedMult > 0 ) ? true : false;
            passChargedEmFrac = ( cef < chargedEmFracCut ) ? true : false;
        } // if( TMath::Abs( eta ) <= 2.4 )

    } // if ( TMath::Abs( eta ) <= 2.7 )
    else if ( TMath::Abs( eta ) <= 3.0) {

        passNEF = ( nef > 0.01 ) ? true : false;
        passNHF = ( nhf < 0.98 ) ? true : false;
        passNeutralMult = ( neutralMult > 2 ) ? true : false;

    } // else if ( TMath::Abs( eta ) <= 3.0)
    else  {
        passNEF = ( nef < 0.9 ) ? true : false;
        passNeutralMult = (neutralMult > 10 ) ? true : false; // CAUTION: JET MET says it should be >10
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
double DiJetAnalysis::boostEta2CM(const double &eta) {
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::boostEta2CM -- begin" << std::endl;
    }
    double etaCM = eta;

    // Apply lab frame boost to CM
    if ( fCollisionSystem == 0 ) { // pp
        // For pp do nothing. Already in the CM frame
    }

    else if ( fCollisionSystem == 1 ) { // pPb
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
    } 
    else if ( fCollisionSystem == 2 ) { // PbPb
        // For PbPb do nothing. Already in the CM frame
    }
    else {
        // Unknown collision system
        // Do nothing
    }

    if ( fVerbose ) {
        std::cout << Form("eta: %5.2f  ->  etaCM: %5.2f", eta, etaCM) << std::endl;
        std::cout << "DiJetAnalysis::boostEta2CM -- end" << std::endl;
    }
    return etaCM;
}

//________________
double DiJetAnalysis::etaLab(const double &eta) {
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::etaLab -- begin" << std::endl;
    }

    double etaL = eta;
    // Check collision system
    if ( fCollisionSystem == 0 ) { 
        // For pp apply eta shift (to move from CM to lab frame, to match pPb)
        etaL += fEtaShift;
    }
    else if ( fCollisionSystem == 1 ) { 
        // For pPb we already in the lab frame. Just need to properly address
        // beam direction
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
    else if ( fCollisionSystem == 2 ) { 
        // For PbPb apply eta shift (to move from CM to lab frame, to match pPb)
        etaL += fEtaShift;
    }
    else {
        // Unknown collision system
        // Do nothing
    }

    if ( fVerbose ) {
        std::cout << Form("eta: %5.2f  ->  etaLab: %5.2f", eta, etaL) << std::endl;
        std::cout << "DiJetAnalysis::etaLab -- end" << std::endl;
    }
    return etaL;
}
    
//________________
bool DiJetAnalysis::isGoodRecoJet(const RecoJet* jet) {
    bool goodJet{false};
    bool goodKine{false};
    bool hasMatching{false};

    double etaCut[2] {fJetEtaLab[0], fJetEtaLab[1]};

    double eta = jet->eta();

    if ( fSelectJetsInCMFrame ) {

        eta = boostEta2CM( eta );

        etaCut[0] = fJetEtaCM[0]; 
        etaCut[1] = fJetEtaCM[1];
    }
    else {
        eta = etaLab( eta );
    }

    if ( jet->ptJECCorr() > 20 && etaCut[0] < eta && eta < etaCut[1] ) {
        goodKine = {true};
    }

    if ( fIsMc ) {
        if ( jet->hasMatching() ) {
            hasMatching = {true};
        }
    }
    else {
        hasMatching = {true};
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
void DiJetAnalysis::processGenJets(const Event* event, const double &weight) {

    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::processGenJets -- begin" << std::endl;
    }

    fMcReweight = {1.};
    fIsGenDijetLabFound = {false}; 
    fIsGenDijetCMFound = {false};

    // Loop over generated jets and search for leading and subleading jets
    double ptLead{-1.}, ptSubLead{-1.};
    int idLead{-1}, idSubLead{-1};
    GenJetIterator genJetIter;
    int counter{0};
    for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

        double pt = (*genJetIter)->pt();
        double eta = etaLab( (*genJetIter)->eta() );
        double phi = (*genJetIter)->phi();

        if ( fVerbose ) {
            std::cout << "Gen jet #" << counter << " ";
            (*genJetIter)->print();
        }

        counter++;

        // Apply single-jet selection to gen jets
        //if ( !isGoodGenJet( *genJetIter ) ) continue;

        // Find leading and subleading jets
        findLeadSubleadJets( pt, (counter-1), ptLead, ptSubLead, idLead, idSubLead );
        
        // Fill inclusive jet pt
        fHM->hGenInclusiveJetPt->Fill(pt, weight);
        fHM->hGenInclusiveJetPtEta->Fill(eta, pt, weight);
        fHM->hGenInclusiveJetPtEtaPtHat->Fill(eta, pt, event->ptHat(), weight);

        // if ( fVerbose ) {
        //     std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idLead: %d idSubLead: %d\n", 
        //                       ptLead, ptSubLead, idLead, idSubLead);
        // }

        if ( pt > 30. ) {
            fHM->hGenGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), weight);
            fHM->hGenGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), weight );
        }
    } // for ( genJetIter = event->genJetCollection()->begin();

    //
    // Check for gen dijet
    //
    if ( idLead>=0 && idSubLead>=0 ) {

        if ( fVerbose ) {
            std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f idLead: %d idSubLead: %d\n", 
                              ptLead, ptSubLead, idLead, idSubLead);
        }

        GenJet* leadJet = event->genJetCollection()->at( idLead );
        ptLead = leadJet->pt();
        double phiLead = leadJet->phi();

        GenJet* subLeadJet = event->genJetCollection()->at( idSubLead );
        ptSubLead = subLeadJet->pt();
        double phiSubLead = subLeadJet->phi();

        double dijetPt = 0.5 * (ptLead + ptSubLead);
        double dijetDphi = deltaPhi(phiLead, phiSubLead);


        // Specifically for pPb
        double etaLead = etaLab( leadJet->eta() );
        double etaSubLead = etaLab( subLeadJet->eta() );    
        double dijetEta = 0.5 * (etaLead + etaSubLead);
        double dijetDetaCM = 0.5 * ( etaLead - etaSubLead );

        double x_Pb = 2. * dijetPt / fCollisionEnergy * TMath::Exp( -1. * dijetDetaCM ) * TMath::CosH( dijetDetaCM );
        double x_p = 2. * dijetPt / fCollisionEnergy * TMath::Exp( dijetDetaCM ) * TMath::CosH( dijetDetaCM );
        double xPbOverXp = x_Pb / x_p;

        if ( fVerbose ) {
            std::cout << Form("Gen dijet in lab frame: ptLead: %5.2f etaLead: %5.2f phiLead: %5.2f\n", ptLead, etaLead, phiLead);
            std::cout << Form("Gen dijet in lab frame: ptSubLead: %5.2f etaSubLead: %5.2f phiSubLead: %5.2f\n", ptSubLead, etaSubLead, phiSubLead);
            std::cout << Form("Gen dijet in lab frame: dijetPt: %5.2f dijetDphi: %5.2f\n", dijetPt, dijetDphi);
            std::cout << Form("Gen dijet in lab frame: dijetEta: %5.2f dijetDetaCM: %5.2f\n", dijetEta, dijetDetaCM);
            std::cout << Form("Gen dijet in lab frame: x_Pb: %5.2f x_p: %5.2f xPbOverXp: %5.2f\n", x_Pb, x_p, xPbOverXp);
        }

        fHM->hGenInclusiveDijetDetaCM->Fill( dijetDetaCM, 1. );
        fHM->hGenInclusiveDijetDetaCMWeighted->Fill( dijetDetaCM, weight );
        fHM->hGenInclusiveDijetDetaCMPt->Fill( dijetDetaCM, dijetPt, 1. );
        fHM->hGenInclusiveDijetDetaCMPtWeighted->Fill( dijetDetaCM, dijetPt, weight );
        fHM->hGenInclusiveDijetEtaDetaCMPt->Fill( dijetEta, dijetDetaCM, dijetPt, 1. );
        fHM->hGenInclusiveDijetEtaDetaCMPtWeighted->Fill( dijetEta, dijetDetaCM, dijetPt, weight );
        fHM->hGenInclusiveDijetXPb->Fill( x_Pb, 1. );
        fHM->hGenInclusiveDijetXPbWeighted->Fill( x_Pb, weight );
        fHM->hGenInclusiveDijetXp->Fill( x_p, 1. );
        fHM->hGenInclusiveDijetXpWeighted->Fill( x_p, weight );
        fHM->hGenInclusiveDijetXPbOverXp->Fill( xPbOverXp, 1. );
        fHM->hGenInclusiveDijetXPbOverXpWeighted->Fill( xPbOverXp, weight );
        fHM->hGenInclusiveDijetXPbOverXpEta->Fill( xPbOverXp, dijetEta, 1. );
        fHM->hGenInclusiveDijetXPbOverXpEtaWeighted->Fill( xPbOverXp, dijetEta, weight );

        //
        // Lab frame
        //

        bool goodDijetLab = isGoodDijet(ptLead, leadJet->eta(), ptSubLead, subLeadJet->eta(), dijetDphi, false);
        if ( fVerbose ) {
            std::cout << Form("Gen dijet in lab frame is %s\n", ((goodDijetLab) ? "[good]" : "[bad]") ); 
        }

        // Analyze gen dijets in lab frame
        if ( goodDijetLab ) {

            fIsGenDijetLabFound = {true};

            // Flush the eta values to reflect the frame
            double etaLead = etaLab( leadJet->eta() );
            double etaSubLead = etaLab( subLeadJet->eta() );
            double dijetEta = 0.5 * (etaLead + etaSubLead);
            double dijetDetaCM = 0.5 * ( etaLead - etaSubLead );

            double x_Pb = 2. * dijetPt / fCollisionEnergy * TMath::Exp( -1. * dijetDetaCM ) * TMath::CosH( dijetDetaCM );
            double x_p = 2. * dijetPt / fCollisionEnergy * TMath::Exp( dijetDetaCM ) * TMath::CosH( dijetDetaCM );
            double xPbOverXp = x_Pb / x_p;

            fHM->hGenSelectedDijetDetaCM->Fill( dijetDetaCM, 1. );
            fHM->hGenSelectedDijetDetaCMWeighted->Fill( dijetDetaCM, weight );
            fHM->hGenSelectedDijetDetaCMPt->Fill( dijetDetaCM, dijetPt, 1. );
            fHM->hGenSelectedDijetDetaCMPtWeighted->Fill( dijetDetaCM, dijetPt, weight );
            fHM->hGenSelectedDijetEtaDetaCMPt->Fill( dijetEta, dijetDetaCM, dijetPt, 1. );
            fHM->hGenSelectedDijetEtaDetaCMPtWeighted->Fill( dijetEta, dijetDetaCM, dijetPt, weight );
            fHM->hGenSelectedDijetXPb->Fill( x_Pb, 1. );
            fHM->hGenSelectedDijetXPbWeighted->Fill( x_Pb, weight );
            fHM->hGenSelectedDijetXp->Fill( x_p, 1. );
            fHM->hGenSelectedDijetXpWeighted->Fill( x_p, weight );
            fHM->hGenSelectedDijetXPbOverXp->Fill( xPbOverXp, 1. );
            fHM->hGenSelectedDijetXPbOverXpWeighted->Fill( xPbOverXp, weight );
            fHM->hGenSelectedDijetXPbOverXpEta->Fill( xPbOverXp, dijetEta, 1. );
            fHM->hGenSelectedDijetXPbOverXpEtaWeighted->Fill( xPbOverXp, dijetEta, weight );
            
            fHM->hGenPtLeadPtSublead->Fill( ptLead, ptSubLead, weight );
            fHM->hGenEtaLeadEtaSublead->Fill( etaLead, etaSubLead, weight );
            fHM->hGenPtLeadPtSubleadMcReweight->Fill( ptLead, ptSubLead, weight * fMcReweight );
            fHM->hGenEtaLeadEtaSubleadMcReweight->Fill( etaLead, etaSubLead, weight * fMcReweight );

            double genDijetLeadSublead[9] {dijetPt, dijetEta, dijetDphi, 
                                           ptLead, etaLead, phiLead, 
                                           ptSubLead, etaSubLead, phiSubLead };

            fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(genDijetLeadSublead);
            fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(genDijetLeadSublead, weight * fMcReweight );
            fHM->hGenDijetEta->Fill(dijetEta, weight * fMcReweight );
            fHM->hGenDijetPtEtaDphi->Fill(dijetPt, dijetEta, dijetDphi, 1.);
            fHM->hGenDijetPtEtaDphiWeighted->Fill(dijetPt, dijetEta, dijetDphi, weight * fMcReweight );
            (dijetEta >= 0) ? fHM->hGenDijetPtEtaForward->Fill(dijetPt, dijetEta) : fHM->hGenDijetPtEtaBackward->Fill(dijetPt, TMath::Abs(dijetEta));
            (dijetEta >= 0) ? fHM->hGenDijetPtEtaForwardWeighted->Fill(dijetPt, dijetEta, weight * fMcReweight) : fHM->hGenDijetPtEtaBackwardWeighted->Fill(dijetPt, TMath::Abs(dijetEta), weight * fMcReweight);

            // Find exact ptAve bin
            int ptAveBin = findDijetPtAveBin( dijetPt );
            int ptAveOldBin = findDijetPtAveOldBin( dijetPt );

            // New ptAve and eta binning
            if ( ptAveBin >=0 ) {
                fHM->hGenDijetEta1D[ptAveBin]->Fill( dijetEta, 1. );
                fHM->hGenDijetEta1DWeighted[ptAveBin]->Fill( dijetEta, weight * fMcReweight );
                fHM->hGenDijetEtaLeadVsEtaSubLead2D[ptAveBin]->Fill( etaLead, etaSubLead, 1. );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DWeighted[ptAveBin]->Fill( etaLead, etaSubLead, weight * fMcReweight );

                (dijetEta >= 0) ? fHM->hGenDijetEtaForward1D[ptAveBin]->Fill(dijetEta, 1.) : fHM->hGenDijetEtaForward1D[ptAveBin]->Fill(TMath::Abs(dijetEta), 1.);
                (dijetEta >= 0) ? fHM->hGenDijetEtaForward1DWeighted[ptAveBin]->Fill(dijetEta, weight * fMcReweight) : fHM->hGenDijetEtaForward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetEta), weight * fMcReweight);
            }

            // Old ptAve binning
            if ( ptAveOldBin >=0 ) {
                // New eta binning
                fHM->hGenDijetEta1DOldPt[ptAveOldBin]->Fill( dijetEta, 1. );
                fHM->hGenDijetEta1DOldPtWeighted[ptAveOldBin]->Fill( dijetEta, weight * fMcReweight );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPt[ptAveOldBin]->Fill( etaLead, etaSubLead, 1. );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[ptAveOldBin]->Fill( etaLead, etaSubLead, weight * fMcReweight );

                (dijetEta >= 0) ? fHM->hGenDijetEtaForward1DOldPt[ptAveOldBin]->Fill(dijetEta, 1.) : fHM->hGenDijetEtaForward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetEta), 1.);
                (dijetEta >= 0) ? fHM->hGenDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetEta, weight * fMcReweight) : fHM->hGenDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetEta), weight * fMcReweight);

                // Old eta binning
                fHM->hGenDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetEta, 1. );
                fHM->hGenDijetEta1DOldPtBinningWeighted[ptAveOldBin]->Fill( dijetEta, weight * fMcReweight );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[ptAveOldBin]->Fill( etaLead, etaSubLead, 1. );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[ptAveOldBin]->Fill( etaLead, etaSubLead, weight * fMcReweight );
                (dijetEta >= 0) ? fHM->hGenDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(dijetEta, 1.) : fHM->hGenDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetEta), 1.);
                (dijetEta >= 0) ? fHM->hGenDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetEta, weight * fMcReweight) : fHM->hGenDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetEta), weight * fMcReweight);
            }
        } // if ( goodDijetLab )

        //
        // CM frame
        //

        bool goodDijetCM = isGoodDijet(ptLead, leadJet->eta(), ptSubLead, subLeadJet->eta(), dijetDphi, true);
        if ( fVerbose ) {
            std::cout << Form("Gen dijet in CM frame is %s\n", ((goodDijetCM) ? "[good]" : "[bad]") ); 
        }

        // Analyze gen dijets in CM frame
        if ( goodDijetCM ) {

            fIsGenDijetCMFound = {true};

            double etaLead = boostEta2CM( leadJet->eta() );
            double etaSubLead = boostEta2CM( subLeadJet->eta() );
            double dijetEtaCM = 0.5 * (etaLead + etaSubLead);

            fHM->hGenEtaCMLeadEtaCMSublead->Fill( etaLead, etaSubLead, weight );

            fHM->hGenDijetEtaCM->Fill(dijetEtaCM, weight * fMcReweight );
            fHM->hGenDijetPtEtaDphiCM->Fill(dijetPt, dijetEtaCM, dijetDphi, 1.);
            fHM->hGenDijetPtEtaDphiCMWeighted->Fill(dijetPt, dijetEtaCM, dijetDphi, weight * fMcReweight );
            (dijetEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForward->Fill(dijetPt, dijetEtaCM) : fHM->hGenDijetPtEtaCMBackward->Fill(dijetPt, TMath::Abs(dijetEtaCM));
            (dijetEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForwardWeighted->Fill(dijetPt, dijetEtaCM, weight * fMcReweight) : fHM->hGenDijetPtEtaCMBackwardWeighted->Fill(dijetPt, TMath::Abs(dijetEtaCM), weight * fMcReweight);

            // Find exact ptAve bin
            int ptAveBin = findDijetPtAveBin( dijetPt );
            int ptAveOldBin = findDijetPtAveOldBin( dijetPt );

            // New ptAve and eta binning
            if ( ptAveBin >=0 ) {
                fHM->hGenDijetEta1DCM[ptAveBin]->Fill( dijetEtaCM, 1. );
                fHM->hGenDijetEta1DCMWeighted[ptAveBin]->Fill( dijetEtaCM, weight * fMcReweight );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DCM[ptAveBin]->Fill( etaLead, etaSubLead, 1. );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[ptAveBin]->Fill( etaLead, etaSubLead, weight * fMcReweight );
                (dijetEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1D[ptAveBin]->Fill(dijetEtaCM, 1.) : fHM->hGenDijetEtaCMBackward1D[ptAveBin]->Fill(TMath::Abs(dijetEtaCM), 1.);
                (dijetEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DWeighted[ptAveBin]->Fill(dijetEtaCM, weight * fMcReweight) : fHM->hGenDijetEtaCMBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetEtaCM), weight * fMcReweight);
            } // if ( ptAveBin >=0 )

            // Old ptAve binning
            if ( ptAveOldBin >=0 ) {
                fHM->hGenDijetEta1DOldPtCM[ptAveOldBin]->Fill( dijetEtaCM, 1. );
                fHM->hGenDijetEta1DOldPtCMWeighted[ptAveOldBin]->Fill( dijetEtaCM, weight * fMcReweight );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[ptAveOldBin]->Fill( etaLead, etaSubLead, 1. );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[ptAveOldBin]->Fill( etaLead, etaSubLead, weight * fMcReweight );
                (dijetEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPt[ptAveOldBin]->Fill(dijetEtaCM, 1.) : fHM->hGenDijetEtaCMBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetEtaCM), 1.);
                (dijetEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetEtaCM, weight * fMcReweight) : fHM->hGenDijetEtaCMBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetEtaCM), weight * fMcReweight);

                fHM->hGenDijetEta1DOldPtBinningCM[ptAveOldBin]->Fill( dijetEtaCM, 1. );
                fHM->hGenDijetEta1DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetEtaCM, weight * fMcReweight );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[ptAveOldBin]->Fill( etaLead, etaSubLead, 1. );
                fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaLead, etaSubLead, weight * fMcReweight );
                (dijetEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPtBinning[ptAveOldBin]->Fill(dijetEtaCM, 1.) : fHM->hGenDijetEtaCMBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetEtaCM), 1.);
                (dijetEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetEtaCM, weight * fMcReweight) : fHM->hGenDijetEtaCMBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetEtaCM), weight * fMcReweight);
            } // if ( ptAveOldBin >=0 )
        } // if ( goodDijetCM )

    } // if ( idLead>=0 && idSubLead>=0 )

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processGenJets -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::findMcWeight(const double& ptLead, const double& ptSublead) {

    if ( fUseMcReweighting !=0 ) {
        if ( fVerbose ) {
            std::cout << Form("DiJetAnalysis::findMcWeight - ptLead: %5.1f ptSublead: %5.1f\n", ptLead, ptSublead);
        }

        int ptLeadBin{-1}; 
        if ( ptLead >= fJetPtLow && ptLead <= fJetPtHi ) {
            ptLeadBin = ( ptLead - fJetPtLow ) / fJetPtStep;
        }
        int ptSubleadBin{-1};
        if ( ptSublead >= fJetPtLow && ptSublead <= fJetPtHi ) {
            ptSubleadBin = ( ptSublead - fJetPtLow ) / fJetPtStep;
        }
        double val = ( ptLeadBin >=0 && ptSubleadBin >= 0 ) ? 
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
bool DiJetAnalysis::isOverweightedEvent(const double& ptLead, const double& ptHat) {
    return (  ( ptLead / ptHat ) > 1.5);
}

//________________
void DiJetAnalysis::processRecoJets(const Event* event, const double &weight) {

    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::processRecoJets -- begin" << std::endl;
    }

    fMcReweight = {1.};
    fIsRecoDijetLabFound = {false};
    fIsRecoDijetCMFound = {false};

    // ptHat value
    double ptHat = event->ptHat();

    // Loop over reconstructed jets and search for leading and subleading jets
    double ptRecoLead{-1.}, ptRecoSubLead{-1.};
    int idRecoLead{-1}, idRecoSubLead{-1};
    RecoJetIterator recoJetIter;
    int counter{0};

    //
    // Loop over reco jets to search for the leading and subleading jets
    //
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        double pt = (*recoJetIter)->ptJECCorr();
        counter++;

        // Check selection criteria
        bool passTrkMax = isGoodTrkMax( (*recoJetIter) );
        bool passJetId = isGoodJetId( (*recoJetIter) );

        if ( fUseJetIdSelection && !passJetId ) { 
            // Do not forget to increment the counter
            if ( fVerbose ) { 
                std::cout << "JetId selection failed. Skip jet" << std::endl; 
            }
            continue; 
        }
        if ( !fUseJetIdSelection && !passTrkMax ) {
            // Do not forget to increment the counter
            if ( fVerbose ) { 
                std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
            }
            continue; 
        }

        // Find leading and subleading reco jets
        findLeadSubleadJets( pt, (counter-1), ptRecoLead, ptRecoSubLead, idRecoLead, idRecoSubLead );
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    // Checking overweighting
    if ( fIsMc ) {
        // Leading jet pt over ptHat vs leading jet pt
        fHM->hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Fill( ptRecoLead/ptHat, ptRecoLead, 1. );
        fHM->hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Fill( ptRecoLead/ptHat, ptRecoLead, weight );

        if ( isOverweightedEvent( ptRecoLead, ptHat ) ) {
            if ( fVerbose ) {
                std::cout << Form("Overweighted event. ptLead/ptHat = %3.2f", ptRecoLead/ptHat) << std::endl;
            }
            fIsOverweightedEvent = {true};
            return;
        } // if ( isOverweightedEvent( ptRecoLead, ptHat ) )
    } // if ( fIsMc )


    //
    // Inclusive jet analysis
    //
    counter = 0;
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        double pt = (*recoJetIter)->ptJECCorr();
        double eta = etaLab( (*recoJetIter)->eta() );
        double phi = (*recoJetIter)->phi();
        double ptRaw = (*recoJetIter)->pt();

        if ( fVerbose ) {
            std::cout << "Reco jet #" << counter << " ";
            (*recoJetIter)->print();
        }

        counter++;

        // JetId parameters
        int chargedMult = (*recoJetIter)->jtPfCHM() + (*recoJetIter)->jtPfCEM() + (*recoJetIter)->jtPfMUM();
        int neutralMult = (*recoJetIter)->jtPfNHM() + (*recoJetIter)->jtPfNEM();
        int numberOfConstituents = chargedMult + neutralMult;

        int dummyIter{0};
        if ( TMath::Abs( eta ) <= 2.4 ) { dummyIter = {0}; }
        else if ( TMath::Abs( eta ) <= 2.7 ) { dummyIter = {1}; }
        else if ( TMath::Abs( eta ) <= 3.0 ) { dummyIter = {2}; }
        else { dummyIter = {3}; }

        // JetId histograms
        fHM->hRecoInclusiveJetNHF[dummyIter]->Fill( (*recoJetIter)->jtPfNHF(), weight );
        fHM->hRecoInclusiveJetNEmF[dummyIter]->Fill( (*recoJetIter)->jtPfNEF(), weight );
        fHM->hRecoInclusiveJetNumOfConst[dummyIter]->Fill( numberOfConstituents, weight );
        fHM->hRecoInclusiveJetMUF[dummyIter]->Fill( (*recoJetIter)->jtPfMUF(), weight );
        fHM->hRecoInclusiveJetCHF[dummyIter]->Fill( (*recoJetIter)->jtPfCHF(), weight );
        fHM->hRecoInclusiveJetChargedMult[dummyIter]->Fill( chargedMult, weight );
        fHM->hRecoInclusiveJetCEmF[dummyIter]->Fill( (*recoJetIter)->jtPfCEF(), weight );
        fHM->hRecoInclusiveJetNumOfNeutPart[dummyIter]->Fill( neutralMult, weight );

        // Check selection criteria
        bool passTrkMax = isGoodTrkMax( (*recoJetIter) );
        bool passJetId = isGoodJetId( (*recoJetIter) );

        if ( fUseJetIdSelection && !passJetId ) { 
            // Do not forget to increment the counter
            if ( fVerbose ) { 
                std::cout << "JetId selection failed. Skip jet" << std::endl; 
            }
            continue; 
        }
        if ( !fUseJetIdSelection && !passTrkMax ) {
            // Do not forget to increment the counter
            if ( fVerbose ) { 
                std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
            }
            continue; 
        }

        fHM->hRecoInclusiveJetPt->Fill(pt, weight);
        fHM->hRecoInclusiveAllJetPtVsEta->Fill(eta, pt, weight);
        fHM->hRecoInclusiveJetPtRawVsEta->Fill(eta, ptRaw, weight);
        fHM->hRecoInclusiveJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);


        if ( (counter-1) == idRecoLead ) {
            fHM->hRecoLeadJetAllPtVsEta->Fill( eta, pt, weight );
        }

        if (  (counter-1) == idRecoSubLead ) {
            fHM->hRecoSubLeadJetAllPtVsEta->Fill( eta, pt, weight );
        }

        // On MC check reco jet matching to gen
        if ( fIsMc ) {

            // If reco jet has matching to gen jet
            if ( (*recoJetIter)->hasMatching() ) {

                GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
                double genPt = matchedJet->pt();
                double genEta = etaLab( matchedJet->eta() );
                double genPhi = matchedJet->phi();

                double JES = pt/genPt;
                double res[4] { JES, genPt, genEta, genPhi };
                res[0] = JES;
                res[1] = genPt; 
                res[2] = genEta;
                res[3] = genPhi;

                double res1[4] = { JES, genPt, genEta, ptHat };
                double res2[4] = { JES, pt, eta, ptHat };

                fHM->hRefInclusiveJetPt->Fill( genPt, weight );
                fHM->hRefInclusiveJetPtEta->Fill( genEta, genPt, weight );
                fHM->hRefInclusiveJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );
                fHM->hRecoInclusiveMatchedJetPt->Fill( pt, weight );
                fHM->hRecoInclusiveMatchedJetPtVsEta->Fill( eta, pt, weight );

                // Jet energy scale
                fHM->hJESInclusiveJetPtEtaPhi->Fill(res, 1.);
                fHM->hJESInclusiveJetPtEtaPhiWeighted->Fill( res, weight );

                fHM->hRecoInclusiveJetJECFactorVsPtEta->Fill( pt / ptRaw, genPt, genEta, weight );
                // fHM->hRecoInclusiveJetJEC2FactorVsPtGen->Fill( JES2, genPt, genEta, weight );
                fHM->hRecoInclusiveJetPtRawOverPtRefVsPtEta->Fill( ptRaw/genPt, genPt, genEta, weight );
                fHM->hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->Fill( ptRaw/genPt, genPt, genEta, weight );
                fHM->hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning->Fill( ptRaw/genPt, ptRaw, eta, weight );

                // Fill JES vs pt for |eta| < 1.4 (midrapidity)
                if ( TMath::Abs( genEta ) < 1.4 ) {
                    fHM->hInclusiveJetJESVsPtGen->Fill( genPt, JES, weight );
                }
                fHM->hInclusiveJetJESGenPtGenEtaPtHatWeighted->Fill( res1, weight );
                fHM->hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Fill( res2, weight );

                double correl[5] { pt, ptRaw, genPt, eta, genEta };
                fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill( correl );
                fHM->hRecoInclusiveJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill( correl, weight );

                // Leading jet
                if ( (counter-1) == idRecoLead ) {
                    fHM->hRecoLeadJetMatchedPtVsEta->Fill( eta, pt, weight);
                    fHM->hLeadingJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                }

                // Subleading jet
                if ( (counter-1) == idRecoSubLead ) {
                    fHM->hRecoSubLeadJetMatchedPtVsEta->Fill( eta, pt, weight);
                    fHM->hSubleadingJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                }

            } // if ( (*recoJetIter)->hasMatching() )
            else {
                // Fill unmatched jets
                fHM->hRecoInclusiveUnmatchedJetPtVsEta->Fill(eta, pt, weight);

                // Leading jet
                if ( (counter-1) == idRecoLead ) {
                    fHM->hRecoLeadJetUnmatchedPtVsEta->Fill( eta, pt, weight);
                }

                // Subleading jet
                if ( (counter-1) == idRecoSubLead ) {
                    fHM->hRecoSubLeadJetUnmatchedPtVsEta->Fill( eta, pt, weight);
                }
            } // else
        } // if ( fIsMc )

        if ( pt > 30. ) {
            fHM->hRecoGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), weight );
            fHM->hRecoGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), weight );
        }

    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    // If leading and subleading jets are found
    if ( idRecoLead>=0 && idRecoSubLead>=0 ) {

        RecoJet* recoLeadJet = event->recoJetCollection()->at( idRecoLead );
        RecoJet* recoSubLeadJet = event->recoJetCollection()->at( idRecoSubLead );
        GenJet* refLeadJet = {nullptr};
        GenJet* refSubLeadJet = {nullptr};

        double ptRecoLead = recoLeadJet->ptJECCorr();
        double ptRawRecoLead = recoLeadJet->pt();
        double etaRecoLead = recoLeadJet->eta();
        double phiRecoLead = recoLeadJet->phi();

        double ptRecoSubLead = recoSubLeadJet->ptJECCorr();
        double ptRawRecoSubLead = recoSubLeadJet->pt();
        double etaRecoSubLead = recoSubLeadJet->eta();
        double phiRecoSubLead = recoSubLeadJet->phi();

        double dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
        double dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);

        // For dijet analysis (in Monte Carlo) reco leading and subleading jets must have matching MC partners
        if ( fIsMc ) {
            if ( !recoLeadJet->hasMatching() || !recoSubLeadJet->hasMatching() ) {
                if ( fVerbose ) {
                    std::cout << "Reco dijet has unmatched jets. Skip dijet analysis\n";
                }
                return;
            }
        } // if ( fIsMc )


        //
        // Lab frame
        //

        bool goodDijetLab = isGoodDijet(ptRecoLead, recoLeadJet->eta(), ptRecoSubLead, recoSubLeadJet->eta(), TMath::Abs( dijetRecoDphi ), false);
        if ( fVerbose ) {
            std::cout << Form("Reco dijet in lab frame is %s\n", ((goodDijetLab) ? "[good]" : "[bad]") ); 
        }

        // Analyze reco dijets in lab frame
        if ( goodDijetLab ) {

            fIsRecoDijetLabFound = {true};

            // Flush the eta values to reflect the frame
            etaRecoLead = etaLab( recoLeadJet->eta() );
            etaRecoSubLead = etaLab( recoSubLeadJet->eta() );
            double dijetRecoEta = 0.5 * (etaRecoLead + etaRecoSubLead);

            if ( fVerbose ) {
                std::cout << "Reco dijet parameters in the lab frame: " << std::endl;
                std::cout << Form("ptLead: %5.2f ptSubLead: %5.2f phiLead: %5.2f phiSubLead: %5.2f\n", ptRecoLead, ptRecoSubLead, phiRecoLead, phiRecoSubLead);
                std::cout << Form("etaLead: %5.2f -> (lab) %5.2f etaSubLead: %5.2f -> (lab) %5.2f\n", recoLeadJet->eta(), etaRecoLead, recoSubLeadJet->eta(), etaRecoSubLead);
                std::cout << Form("dijet ptAve: %5.2f dijet eta: %5.2f dijet dphi: %5.2f\n", dijetRecoPt, dijetRecoEta, dijetRecoDphi);
            }

            // Correlation between leading and subleading
            fHM->hRecoPtLeadPtSublead->Fill( ptRecoLead, ptRecoSubLead, weight );
            fHM->hRecoEtaLeadEtaSublead->Fill( etaRecoLead, etaRecoSubLead, weight );
            fHM->hRecoPtLeadPtSubleadMcReweight->Fill( ptRecoLead, ptRecoSubLead, weight * fMcReweight );
            fHM->hRecoEtaLeadEtaSubleadMcReweight->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

            double dijetRecoInfo[9] { dijetRecoPt, dijetRecoEta, dijetRecoDphi,
                                      ptRecoLead, etaRecoLead, phiRecoLead,
                                      ptRecoSubLead, etaRecoSubLead, phiRecoSubLead };
            fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(dijetRecoInfo);
            fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(dijetRecoInfo, weight * fMcReweight);
            fHM->hRecoDijetEta->Fill( dijetRecoEta, weight * fMcReweight);
            fHM->hRecoDijetPtEta->Fill( dijetRecoPt, dijetRecoEta, weight * fMcReweight);
            fHM->hRecoDijetPtEtaDphi->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, 1. );
            fHM->hRecoDijetPtEtaDphiWeighted->Fill( dijetRecoPt, dijetRecoEta, dijetRecoDphi, weight * fMcReweight);
            if ( dijetRecoEta >= 0 ) {
                fHM->hRecoDijetPtEtaForward->Fill(dijetRecoPt, dijetRecoEta, 1.);
                fHM->hRecoDijetPtEtaForwardWeighted->Fill(dijetRecoPt, dijetRecoEta, weight * fMcReweight);
            }
            else {
                fHM->hRecoDijetPtEtaBackward->Fill(dijetRecoPt, TMath::Abs(dijetRecoEta), 1.);
                fHM->hRecoDijetPtEtaBackwardWeighted->Fill(dijetRecoPt, TMath::Abs(dijetRecoEta), weight * fMcReweight);
            }

            // Find dijet exact dijet pT bins
            int ptAveBin = findDijetPtAveBin( dijetRecoPt );
            int ptAveOldBin = findDijetPtAveOldBin( dijetRecoPt );

            // New ptAve and eta binning
            if ( ptAveBin >=0 ) {
                fHM->hRecoDijetEta1D[ptAveBin]->Fill( dijetRecoEta, 1. );
                fHM->hRecoDijetEta1DWeighted[ptAveBin]->Fill( dijetRecoEta, weight * fMcReweight );

                fHM->hRecoDijetEtaLeadVsEtaSubLead2D[ptAveBin]->Fill( etaRecoLead, etaRecoSubLead, 1. );
                fHM->hRecoDijetEtaLeadVsEtaSubLead2DWeighted[ptAveBin]->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

                (dijetRecoEta >= 0) ? fHM->hRecoDijetEtaForward1D[ptAveBin]->Fill(dijetRecoEta, 1.) : fHM->hRecoDijetEtaBackward1D[ptAveBin]->Fill(TMath::Abs(dijetRecoEta), 1.);
                (dijetRecoEta >= 0) ? fHM->hRecoDijetEtaForward1DWeighted[ptAveBin]->Fill(dijetRecoEta, weight * fMcReweight) : fHM->hRecoDijetEtaBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetRecoEta), weight * fMcReweight);
            } // if ( ptAveBin >=0 )

            // Old ptAve binning
            if ( ptAveOldBin >=0 ) {

                // New eta binning
                fHM->hRecoDijetEta1DOldPt[ptAveOldBin]->Fill( dijetRecoEta, 1. );
                fHM->hRecoDijetEta1DOldPtWeighted[ptAveOldBin]->Fill( dijetRecoEta, weight * fMcReweight );

                fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPt[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, 1. );
                fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

                (dijetRecoEta >= 0) ? fHM->hRecoDijetEtaForward1DOldPt[ptAveOldBin]->Fill(dijetRecoEta, 1.) : fHM->hRecoDijetEtaBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEta), 1.);
                (dijetRecoEta >= 0) ? fHM->hRecoDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetRecoEta, weight * fMcReweight) : fHM->hRecoDijetEtaBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEta), weight * fMcReweight);

                // Old eta binning
                fHM->hRecoDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetRecoEta, 1. );
                fHM->hRecoDijetEta1DOldPtBinningWeighted[ptAveOldBin]->Fill( dijetRecoEta, weight * fMcReweight );

                fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, 1. );
                fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

                (dijetRecoEta >= 0) ? fHM->hRecoDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(dijetRecoEta, 1.) : fHM->hRecoDijetEtaBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEta), 1.);
                (dijetRecoEta >= 0) ? fHM->hRecoDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetRecoEta, weight * fMcReweight) : fHM->hRecoDijetEtaBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEta), weight * fMcReweight);
            } // if ( ptAveOldBin >=0 )

            // In case of MC
            if ( fIsMc ) {
                refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
                if ( !refLeadJet ) {
                    std::cerr << "Error: Leading jet has no matching gen jet\n";
                    return;
                }
                refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );
                if ( !refSubLeadJet ) {
                    std::cerr << "Error: Subleading jet has no matching gen jet\n";
                    return;
                }

                double ptRefLead = refLeadJet->pt();
                double etaRefLead = etaLab( refLeadJet->eta() );
                double phiRefLead = refLeadJet->phi();

                double ptRefSubLead = refSubLeadJet->pt();
                double etaRefSubLead = etaLab( refSubLeadJet->eta() );
                double phiRefSubLead = refSubLeadJet->phi();

                double dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
                double dijetRefEta = 0.5 * ( etaRefLead + etaRefSubLead );
                double dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);

                if ( fVerbose ) {
                    std::cout << Form("Ref dijet parameters: pt: %5.1f eta: %5.2f dphi: %5.2f\n", dijetRefPt, dijetRefEta, dijetRefDphi);
                }

                // Leading jet information
                double correl[5] { ptRecoLead, ptRawRecoLead, ptRefLead, etaRecoLead, etaRefLead };
                fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, weight * fMcReweight );

                // Subleading jet information
                correl[0] = ptRecoSubLead;
                correl[1] = ptRawRecoSubLead;
                correl[2] = ptRefSubLead;
                correl[3] = etaRecoSubLead; 
                correl[4] = etaRefSubLead;
                fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
                fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, weight * fMcReweight);

                fHM->hRefPtLeadPtSublead->Fill( ptRefLead, ptRefSubLead, weight );
                fHM->hRefEtaLeadEtaSublead->Fill( ptRefLead, ptRefSubLead, weight );
                fHM->hRefPtLeadPtSubleadMcReweight->Fill( ptRefLead, ptRefSubLead, weight * fMcReweight );
                fHM->hRefEtaLeadEtaSubleadMcReweight->Fill( ptRefLead, ptRefSubLead, weight * fMcReweight );

                double dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                               ptRecoLead, etaRecoLead,
                                               ptRecoSubLead, etaRecoSubLead,
                                               dijetRefPt, dijetRefEta,
                                               ptRefLead, etaRefLead,
                                               ptRefSubLead, etaRefSubLead };

                double dijetUnfold[4] = { dijetRecoPt, dijetRecoEta, dijetRefPt, dijetRefEta };

                fHM->hRecoDijetPtEtaRefDijetPtEta->Fill(dijetUnfold, 1.);
                fHM->hRecoDijetPtEtaRefDijetPtEtaWeighted->Fill(dijetUnfold, weight * fMcReweight);

                fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Fill(dijetRecoUnfold);
                fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, weight * fMcReweight );
                fHM->hRefDijetEta->Fill( dijetRefEta, weight * fMcReweight );
                fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEta, dijetRefEta, weight * fMcReweight );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, 1.);
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Fill( dijetRecoEta, dijetRefEta, dijetRecoPt, weight * fMcReweight );
                fHM->hRefDijetPtEtaDphi->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, 1. );
                fHM->hRefDijetPtEtaDphiWeighted->Fill( dijetRefPt, dijetRefEta, dijetRefDphi, weight * fMcReweight );
                (dijetRefEta >= 0) ? fHM->hRefDijetPtEtaForward->Fill(dijetRefPt, dijetRefEta) : fHM->hRefDijetPtEtaBackward->Fill(dijetRefPt, TMath::Abs(dijetRefEta));
                (dijetRefEta >= 0) ? fHM->hRefDijetPtEtaForwardWeighted->Fill(dijetRefPt, dijetRefEta, weight * fMcReweight) : fHM->hRefDijetPtEtaBackwardWeighted->Fill(dijetRefPt, TMath::Abs(dijetRefEta), weight * fMcReweight);         

                // Find exact dijet pT bins
                int ptAveBin = findDijetPtAveBin( dijetRefPt );
                int ptAveOldBin = findDijetPtAveOldBin( dijetRefPt );

                // New ptAve and eta binning
                if ( ptAveBin >=0 ) {
                    fHM->hRefDijetEta1D[ptAveBin]->Fill( dijetRefEta, 1. );
                    fHM->hRefDijetEta1DWeighted[ptAveBin]->Fill( dijetRefEta, weight * fMcReweight );

                    fHM->hRefEtaLeadVsEtaSubLead2D[ptAveBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                    fHM->hRefEtaLeadVsEtaSubLead2DWeighted[ptAveBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                    fHM->hRecoVsRefDijetEta2D[ptAveBin]->Fill( dijetRecoEta, dijetRefEta, 1. );
                    fHM->hRecoVsRefDijetEta2DWeighted[ptAveBin]->Fill( dijetRecoEta, dijetRefEta, weight * fMcReweight );
                    fHM->hRecoVsRefLeadJetEta2D[ptAveBin]->Fill( etaRecoLead, etaRefLead, 1. );
                    fHM->hRecoVsRefLeadJetEta2DWeighted[ptAveBin]->Fill( etaRecoLead, etaRefLead, weight * fMcReweight );
                    fHM->hRecoVsRefSubLeadJetEta2D[ptAveBin]->Fill( etaRecoSubLead, etaRefSubLead, 1. );
                    fHM->hRecoVsRefSubLeadJetEta2DWeighted[ptAveBin]->Fill( etaRecoSubLead, etaRefSubLead, weight * fMcReweight );

                    (dijetRefEta >= 0) ? fHM->hRefDijetEtaForward1D[ptAveBin]->Fill(dijetRefEta, 1.) : fHM->hRefDijetEtaBackward1D[ptAveBin]->Fill(TMath::Abs(dijetRefEta), 1.);
                    (dijetRefEta >= 0) ? fHM->hRefDijetEtaForward1DWeighted[ptAveBin]->Fill(dijetRefEta, weight * fMcReweight) : fHM->hRefDijetEtaBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetRefEta), weight * fMcReweight);
                } // if ( ptAveBin >=0 )

                // Old ptAve binning
                if ( ptAveOldBin >=0 ) {

                    // New eta binning
                    fHM->hRefDijetEta1DOldPt[ptAveOldBin]->Fill( dijetRefEta, 1. );
                    fHM->hRefDijetEta1DOldPt[ptAveOldBin]->Fill( dijetRefEta, weight * fMcReweight );

                    fHM->hRefEtaLeadVsEtaSubLead2DOldPt[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                    fHM->hRefEtaLeadVsEtaSubLead2DOldPt[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                    fHM->hRecoVsRefDijetEta2DOldPt[ptAveOldBin]->Fill( dijetRecoEta, dijetRefEta, 1. );
                    fHM->hRecoVsRefDijetEta2DOldPt[ptAveOldBin]->Fill( dijetRecoEta, dijetRefEta, weight * fMcReweight );
                    fHM->hRecoVsRefLeadJetEta2DOldPt[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, 1. );
                    fHM->hRecoVsRefLeadJetEta2DOldPt[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, weight * fMcReweight );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPt[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, 1. );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPt[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, weight * fMcReweight );

                    (dijetRefEta >= 0) ? fHM->hRefDijetEtaForward1DOldPt[ptAveOldBin]->Fill(dijetRefEta, 1.) : fHM->hRefDijetEtaBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), 1.);
                    (dijetRefEta >= 0) ? fHM->hRefDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetRefEta, weight * fMcReweight) : fHM->hRefDijetEtaBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), weight * fMcReweight);

                    // Old eta binning
                    fHM->hRefDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetRefEta, 1. );
                    fHM->hRefDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetRefEta, weight * fMcReweight );

                    fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinning[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                    fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinning[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                    fHM->hRecoVsRefDijetEta2DOldPtBinning[ptAveOldBin]->Fill( dijetRecoEta, dijetRefEta, 1. );
                    fHM->hRecoVsRefDijetEta2DOldPtBinning[ptAveOldBin]->Fill( dijetRecoEta, dijetRefEta, weight * fMcReweight );
                    fHM->hRecoVsRefLeadJetEta2DOldPtBinning[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, 1. );
                    fHM->hRecoVsRefLeadJetEta2DOldPtBinning[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, weight * fMcReweight );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPtBinning[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, 1. );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPtBinning[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, weight * fMcReweight );

                    (dijetRefEta >= 0) ? fHM->hRefDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(dijetRefEta, 1.) : fHM->hRefDijetEtaBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), 1.);
                    (dijetRefEta >= 0) ? fHM->hRefDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetRefEta, weight * fMcReweight) : fHM->hRefDijetEtaBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), weight * fMcReweight);
                } // if ( ptAveOldBin >=0 )
            } // if ( fIsMc )

            if ( fVerbose ) {
                std::cout << "\t[DONE]\n";
            }
        } // if ( fIsDijetFound )

        //
        // CM frame
        // 

        bool goodDijetCM = isGoodDijet(ptRecoLead, recoLeadJet->eta(), ptRecoSubLead, recoSubLeadJet->eta(), TMath::Abs( dijetRecoDphi ), true);
        if ( fVerbose ) {
            std::cout << Form("Reco dijet in CM frame is %s\n", ((goodDijetCM) ? "[good]" : "[bad]") ); 
        }

        // Analyze reco dijets in CM frame
        if ( goodDijetCM ) {


            fIsRecoDijetCMFound = {true};

            // Flush the eta values to reflect the frame
            double etaRecoLead = boostEta2CM( recoLeadJet->eta() );
            double etaRecoSubLead = boostEta2CM( recoSubLeadJet->eta() );
            // double dijetRecoEta = dijetEtaInFrame(recoLeadJet->eta(), recoSubLeadJet->eta(), false);
            double dijetRecoEtaCM = 0.5 * (etaRecoLead + etaRecoSubLead);

    
            if ( fVerbose ) {
                std::cout << "Reco dijet parameters in the c.m. frame: " << std::endl;
                std::cout << Form("ptLead: %5.2f ptSubLead: %5.2f phiLead: %5.2f phiSubLead: %5.2f\n", ptRecoLead, ptRecoSubLead, phiRecoLead, phiRecoSubLead);
                std::cout << Form("etaLead: %5.2f -> (c.m.) %5.2f etaSubLead: %5.2f -> (c.m.) %5.2f\n", recoLeadJet->eta(), etaRecoLead, recoSubLeadJet->eta(), etaRecoSubLead);
                std::cout << Form("dijet ptAve: %5.2f dijet eta (CM): %5.2f dijet dphi: %5.2f\n", dijetRecoPt, dijetRecoEtaCM, dijetRecoDphi);
            }

            fHM->hRecoEtaCMLeadEtaCMSublead->Fill( etaRecoLead, etaRecoSubLead, weight );
            fHM->hRecoDijetEtaCM->Fill( dijetRecoEtaCM, weight * fMcReweight);
            fHM->hRecoDijetPtEtaDphiCM->Fill( dijetRecoPt, dijetRecoEtaCM, dijetRecoDphi, 1. );
            fHM->hRecoDijetPtEtaDphiCMWeighted->Fill( dijetRecoPt, dijetRecoEtaCM, dijetRecoDphi, weight * fMcReweight);
            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaCMForward->Fill(dijetRecoPt, dijetRecoEtaCM, 1.) : fHM->hRecoDijetPtEtaCMBackward->Fill(dijetRecoPt, TMath::Abs(dijetRecoEtaCM), 1.);
            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaCMForwardWeighted->Fill(dijetRecoPt, dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetPtEtaCMBackwardWeighted->Fill(dijetRecoPt, TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);

            // Find exact dijet pT bin
            int ptAveBin = findDijetPtAveBin( dijetRecoPt );
            int ptAveOldBin = findDijetPtAveOldBin( dijetRecoPt );

            // New ptAve and eta binning
            if ( ptAveBin >=0 ) {
                fHM->hRecoDijetEta1DCM[ptAveBin]->Fill( dijetRecoEtaCM, 1. );
                fHM->hRecoDijetEta1DCMWeighted[ptAveBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

                fHM->hRecoEtaLeadVsEtaSubLead2DCM[ptAveBin]->Fill( etaRecoLead, etaRecoSubLead, 1. );
                fHM->hRecoEtaLeadVsEtaSubLead2DCMWeighted[ptAveBin]->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1D[ptAveBin]->Fill(dijetRecoEtaCM, 1.) : fHM->hRecoDijetEtaCMBackward1D[ptAveBin]->Fill(TMath::Abs(dijetRecoEtaCM), 1.);
                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DWeighted[ptAveBin]->Fill(dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetEtaCMBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);
            } // if ( ptAveBin >=0 )

            // Old ptAve binning
            if ( ptAveOldBin >=0 ) {
                // New eta binning
                fHM->hRecoDijetEta1DOldPtCM[ptAveOldBin]->Fill( dijetRecoEtaCM, 1. );
                fHM->hRecoDijetEta1DOldPtCMWeighted[ptAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

                fHM->hRecoEtaLeadVsEtaSubLead2DOldPtCM[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, 1. );
                fHM->hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPt[ptAveOldBin]->Fill(dijetRecoEtaCM, 1.) : fHM->hRecoDijetEtaCMBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), 1.);
                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetEtaCMBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);

                // Old eta binning
                fHM->hRecoDijetEta1DOldPtBinningCM[ptAveOldBin]->Fill( dijetRecoEtaCM, 1. );
                fHM->hRecoDijetEta1DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

                fHM->hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, 1. );
                fHM->hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaRecoLead, etaRecoSubLead, weight * fMcReweight );

                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPtBinning[ptAveOldBin]->Fill(dijetRecoEtaCM, 1.) : fHM->hRecoDijetEtaCMBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), 1.);
                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetEtaCMBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);
            } // if ( ptAveOldBin >=0 )

            // In case of MC
            if ( fIsMc ) {

                refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
                if ( !refLeadJet ) {
                    std::cerr << "Error: Leading jet has no matching gen jet\n";
                    return;
                }
                refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );
                if ( !refSubLeadJet ) {
                    std::cerr << "Error: Subleading jet has no matching gen jet\n";
                    return;
                }

                double ptRefLead = refLeadJet->pt();
                double etaRefLead = boostEta2CM( refLeadJet->eta() );
                double phiRefLead = refLeadJet->phi();

                double ptRefSubLead = refSubLeadJet->pt();
                double etaRefSubLead = boostEta2CM( refSubLeadJet->eta() );
                double phiRefSubLead = refSubLeadJet->phi();

                double dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
                double dijetRefEtaCM= 0.5 * (etaRefLead + etaRefSubLead);
                double dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);

                fHM->hRefEtaCMLeadEtaCMSublead->Fill( etaRefLead, etaRefSubLead, weight );

                fHM->hRefDijetEtaCM->Fill( dijetRefEtaCM, weight );
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPt, 1.);
                fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPt, weight * fMcReweight );
                fHM->hRefDijetPtEtaDphiCM->Fill( dijetRefPt, dijetRefEtaCM, dijetRefDphi, 1. );
                fHM->hRefDijetPtEtaDphiCMWeighted->Fill( dijetRefPt, dijetRefEtaCM, dijetRefDphi, weight * fMcReweight );
            
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForward->Fill(dijetRefPt, dijetRefEtaCM, 1.) : fHM->hRefDijetPtEtaCMBackward->Fill(dijetRefPt, TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForwardWeighted->Fill(dijetRefPt, dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetPtEtaCMBackwardWeighted->Fill(dijetRefPt, TMath::Abs(dijetRefEtaCM), weight * fMcReweight);

                // Find exact dijet ptAve bin
                int ptAveBin = findDijetPtAveBin( dijetRefPt );
                int ptAveOldBin = findDijetPtAveOldBin( dijetRefPt );

                // New ptAve and eta binning
                if ( ptAveBin >=0 ) {
                    fHM->hRefDijetEta1DCM[ptAveBin]->Fill( dijetRefEtaCM, 1. );
                    fHM->hRefDijetEta1DCMWeighted[ptAveBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                    fHM->hRefEtaLeadVsEtaSubLead2DCM[ptAveBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                    fHM->hRefEtaLeadVsEtaSubLead2DCMWeighted[ptAveBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                    fHM->hRecoVsRefDijetEta2DCM[ptAveBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, 1. );
                    fHM->hRecoVsRefDijetEta2DCMWeighted[ptAveBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, weight * fMcReweight );
                    fHM->hRecoVsRefLeadJetEta2DCM[ptAveBin]->Fill( etaRecoLead, etaRefLead, 1. );
                    fHM->hRecoVsRefLeadJetEta2DCMWeighted[ptAveBin]->Fill( etaRecoLead, etaRefLead, weight * fMcReweight );
                    fHM->hRecoVsRefSubLeadJetEta2DCM[ptAveBin]->Fill( etaRecoSubLead, etaRefSubLead, 1. );
                    fHM->hRecoVsRefSubLeadJetEta2DCMWeighted[ptAveBin]->Fill( etaRecoSubLead, etaRefSubLead, weight * fMcReweight );

                    (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1D[ptAveBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefDijetEtaCMBackward1D[ptAveBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                    (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DWeighted[ptAveBin]->Fill(dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetEtaCMBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight);
                }

                // Old ptAve binning
                if ( ptAveOldBin >=0 ) {
                    // New eta binning
                    fHM->hRefDijetEta1DOldPtCM[ptAveOldBin]->Fill( dijetRefEtaCM, 1. );
                    fHM->hRefDijetEta1DOldPtCMWeighted[ptAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                    fHM->hRefEtaLeadVsEtaSubLead2DOldPtCM[ptAveBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                    fHM->hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[ptAveBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                    fHM->hRecoVsRefDijetEta2DOldPtCM[ptAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, 1. );
                    fHM->hRecoVsRefDijetEta2DOldPtCMWeighted[ptAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, weight * fMcReweight );
                    fHM->hRecoVsRefLeadJetEta2DOldPtCM[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, 1. );
                    fHM->hRecoVsRefLeadJetEta2DOldPtCMWeighted[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, weight * fMcReweight );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPtCM[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, 1. );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, weight * fMcReweight );

                    (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPt[ptAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefDijetEtaCMBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                    (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetEtaCMBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight);

                    // Old eta binning
                    fHM->hRefDijetEta1DOldPtBinningCM[ptAveOldBin]->Fill( dijetRefEtaCM, 1. );
                    fHM->hRefDijetEta1DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                    fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                    fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                    fHM->hRecoVsRefDijetEta2DOldPtBinningCM[ptAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, 1. );
                    fHM->hRecoVsRefDijetEta2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, weight * fMcReweight );
                    fHM->hRecoVsRefLeadJetEta2DOldPtBinningCM[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, 1. );
                    fHM->hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaRecoLead, etaRefLead, weight * fMcReweight );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPtBinningCM[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, 1. );
                    fHM->hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaRecoSubLead, etaRefSubLead, weight * fMcReweight );

                    (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPtBinning[ptAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefDijetEtaCMBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                    (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetEtaCMBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight);
                } // if ( ptAveOldBin >=0 )
            } // if ( fIsMc )
        } // if ( goodDijetCM )
    } // if ( idRecoLead>=0 && idRecoSubLead>=0 )

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processRecoJets -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRefJets(const Event* event, const double &weight) {
    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::processRefJets -- begin" << std::endl;
    }

    fMcReweight = {1.};
    fIsRefSelDijetLabFound = {false};
    fIsRefSelDijetCMFound = {false};

    // ptHat value
    double ptHat = event->ptHat();

    // Loop over reconstructed jets and select leading and subleading jets
    // using reference jets
    RecoJetIterator recoJetIter;
    int counter{0};
    double ptRefLead{-100.}, ptRefSubLead{-100.};
    int idRecoLead{-1}, idRecoSubLead{-1};
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        counter++;
        // Analyze only jets that have matched partner
        if ( !(*recoJetIter)->hasMatching() ) {
            // Do not forget to increment the counter
            continue; 
        }

        // Check selection criteria
        bool passTrkMax = isGoodTrkMax( (*recoJetIter) );
        bool passJetId = isGoodJetId( (*recoJetIter) );

        if ( fUseJetIdSelection && !passJetId ) { 
            // Do not forget to increment the counter
            if ( fVerbose ) { 
                std::cout << "JetId selection failed. Skip jet" << std::endl; 
            }
            continue; 
        }
        if ( !fUseJetIdSelection && !passTrkMax ) {
            // Do not forget to increment the counter
            if ( fVerbose ) { 
                std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
            }
            continue; 
        }

        // Retrieve matched gen jet
        GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
        double genPt = matchedJet->pt();
        double genEta = etaLab( matchedJet->eta() );
        double genPhi = matchedJet->phi();

        if ( fVerbose ) {
            std::cout << "Ref jet info for reco jet #" << counter-1;
            matchedJet->print();
            std::cout << "Reco jet #" << counter-1 << " ";
            (*recoJetIter)->print();
        }

        fHM->hRefSelInclusiveJetPt->Fill( genPt, weight * fMcReweight );
        fHM->hRefSelInclusiveJetPtEta->Fill(genEta, genPt, weight * fMcReweight);
        fHM->hRefSelInclusiveJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight * fMcReweight);

        // Find leading and subleading ref jets, but store id of reco jets
        findLeadSubleadJets( genPt, (counter-1), ptRefLead, ptRefSubLead, idRecoLead, idRecoSubLead );
        
        if ( fVerbose ) {
            std::cout << Form("ref lead pT: %5.2f refSubLead pT: %5.2f idRecoLead: %d idRecoSubLead: %d genId: %d \n", 
                              ptRefLead, ptRefSubLead, idRecoLead, idRecoSubLead, (*recoJetIter)->genJetId());
        }

    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    //
    // Check if leading and subleading jets were found
    //
    if (idRecoLead>=0 && idRecoSubLead>=0) {

        // Retrieve leading and subleading jets
        RecoJet *recoLeadJet = event->recoJetCollection()->at( idRecoLead );
        RecoJet *recoSubLeadJet = event->recoJetCollection()->at( idRecoSubLead );
        GenJet* refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
        GenJet* refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );

        if ( fVerbose ) {
            std::cout << "Leading ref-selected jet parameters: ";
            refLeadJet->print();
            std::cout << "Corresponding leading reco jet parameters: ";
            recoLeadJet->print();
            std::cout << "Subleading ref-selected jet parameters: ";
            refSubLeadJet->print();
            std::cout << "Corresponding subleading reco jet parameters: ";
            recoSubLeadJet->print();
        }

        // Retrieve kinematic information for reference leading and subleading jets, dijets
        double ptRefLead = refLeadJet->pt();
        double etaRefLead = refLeadJet->eta();
        double phiRefLead = refLeadJet->phi();

        double ptRefSubLead = refSubLeadJet->pt();
        double etaRefSubLead = refSubLeadJet->eta();
        double phiRefSubLead = refSubLeadJet->phi();

        double dijetRefPt = 0.5 * (ptRefLead + ptRefSubLead);
        double dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);

        //
        // Lab frame
        //

        bool goodDijetLab = isGoodDijet(ptRefLead, refLeadJet->eta(), ptRefSubLead, refSubLeadJet->eta(), TMath::Abs( dijetRefDphi ), false);
        if ( fVerbose ) {
            std::cout << Form("Ref dijet in lab frame is %s\n", ((goodDijetLab) ? "[good]" : "[bad]") ); 
        }

        // Analyze ref-selected dijets in lab frame
        if ( goodDijetLab ) {

            fIsRefSelDijetLabFound = {true};

            // Flush the eta values of reference jets to reflect the frame
            double etaRefLead = etaLab( refLeadJet->eta() );
            double etaRefSubLead = etaLab( refSubLeadJet->eta() );
            double dijetRefEta = 0.5 * (etaRefLead + etaRefSubLead);

            // Set parameters for the reco partners
            double ptRecoLead = recoLeadJet->ptJECCorr();
            double etaRecoLead = etaLab( recoLeadJet->eta() );
            double phiRecoLead = recoLeadJet->phi();

            double ptRecoSubLead = recoSubLeadJet->ptJECCorr();
            double etaRecoSubLead = etaLab( recoSubLeadJet->eta() );
            double phiRecoSubLead = recoSubLeadJet->phi();

            double dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
            double dijetRecoEta = 0.5 * (etaRecoLead + etaRecoSubLead);
            double dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);

            // Dijet reco vs ref for unfolding
            double dijetRecoUnfold[12] = { dijetRecoPt, dijetRecoEta,
                                           ptRecoLead, etaRecoLead,
                                           ptRecoSubLead, etaRecoSubLead,
                                           dijetRefPt, dijetRefEta,
                                           ptRefLead, etaRefLead,
                                           ptRefSubLead, etaRefSubLead };    

            fHM->hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, weight * fMcReweight );
            fHM->hRefSelDijetEta->Fill(dijetRefEta, weight * fMcReweight );
            fHM->hRefSelDijetPtEtaDphi->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, 1.);
            fHM->hRefSelDijetPtEtaDphiWeighted->Fill(dijetRefPt, dijetRefEta, dijetRefDphi, weight * fMcReweight );

            // Find exact dijet pT bins
            int ptAveBin = findDijetPtAveBin( dijetRefPt );
            int ptAveOldBin = findDijetPtAveOldBin( dijetRefPt );

            // New ptAve and eta binning
            if ( ptAveBin >= 0 ) {
                fHM->hRefSelDijetEta1D[ptAveBin]->Fill( dijetRefEta, 1. );
                fHM->hRefSelDijetEta1DWeighted[ptAveBin]->Fill( dijetRefEta, weight * fMcReweight );

                fHM->hRefSelRecoDijetEta1D[ptAveBin]->Fill( dijetRecoEta, 1. );
                fHM->hRefSelRecoDijetEta1DWeighted[ptAveBin]->Fill( dijetRecoEta, weight * fMcReweight );

                fHM->hRefSelEtaLeadVsEtaSubLead2D[ptAveBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                fHM->hRefSelEtaLeadVsEtaSubLead2DWeighted[ptAveBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                (dijetRefEta >= 0) ? fHM->hRefSelDijetEtaForward1D[ptAveBin]->Fill(dijetRefEta, 1.) : fHM->hRefSelDijetEtaBackward1D[ptAveBin]->Fill(TMath::Abs(dijetRefEta), 1.);
                (dijetRefEta >= 0) ? fHM->hRefSelDijetEtaForward1DWeighted[ptAveBin]->Fill(dijetRefEta, weight * fMcReweight ) : fHM->hRefSelDijetEtaBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetRefEta), weight * fMcReweight );
            } // if ( ptAveBin >= 0 )

            // Old ptAve binning
            if ( ptAveOldBin >= 0 ) {
                // New eta binning
                fHM->hRefSelDijetEta1DOldPt[ptAveOldBin]->Fill( dijetRefEta, 1. );
                fHM->hRefSelDijetEta1DOldPtWeighted[ptAveOldBin]->Fill( dijetRefEta, weight * fMcReweight );

                fHM->hRefSelRecoDijetEta1DOldPt[ptAveOldBin]->Fill( dijetRecoEta, 1. );
                fHM->hRefSelRecoDijetEta1DOldPtWeighted[ptAveOldBin]->Fill( dijetRecoEta, weight * fMcReweight );

                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPt[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                (dijetRefEta >= 0) ? fHM->hRefSelDijetEtaForward1DOldPt[ptAveOldBin]->Fill(dijetRefEta, 1.) : fHM->hRefSelDijetEtaBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), 1.);
                (dijetRefEta >= 0) ? fHM->hRefSelDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetRefEta, weight * fMcReweight ) : fHM->hRefSelDijetEtaBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), weight * fMcReweight );

                // Old eta binning
                fHM->hRefSelDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetRefEta, 1. );
                fHM->hRefSelDijetEta1DOldPtBinningWeighted[ptAveOldBin]->Fill( dijetRefEta, weight * fMcReweight );

                fHM->hRefSelRecoDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetRecoEta, 1. );
                fHM->hRefSelRecoDijetEta1DOldPtBinningWeighted[ptAveOldBin]->Fill( dijetRecoEta, weight * fMcReweight );

                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                (dijetRefEta >= 0) ? fHM->hRefSelDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(dijetRefEta, 1.) : fHM->hRefSelDijetEtaBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), 1.);
                (dijetRefEta >= 0) ? fHM->hRefSelDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetRefEta, weight * fMcReweight ) : fHM->hRefSelDijetEtaBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEta), weight * fMcReweight );
            } // if ( ptAveOldBin >= 0 )

        } // if ( goodDijetLab )

        //
        // CM frame
        //

        bool goodDijetCM = isGoodDijet(ptRefLead, refLeadJet->eta(), ptRefSubLead, refSubLeadJet->eta(), TMath::Abs( dijetRefDphi ), true);
        if ( fVerbose ) {
            std::cout << Form("Ref dijet in CM frame is %s\n", ((goodDijetCM) ? "[good]" : "[bad]") ); 
        }

        // Analyze ref-selected dijets in CM frame
        if ( goodDijetCM ) {

            fIsRefSelDijetCMFound = {true};

            // Flush the eta values of reference jets to reflect the frame
            double etaRefLead = boostEta2CM( refLeadJet->eta() );
            double etaRefSubLead = boostEta2CM( refSubLeadJet->eta() );
            double dijetRefEtaCM = 0.5 * (etaRefLead + etaRefSubLead);

            // Set parameters for the reco partners
            double ptRecoLead = recoLeadJet->ptJECCorr();
            double etaRecoLead = boostEta2CM( recoLeadJet->eta() );
            double phiRecoLead = recoLeadJet->phi();

            double ptRecoSubLead = recoSubLeadJet->ptJECCorr();
            double etaRecoSubLead = boostEta2CM( recoSubLeadJet->eta() );
            double phiRecoSubLead = recoSubLeadJet->phi();

            double dijetRecoPt = 0.5 * (ptRecoLead + ptRecoSubLead);
            double dijetRecoEtaCM = 0.5 * (etaRecoLead + etaRecoSubLead);
            double dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);

            fHM->hRefSelDijetEtaCM->Fill(dijetRefEtaCM, weight * fMcReweight );
            fHM->hRefSelDijetPtEtaDphiCM->Fill(dijetRefPt, dijetRefEtaCM, dijetRefDphi, 1.);
            fHM->hRefSelDijetPtEtaDphiCMWeighted->Fill(dijetRefPt, dijetRefEtaCM, dijetRefDphi, weight * fMcReweight );

            // Find exact dijet pT bin
            int ptAveBin = findDijetPtAveBin( dijetRefPt );
            int ptAveOldBin = findDijetPtAveOldBin( dijetRefPt );

            // New ptAve and eta binning
            if ( ptAveBin >= 0 ) {
                fHM->hRefSelDijetEta1DCM[ptAveBin]->Fill( dijetRefEtaCM, 1. );
                fHM->hRefSelDijetEta1DCMWeighted[ptAveBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                fHM->hRefSelRecoDijetEta1DCM[ptAveBin]->Fill( dijetRecoEtaCM, 1. );
                fHM->hRefSelRecoDijetEta1DCMWeighted[ptAveBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

                fHM->hRefSelEtaLeadVsEtaSubLead2DCM[ptAveBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                fHM->hRefSelEtaLeadVsEtaSubLead2DCMWeighted[ptAveBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1D[ptAveBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefSelDijetEtaCMBackward1D[ptAveBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DWeighted[ptAveBin]->Fill(dijetRefEtaCM, weight * fMcReweight ) : fHM->hRefSelDijetEtaCMBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight );
            } // if ( ptAveBin >= 0 )

            // Old ptAve binning
            if ( ptAveOldBin >= 0 ) {
                // New eta binning
                fHM->hRefSelDijetEta1DOldPtCM[ptAveOldBin]->Fill( dijetRefEtaCM, 1. );
                fHM->hRefSelDijetEta1DOldPtCMWeighted[ptAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                fHM->hRefSelRecoDijetEta1DOldPtCM[ptAveOldBin]->Fill( dijetRecoEtaCM, 1. );
                fHM->hRefSelRecoDijetEta1DOldPtCMWeighted[ptAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtCM[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPt[ptAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefSelDijetEtaCMBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight ) : fHM->hRefSelDijetEtaCMBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight );

                // Old eta binning
                fHM->hRefSelDijetEta1DOldPtBinningCM[ptAveOldBin]->Fill( dijetRefEtaCM, 1. );
                fHM->hRefSelDijetEta1DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                fHM->hRefSelRecoDijetEta1DOldPtBinningCM[ptAveOldBin]->Fill( dijetRecoEtaCM, 1. );
                fHM->hRefSelRecoDijetEta1DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, 1. );
                fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaRefLead, etaRefSubLead, weight * fMcReweight );

                (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPtBinning[ptAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefSelDijetEtaCMBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight ) : fHM->hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight );
            } // if ( ptAveOldBin >= 0 )
        }

    } // if (idRecoLead>=0 && idRecoSubLead>=0)

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processRefJets -- end" << std::endl;
    }
}

//________________
bool DiJetAnalysis::isGoodDijet(const double& ptLead, const double& etaLead, 
                                const double& ptSubLead, const double& etaSubLead, 
                                const double& dphi, const bool& isCM) {

    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::isGoodDijet -- begin" << std::endl;
    }

    double eta1 = ( isCM ) ? boostEta2CM(etaLead) : etaLab(etaLead);
    double eta2 = ( isCM ) ? boostEta2CM(etaSubLead) : etaLab(etaSubLead);
    double etaCut[2] = { fJetEtaLab[0], fJetEtaLab[1] };
    if ( isCM ) {
        etaCut[0] = fJetEtaCM[0];
        etaCut[1] = fJetEtaCM[1];
    }
    bool isGood = ( ptLead > fLeadJetPtLow &&
                    etaCut[0] <= eta1 && eta1 < etaCut[1] &&
                    ptSubLead > fSubleadJetPtLow && 
                    etaCut[0] <= eta2 && eta2 < etaCut[1] &&
                    TMath::Abs( dphi ) > fDijetPhiCut );

    // Check reweight
    if ( fIsMc && fUseMcReweighting != 0 ) {
        findMcWeight(ptLead, ptSubLead);
    }
    else {
        fMcReweight = {1.};
    }

    if ( fVerbose ) {
        std::cout << Form("Dijet status: %s\n", ( (isGood) ? "[good]" : "[bad]" ) );
        std::cout << Form("Leading jet pT %5.2f > %5.2f GeV: \t%s\n", ptLead, fLeadJetPtLow, ( (ptLead > fLeadJetPtLow) ? "[good]" : "[bad]" ) );
        std::cout << Form("Subleading jet pT %5.2f > %5.2f GeV: \t%s\n", ptSubLead, fSubleadJetPtLow, ( (ptSubLead > fSubleadJetPtLow) ? "[good]" : "[bad]" ) );
        std::cout << Form("Leading jet eta %3.2f <= %3.2f < %3.2f: \t%s\n", etaCut[0], eta1, etaCut[1], ( (etaCut[0] <= eta1 && eta1 < etaCut[1]) ? "[good]" : "[bad]" ) );
        std::cout << Form("Subleading jet eta %3.2f <= %3.2f < %3.2f: \t%s\n", etaCut[0], eta2, etaCut[1], ( (etaCut[0] <= eta2 && eta2 < etaCut[1]) ? "[good]" : "[bad]" ) );
        std::cout << Form("Delta phi %3.2f > %3.2f: \t%s\n", TMath::Abs( dphi ), fDijetPhiCut, ( (TMath::Abs( dphi ) > fDijetPhiCut) ? "[good]" : "[bad]" ) );
        std::cout << Form("Reweighting factor: %5.2f\n", fMcReweight);
        std::cout << "DiJetAnalysis::isGoodDijet -- end" << std::endl;
    }
    return isGood;
}


//________________
void DiJetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    if ( fVerbose ) {
        std::cout << "\n\n++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "++++++++++++++++++++++++++++++++++++++++" << std::endl;
        std::cout << "DiJetAnalysis::processEvent -- begin" << std::endl;
    }

    if ( !fHM ) {
        std::cout << "[Warning] No histogram manager connected to the DiJetAnalysis\n";
    }

    // Must be flushed for each event !!!!
    fIsGenDijetLabFound = {false};
    fIsGenDijetCMFound = {false};
    fIsRecoDijetLabFound = {false};
    fIsRecoDijetCMFound = {false};
    fIsRefSelDijetLabFound = {false};
    fIsRefSelDijetCMFound = {false};
    fIsOverweightedEvent = {false};

    //
    // Event quantities
    //

    // ptHat
    double ptHat = event->ptHat();
    // Vertex z position
    double vz = event->vz();
    // ptHat weight 
    double ptHatW = event->ptHatWeight();
    // Centrality
    double centrality = event->centrality();
    // Centrality weight
    double centW = event->centralityWeight();
    if ( fCollisionSystem != 2 ) { // Apply centrality weight only for PbPb
        centW = 1.;
    }
    // Final weight
    double weight{1.};

    // Check correctness of MC sample for pPb 8160
    if ( fCollisionSystem == 1 ) { 
        if ( fIsMc ) {
            // Skip events with ptHat that is outside the ranged embedded
            if ( ptHat <= fPtHatRange[0] || ptHat > fPtHatRange[1] ) {
                if ( fVerbose ) {
                    std::cout << Form("[WARNING] Bad ptHat value: %4.1f < %4.1f <= %4.1f\n", fPtHatRange[0], ptHat, fPtHatRange[1]);
                }
                return;
            }

            // For MC we need to flip the direction of Pb-going in order to properly reweight distributions
            if ( fIsPbGoingDir ) {
                vz = -vz;
            }
        }
        else {
            if ( !fIsPbGoingDir ) {
                vz = -vz;
            }
        }
    } // if ( fCollisionSystem == 1 )

    // Calculate event weight
    weight = eventWeight(ptHat, vz, centW, ptHatW);

    // Process and analyze reco jets
    fHM->hRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), 1. );
    processRecoJets(event, weight);

    if (fIsOverweightedEvent) {
        return;
    }

    // Fill event histograms
    fHM->hVz->Fill( vz,  1. );
    fHM->hVzWeighted->Fill( vz, weight );

    fHM->hPtHat->Fill( ptHat, 1. );
    fHM->hPtHatWeighted->Fill( ptHat, weight );

    fHM->hHiBin->Fill( event->hiBin(), 1. );
    fHM->hHiBinWeighted->Fill( event->hiBin(), weight );

    if ( fIsMc ) {
        fHM->hGenJetCollectionSize->Fill( event->genJetCollection()->size(), 1. );
        fHM->hGenVsRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), event->genJetCollection()->size(), 1. );
        // Process and analyze gen jets
        processGenJets(event, weight);
        processRefJets(event, weight);
    }

    if ( fIsGenDijetLabFound ) {
        fHM->hVzGenDijetLab->Fill( vz, 1. );
        fHM->hVzGenDijetLabWeighted->Fill( vz, weight );
        fHM->hHiBinGenDijetLab->Fill( event->hiBin(), 1. );
        fHM->hHiBinGenDijetLabWeighted->Fill( event->hiBin(), weight );
    }

    if ( fIsGenDijetCMFound ) {
        fHM->hVzGenDijetCM->Fill( vz, 1. );
        fHM->hVzGenDijetCMWeighted->Fill( vz, weight );
        fHM->hHiBinGenDijetCM->Fill( event->hiBin(), 1. );
        fHM->hHiBinGenDijetCMWeighted->Fill( event->hiBin(), weight );
    }

    if ( fIsRecoDijetLabFound ) {
        fHM->hVzRecoDijetLab->Fill( vz, 1. );
        fHM->hVzRecoDijetLabWeighted->Fill( vz, weight );
        fHM->hHiBinRecoDijetLab->Fill( event->hiBin(), 1. );
        fHM->hHiBinRecoDijetLabWeighted->Fill( event->hiBin(), weight );
    }

    if ( fIsRecoDijetCMFound ) {
        fHM->hVzRecoDijetCM->Fill( vz, 1. );
        fHM->hVzRecoDijetCMWeighted->Fill( vz, weight );
        fHM->hHiBinRecoDijetCM->Fill( event->hiBin(), 1. );
        fHM->hHiBinRecoDijetCMWeighted->Fill( event->hiBin(), weight );
    }

    if ( fIsRefSelDijetLabFound ) {
        fHM->hVzRefSelDijetLab->Fill( vz, 1. );
        fHM->hVzRefSelDijetLabWeighted->Fill( vz, weight );
        fHM->hHiBinRefSelDijetLab->Fill( event->hiBin(), 1. );
        fHM->hHiBinRefSelDijetLabWeighted->Fill( event->hiBin(), weight );
    }

    if ( fIsRefSelDijetCMFound ) {
        fHM->hVzRefSelDijetCM->Fill( vz, 1. );
        fHM->hVzRefSelDijetCMWeighted->Fill( vz, weight );
        fHM->hHiBinRefSelDijetCM->Fill( event->hiBin(), 1. );
        fHM->hHiBinRefSelDijetCMWeighted->Fill( event->hiBin(), weight );
    }

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processEvent -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::finish() {
    // Save data and close files
    // fCycleCounter++;
    // std::cout << Form("DiJetAnalysis::processEvent [INFO] Total events processed: %d Sample fraction: %3.2f%%", 
    //                   (fCycleCounter * 50000) + fEventCounter, 
    //                   (double)(fCycleCounter * 50000 + fEventCounter) / fNEventsInSample )
    //           << std::endl;
    // std::cout << Form("DiJetAnalysis::processEvent [INFO]: Total number of events processed: %d", fTotalCounter) << std::endl;
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
