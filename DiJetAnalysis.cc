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
#include <algorithm>

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
    fMcReweight{1}, 
    fRecoIdLead{-1}, fRecoIdSubLead{-1}, fGenIdLead{-1}, fGenIdSubLead{-1}, fRefSelRecoIdLead{-1}, fRefSelRecoIdSubLead{-1},
    fRecoPtSortedJetIds{}, fGenPtSortedJetIds{}, fRefSelRecoPtSortedJetIds{},
    fPtAveBins{}, fPtAveOldBins{} {

    fJetEtaLab[0] = -3.; fJetEtaLab[1] = 3.;
    fJetEtaCM[0] = -2.5; fJetEtaCM[1] = 2.5;
    fPtHatRange[0] = {0};
    fPtHatRange[1] = {100000000};
    for (int i=0; i<fJetPtBins; i++) {
        for (int j=0; j<fJetPtBins; j++) {
            fJetPtLeadPtSubleadReweightMatrix[i][j] = 1;
        }
    } // for (int i=0; i<fJetPtBins; i++)

    float dijetPtVals[17] {  50.,  60.,   70.,  80.,  90.,
                              100., 110.,  120., 130., 140.,
                              150., 160.,  180., 200., 250., 
                              300., 500.};
    int sizeOfPtVals = sizeof(dijetPtVals)/sizeof(dijetPtVals[0]);

    fPtAveBins.assign(dijetPtVals, dijetPtVals + sizeOfPtVals);

    float dijetPtOldVals[7] {25., 55., 75., 95., 115., 150., 400.};
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

            float nCorr[75][75] = {
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
int DiJetAnalysis::findDijetPtAveBin(const float &ptAve) {
    // if ( fVerbose ) {
    //     std::cout << "\nDiJetAnalysis::findDijetPtAveBin -- begin" << std::endl;
    // }
    int bin{-1};
    if ( fPtAveBins[0] < ptAve && ptAve < fPtAveBins.at( fPtAveBins.size()-1 ) ) {
        for (unsigned int i=0; i<fPtAveBins.size()-1; i++) {
            if ( fPtAveBins[i] <= ptAve && ptAve < fPtAveBins[i+1] ) {
                bin = i;
                break;
            }
        }
    }

    // if ( fVerbose ) {
    //     std::cout << Form("ptAve: %5.2f bin: %d\n", ptAve, bin);
    //     std::cout << "DiJetAnalysis::findDijetPtAveBin -- end" << std::endl;
    // }
    return bin;
}

//________________
int DiJetAnalysis::findDijetPtAveOldBin(const float &ptAve) {
    // if ( fVerbose ) {
    //     std::cout << "DiJetAnalysis::findDijetPtAveOldBin -- begin" << std::endl;
    // }
    int bin{-1};
    if ( fPtAveOldBins[0] < ptAve && ptAve < fPtAveOldBins.at( fPtAveOldBins.size()-1 ) ) {
        for (unsigned int i=0; i<fPtAveOldBins.size()-1; i++) {
            if ( fPtAveOldBins[i] <= ptAve && ptAve < fPtAveOldBins[i+1] ) {
                bin = i;
                break;
            }
        }
    }

    // if ( fVerbose ) {
    //     std::cout << Form("ptAve: %5.2f bin: %d\n", ptAve, bin);
    //     std::cout << "DiJetAnalysis::findDijetPtAveOldBin -- end" << std::endl;
    // }
    return bin;
}

//________________
double DiJetAnalysis::eventWeight(const float& ptHat, const float& vz, 
                                  const float& centWeight, const float& ptHatW) {

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
void DiJetAnalysis::findLeadSubleadJets(const float &pt, const int &counter,
                                        float &ptLead, float &ptSublead,
                                        int &idLead, int &idSubLead) {
    // Find leading and subleading jets
    // if ( fVerbose ) {
    //     std::cout << "DiJetAnalysis::findLeadSubleadJets -- begin" << std::endl;
    // }

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

    // if ( fVerbose ) {
    //     std::cout << Form("Lead pT: %5.2f SubLead pT: %5.2f Lead id: %d SubLead id: %d\n", 
    //                       ptLead, ptSublead, idLead, idSubLead);
    //     std::cout << "DiJetAnalysis::findLeadSubleadJets - end" << std::endl;
    // }
}

//________________
float DiJetAnalysis::deltaPhi(const float& phi1, const float &phi2) {
    float dphi = phi1 - phi2;
    if ( dphi > TMath::Pi() ) dphi -= TMath::TwoPi();
    if ( dphi < -TMath::Pi() ) dphi += TMath::TwoPi();
    return dphi;
}

//________________
bool DiJetAnalysis::isGoodGenJet(const GenJet* jet) {
    bool goodJet{false};
    float etaCut[2] {fJetEtaLab[0], fJetEtaLab[1]}; 
    float eta = jet->eta();

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
float DiJetAnalysis::dijetEtaInFrame(const float& eta1, const float& eta2, bool isCM) {
    float etaDijet = 0.5 * (eta1 + eta2);
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
    float rawPt = jet->rawPt();
    float trackMaxPt = jet->trackMaxPt();
    if ( TMath::Abs( jet->eta() ) < 2.4 &&
         ( trackMaxPt/rawPt < 0.01 ||
           trackMaxPt/rawPt > 0.98) ) {
        goodTrackMax = {false};
    }

    // if ( fVerbose ) {
    //     std::cout << "TrackMaxPt/rawPt: " << trackMaxPt/rawPt << ( (goodTrackMax) ? " [passed]" : " [failed]" ) 
    //               << ( (trackMaxPt/rawPt < 0.01) ? " too low value " : "" ) << ( (trackMaxPt/rawPt > 0.98) ? " too large value " : "" )
    //               << std::endl;
    // }

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

    // if ( fVerbose ) {
    //     std::cout << "JetId selection results: " << ( (passJetId) ? "[passed]" : "[failed]" ) << " Reasons ";
    //     std::cout << Form("passNHF: %d \tpassNEF: %d \tpassNumConst: %d \tpassMuonFrac: %d \tpassChFrac: %d \tpassChMult: %d \tpassChEmFrac: %d \tpassNeutMult: %d\n", 
    //                       passNHF, passNEF, passNumOfConstituents, passMuonFrac, passChargedFrac, 
    //                       passChargedMult , passChargedEmFrac , passNeutralMult);
    // }
		
	return passJetId;
}

//________________
float DiJetAnalysis::boostEta2CM(const float &eta) {
    // if ( fVerbose ) {
    //     std::cout << "DiJetAnalysis::boostEta2CM -- begin" << std::endl;
    // }
    float etaCM = eta;

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

    // if ( fVerbose ) {
    //     std::cout << Form("eta: %5.2f  ->  etaCM: %5.2f", eta, etaCM) << std::endl;
    //     std::cout << "DiJetAnalysis::boostEta2CM -- end" << std::endl;
    // }
    return etaCM;
}

//________________
float DiJetAnalysis::etaLab(const float &eta) {
    // if ( fVerbose ) {
    //     std::cout << "DiJetAnalysis::etaLab -- begin" << std::endl;
    // }

    float etaL = eta;
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

    // if ( fVerbose ) {
    //     std::cout << Form("eta: %5.2f  ->  etaLab: %5.2f", eta, etaL) << std::endl;
    //     std::cout << "DiJetAnalysis::etaLab -- end" << std::endl;
    // }
    return etaL;
}
    
//________________
bool DiJetAnalysis::isGoodRecoJet(const RecoJet* jet) {
    bool goodJet{false};
    bool goodKine{false};
    bool hasMatching{false};

    float etaCut[2] {fJetEtaLab[0], fJetEtaLab[1]};

    float eta = jet->eta();

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
void DiJetAnalysis::findMcWeight(const float& ptLead, const float& ptSublead) {

    fMcReweight = {1};
    // if ( fUseMcReweighting !=0 ) {
    //     if ( fVerbose ) {
    //         std::cout << Form("DiJetAnalysis::findMcWeight - ptLead: %5.1f ptSublead: %5.1f\n", ptLead, ptSublead);
    //     }

    //     int ptLeadBin{-1}; 
    //     if ( ptLead >= fJetPtLow && ptLead <= fJetPtHi ) {
    //         ptLeadBin = ( ptLead - fJetPtLow ) / fJetPtStep;
    //     }
    //     int ptSubleadBin{-1};
    //     if ( ptSublead >= fJetPtLow && ptSublead <= fJetPtHi ) {
    //         ptSubleadBin = ( ptSublead - fJetPtLow ) / fJetPtStep;
    //     }
    //     float val = ( ptLeadBin >=0 && ptSubleadBin >= 0 ) ? 
    //                 fJetPtLeadPtSubleadReweightMatrix[ptLeadBin][ptSubleadBin] : 1.;

    //     if ( fVerbose ) {
    //         std::cout << Form("\t ptLeadBin: %d ptSubleadBin: %d weight: %6.3f\n", ptLeadBin, ptSubleadBin, val);
    //     }    

    //     fMcReweight = val;
    // }
    // else {
    //     fMcReweight = {1};
    // }
}

//________________
void DiJetAnalysis::makePtSortedJetVectors(const Event* event) {
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::makePtSortedJetVectors -- begin" << std::endl;
    }

    // Clear vectors
    fRecoPtSortedJetIds.clear();
    if ( fIsMc ) {
        fGenPtSortedJetIds.clear();
        fRefSelRecoPtSortedJetIds.clear();
    }

    // Reco jet iterators
    RecoJetIterator recoJetIter;
    // Gen jet iterators
    GenJetIterator genJetIter;

    //
    // Reco jets
    //

    // Jet counter
    int recoJetCounter{0};
    // Loop over reconstructed jets and store indices of good jets
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        recoJetCounter++;

        // Check selection criteria
        bool passTrkMax = isGoodTrkMax( (*recoJetIter) );
        bool passJetId = isGoodJetId( (*recoJetIter) );

        if ( fUseJetIdSelection && !passJetId ) { 
            // if ( fVerbose ) { 
            //     std::cout << "JetId selection failed. Skip jet" << std::endl; 
            // }
            continue; 
        }
        if ( !fUseJetIdSelection && !passTrkMax ) {
            // if ( fVerbose ) { 
            //     // std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
            // }
            continue; 
        }

        fRecoPtSortedJetIds.push_back( recoJetCounter-1 );
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    // Sort indices based on the jet corrected pT (from high to low)
    std::sort( fRecoPtSortedJetIds.begin(), fRecoPtSortedJetIds.end(), [&](int i, int j) { 
        return event->recoJetCollection()->at(i)->ptJECCorr() > event->recoJetCollection()->at(j)->ptJECCorr(); 
    } );

    if ( fVerbose ) {
        // Print sorted indices and corresponding jet pT    
        for (const auto& id : fRecoPtSortedJetIds) {
            std::cout << Form("Sorted jet index: %d | pT: %5.1f\n", id, event->recoJetCollection()->at(id)->ptJECCorr());
        }
    }

    if ( fRecoPtSortedJetIds.size() >= 2) {
        // Indices of leading and subleading reco jets
        fRecoIdLead = fRecoPtSortedJetIds.at(0);
        fRecoIdSubLead = fRecoPtSortedJetIds.at(1);
    }
    else if ( fRecoPtSortedJetIds.size() == 1 ) {
        fRecoIdLead = fRecoPtSortedJetIds.at(0);
        fRecoIdSubLead = -1;
    }

    //
    // Pure Monte Carlo part
    //

    if ( fIsMc ) {

        //
        // Gen jets
        //

        int genJetCounter{0};
        // Loop over generated jets and store indices of good jets
        for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {
            genJetCounter++;
            fGenPtSortedJetIds.push_back( genJetCounter-1 );
        } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )

        // Sort indices based on the jet pT (from high to low)
        std::sort( fGenPtSortedJetIds.begin(), fGenPtSortedJetIds.end(), [&](int i, int j) { 
            return event->genJetCollection()->at(i)->pt() > event->genJetCollection()->at(j)->pt(); 
        } );

        if ( fVerbose ) {
            // Print sorted indices and corresponding jet pT
            for (const auto& id : fGenPtSortedJetIds) {
                std::cout << Form("Sorted gen jet index: %d | pT: %5.1f\n", id, event->genJetCollection()->at(id)->pt());
            }
        }

        if ( fGenPtSortedJetIds.size() >= 2) {
            // Indices of leading and subleading gen jets
            fGenIdLead = fGenPtSortedJetIds.at(0);
            fGenIdSubLead = fGenPtSortedJetIds.at(1);
        }
        else if ( fGenPtSortedJetIds.size() == 1 ) {
            fGenIdLead = fGenPtSortedJetIds.at(0);
            fGenIdSubLead = -1;
        }

        //
        // Ref-selected reco jets
        //

        int refSelRecoJetCounter{0};
        // Loop over reconstructed jets and select those only that have matching gen jets
        for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {
            refSelRecoJetCounter++;
            if ( !(*recoJetIter)->hasMatching() ) continue;
            // Check selection criteria
            bool passTrkMax = isGoodTrkMax( (*recoJetIter) );
            bool passJetId = isGoodJetId( (*recoJetIter) );
    
            if ( fUseJetIdSelection && !passJetId ) { 
                // Do not forget to increment the counter
                // if ( fVerbose ) { 
                //     std::cout << "JetId selection failed. Skip jet" << std::endl; 
                // }
                continue; 
            }
            if ( !fUseJetIdSelection && !passTrkMax ) {
                // Do not forget to increment the counter
                // if ( fVerbose ) { 
                //     std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
                // }
                continue; 
            }
            fRefSelRecoPtSortedJetIds.push_back( refSelRecoJetCounter-1 );
        } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

        // Sort indices based on the gen jet pT from the reco matched jets (from high to low)
        std::sort( fRefSelRecoPtSortedJetIds.begin(), fRefSelRecoPtSortedJetIds.end(), [&](int i, int j) { 
            return event->genJetCollection()->at( event->recoJetCollection()->at(i)->genJetId() )->pt() > event->genJetCollection()->at( event->recoJetCollection()->at(j)->genJetId() )->pt(); 
        } );

        if ( fVerbose ) {
            // Print sorted indices and corresponding gen jet pT
            for (const auto& id : fRefSelRecoPtSortedJetIds) {
                std::cout << Form("Sorted ref-selected reco jet index: %d | gen pT: %5.1f | reco pT: %5.1f\n", id, event->genJetCollection()->at( event->recoJetCollection()->at(id)->genJetId() )->pt(), event->recoJetCollection()->at(id)->ptJECCorr());
            }
        }

        if ( fRefSelRecoPtSortedJetIds.size() >= 2) {
            // Indices of leading and subleading ref-selected reco jets
            fRefSelRecoIdLead = fRefSelRecoPtSortedJetIds.at(0);
            fRefSelRecoIdSubLead = fRefSelRecoPtSortedJetIds.at(1);
        }
        else if ( fRefSelRecoPtSortedJetIds.size() == 1 ) {
            fRefSelRecoIdLead = fRefSelRecoPtSortedJetIds.at(0);
            fRefSelRecoIdSubLead = -1;
        }

    } // if ( fIsMc )

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::makePtSortedJetVectors -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::processInclusiveJets(const Event* event, const double& weight) {
    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processInclusiveJets -- begin" << std::endl;
    }

    fHM->hRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), 1. );
    processRecoJets( event, weight );

    if ( fIsMc ) {
        fHM->hGenJetCollectionSize->Fill( event->genJetCollection()->size(), 1. );
        fHM->hGenVsRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), event->genJetCollection()->size(), 1. );
        processGenJets( event, weight );
        processRefJets( event, weight );
    } // if ( fIsMc )

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processInclusiveJets -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRecoJets(const Event* event, const double &weight) {

    if ( weight <= 0. ) {
        std::cout << "DiJetAnalysis::processRecoJets -- weight is zero or negative. Skip processing." << std::endl;
        return;
    }

    // ptHat value
    float ptHat = event->ptHat();

    // Jet iterators
    RecoJetIterator recoJetIter;

    // Jet counter
    int recoJetCounter{0};
    if ( event->recoJetCollection()->size() > 0 ) {

        // Loop over reconstructed jets and search for leading and subleading jets
        for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

            float pt = (*recoJetIter)->ptJECCorr();
            float eta = etaLab( (*recoJetIter)->eta() );
            float phi = (*recoJetIter)->phi();
            float ptRaw = (*recoJetIter)->pt();
    
            if ( fVerbose ) {
                if ( !fIsMc ) {
                    std::cout << Form("Reco jet #%d pt: %5.1f eta: %3.2f phi: %3.2f ptRaw: %5.1f", recoJetCounter, pt, eta, phi, ptRaw) << std::endl;
                }
                else {
                    std::cout << Form("Reco jet #%d pt: %5.1f eta: %3.2f phi: %3.2f ptRaw: %5.1f genJetId: %d", recoJetCounter, pt, eta, phi, ptRaw, (*recoJetIter)->genJetId() ) << std::endl;
                }
            }
    
            recoJetCounter++;

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
                // if ( fVerbose ) { 
                //     std::cout << "JetId selection failed. Skip jet" << std::endl; 
                // }
                continue; 
            }
            if ( !fUseJetIdSelection && !passTrkMax ) {
                // if ( fVerbose ) { 
                //     // std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
                // }
                continue; 
            }

            fHM->hRecoInclusiveJetPt->Fill(pt, weight);
            fHM->hRecoInclusiveAllJetPtVsEta->Fill(eta, pt, weight);
            fHM->hRecoInclusiveJetPtRawVsEta->Fill(eta, ptRaw, weight);
            fHM->hRecoInclusiveJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

            // Inclusive leading jet
            if ( (recoJetCounter - 1) == fRecoIdLead && (fRecoIdLead >= 0) ) {
                fHM->hRecoLeadJetInclusivePtVsEta->Fill( eta, pt, weight );
                fHM->hRecoLeadJetInclusivePtEtaPtHat->Fill( eta, pt, ptHat, weight );
            }
            // Inclusive subleading jet
            if ( (recoJetCounter - 1) == fRecoIdSubLead && (fRecoIdSubLead >= 0) ) {
                fHM->hRecoSubLeadJetInclusivePtVsEta->Fill( eta, pt, weight );
                fHM->hRecoSubLeadJetInclusivePtEtaPtHat->Fill( eta, pt, ptHat, weight );
            }

            if ( pt > 30. ) {
                fHM->hRecoGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), weight );
                fHM->hRecoGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), weight );
            }

            // On MC check reco jet matching to gen
            if ( fIsMc ) {

                // If reco jet has matching to gen jet
                if ( (*recoJetIter)->hasMatching() ) {

                    GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
                    if ( !matchedJet ) {
                        std::cout << Form("Cannot retrieve gen jet with id: %d", (*recoJetIter)->genJetId() ) << std::endl;
                        continue;
                    }
                    fHM->hRecoMatchedJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
                    
                    float genPt = matchedJet->pt();
                    float genEta = etaLab( matchedJet->eta() );
                    float genPhi = matchedJet->phi();

                    float JES = pt/genPt;
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
                    fHM->hRecoInclusiveMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );

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
                    if ( (recoJetCounter-1) == fRecoIdLead && (fRecoIdLead >= 0) ) {
                        fHM->hRecoLeadJetMatchedPtVsEta->Fill( eta, pt, weight);
                        fHM->hLeadingJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                        fHM->hRecoLeadJetMatchedPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                        fHM->hRefLeadJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );
                    }

                    // Subleading jet
                    if ( (recoJetCounter-1) == fRecoIdSubLead && (fRecoIdSubLead >= 0) ) {
                        fHM->hRecoSubLeadJetMatchedPtVsEta->Fill( eta, pt, weight);
                        fHM->hSubleadingJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                        fHM->hRecoSubLeadJetUnmatchedPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                        fHM->hRefSubLeadJetPtEtaPtHat->Fill( genEta, genPt, ptHat, weight );
                    }

                } // if ( (*recoJetIter)->hasMatching() )
                else {
                    // Fill unmatched jets
                    fHM->hRecoInclusiveUnmatchedJetPtVsEta->Fill(eta, pt, weight);
                    fHM->hRecoUnmatchedJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
                    fHM->hRecoInclusiveUnmatchedJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

                    // Leading jet
                    if ( (recoJetCounter-1) == fRecoIdLead && (fRecoIdLead >= 0) ) {
                        fHM->hRecoLeadJetUnmatchedPtVsEta->Fill( eta, pt, weight);
                        fHM->hRecoLeadJetUnmatchedPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                    }

                    // Subleading jet
                    if ( (recoJetCounter-1) == fRecoIdSubLead && (fRecoIdSubLead >= 0) ) {
                        fHM->hRecoSubLeadJetUnmatchedPtVsEta->Fill( eta, pt, weight);
                        fHM->hRecoSubLeadJetUnmatchedPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                        
                    }
                } // else
            } // if ( fIsMc )
        } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )
    } // if ( event->recoJetCollection()->size() > 0 )

    if ( fVerbose ) {
        std::cout << Form("Reco jet idLead: %d  idSubLead: %d\n", fRecoIdLead, fRecoIdSubLead);
    }
}

//________________
void DiJetAnalysis::processGenJets(const Event* event, const double &weight) {

    if ( weight <= 0. ) {
        std::cout << "DiJetAnalysis::processGenJets -- weight is zero or negative. Skip processing." << std::endl;
        return;
    }
    // ptHat value
    float ptHat = event->ptHat();

    // Jet iterators
    GenJetIterator genJetIter;

    // Jet counter
    int genJetCounter{0};
    if ( event->genJetCollection()->size() > 0 ) {

        // Loop over generated jets and search for leading and subleading jets
        for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

            float pt = (*genJetIter)->pt();
            float eta = etaLab( (*genJetIter)->eta() );
            float phi = (*genJetIter)->phi();
    
            if ( fVerbose ) {
                std::cout << Form("Gen jet #%d: pt: %5.2f eta: %5.2f phi: %5.2f\n", genJetCounter, pt, eta, phi);
                // std::cout << "Gen jet #" << counter << " ";
                // (*genJetIter)->print();
            }
    
            genJetCounter++;
    
            // Apply single-jet selection to gen jets
            //if ( !isGoodGenJet( *genJetIter ) ) continue;

            // Fill inclusive jet pt
            fHM->hGenInclusiveJetPt->Fill(pt, weight);
            fHM->hGenInclusiveJetPtEta->Fill(eta, pt, weight);
            fHM->hGenInclusiveJetPtEtaPtHat->Fill(eta, pt, event->ptHat(), weight);

            if ( (genJetCounter - 1) ==  fGenIdLead && (fGenIdLead >= 0) ) {
                fHM->hGenLeadJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
            }
            if ( (genJetCounter - 1) ==  fGenIdSubLead && (fGenIdSubLead >= 0) ) {
                fHM->hGenSubLeadJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
            }
    
            if ( pt > 30. ) {
                fHM->hGenGoodInclusiveJetEtaLabFrame->Fill( etaLab(eta), weight);
                fHM->hGenGoodInclusiveJetEtaCMFrame->Fill( boostEta2CM(eta), weight );
            }

        } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )
    } // if ( event->genJetCollection()->size() > 0 )

    if ( fVerbose ) {
        std::cout << Form("Gen jet idLead: %d  idSubLead: %d\n", fGenIdLead, fGenIdSubLead);
    }
}

//________________
void DiJetAnalysis::processRefJets(const Event* event, const double &weight) {

    if ( weight <= 0. ) {
        std::cout << "DiJetAnalysis::processRefJets -- weight is zero or negative. Skip processing." << std::endl;
        return;
    }

    // ptHat value
    float ptHat = event->ptHat();

    // Jet iterators
    RecoJetIterator recoJetIter;

    // Jet counter
    int refSelJetCounter{0};
    if ( event->recoJetCollection()->size() > 0 ) {

        // Loop over reconstructed jets and search for leading and subleading gen-matched jets
        for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

            refSelJetCounter++;
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
                // if ( fVerbose ) { 
                //     std::cout << "JetId selection failed. Skip jet" << std::endl; 
                // }
                continue; 
            }
            if ( !fUseJetIdSelection && !passTrkMax ) {
                // Do not forget to increment the counter
                // if ( fVerbose ) { 
                //     std::cout << "TrackMaxPt/rawPt selection failed. Skip jet" << std::endl; 
                // }
                continue; 
            }

            // Retrieve matched gen jet
            GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
            float genPt = matchedJet->pt();
            float genEta = etaLab( matchedJet->eta() );
            float genPhi = matchedJet->phi();

            if ( fVerbose ) {
                std::cout << Form("Ref jet #%d pt: %5.2f eta: %5.2f phi: %5.2f --> Reco jet #%d pt: %5.2f eta: %5.2f phi: %5.2f\n", (*recoJetIter)->genJetId(), genPt, genEta, genPhi, refSelJetCounter-1, (*recoJetIter)->ptJECCorr(), (*recoJetIter)->eta(), (*recoJetIter)->phi());
                // matchedJet->print();
                // std::cout << "Reco jet #" << counter-1 << " ";
                // (*recoJetIter)->print();
            }

            fHM->hRefSelInclusiveJetPt->Fill( genPt, weight * fMcReweight );
            fHM->hRefSelInclusiveJetPtEta->Fill(genEta, genPt, weight * fMcReweight);
            fHM->hRefSelInclusiveJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight * fMcReweight);

            if ( (refSelJetCounter - 1) == fRefSelRecoIdLead && (fRefSelRecoIdLead >= 0) ) {
                fHM->hRefSelLeadJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight * fMcReweight);
            }
            if ( (refSelJetCounter - 1) == fRefSelRecoIdSubLead && (fRefSelRecoIdSubLead >= 0) ) {
                fHM->hRefSelSubLeadJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight * fMcReweight);
            }
        } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )
    } // if ( event->recoJetCollection()->size() > 0 )
    
    if ( fVerbose ) {
        std::cout << Form("Ref-selected reco jet idLead: %d  idSubLead: %d\n", fRefSelRecoIdLead, fRefSelRecoIdSubLead);
    }
}

//________________
bool DiJetAnalysis::isOverweightedEvent(const Event* event, const double& weight) {

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::isOverweightedEvent -- begin" << std::endl;
    }

    // Clear the state
    bool isRecoOverweightedEvent = {false};
    bool isGenOverweightedEvent = {false};

    // Event ptHat value
    float ptHat = event->ptHat();

    //
    // Process reco jets in order to remove x-jets
    //
    if ( fRecoIdLead >= 0 && fRecoIdSubLead >= 0 ) {
        RecoJet *recoLeadJet = event->recoJetCollection()->at( fRecoIdLead );
        RecoJet *recoSubLeadJet = event->recoJetCollection()->at( fRecoIdSubLead );
        
        float ptRecoLead = recoLeadJet->pt();
        float ptRecoSubLead = recoSubLeadJet->pt();
        float dijetRecoPt = ptRecoLead + ptRecoSubLead;
        float dijetRecoPtAve = 0.5 * dijetRecoPt;

        fHM->hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Fill( ptRecoLead, ptRecoLead/ptHat, 1. );
        fHM->hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Fill( ptRecoLead, ptRecoLead/ptHat, weight );
        
        fHM->hRecoDijetPtOverPtHatVsDijetPt->Fill(dijetRecoPt, dijetRecoPt/ptHat, 1.);
        fHM->hRecoDijetPtOverPtHatVsDijetPtWeighted->Fill(dijetRecoPt, dijetRecoPt/ptHat, weight);
        fHM->hRecoDijetPtAveOverPtHatVsDijetPtAve->Fill(dijetRecoPtAve, dijetRecoPtAve/ptHat, 1.);
        fHM->hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Fill(dijetRecoPtAve, dijetRecoPtAve/ptHat, weight);

        if ( isOverweighted( ptRecoLead, dijetRecoPtAve, ptHat ) ) {
            if ( fVerbose ) {
                std::cout << Form("Overweighted event [reco]: ptLead/ptHat = %3.2f ptAve/ptHat = %3.2f", ptRecoLead/ptHat, dijetRecoPtAve/ptHat) << std::endl;
            }
            isRecoOverweightedEvent = {true};
        } // if ( isOverweightedEvent( ptRecoLead, ptHat ) )
    } // if ( fRecoIdLead >= 0 && fRecoIdSubLead >= 0 )
    else {
        // Skip events with less than 2 jets
        isRecoOverweightedEvent = {true};
    }

    //
    // Process gen jets in order to remove overweighted gen events
    //

    if ( fGenIdLead >= 0 && fGenIdSubLead >= 0 ) {

        GenJet* leadJet = event->genJetCollection()->at( fGenIdLead );
        GenJet* subLeadJet = event->genJetCollection()->at( fGenIdSubLead );
        float ptGenLead = leadJet->pt();
        float ptGenSubLead = subLeadJet->pt();
        float dijetGenPt = ptGenLead + ptGenSubLead;
        float dijetGenPtAve = 0.5 * dijetGenPt;

        fHM->hGenLeadingJetPtOverPtHatVsLeadingJetPt->Fill( ptGenLead, ptGenLead/ptHat, 1. );
        fHM->hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Fill( ptGenLead, ptGenLead/ptHat, weight );

        fHM->hGenDijetPtOverPtHatVsDijetPt->Fill(dijetGenPt, dijetGenPt/ptHat, 1.);
        fHM->hGenDijetPtOverPtHatVsDijetPtWeighted->Fill(dijetGenPt, dijetGenPt/ptHat, weight);
        fHM->hGenDijetPtAveOverPtHatVsDijetPtAve->Fill(dijetGenPtAve, dijetGenPtAve/ptHat, 1.);
        fHM->hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Fill(dijetGenPtAve, dijetGenPtAve/ptHat, weight);

        if ( isOverweighted( ptGenLead, dijetGenPtAve, ptHat ) ) {
            if ( fVerbose ) {
                std::cout << Form("Overweighted event [Gen]: ptLead/ptHat = %3.2f ptAve/ptHat = %3.2f", ptGenLead/ptHat, dijetGenPtAve/ptHat) << std::endl;
            }
            isGenOverweightedEvent = {true};
        } // if ( isOverweightedEvent( ptGenLead, ptHat ) )
    } // if ( fGenIdLead >= 0 && fGenIdSubLead >= 0 )
    else {
        // Skip events with less than 2 jets
        isGenOverweightedEvent = {true};
    }


    bool overweightedEvent = isRecoOverweightedEvent || isGenOverweightedEvent;

    if ( fVerbose ) {
        std::cout << Form("Event overweighted: %s Reco overweighted: %s Gen overweighted: %s", 
                          ((overweightedEvent) ? "[true]" : "[false]"),
                          ((isRecoOverweightedEvent) ? "[true]" : "[false]"), 
                          ((isGenOverweightedEvent) ? "[true]" : "[false]")) << std::endl;
        std::cout << "DiJetAnalysis::isOverweightedEvent -- end" << std::endl;
    }
    return overweightedEvent;
}

//________________
bool DiJetAnalysis::isOverweighted(const float& ptLead, const float& dijetPtAve, const float& ptHat) {
    return (  ( ( ptLead / ptHat ) > 2.5) || ( ( dijetPtAve / ptHat ) > 1.7) );
}


//________________
void DiJetAnalysis::processGenDijets(const Event* event, const double &weight) {

    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::processGenDijets -- begin" << std::endl;
    }

    if ( weight <= 0. ) {
        std::cerr << "Error: weight is not positive: " << weight << std::endl;
        return;
    }

    fMcReweight = {1.};
    fIsGenDijetLabFound = {false}; 
    fIsGenDijetCMFound = {false};

    //
    // Check for gen dijet
    //
    if ( fGenIdLead < 0 || fGenIdSubLead < 0 ) {
        if ( fVerbose ) {
            std::cout << "Gen dijet not found" << std::endl;
            std::cout << "DiJetAnalysis::processGenDijets -- end" << std::endl;
        }
        return;
    } 

    GenJet* leadJet = event->genJetCollection()->at( fGenIdLead );
    float ptGenLead = leadJet->pt();
    float phiGenLead = leadJet->phi();

    GenJet* subLeadJet = event->genJetCollection()->at( fGenIdSubLead );
    float ptGenSubLead = subLeadJet->pt();
    float phiGenSubLead = subLeadJet->phi();
    

    if ( fVerbose ) {
        std::cout << Form("Gen dijet idLead: %d ptLead: %5.2f idSubLead: %d ptSubLead: %5.2f\n", fGenIdLead, ptGenLead, fGenIdSubLead, ptGenSubLead);
    }

    float dijetGenPtAve = 0.5 * (ptGenLead + ptGenSubLead);
    float dijetGenDphi = deltaPhi(phiGenLead, phiGenSubLead);


    // Specifically for pPb
    float etaLead = etaLab( leadJet->eta() );
    float etaSubLead = etaLab( subLeadJet->eta() );    
    float dijetEta = 0.5 * (etaLead + etaSubLead);
    float dijetGenDetaCM = 0.5 * ( etaLead - etaSubLead );

    float x_Pb = 2. * dijetGenPtAve / fCollisionEnergy * TMath::Exp( -1. * dijetGenDetaCM ) * TMath::CosH( dijetGenDetaCM );
    float x_p = 2. * dijetGenPtAve / fCollisionEnergy * TMath::Exp( dijetGenDetaCM ) * TMath::CosH( dijetGenDetaCM );
    float xPbOverXp = x_Pb / x_p;

    // if ( fVerbose ) {
    //     std::cout << "Inclusive gen dijet\n";
    //     std::cout << Form("Gen lead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", ptGenLead, etaLead, phiGenLead);
    //     std::cout << Form("Gen sublead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", ptGenSubLead, etaSubLead, phiGenSubLead);
    //     std::cout << Form("Gen dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f delta eta: %5.2f x_Pb: %5.2f x_p: %5.2f xPbOverXp: %5.2f\n", dijetGenPtAve, dijetEta, dijetGenDphi, dijetGenDetaCM, x_Pb, x_p, xPbOverXp);
    // }

    fHM->hGenInclusiveDijetDetaCM->Fill( dijetGenDetaCM, 1. );
    fHM->hGenInclusiveDijetDetaCMWeighted->Fill( dijetGenDetaCM, weight );
    fHM->hGenInclusiveDijetDetaCMPt->Fill( dijetGenDetaCM, dijetGenPtAve, 1. );
    fHM->hGenInclusiveDijetDetaCMPtWeighted->Fill( dijetGenDetaCM, dijetGenPtAve, weight );
    fHM->hGenInclusiveDijetEtaDetaCMPt->Fill( dijetEta, dijetGenDetaCM, dijetGenPtAve, 1. );
    fHM->hGenInclusiveDijetEtaDetaCMPtWeighted->Fill( dijetEta, dijetGenDetaCM, dijetGenPtAve, weight );
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

    bool fIsGenDijetLabFound = isGoodDijet(ptGenLead, leadJet->eta(), ptGenSubLead, subLeadJet->eta(), dijetGenDphi, false);
    if ( fVerbose ) {
        std::cout << Form("Gen dijet in lab frame is %s\n", ((fIsGenDijetLabFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze gen dijets in lab frame
    if ( fIsGenDijetLabFound ) {

        // Flush the eta values to reflect the frame
        float etaGenLeadLab = etaLab( leadJet->eta() );
        float etaGenSubLeadLab = etaLab( subLeadJet->eta() );
        float dijetGenEtaLab = 0.5 * (etaGenLeadLab + etaGenSubLeadLab);
        float dijetGenDetaCM = 0.5 * ( etaGenLeadLab - etaGenSubLeadLab );

        float x_Pb = 2. * dijetGenPtAve / fCollisionEnergy * TMath::Exp( -1. * dijetGenDetaCM ) * TMath::CosH( dijetGenDetaCM );
        float x_p = 2. * dijetGenPtAve / fCollisionEnergy * TMath::Exp( dijetGenDetaCM ) * TMath::CosH( dijetGenDetaCM );
        float xPbOverXp = x_Pb / x_p;

        if ( fVerbose ) {
            std::cout << "Lab frame gen dijet\n";
            std::cout << Form("Gen lead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", ptGenLead, etaGenLeadLab, phiGenLead);
            std::cout << Form("Gen sublead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", ptGenSubLead, etaGenSubLeadLab, phiGenSubLead);
            std::cout << Form("Gen dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f delta eta: %5.2f x_Pb: %5.2f x_p: %5.2f xPbOverXp: %5.2f\n", dijetGenPtAve, dijetGenEtaLab, dijetGenDphi, dijetGenDetaCM, x_Pb, x_p, xPbOverXp);
        }

        fHM->hGenSelectedDijetDetaCM->Fill( dijetGenDetaCM, 1. );
        fHM->hGenSelectedDijetDetaCMWeighted->Fill( dijetGenDetaCM, weight );
        fHM->hGenSelectedDijetDetaCMPt->Fill( dijetGenDetaCM, dijetGenPtAve, 1. );
        fHM->hGenSelectedDijetDetaCMPtWeighted->Fill( dijetGenDetaCM, dijetGenPtAve, weight );
        fHM->hGenSelectedDijetEtaDetaCMPt->Fill( dijetGenEtaLab, dijetGenDetaCM, dijetGenPtAve, 1. );
        fHM->hGenSelectedDijetEtaDetaCMPtWeighted->Fill( dijetGenEtaLab, dijetGenDetaCM, dijetGenPtAve, weight );
        fHM->hGenSelectedDijetXPb->Fill( x_Pb, 1. );
        fHM->hGenSelectedDijetXPbWeighted->Fill( x_Pb, weight );
        fHM->hGenSelectedDijetXp->Fill( x_p, 1. );
        fHM->hGenSelectedDijetXpWeighted->Fill( x_p, weight );
        fHM->hGenSelectedDijetXPbOverXp->Fill( xPbOverXp, 1. );
        fHM->hGenSelectedDijetXPbOverXpWeighted->Fill( xPbOverXp, weight );
        fHM->hGenSelectedDijetXPbOverXpEta->Fill( xPbOverXp, dijetGenEtaLab, 1. );
        fHM->hGenSelectedDijetXPbOverXpEtaWeighted->Fill( xPbOverXp, dijetGenEtaLab, weight );
        
        fHM->hGenPtLeadPtSublead->Fill( ptGenLead, ptGenSubLead, weight );
        fHM->hGenEtaLeadEtaSublead->Fill( etaGenLeadLab, etaGenSubLeadLab, weight );
        fHM->hGenPtLeadPtSubleadMcReweight->Fill( ptGenLead, ptGenSubLead, weight * fMcReweight );
        fHM->hGenEtaLeadEtaSubleadMcReweight->Fill( etaGenLeadLab, etaGenSubLeadLab, weight * fMcReweight );

        double genDijetLeadSublead[9] {dijetGenPtAve, dijetGenEtaLab, dijetGenDphi, 
                                        ptGenLead, etaGenLeadLab, phiGenLead, 
                                        ptGenSubLead, etaGenSubLeadLab, phiGenSubLead };

        fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(genDijetLeadSublead);
        fHM->hGenDijetPtEtaPhiDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(genDijetLeadSublead, weight * fMcReweight );
        fHM->hGenDijetEta->Fill(dijetGenEtaLab, weight * fMcReweight );
        fHM->hGenDijetPtEtaDphi->Fill(dijetGenPtAve, dijetGenEtaLab, dijetGenDphi, 1.);
        fHM->hGenDijetPtEtaDphiWeighted->Fill(dijetGenPtAve, dijetGenEtaLab, dijetGenDphi, weight * fMcReweight );
        (dijetGenEtaLab >= 0) ? fHM->hGenDijetPtEtaForward->Fill(dijetGenPtAve, dijetGenEtaLab) : fHM->hGenDijetPtEtaBackward->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaLab));
        (dijetGenEtaLab >= 0) ? fHM->hGenDijetPtEtaForwardWeighted->Fill(dijetGenPtAve, dijetGenEtaLab, weight * fMcReweight) : fHM->hGenDijetPtEtaBackwardWeighted->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaLab), weight * fMcReweight);

        // Find exact ptAve bin
        int ptAveBin = findDijetPtAveBin( dijetGenPtAve );
        int ptAveOldBin = findDijetPtAveOldBin( dijetGenPtAve );

        // New ptAve and eta binning
        if ( ptAveBin >= 0 ) {
            fHM->hGenDijetEta1D[ptAveBin]->Fill( dijetGenEtaLab, 1. );
            fHM->hGenDijetEta1DWeighted[ptAveBin]->Fill( dijetGenEtaLab, weight * fMcReweight );
            fHM->hGenDijetEtaLeadVsEtaSubLead2D[ptAveBin]->Fill( etaGenLeadLab, etaGenSubLeadLab, 1. );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DWeighted[ptAveBin]->Fill( etaGenLeadLab, etaGenSubLeadLab, weight * fMcReweight );

            (dijetGenEtaLab >= 0) ? fHM->hGenDijetEtaForward1D[ptAveBin]->Fill(dijetGenEtaLab, 1.) : fHM->hGenDijetEtaForward1D[ptAveBin]->Fill(TMath::Abs(dijetGenEtaLab), 1.);
            (dijetGenEtaLab >= 0) ? fHM->hGenDijetEtaForward1DWeighted[ptAveBin]->Fill(dijetGenEtaLab, weight * fMcReweight) : fHM->hGenDijetEtaForward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetGenEtaLab), weight * fMcReweight);
        }

        // Old ptAve binning
        if ( ptAveOldBin >= 0 ) {
            // New eta binning
            fHM->hGenDijetEta1DOldPt[ptAveOldBin]->Fill( dijetGenEtaLab, 1. );
            fHM->hGenDijetEta1DOldPtWeighted[ptAveOldBin]->Fill( dijetGenEtaLab, weight * fMcReweight );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPt[ptAveOldBin]->Fill( etaGenLeadLab, etaGenSubLeadLab, 1. );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtWeighted[ptAveOldBin]->Fill( etaGenLeadLab, etaGenSubLeadLab, weight * fMcReweight );

            (dijetGenEtaLab >= 0) ? fHM->hGenDijetEtaForward1DOldPt[ptAveOldBin]->Fill(dijetGenEtaLab, 1.) : fHM->hGenDijetEtaForward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaLab), 1.);
            (dijetGenEtaLab >= 0) ? fHM->hGenDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetGenEtaLab, weight * fMcReweight) : fHM->hGenDijetEtaForward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaLab), weight * fMcReweight);

            // Old eta binning
            fHM->hGenDijetEta1DOldPtBinning[ptAveOldBin]->Fill( dijetGenEtaLab, 1. );
            fHM->hGenDijetEta1DOldPtBinningWeighted[ptAveOldBin]->Fill( dijetGenEtaLab, weight * fMcReweight );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinning[ptAveOldBin]->Fill( etaGenLeadLab, etaGenSubLeadLab, 1. );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[ptAveOldBin]->Fill( etaGenLeadLab, etaGenSubLeadLab, weight * fMcReweight );
            (dijetGenEtaLab >= 0) ? fHM->hGenDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(dijetGenEtaLab, 1.) : fHM->hGenDijetEtaForward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaLab), 1.);
            (dijetGenEtaLab >= 0) ? fHM->hGenDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetGenEtaLab, weight * fMcReweight) : fHM->hGenDijetEtaForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaLab), weight * fMcReweight);
        }
    } // if ( goodDijetLab )

    //
    // CM frame
    //

    bool fIsGenDijetCMFound = isGoodDijet(ptGenLead, leadJet->eta(), ptGenSubLead, subLeadJet->eta(), dijetGenDphi, true);
    if ( fVerbose ) {
        std::cout << Form("Gen dijet in CM frame is %s\n", ((fIsGenDijetCMFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze gen dijets in CM frame
    if ( fIsGenDijetCMFound ) {

        float etaGenLeadCM = boostEta2CM( leadJet->eta() );
        float etaGenSubLeadCM = boostEta2CM( subLeadJet->eta() );
        float dijetGenEtaCM = 0.5 * (etaGenLeadCM + etaGenSubLeadCM);
        float dijetGenDetaCM = 0.5 * ( etaGenLeadCM - etaGenSubLeadCM );

        float x_Pb = 2. * dijetGenPtAve / fCollisionEnergy * TMath::Exp( -1. * dijetGenDetaCM ) * TMath::CosH( dijetGenDetaCM );
        float x_p = 2. * dijetGenPtAve / fCollisionEnergy * TMath::Exp( dijetGenDetaCM ) * TMath::CosH( dijetGenDetaCM );
        float xPbOverXp = x_Pb / x_p;

        if ( fVerbose ) {
            std::cout << "Lab frame gen dijet\n";
            std::cout << Form("Gen lead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", ptGenLead, etaGenLeadCM, phiGenLead);
            std::cout << Form("Gen sublead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", ptGenSubLead, etaGenSubLeadCM, phiGenSubLead);
            std::cout << Form("Gen dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f delta eta: %5.2f x_Pb: %5.2f x_p: %5.2f xPbOverXp: %5.2f\n", dijetGenPtAve, dijetGenEtaCM, dijetGenDphi, dijetGenDetaCM, x_Pb, x_p, xPbOverXp);
        }

        fHM->hGenEtaCMLeadEtaCMSublead->Fill( etaGenLeadCM, etaGenSubLeadCM, weight );

        fHM->hGenDijetEtaCM->Fill(dijetGenEtaCM, weight * fMcReweight );
        fHM->hGenDijetPtEtaDphiCM->Fill(dijetGenPtAve, dijetGenEtaCM, dijetGenDphi, 1.);
        fHM->hGenDijetPtEtaDphiCMWeighted->Fill(dijetGenPtAve, dijetGenEtaCM, dijetGenDphi, weight * fMcReweight );
        (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForward->Fill(dijetGenPtAve, dijetGenEtaCM) : fHM->hGenDijetPtEtaCMBackward->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaCM));
        (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForwardWeighted->Fill(dijetGenPtAve, dijetGenEtaCM, weight * fMcReweight) : fHM->hGenDijetPtEtaCMBackwardWeighted->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaCM), weight * fMcReweight);

        // Find exact ptAve bin
        int ptAveBin = findDijetPtAveBin( dijetGenPtAve );
        int ptAveOldBin = findDijetPtAveOldBin( dijetGenPtAve );

        // New ptAve and eta binning
        if ( ptAveBin >= 0 ) {
            fHM->hGenDijetEta1DCM[ptAveBin]->Fill( dijetGenEtaCM, 1. );
            fHM->hGenDijetEta1DCMWeighted[ptAveBin]->Fill( dijetGenEtaCM, weight * fMcReweight );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DCM[ptAveBin]->Fill( etaGenLeadCM, etaGenSubLeadCM, 1. );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DCMWeighted[ptAveBin]->Fill( etaGenLeadCM, etaGenSubLeadCM, weight * fMcReweight );
            (dijetGenEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1D[ptAveBin]->Fill(dijetGenEtaCM, 1.) : fHM->hGenDijetEtaCMBackward1D[ptAveBin]->Fill(TMath::Abs(dijetGenEtaCM), 1.);
            (dijetGenEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DWeighted[ptAveBin]->Fill(dijetGenEtaCM, weight * fMcReweight) : fHM->hGenDijetEtaCMBackward1DWeighted[ptAveBin]->Fill(TMath::Abs(dijetGenEtaCM), weight * fMcReweight);
        } // if ( ptAveBin >=0 )

        // Old ptAve binning
        if ( ptAveOldBin >= 0 ) {
            fHM->hGenDijetEta1DOldPtCM[ptAveOldBin]->Fill( dijetGenEtaCM, 1. );
            fHM->hGenDijetEta1DOldPtCMWeighted[ptAveOldBin]->Fill( dijetGenEtaCM, weight * fMcReweight );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtCM[ptAveOldBin]->Fill( etaGenLeadCM, etaGenSubLeadCM, 1. );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtCMWeighted[ptAveOldBin]->Fill( etaGenLeadCM, etaGenSubLeadCM, weight * fMcReweight );
            (dijetGenEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPt[ptAveOldBin]->Fill(dijetGenEtaCM, 1.) : fHM->hGenDijetEtaCMBackward1DOldPt[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaCM), 1.);
            (dijetGenEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPtWeighted[ptAveOldBin]->Fill(dijetGenEtaCM, weight * fMcReweight) : fHM->hGenDijetEtaCMBackward1DOldPtWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaCM), weight * fMcReweight);

            fHM->hGenDijetEta1DOldPtBinningCM[ptAveOldBin]->Fill( dijetGenEtaCM, 1. );
            fHM->hGenDijetEta1DOldPtBinningCMWeighted[ptAveOldBin]->Fill( dijetGenEtaCM, weight * fMcReweight );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCM[ptAveOldBin]->Fill( etaGenLeadCM, etaGenSubLeadCM, 1. );
            fHM->hGenDijetEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[ptAveOldBin]->Fill( etaGenLeadCM, etaGenSubLeadCM, weight * fMcReweight );
            (dijetGenEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPtBinning[ptAveOldBin]->Fill(dijetGenEtaCM, 1.) : fHM->hGenDijetEtaCMBackward1DOldPtBinning[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaCM), 1.);
            (dijetGenEtaCM >= 0) ? fHM->hGenDijetEtaCMForward1DOldPtBinningWeighted[ptAveOldBin]->Fill(dijetGenEtaCM, weight * fMcReweight) : fHM->hGenDijetEtaCMBackward1DOldPtBinningWeighted[ptAveOldBin]->Fill(TMath::Abs(dijetGenEtaCM), weight * fMcReweight);
        } // if ( ptAveOldBin >=0 )
    } // if ( goodDijetCM )

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processGenDijets -- end" << std::endl;
    }
}


//________________
void DiJetAnalysis::processRecoDijets(const Event* event, const double &weight) {

    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::processRecoDijets -- begin" << std::endl;
    }

    if ( weight <= 0. ) {
        std::cerr << "Error: weight is not positive: " << weight << std::endl;
        return;
    }

    fMcReweight = {1.};
    fIsRecoDijetLabFound = {false};
    fIsRecoDijetCMFound = {false};

    // ptHat value
    float ptHat = event->ptHat();

    if ( fRecoIdLead < 0 || fRecoIdSubLead < 0 ) {
        if ( fVerbose ) {
            std::cout << "Reco dijet not found" << std::endl;
            std::cout << "DiJetAnalysis::processRecoDijets -- end" << std::endl;
        }
        return;
    }

    // Leading jet
    RecoJet* recoLeadJet = event->recoJetCollection()->at( fRecoIdLead );
    float ptRecoLead = recoLeadJet->ptJECCorr();
    float ptRawRecoLead = recoLeadJet->pt();
    float etaRecoLead = recoLeadJet->eta();
    float phiRecoLead = recoLeadJet->phi();

    // Subleading jet
    RecoJet* recoSubLeadJet = event->recoJetCollection()->at( fRecoIdSubLead );
    float ptRecoSubLead = recoSubLeadJet->ptJECCorr();
    float ptRawRecoSubLead = recoSubLeadJet->pt();
    float etaRecoSubLead = recoSubLeadJet->eta();
    float phiRecoSubLead = recoSubLeadJet->phi();

    if ( fVerbose ) {
        std::cout << Form("Reco dijet: idLead: %d ptLead: %5.2f idSubLead: %d ptSubLead: %5.2f", fRecoIdLead, ptRecoSubLead, fRecoIdSubLead, ptRecoSubLead) << std::endl;
    }

    GenJet* refLeadJet = {nullptr};
    GenJet* refSubLeadJet = {nullptr};

    // For Monte Carlo check if reco jets have matching gen jets. Skip reco dijet if not
    if ( fIsMc ) {
        if ( !recoLeadJet->hasMatching() || !recoSubLeadJet->hasMatching() ) {
            if ( fVerbose ) {
                std::cout << Form("Reco dijet has unmatched jets. idRecoLead: %d idRecoSubLead: %d Lead has matching: %s SubLead has matching: %s", 
                                fRecoIdLead, fRecoIdSubLead, 
                                (recoLeadJet->hasMatching() ? "[true]" : "[false]"), 
                                (recoSubLeadJet->hasMatching() ? "[true]" : "[false]"));
                std::cout << "\tSkip reco dijet\n";
            }
            return;
        }

        // Matching gen jet for leading reco jet
        refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
        if ( !refLeadJet ) {
            std::cerr << "Error: Leading jet has no matching gen jet\n";
            return;
        }

        // Matching gen jet for subleading reco jet
        refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );
        if ( !refSubLeadJet ) {
            std::cerr << "Error: Subleading jet has no matching gen jet\n";
            return;
        }
    } // if ( fIsMc )

    // Dijet parameters
    float dijetRecoPtAve = 0.5 * (ptRecoLead + ptRecoSubLead);
    float dijetRecoDphi = deltaPhi(phiRecoLead, phiRecoSubLead);

    //
    // Lab frame
    //

    bool fIsRecoDijetLabFound = isGoodDijet(ptRecoLead, recoLeadJet->eta(), 
                                            ptRecoSubLead, recoSubLeadJet->eta(), 
                                            TMath::Abs( dijetRecoDphi ), 
                                            false);
    if ( fVerbose ) {
        std::cout << Form("Reco dijet in lab frame is %s\n", ((fIsRecoDijetLabFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze reco dijets in lab frame
    if ( fIsRecoDijetLabFound ) {

        // Flush the eta values to reflect the frame
        float etaRecoLeadLab = etaLab( recoLeadJet->eta() );
        float etaRecoSubLeadLab = etaLab( recoSubLeadJet->eta() );
        float dijetRecoEtaLab = 0.5 * (etaRecoLeadLab + etaRecoSubLeadLab);

        if ( fVerbose ) {
            std::cout << "Reco dijet parameters in the lab frame: " << std::endl;
            std::cout << Form("Reco lead pt: %5.2f eta: %5.2f phi: %5.2f", ptRecoLead, etaRecoLeadLab, phiRecoLead) << std::endl;
            std::cout << Form("Reco sublead pt: %5.2f eta: %5.2f phi: %5.2f", ptRecoSubLead, etaRecoSubLeadLab, phiRecoSubLead) << std::endl;
            std::cout << Form("Reco dijet ptAve: %5.2f dijet eta: %5.2f dijet dphi: %5.2f\n", dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDphi);
        }

        // Correlation between leading and subleading
        fHM->hRecoPtLeadPtSublead->Fill( ptRecoLead, ptRecoSubLead, weight );
        fHM->hRecoEtaLeadEtaSublead->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight );
        fHM->hRecoPtLeadPtSubleadMcReweight->Fill( ptRecoLead, ptRecoSubLead, weight * fMcReweight );
        fHM->hRecoEtaLeadEtaSubleadMcReweight->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight * fMcReweight );

        double dijetRecoInfo[9] { dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDphi,
                                    ptRecoLead, etaRecoLeadLab, phiRecoLead,
                                    ptRecoSubLead, etaRecoSubLeadLab, phiRecoSubLead };
        fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhi->Fill(dijetRecoInfo);
        fHM->hRecoDijetPtEtaDeltaPhiLeadJetPtEtaPhiSubleadJetPtEtaPhiWeighted->Fill(dijetRecoInfo, weight * fMcReweight);
        fHM->hRecoDijetEta->Fill( dijetRecoEtaLab, weight * fMcReweight);
        fHM->hRecoDijetPtEta->Fill( dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight);
        fHM->hRecoDijetPtEtaDphi->Fill( dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDphi, 1. );
        fHM->hRecoDijetPtEtaDphiWeighted->Fill( dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDphi, weight * fMcReweight);
        if ( dijetRecoEtaLab >= 0 ) {
            fHM->hRecoDijetPtEtaForward->Fill(dijetRecoPtAve, dijetRecoEtaLab, 1.);
            fHM->hRecoDijetPtEtaForwardWeighted->Fill(dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight);
        }
        else {
            fHM->hRecoDijetPtEtaBackward->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaLab), 1.);
            fHM->hRecoDijetPtEtaBackwardWeighted->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaLab), weight * fMcReweight);
        }

        // Find dijet exact dijet pT bins
        int recoPtAveBin = findDijetPtAveBin( dijetRecoPtAve );
        int recoPtAveOldBin = findDijetPtAveOldBin( dijetRecoPtAve );

        // New ptAve and eta binning
        if ( recoPtAveBin >=0 ) {
            fHM->hRecoDijetEta1D[recoPtAveBin]->Fill( dijetRecoEtaLab, 1. );
            fHM->hRecoDijetEta1DWeighted[recoPtAveBin]->Fill( dijetRecoEtaLab, weight * fMcReweight );

            fHM->hRecoDijetEtaLeadVsEtaSubLead2D[recoPtAveBin]->Fill( etaRecoLeadLab, etaRecoSubLeadLab, 1. );
            fHM->hRecoDijetEtaLeadVsEtaSubLead2DWeighted[recoPtAveBin]->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight * fMcReweight );

            (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetEtaForward1D[recoPtAveBin]->Fill(dijetRecoEtaLab, 1.) : fHM->hRecoDijetEtaBackward1D[recoPtAveBin]->Fill(TMath::Abs(dijetRecoEtaLab), 1.);
            (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetEtaForward1DWeighted[recoPtAveBin]->Fill(dijetRecoEtaLab, weight * fMcReweight) : fHM->hRecoDijetEtaBackward1DWeighted[recoPtAveBin]->Fill(TMath::Abs(dijetRecoEtaLab), weight * fMcReweight);
        } // if ( recoPtAveBin >=0 )

        // Old ptAve binning
        if ( recoPtAveOldBin >=0 ) {

            // New eta binning
            fHM->hRecoDijetEta1DOldPt[recoPtAveOldBin]->Fill( dijetRecoEtaLab, 1. );
            fHM->hRecoDijetEta1DOldPtWeighted[recoPtAveOldBin]->Fill( dijetRecoEtaLab, weight * fMcReweight );

            fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPt[recoPtAveOldBin]->Fill( etaRecoLeadLab, etaRecoSubLeadLab, 1. );
            fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPtWeighted[recoPtAveOldBin]->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight * fMcReweight );

            (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetEtaForward1DOldPt[recoPtAveOldBin]->Fill(dijetRecoEtaLab, 1.) : fHM->hRecoDijetEtaBackward1DOldPt[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaLab), 1.);
            (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetEtaForward1DOldPtWeighted[recoPtAveOldBin]->Fill(dijetRecoEtaLab, weight * fMcReweight) : fHM->hRecoDijetEtaBackward1DOldPtWeighted[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaLab), weight * fMcReweight);

            // Old eta binning
            fHM->hRecoDijetEta1DOldPtBinning[recoPtAveOldBin]->Fill( dijetRecoEtaLab, 1. );
            fHM->hRecoDijetEta1DOldPtBinningWeighted[recoPtAveOldBin]->Fill( dijetRecoEtaLab, weight * fMcReweight );

            fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinning[recoPtAveOldBin]->Fill( etaRecoLeadLab, etaRecoSubLeadLab, 1. );
            fHM->hRecoDijetEtaLeadVsEtaSubLead2DOldPtBinningWeighted[recoPtAveOldBin]->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight * fMcReweight );

            (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetEtaForward1DOldPtBinning[recoPtAveOldBin]->Fill(dijetRecoEtaLab, 1.) : fHM->hRecoDijetEtaBackward1DOldPtBinning[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaLab), 1.);
            (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetEtaForward1DOldPtBinningWeighted[recoPtAveOldBin]->Fill(dijetRecoEtaLab, weight * fMcReweight) : fHM->hRecoDijetEtaBackward1DOldPtBinningWeighted[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaLab), weight * fMcReweight);
        } // if ( recoPtAveOldBin >=0 )

        // In case of MC
        if ( fIsMc ) {
            float ptRefLeadLab = refLeadJet->pt();
            float etaRefLeadLab = etaLab( refLeadJet->eta() );
            float phiRefLeadLab = refLeadJet->phi();

            float ptRefSubLeadLab = refSubLeadJet->pt();
            float etaRefSubLeadLab = etaLab( refSubLeadJet->eta() );
            float phiRefSubLeadLab = refSubLeadJet->phi();

            float dijetRefPtAveLab = 0.5 * (ptRefLeadLab + ptRefSubLeadLab);
            float dijetRefEtaLab = 0.5 * ( etaRefLeadLab + etaRefSubLeadLab );
            float dijetRefDphiLab = deltaPhi(phiRefLeadLab, phiRefSubLeadLab);

            if ( fVerbose ) {
                std::cout << Form("<-- Ref lead pT: %5.1f eta: %5.2f phi: %5.2f\n", ptRefLeadLab, etaRefLeadLab, phiRefLeadLab);
                std::cout << Form("<-- Ref sublead pT: %5.1f eta: %5.2f phi: %5.2f\n", ptRefSubLeadLab, etaRefSubLeadLab, phiRefSubLeadLab);
                std::cout << Form("<-- Ref dijet ptAve: %5.1f eta: %5.2f dphi: %5.2f\n", dijetRefPtAveLab, dijetRefEtaLab, dijetRefDphiLab);
            }

            // Leading jet information
            double correl[5] { ptRecoLead, ptRawRecoLead, ptRefLeadLab, etaRecoLeadLab, etaRefLeadLab };
            fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
            fHM->hRecoLeadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, weight * fMcReweight );

            // Subleading jet information
            correl[0] = ptRecoSubLead;
            correl[1] = ptRawRecoSubLead;
            correl[2] = ptRefSubLeadLab;
            correl[3] = etaRecoSubLeadLab; 
            correl[4] = etaRefSubLeadLab;
            fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGen->Fill(correl);
            fHM->hRecoSubleadingJetPtCorrPtRawPtRefEtaCorrEtaGenWeighted->Fill(correl, weight * fMcReweight);

            fHM->hRefPtLeadPtSublead->Fill( ptRefLeadLab, ptRefSubLeadLab, weight );
            fHM->hRefEtaLeadEtaSublead->Fill( ptRefLeadLab, ptRefSubLeadLab, weight );
            fHM->hRefPtLeadPtSubleadMcReweight->Fill( ptRefLeadLab, ptRefSubLeadLab, weight * fMcReweight );
            fHM->hRefEtaLeadEtaSubleadMcReweight->Fill( ptRefLeadLab, ptRefSubLeadLab, weight * fMcReweight );

            double dijetRecoUnfold[12] = { dijetRecoPtAve, dijetRecoEtaLab,
                                            ptRecoLead, etaRecoLeadLab,
                                            ptRecoSubLead, etaRecoSubLeadLab,
                                            dijetRefPtAveLab, dijetRefEtaLab,
                                            ptRefLeadLab, etaRefLeadLab,
                                            ptRefSubLeadLab, etaRefSubLeadLab };

            double dijetUnfold[4] = { dijetRecoPtAve, dijetRecoEtaLab, dijetRefPtAveLab, dijetRefEtaLab };

            fHM->hRecoDijetPtEtaRefDijetPtEta->Fill(dijetUnfold, 1.);
            fHM->hRecoDijetPtEtaRefDijetPtEtaWeighted->Fill(dijetUnfold, weight * fMcReweight);

            fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEta->Fill(dijetRecoUnfold);
            fHM->hRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, weight * fMcReweight );
            fHM->hRefDijetEta->Fill( dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEtaLab, dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Fill( dijetRecoEtaLab, dijetRefEtaLab, dijetRecoPtAve, 1.);
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Fill( dijetRecoEtaLab, dijetRefEtaLab, dijetRecoPtAve, weight * fMcReweight );
            fHM->hRefDijetPtEtaDphi->Fill( dijetRefPtAveLab, dijetRefEtaLab, dijetRefDphiLab, 1. );
            fHM->hRefDijetPtEtaDphiWeighted->Fill( dijetRefPtAveLab, dijetRefEtaLab, dijetRefDphiLab, weight * fMcReweight );
            (dijetRefEtaLab >= 0) ? fHM->hRefDijetPtEtaForward->Fill(dijetRefPtAveLab, dijetRefEtaLab) : fHM->hRefDijetPtEtaBackward->Fill(dijetRefPtAveLab, TMath::Abs(dijetRefEtaLab));
            (dijetRefEtaLab >= 0) ? fHM->hRefDijetPtEtaForwardWeighted->Fill(dijetRefPtAveLab, dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetPtEtaBackwardWeighted->Fill(dijetRefPtAveLab, TMath::Abs(dijetRefEtaLab), weight * fMcReweight);         

            // Find exact dijet pT bins
            int refPtAveBin = findDijetPtAveBin( dijetRefPtAveLab );
            int refPtAveOldBin = findDijetPtAveOldBin( dijetRefPtAveLab );

            // New ptAve and eta binning
            if ( refPtAveBin >=0 ) {
                fHM->hRefDijetEta1D[refPtAveBin]->Fill( dijetRefEtaLab, 1. );
                fHM->hRefDijetEta1DWeighted[refPtAveBin]->Fill( dijetRefEtaLab, weight * fMcReweight );

                fHM->hRefEtaLeadVsEtaSubLead2D[refPtAveBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, 1. );
                fHM->hRefEtaLeadVsEtaSubLead2DWeighted[refPtAveBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

                fHM->hRecoVsRefDijetEta2D[refPtAveBin]->Fill( dijetRecoEtaLab, dijetRefEtaLab, 1. );
                fHM->hRecoVsRefDijetEta2DWeighted[refPtAveBin]->Fill( dijetRecoEtaLab, dijetRefEtaLab, weight * fMcReweight );
                fHM->hRecoVsRefLeadJetEta2D[refPtAveBin]->Fill( etaRecoLeadLab, etaRefLeadLab, 1. );
                fHM->hRecoVsRefLeadJetEta2DWeighted[refPtAveBin]->Fill( etaRecoLeadLab, etaRefLeadLab, weight * fMcReweight );
                fHM->hRecoVsRefSubLeadJetEta2D[refPtAveBin]->Fill( etaRecoSubLeadLab, etaRefSubLeadLab, 1. );
                fHM->hRecoVsRefSubLeadJetEta2DWeighted[refPtAveBin]->Fill( etaRecoSubLeadLab, etaRefSubLeadLab, weight * fMcReweight );

                (dijetRefEtaLab >= 0) ? fHM->hRefDijetEtaForward1D[refPtAveBin]->Fill(dijetRefEtaLab, 1.) : fHM->hRefDijetEtaBackward1D[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaLab), 1.);
                (dijetRefEtaLab >= 0) ? fHM->hRefDijetEtaForward1DWeighted[refPtAveBin]->Fill(dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetEtaBackward1DWeighted[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaLab), weight * fMcReweight);
            } // if ( refPtAveBin >=0 )

            // Old ptAve binning
            if ( refPtAveOldBin >=0 ) {

                // New eta binning
                fHM->hRefDijetEta1DOldPt[refPtAveOldBin]->Fill( dijetRefEtaLab, 1. );
                fHM->hRefDijetEta1DOldPt[refPtAveOldBin]->Fill( dijetRefEtaLab, weight * fMcReweight );

                fHM->hRefEtaLeadVsEtaSubLead2DOldPt[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, 1. );
                fHM->hRefEtaLeadVsEtaSubLead2DOldPt[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

                fHM->hRecoVsRefDijetEta2DOldPt[refPtAveOldBin]->Fill( dijetRecoEtaLab, dijetRefEtaLab, 1. );
                fHM->hRecoVsRefDijetEta2DOldPt[refPtAveOldBin]->Fill( dijetRecoEtaLab, dijetRefEtaLab, weight * fMcReweight );
                fHM->hRecoVsRefLeadJetEta2DOldPt[refPtAveOldBin]->Fill( etaRecoLeadLab, etaRefLeadLab, 1. );
                fHM->hRecoVsRefLeadJetEta2DOldPt[refPtAveOldBin]->Fill( etaRecoLeadLab, etaRefLeadLab, weight * fMcReweight );
                fHM->hRecoVsRefSubLeadJetEta2DOldPt[refPtAveOldBin]->Fill( etaRecoSubLeadLab, etaRefSubLeadLab, 1. );
                fHM->hRecoVsRefSubLeadJetEta2DOldPt[refPtAveOldBin]->Fill( etaRecoSubLeadLab, etaRefSubLeadLab, weight * fMcReweight );

                (dijetRefEtaLab >= 0) ? fHM->hRefDijetEtaForward1DOldPt[refPtAveOldBin]->Fill(dijetRefEtaLab, 1.) : fHM->hRefDijetEtaBackward1DOldPt[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), 1.);
                (dijetRefEtaLab >= 0) ? fHM->hRefDijetEtaForward1DOldPtWeighted[refPtAveOldBin]->Fill(dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetEtaBackward1DOldPtWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), weight * fMcReweight);

                // Old eta binning
                fHM->hRefDijetEta1DOldPtBinning[refPtAveOldBin]->Fill( dijetRefEtaLab, 1. );
                fHM->hRefDijetEta1DOldPtBinning[refPtAveOldBin]->Fill( dijetRefEtaLab, weight * fMcReweight );

                fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinning[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, 1. );
                fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinning[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

                fHM->hRecoVsRefDijetEta2DOldPtBinning[refPtAveOldBin]->Fill( dijetRecoEtaLab, dijetRefEtaLab, 1. );
                fHM->hRecoVsRefDijetEta2DOldPtBinning[refPtAveOldBin]->Fill( dijetRecoEtaLab, dijetRefEtaLab, weight * fMcReweight );
                fHM->hRecoVsRefLeadJetEta2DOldPtBinning[refPtAveOldBin]->Fill( etaRecoLeadLab, etaRefLeadLab, 1. );
                fHM->hRecoVsRefLeadJetEta2DOldPtBinning[refPtAveOldBin]->Fill( etaRecoLeadLab, etaRefLeadLab, weight * fMcReweight );
                fHM->hRecoVsRefSubLeadJetEta2DOldPtBinning[refPtAveOldBin]->Fill( etaRecoSubLeadLab, etaRefSubLeadLab, 1. );
                fHM->hRecoVsRefSubLeadJetEta2DOldPtBinning[refPtAveOldBin]->Fill( etaRecoSubLeadLab, etaRefSubLeadLab, weight * fMcReweight );

                (dijetRefEtaLab >= 0) ? fHM->hRefDijetEtaForward1DOldPtBinning[refPtAveOldBin]->Fill(dijetRefEtaLab, 1.) : fHM->hRefDijetEtaBackward1DOldPtBinning[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), 1.);
                (dijetRefEtaLab >= 0) ? fHM->hRefDijetEtaForward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetEtaBackward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), weight * fMcReweight);
            } // if ( refPtAveOldBin >=0 )
        } // if ( fIsMc )
    } // if ( fIsRecoDijetLabFound )

    //
    // CM frame
    // 

    bool fIsRecoDijetCMFound = isGoodDijet(ptRecoLead, recoLeadJet->eta(), ptRecoSubLead, recoSubLeadJet->eta(), TMath::Abs( dijetRecoDphi ), true);
    if ( fVerbose ) {
        std::cout << Form("Reco dijet in CM frame is %s\n", ((fIsRecoDijetCMFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze reco dijets in CM frame
    if ( fIsRecoDijetCMFound ) {

        // Flush the eta values to reflect the frame
        float etaRecoLeadCM = boostEta2CM( recoLeadJet->eta() );
        float etaRecoSubLeadCM = boostEta2CM( recoSubLeadJet->eta() );
        float dijetRecoEtaCM = 0.5 * (etaRecoLeadCM + etaRecoSubLeadCM);

        if ( fVerbose ) {
            std::cout << "Reco dijet parameters in the C.M. frame: " << std::endl;
            std::cout << Form("Reco lead pt: %5.2f eta: %5.2f phi: %5.2f", ptRecoLead, etaRecoLeadCM, phiRecoLead) << std::endl;
            std::cout << Form("Reco sublead pt: %5.2f eta: %5.2f phi: %5.2f", ptRecoSubLead, etaRecoSubLeadCM, phiRecoSubLead) << std::endl;
            std::cout << Form("dijet ptAve: %5.2f dijet eta (CM): %5.2f dijet dphi: %5.2f\n", dijetRecoPtAve, dijetRecoEtaCM, dijetRecoDphi);
        }

        fHM->hRecoEtaCMLeadEtaCMSublead->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight );
        fHM->hRecoDijetEtaCM->Fill( dijetRecoEtaCM, weight * fMcReweight);
        fHM->hRecoDijetPtEtaDphiCM->Fill( dijetRecoPtAve, dijetRecoEtaCM, dijetRecoDphi, 1. );
        fHM->hRecoDijetPtEtaDphiCMWeighted->Fill( dijetRecoPtAve, dijetRecoEtaCM, dijetRecoDphi, weight * fMcReweight);
        (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaCMForward->Fill(dijetRecoPtAve, dijetRecoEtaCM, 1.) : fHM->hRecoDijetPtEtaCMBackward->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaCM), 1.);
        (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaCMForwardWeighted->Fill(dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetPtEtaCMBackwardWeighted->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);

        // Find exact dijet pT bin
        int recoPtAveBin = findDijetPtAveBin( dijetRecoPtAve );
        int recoPtAveOldBin = findDijetPtAveOldBin( dijetRecoPtAve );

        // New ptAve and eta binning
        if ( recoPtAveBin >=0 ) {
            fHM->hRecoDijetEta1DCM[recoPtAveBin]->Fill( dijetRecoEtaCM, 1. );
            fHM->hRecoDijetEta1DCMWeighted[recoPtAveBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

            fHM->hRecoEtaLeadVsEtaSubLead2DCM[recoPtAveBin]->Fill( etaRecoLeadCM, etaRecoSubLeadCM, 1. );
            fHM->hRecoEtaLeadVsEtaSubLead2DCMWeighted[recoPtAveBin]->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight * fMcReweight );

            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1D[recoPtAveBin]->Fill(dijetRecoEtaCM, 1.) : fHM->hRecoDijetEtaCMBackward1D[recoPtAveBin]->Fill(TMath::Abs(dijetRecoEtaCM), 1.);
            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DWeighted[recoPtAveBin]->Fill(dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetEtaCMBackward1DWeighted[recoPtAveBin]->Fill(TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);
        } // if ( recoPtAveBin >=0 )

        // Old ptAve binning
        if ( recoPtAveOldBin >=0 ) {
            // New eta binning
            fHM->hRecoDijetEta1DOldPtCM[recoPtAveOldBin]->Fill( dijetRecoEtaCM, 1. );
            fHM->hRecoDijetEta1DOldPtCMWeighted[recoPtAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

            fHM->hRecoEtaLeadVsEtaSubLead2DOldPtCM[recoPtAveOldBin]->Fill( etaRecoLeadCM, etaRecoSubLeadCM, 1. );
            fHM->hRecoEtaLeadVsEtaSubLead2DOldPtCMWeighted[recoPtAveOldBin]->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight * fMcReweight );

            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPt[recoPtAveOldBin]->Fill(dijetRecoEtaCM, 1.) : fHM->hRecoDijetEtaCMBackward1DOldPt[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), 1.);
            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPtWeighted[recoPtAveOldBin]->Fill(dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetEtaCMBackward1DOldPtWeighted[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);

            // Old eta binning
            fHM->hRecoDijetEta1DOldPtBinningCM[recoPtAveOldBin]->Fill( dijetRecoEtaCM, 1. );
            fHM->hRecoDijetEta1DOldPtBinningCMWeighted[recoPtAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

            fHM->hRecoEtaLeadVsEtaSubLead2DOldPtBinningCM[recoPtAveOldBin]->Fill( etaRecoLeadCM, etaRecoSubLeadCM, 1. );
            fHM->hRecoEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[recoPtAveOldBin]->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight * fMcReweight );

            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPtBinning[recoPtAveOldBin]->Fill(dijetRecoEtaCM, 1.) : fHM->hRecoDijetEtaCMBackward1DOldPtBinning[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), 1.);
            (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetEtaCMForward1DOldPtBinningWeighted[recoPtAveOldBin]->Fill(dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetEtaCMBackward1DOldPtBinningWeighted[recoPtAveOldBin]->Fill(TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);
        } // if ( recoPtAveOldBin >=0 )

        // In case of MC
        if ( fIsMc ) {

            float ptRefLeadCM = refLeadJet->pt();
            float etaRefLeadCM = boostEta2CM( refLeadJet->eta() );
            float phiRefLeadCM = refLeadJet->phi();

            float ptRefSubLeadCM = refSubLeadJet->pt();
            float etaRefSubLeadCM = boostEta2CM( refSubLeadJet->eta() );
            float phiRefSubLeadCM = refSubLeadJet->phi();

            float dijetRefPtAveCM = 0.5 * (ptRefLeadCM + ptRefSubLeadCM);
            float dijetRefEtaCM= 0.5 * (etaRefLeadCM + etaRefSubLeadCM);
            float dijetRefDphiCM = deltaPhi(phiRefLeadCM, phiRefSubLeadCM);

            if ( fVerbose ) {
                std::cout << Form("<-- Ref lead pT: %5.1f eta: %5.2f phi: %5.2f\n", ptRefLeadCM, etaRefLeadCM, phiRefLeadCM);
                std::cout << Form("<-- Ref sublead pT: %5.1f eta: %5.2f phi: %5.2f\n", ptRefSubLeadCM, etaRefSubLeadCM, phiRefSubLeadCM);
                std::cout << Form("<-- Ref dijet ptAve: %5.1f eta: %5.2f dphi: %5.2f\n", dijetRefPtAveCM, dijetRefEtaCM, dijetRefDphiCM);
            }

            fHM->hRefEtaCMLeadEtaCMSublead->Fill( etaRefLeadCM, etaRefSubLeadCM, weight );

            fHM->hRefDijetEtaCM->Fill( dijetRefEtaCM, weight );
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPtAve, 1.);
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPtAve, weight * fMcReweight );
            fHM->hRefDijetPtEtaDphiCM->Fill( dijetRefPtAveCM, dijetRefEtaCM, dijetRefDphiCM, 1. );
            fHM->hRefDijetPtEtaDphiCMWeighted->Fill( dijetRefPtAveCM, dijetRefEtaCM, dijetRefDphiCM, weight * fMcReweight );
        
            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForward->Fill(dijetRefPtAveCM, dijetRefEtaCM, 1.) : fHM->hRefDijetPtEtaCMBackward->Fill(dijetRefPtAveCM, TMath::Abs(dijetRefEtaCM), 1.);
            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForwardWeighted->Fill(dijetRefPtAveCM, dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetPtEtaCMBackwardWeighted->Fill(dijetRefPtAveCM, TMath::Abs(dijetRefEtaCM), weight * fMcReweight);

            // Find exact dijet ptAve bin
            int refPtAveBin = findDijetPtAveBin( dijetRefPtAveCM );
            int refPtAveOldBin = findDijetPtAveOldBin( dijetRefPtAveCM );

            // New ptAve and eta binning
            if ( refPtAveBin >=0 ) {
                fHM->hRefDijetEta1DCM[refPtAveBin]->Fill( dijetRefEtaCM, 1. );
                fHM->hRefDijetEta1DCMWeighted[refPtAveBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                fHM->hRefEtaLeadVsEtaSubLead2DCM[refPtAveBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, 1. );
                fHM->hRefEtaLeadVsEtaSubLead2DCMWeighted[refPtAveBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );

                fHM->hRecoVsRefDijetEta2DCM[refPtAveBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, 1. );
                fHM->hRecoVsRefDijetEta2DCMWeighted[refPtAveBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, weight * fMcReweight );
                fHM->hRecoVsRefLeadJetEta2DCM[refPtAveBin]->Fill( etaRecoLeadCM, etaRefLeadCM, 1. );
                fHM->hRecoVsRefLeadJetEta2DCMWeighted[refPtAveBin]->Fill( etaRecoLeadCM, etaRefLeadCM, weight * fMcReweight );
                fHM->hRecoVsRefSubLeadJetEta2DCM[refPtAveBin]->Fill( etaRecoSubLeadCM, etaRefSubLeadCM, 1. );
                fHM->hRecoVsRefSubLeadJetEta2DCMWeighted[refPtAveBin]->Fill( etaRecoSubLeadCM, etaRefSubLeadCM, weight * fMcReweight );

                (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1D[refPtAveBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefDijetEtaCMBackward1D[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DWeighted[refPtAveBin]->Fill(dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetEtaCMBackward1DWeighted[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight);
            }

            // Old ptAve binning
            if ( refPtAveOldBin >=0 ) {
                // New eta binning
                fHM->hRefDijetEta1DOldPtCM[refPtAveOldBin]->Fill( dijetRefEtaCM, 1. );
                fHM->hRefDijetEta1DOldPtCMWeighted[refPtAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                fHM->hRefEtaLeadVsEtaSubLead2DOldPtCM[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, 1. );
                fHM->hRefEtaLeadVsEtaSubLead2DOldPtCMWeighted[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );

                fHM->hRecoVsRefDijetEta2DOldPtCM[refPtAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, 1. );
                fHM->hRecoVsRefDijetEta2DOldPtCMWeighted[refPtAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, weight * fMcReweight );
                fHM->hRecoVsRefLeadJetEta2DOldPtCM[refPtAveOldBin]->Fill( etaRecoLeadCM, etaRefLeadCM, 1. );
                fHM->hRecoVsRefLeadJetEta2DOldPtCMWeighted[refPtAveOldBin]->Fill( etaRecoLeadCM, etaRefLeadCM, weight * fMcReweight );
                fHM->hRecoVsRefSubLeadJetEta2DOldPtCM[refPtAveOldBin]->Fill( etaRecoSubLeadCM, etaRefSubLeadCM, 1. );
                fHM->hRecoVsRefSubLeadJetEta2DOldPtCMWeighted[refPtAveOldBin]->Fill( etaRecoSubLeadCM, etaRefSubLeadCM, weight * fMcReweight );

                (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPt[refPtAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefDijetEtaCMBackward1DOldPt[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPtWeighted[refPtAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetEtaCMBackward1DOldPtWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight);

                // Old eta binning
                fHM->hRefDijetEta1DOldPtBinningCM[refPtAveOldBin]->Fill( dijetRefEtaCM, 1. );
                fHM->hRefDijetEta1DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

                fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinningCM[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, 1. );
                fHM->hRefEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );

                fHM->hRecoVsRefDijetEta2DOldPtBinningCM[refPtAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, 1. );
                fHM->hRecoVsRefDijetEta2DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( dijetRecoEtaCM, dijetRefEtaCM, weight * fMcReweight );
                fHM->hRecoVsRefLeadJetEta2DOldPtBinningCM[refPtAveOldBin]->Fill( etaRecoLeadCM, etaRefLeadCM, 1. );
                fHM->hRecoVsRefLeadJetEta2DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( etaRecoLeadCM, etaRefLeadCM, weight * fMcReweight );
                fHM->hRecoVsRefSubLeadJetEta2DOldPtBinningCM[refPtAveOldBin]->Fill( etaRecoSubLeadCM, etaRefSubLeadCM, 1. );
                fHM->hRecoVsRefSubLeadJetEta2DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( etaRecoSubLeadCM, etaRefSubLeadCM, weight * fMcReweight );

                (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPtBinning[refPtAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefDijetEtaCMBackward1DOldPtBinning[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
                (dijetRefEtaCM >= 0) ? fHM->hRefDijetEtaCMForward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetEtaCMBackward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight);
            } // if ( refPtAveOldBin >=0 )
        } // if ( fIsMc )
    } // if ( fIsRecoDijetCMFound )

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processRecoDijets -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::processRefDijets(const Event* event, const double &weight) {
    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::processRefDijets -- begin" << std::endl;
    }

    if ( weight <= 0. ) {
        std::cerr << "Error: weight is not positive: " << weight << std::endl;
        return;
    }

    fMcReweight = {1.};
    fIsRefSelDijetLabFound = {false};
    fIsRefSelDijetCMFound = {false};

    // ptHat value
    float ptHat = event->ptHat();

    if ( fRefSelRecoIdLead < 0 || fRefSelRecoIdSubLead < 0 ) {
        if ( fVerbose ) {
            std::cout << "Ref-selected reco dijet not found\n";
            std::cout << "DiJetAnalysis::processRefDijets -- end" << std::endl;
        }
        return;
    }

    // Retrieve leading and subleading jets
    RecoJet *recoLeadJet = event->recoJetCollection()->at( fRefSelRecoIdLead );
    RecoJet *recoSubLeadJet = event->recoJetCollection()->at( fRefSelRecoIdSubLead );
    GenJet *refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
    GenJet *refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );

    // Retrieve kinematic information for reference leading and subleading jets, dijets
    float ptRefLead = refLeadJet->pt();
    float etaRefLead = refLeadJet->eta();
    float phiRefLead = refLeadJet->phi();

    float ptRefSubLead = refSubLeadJet->pt();
    float etaRefSubLead = refSubLeadJet->eta();
    float phiRefSubLead = refSubLeadJet->phi();

    float dijetRefPtAve = 0.5 * (ptRefLead + ptRefSubLead);
    float dijetRefDphi = deltaPhi(phiRefLead, phiRefSubLead);

    // if ( fVerbose ) {
    //     std::cout << Form("Ref leading jet pt: %5.2f eta: %5.2f phi: %5.2f --> Reco pt: %5.2f eta: %5.2f phi: %5.2f\n", ptRefLead, etaRefLead, phiRefLead, recoLeadJet->ptJECCorr(), recoLeadJet->eta(), recoLeadJet->phi());
    //     std::cout << Form("Ref subleading jet pt: %5.2f eta: %5.2f phi: %5.2f --> Reco pt: %5.2f eta: %5.2f phi: %5.2f\n", ptRefSubLead, etaRefSubLead, phiRefSubLead, recoSubLeadJet->ptJECCorr(), recoSubLeadJet->eta(), recoSubLeadJet->phi());
    //     std::cout << Form("Ref dijet ptAve: %5.2f dphi: %5.2f\n", dijetRefPtAve, dijetRefDphi);
    // }

    //
    // Lab frame
    //

    bool fIsRefSelDijetLabFound = isGoodDijet(ptRefLead, refLeadJet->eta(), ptRefSubLead, refSubLeadJet->eta(), TMath::Abs( dijetRefDphi ), false);
    if ( fVerbose ) {
        std::cout << Form("Ref dijet in lab frame is %s\n", ((fIsRefSelDijetLabFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze ref-selected dijets in lab frame
    if ( fIsRefSelDijetLabFound ) {

        // Flush the eta values of reference jets to reflect the frame
        float etaRefLeadLab = etaLab( refLeadJet->eta() );
        float etaRefSubLeadLab = etaLab( refSubLeadJet->eta() );
        float dijetRefEtaLab = 0.5 * (etaRefLeadLab + etaRefSubLeadLab);

        // Set parameters for the reco partners
        float ptRecoLeadLeadLab = recoLeadJet->ptJECCorr();
        float etaRecoLeadLab = etaLab( recoLeadJet->eta() );
        float phiRecoLeadLab = recoLeadJet->phi();

        float ptRecoSubLeadLab = recoSubLeadJet->ptJECCorr();
        float etaRecoSubLead = etaLab( recoSubLeadJet->eta() );
        float phiRecoSubLead = recoSubLeadJet->phi();

        float dijetRecoPtAveLab = 0.5 * (ptRecoLeadLeadLab + ptRecoSubLeadLab);
        float dijetRecoEtaLab = 0.5 * (etaRecoLeadLab + etaRecoSubLead);
        float dijetRecoDphiLab = deltaPhi(phiRecoLeadLab, phiRecoSubLead);

        if ( fVerbose ) {
            std::cout << "Ref dijet parameters in the lab frame: " << std::endl;
            std::cout << Form("Ref lead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco lead pt: %5.2f eta: %5.2f phi: %5.2f", ptRefLead, etaRefLeadLab, phiRefLead, ptRecoLeadLeadLab, etaRecoLeadLab, phiRecoLeadLab) << std::endl;
            std::cout << Form("Ref sublead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco sublead pt: %5.2f eta: %5.2f phi: %5.2f", ptRefSubLead, etaRefSubLeadLab, phiRefSubLead, ptRecoSubLeadLab, etaRecoSubLead, phiRecoSubLead) << std::endl;
            std::cout << Form("Ref dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f --> Reco dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f\n", dijetRefPtAve, dijetRefEtaLab, dijetRefDphi, dijetRecoPtAveLab, dijetRecoEtaLab, dijetRecoDphiLab) << std::endl;
        }

        // Dijet reco vs ref for unfolding
        double dijetRecoUnfold[12] = { dijetRecoPtAveLab, dijetRecoEtaLab,
                                        ptRecoLeadLeadLab, etaRecoLeadLab,
                                        ptRecoSubLeadLab, etaRecoSubLead,
                                        dijetRefPtAve, dijetRefEtaLab,
                                        ptRefLead, etaRefLeadLab,
                                        ptRefSubLead, etaRefSubLeadLab };    

        fHM->hRefSelRecoDijetPtEtaLeadJetPtEtaSubleadJetPtEtaGenDijetPtEtaLeadPtEtaSubleadPtEtaWeighted->Fill(dijetRecoUnfold, weight * fMcReweight );
        fHM->hRefSelDijetEta->Fill(dijetRefEtaLab, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaDphi->Fill(dijetRefPtAve, dijetRefEtaLab, dijetRefDphi, 1.);
        fHM->hRefSelDijetPtEtaDphiWeighted->Fill(dijetRefPtAve, dijetRefEtaLab, dijetRefDphi, weight * fMcReweight );

        // Find exact dijet pT bins
        int refPtAveBin = findDijetPtAveBin( dijetRefPtAve );
        int refPtAveOldBin = findDijetPtAveOldBin( dijetRefPtAve );

        // New ptAve and eta binning
        if ( refPtAveBin >= 0 ) {
            fHM->hRefSelDijetEta1D[refPtAveBin]->Fill( dijetRefEtaLab, 1. );
            fHM->hRefSelDijetEta1DWeighted[refPtAveBin]->Fill( dijetRefEtaLab, weight * fMcReweight );

            fHM->hRefSelRecoDijetEta1D[refPtAveBin]->Fill( dijetRecoEtaLab, 1. );
            fHM->hRefSelRecoDijetEta1DWeighted[refPtAveBin]->Fill( dijetRecoEtaLab, weight * fMcReweight );

            fHM->hRefSelEtaLeadVsEtaSubLead2D[refPtAveBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, 1. );
            fHM->hRefSelEtaLeadVsEtaSubLead2DWeighted[refPtAveBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

            (dijetRefEtaLab >= 0) ? fHM->hRefSelDijetEtaForward1D[refPtAveBin]->Fill(dijetRefEtaLab, 1.) : fHM->hRefSelDijetEtaBackward1D[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaLab), 1.);
            (dijetRefEtaLab >= 0) ? fHM->hRefSelDijetEtaForward1DWeighted[refPtAveBin]->Fill(dijetRefEtaLab, weight * fMcReweight ) : fHM->hRefSelDijetEtaBackward1DWeighted[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaLab), weight * fMcReweight );
        } // if ( refPtAveBin >= 0 )

        // Old ptAve binning
        if ( refPtAveOldBin >= 0 ) {
            // New eta binning
            fHM->hRefSelDijetEta1DOldPt[refPtAveOldBin]->Fill( dijetRefEtaLab, 1. );
            fHM->hRefSelDijetEta1DOldPtWeighted[refPtAveOldBin]->Fill( dijetRefEtaLab, weight * fMcReweight );

            fHM->hRefSelRecoDijetEta1DOldPt[refPtAveOldBin]->Fill( dijetRecoEtaLab, 1. );
            fHM->hRefSelRecoDijetEta1DOldPtWeighted[refPtAveOldBin]->Fill( dijetRecoEtaLab, weight * fMcReweight );

            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPt[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, 1. );
            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtWeighted[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

            (dijetRefEtaLab >= 0) ? fHM->hRefSelDijetEtaForward1DOldPt[refPtAveOldBin]->Fill(dijetRefEtaLab, 1.) : fHM->hRefSelDijetEtaBackward1DOldPt[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), 1.);
            (dijetRefEtaLab >= 0) ? fHM->hRefSelDijetEtaForward1DOldPtWeighted[refPtAveOldBin]->Fill(dijetRefEtaLab, weight * fMcReweight ) : fHM->hRefSelDijetEtaBackward1DOldPtWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), weight * fMcReweight );

            // Old eta binning
            fHM->hRefSelDijetEta1DOldPtBinning[refPtAveOldBin]->Fill( dijetRefEtaLab, 1. );
            fHM->hRefSelDijetEta1DOldPtBinningWeighted[refPtAveOldBin]->Fill( dijetRefEtaLab, weight * fMcReweight );

            fHM->hRefSelRecoDijetEta1DOldPtBinning[refPtAveOldBin]->Fill( dijetRecoEtaLab, 1. );
            fHM->hRefSelRecoDijetEta1DOldPtBinningWeighted[refPtAveOldBin]->Fill( dijetRecoEtaLab, weight * fMcReweight );

            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinning[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, 1. );
            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinningWeighted[refPtAveOldBin]->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

            (dijetRefEtaLab >= 0) ? fHM->hRefSelDijetEtaForward1DOldPtBinning[refPtAveOldBin]->Fill(dijetRefEtaLab, 1.) : fHM->hRefSelDijetEtaBackward1DOldPtBinning[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), 1.);
            (dijetRefEtaLab >= 0) ? fHM->hRefSelDijetEtaForward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(dijetRefEtaLab, weight * fMcReweight ) : fHM->hRefSelDijetEtaBackward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaLab), weight * fMcReweight );
        } // if ( refPtAveOldBin >= 0 )

    } // if ( fIsRefSelDijetLabFound )

    //
    // CM frame
    //

    bool fIsRefSelDijetCMFound = isGoodDijet(ptRefLead, refLeadJet->eta(), ptRefSubLead, refSubLeadJet->eta(), TMath::Abs( dijetRefDphi ), true);
    if ( fVerbose ) {
        std::cout << Form("Ref dijet in CM frame is %s\n", ((fIsRefSelDijetCMFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze ref-selected dijets in CM frame
    if ( fIsRefSelDijetCMFound ) {

        // Flush the eta values of reference jets to reflect the frame
        float etaRefLeadCM = boostEta2CM( refLeadJet->eta() );
        float etaRefSubLeadCM = boostEta2CM( refSubLeadJet->eta() );
        float dijetRefEtaCM = 0.5 * (etaRefLeadCM + etaRefSubLeadCM);

        // Set parameters for the reco partners
        float ptRecoLeadCM = recoLeadJet->ptJECCorr();
        float etaRecoLeadCM = boostEta2CM( recoLeadJet->eta() );
        float phiRecoLeadCM = recoLeadJet->phi();

        float ptRecoSubLeadCM = recoSubLeadJet->ptJECCorr();
        float etaRecoSubLeadCM = boostEta2CM( recoSubLeadJet->eta() );
        float phiRecoSubLeadCM = recoSubLeadJet->phi();

        float dijetRecoPtAveCM = 0.5 * (ptRecoLeadCM + ptRecoSubLeadCM);
        float dijetRecoEtaCM = 0.5 * (etaRecoLeadCM + etaRecoSubLeadCM);
        float dijetRecoDphiCM = deltaPhi(phiRecoLeadCM, phiRecoSubLeadCM);

        if ( fVerbose ) {
            std::cout << "Ref dijet parameters in the C.M. frame: " << std::endl;
            std::cout << Form("Ref lead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco lead pt: %5.2f eta: %5.2f phi: %5.2f", ptRefLead, etaRefLeadCM, phiRefLead, ptRecoLeadCM, etaRecoLeadCM, phiRecoLeadCM) << std::endl;
            std::cout << Form("Ref sublead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco sublead pt: %5.2f eta: %5.2f phi: %5.2f", ptRefSubLead, etaRefSubLeadCM, phiRefSubLead, ptRecoSubLeadCM, etaRecoSubLeadCM, phiRecoSubLeadCM) << std::endl;
            std::cout << Form("Ref dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f --> Reco dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f\n", dijetRefPtAve, dijetRefEtaCM, dijetRefDphi, dijetRecoPtAveCM, dijetRecoEtaCM, dijetRecoDphiCM) << std::endl;
        }

        fHM->hRefSelDijetEtaCM->Fill(dijetRefEtaCM, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaDphiCM->Fill(dijetRefPtAve, dijetRefEtaCM, dijetRefDphi, 1.);
        fHM->hRefSelDijetPtEtaDphiCMWeighted->Fill(dijetRefPtAve, dijetRefEtaCM, dijetRefDphi, weight * fMcReweight );

        // Find exact dijet pT bin
        int refPtAveBin = findDijetPtAveBin( dijetRefPtAve );
        int refPtAveOldBin = findDijetPtAveOldBin( dijetRefPtAve );

        // New ptAve and eta binning
        if ( refPtAveBin >= 0 ) {
            fHM->hRefSelDijetEta1DCM[refPtAveBin]->Fill( dijetRefEtaCM, 1. );
            fHM->hRefSelDijetEta1DCMWeighted[refPtAveBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

            fHM->hRefSelRecoDijetEta1DCM[refPtAveBin]->Fill( dijetRecoEtaCM, 1. );
            fHM->hRefSelRecoDijetEta1DCMWeighted[refPtAveBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

            fHM->hRefSelEtaLeadVsEtaSubLead2DCM[refPtAveBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, 1. );
            fHM->hRefSelEtaLeadVsEtaSubLead2DCMWeighted[refPtAveBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );

            (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1D[refPtAveBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefSelDijetEtaCMBackward1D[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
            (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DWeighted[refPtAveBin]->Fill(dijetRefEtaCM, weight * fMcReweight ) : fHM->hRefSelDijetEtaCMBackward1DWeighted[refPtAveBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight );
        } // if ( refPtAveBin >= 0 )

        // Old ptAve binning
        if ( refPtAveOldBin >= 0 ) {
            // New eta binning
            fHM->hRefSelDijetEta1DOldPtCM[refPtAveOldBin]->Fill( dijetRefEtaCM, 1. );
            fHM->hRefSelDijetEta1DOldPtCMWeighted[refPtAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

            fHM->hRefSelRecoDijetEta1DOldPtCM[refPtAveOldBin]->Fill( dijetRecoEtaCM, 1. );
            fHM->hRefSelRecoDijetEta1DOldPtCMWeighted[refPtAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtCM[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, 1. );
            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtCMWeighted[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );

            (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPt[refPtAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefSelDijetEtaCMBackward1DOldPt[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
            (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPtWeighted[refPtAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight ) : fHM->hRefSelDijetEtaCMBackward1DOldPtWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight );

            // Old eta binning
            fHM->hRefSelDijetEta1DOldPtBinningCM[refPtAveOldBin]->Fill( dijetRefEtaCM, 1. );
            fHM->hRefSelDijetEta1DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( dijetRefEtaCM, weight * fMcReweight );

            fHM->hRefSelRecoDijetEta1DOldPtBinningCM[refPtAveOldBin]->Fill( dijetRecoEtaCM, 1. );
            fHM->hRefSelRecoDijetEta1DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( dijetRecoEtaCM, weight * fMcReweight );

            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCM[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, 1. );
            fHM->hRefSelEtaLeadVsEtaSubLead2DOldPtBinningCMWeighted[refPtAveOldBin]->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );

            (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPtBinning[refPtAveOldBin]->Fill(dijetRefEtaCM, 1.) : fHM->hRefSelDijetEtaCMBackward1DOldPtBinning[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), 1.);
            (dijetRefEtaCM >= 0) ? fHM->hRefSelDijetEtaCMForward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(dijetRefEtaCM, weight * fMcReweight ) : fHM->hRefSelDijetEtaCMBackward1DOldPtBinningWeighted[refPtAveOldBin]->Fill(TMath::Abs(dijetRefEtaCM), weight * fMcReweight );
        } // if ( refPtAveOldBin >= 0 )
    }

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processRefJets -- end" << std::endl;
    }
}

//________________
bool DiJetAnalysis::isGoodDijet(const float& ptLead, const float& etaLead, 
                                const float& ptSubLead, const float& etaSubLead, 
                                const float& dphi, const bool& isCM) {

    // if ( fVerbose ) {
    //     std::cout << "\nDiJetAnalysis::isGoodDijet -- begin" << std::endl;
    // }

    float eta1 = ( isCM ) ? boostEta2CM(etaLead) : etaLab(etaLead);
    float eta2 = ( isCM ) ? boostEta2CM(etaSubLead) : etaLab(etaSubLead);
    float etaCut[2] = { fJetEtaLab[0], fJetEtaLab[1] };
    if ( isCM ) {
        etaCut[0] = fJetEtaCM[0];
        etaCut[1] = fJetEtaCM[1];
    }
    bool isGood = ( ptLead > fLeadJetPtLow &&
                    etaCut[0] <= eta1 && eta1 < etaCut[1] &&
                    ptSubLead > fSubleadJetPtLow && 
                    etaCut[0] <= eta2 && eta2 < etaCut[1] &&
                    TMath::Abs( dphi ) > fDijetPhiCut );

    fMcReweight = {1.};
    // // Check reweight
    // if ( fIsMc && fUseMcReweighting != 0 ) {
    //     findMcWeight(ptLead, ptSubLead);
    // }
    // else {
    //     fMcReweight = {1.};
    // }

    // if ( fVerbose ) {
    //     std::cout << Form("Dijet status: %s\n", ( (isGood) ? "[good]" : "[bad]" ) );
    //     std::cout << Form("Leading jet pT %5.2f > %5.2f GeV: \t%s\n", ptLead, fLeadJetPtLow, ( (ptLead > fLeadJetPtLow) ? "[good]" : "[bad]" ) );
    //     std::cout << Form("Subleading jet pT %5.2f > %5.2f GeV: \t%s\n", ptSubLead, fSubleadJetPtLow, ( (ptSubLead > fSubleadJetPtLow) ? "[good]" : "[bad]" ) );
    //     std::cout << Form("Leading jet eta %3.2f <= %3.2f < %3.2f: \t%s\n", etaCut[0], eta1, etaCut[1], ( (etaCut[0] <= eta1 && eta1 < etaCut[1]) ? "[good]" : "[bad]" ) );
    //     std::cout << Form("Subleading jet eta %3.2f <= %3.2f < %3.2f: \t%s\n", etaCut[0], eta2, etaCut[1], ( (etaCut[0] <= eta2 && eta2 < etaCut[1]) ? "[good]" : "[bad]" ) );
    //     std::cout << Form("Delta phi %3.2f > %3.2f: \t%s\n", TMath::Abs( dphi ), fDijetPhiCut, ( (TMath::Abs( dphi ) > fDijetPhiCut) ? "[good]" : "[bad]" ) );
    //     std::cout << Form("Reweighting factor: %5.2f\n", fMcReweight);
    //     std::cout << "DiJetAnalysis::isGoodDijet -- end" << std::endl;
    // }
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
    fRecoIdLead = {-1}; 
    fRecoIdSubLead = {-1}; 
    fGenIdLead = {-1}; 
    fGenIdSubLead = {-1}; 
    fRefSelRecoIdLead = {-1}; 
    fRefSelRecoIdSubLead = {-1};

    fIsGenDijetLabFound = {false};
    fIsGenDijetCMFound = {false};
    fIsRecoDijetLabFound = {false};
    fIsRecoDijetCMFound = {false};
    fIsRefSelDijetLabFound = {false};
    fIsRefSelDijetCMFound = {false};
    fRecoPtSortedJetIds.clear();
    fGenPtSortedJetIds.clear();
    fRefSelRecoPtSortedJetIds.clear();

    //
    // Event quantities
    //

    // ptHat
    float ptHat = event->ptHat();
    // Vertex z position
    float vz = event->vz();
    // ptHat weight 
    float ptHatW = event->ptHatWeight();
    // Centrality
    float centrality = event->centrality();
    // Centrality weight
    float centW = event->centralityWeight();
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

    // Fill event histograms
    fHM->hVz->Fill( vz,  1. );
    fHM->hVzWeighted->Fill( vz, weight );

    fHM->hPtHat->Fill( ptHat, 1. );
    fHM->hPtHatWeighted->Fill( ptHat, weight );

    fHM->hHiBin->Fill( event->hiBin(), 1. );
    fHM->hHiBinWeighted->Fill( event->hiBin(), weight );

    //
    // Create pT-sorted jet indices
    //
    makePtSortedJetVectors( event );

    // Check for event overweight in both reco and gen jets 
    // (x-jets and purelly overweighted gen jets) in Monte Carlo
    if ( fIsMc ) {
        bool overweight = isOverweightedEvent( event, weight );
        if ( overweight ) {
            if ( fVerbose ) {
                std::cout << "Overweighted event. Skip it." << std::endl;
            }
            return;
        }
    } // if ( fIsMc )

    //
    // Loop over jet collections (reco, gen, ref)
    //
    processInclusiveJets(event, weight);

    // Process and analyze reco dijets
    processRecoDijets(event, weight);

    // Process and analyze MC jets
    if ( fIsMc ) {
        
        // Process and analyze gen dijets
        processGenDijets(event, weight);
        processRefDijets(event, weight);
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
