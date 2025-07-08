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
#include <unordered_map>
#include <cmath>

// Jet analysis headers
#include "JetCut.h"
#include "DiJet.h"
#include "DiJetCut.h"
#include "DiJetAnalysis.h"

//________________
DiJetAnalysis::DiJetAnalysis() : BaseAnalysis(), 
    fVzWeight{nullptr}, fDijetPtAveWeight{nullptr},
    fUseCentralityWeight{}, fHM{nullptr},
    fEtaShift{0.465}, fIsMc{false}, fCollisionSystem{1}, fCollisionEnergy{8160},
    fIsPbGoingDir{false}, fVerbose{false},
    fNEventsInSample{1000000},
    fIsGenDijetLabFound{false}, fIsGenDijetCMFound{false},
    fIsRecoDijetLabFound{false}, fIsRecoDijetCMFound{false},
    fIsRefSelDijetLabFound{false}, fIsRefSelDijetCMFound{false},
    fUseMcReweighting{0}, fJetPtBins{75},
    fMcReweight{1}, 
    fRecoIdLead{-1}, fRecoIdSubLead{-1}, fGenIdLead{-1}, fGenIdSubLead{-1}, fRefSelRecoIdLead{-1}, fRefSelRecoIdSubLead{-1},
    fRecoPtSortedJetIds{}, fGenPtSortedJetIds{}, fRefSelRecoPtSortedJetIds{},
    fRecoDijet{nullptr}, fGenDijet{nullptr}, fRefDijet{nullptr},
    fRecoJetCut{nullptr}, fGenJetCut{nullptr}, fDiJetCut{nullptr},
    fPtAveBins{}, fPtAveOldBins{} {

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

    fRecoDijet = new DiJet{};
    fGenDijet = new DiJet{};
    fRefDijet = new DiJet{};
}

//________________
DiJetAnalysis::~DiJetAnalysis() {
    if (fHM) { delete fHM; fHM = nullptr; }
    if (fRecoDijet) { delete fRecoDijet; fRecoDijet = nullptr; }
    if (fGenDijet) { delete fGenDijet; fGenDijet = nullptr; }
    if (fRefDijet) { delete fRefDijet; fRefDijet = nullptr; }
    if (fRecoJetCut) { delete fRecoJetCut; fRecoJetCut = nullptr; }
    if (fGenJetCut) { delete fGenJetCut; fGenJetCut = nullptr; }
    if (fDiJetCut) { delete fDiJetCut; fDiJetCut = nullptr; }
    if (fVzWeight) { delete fVzWeight; fVzWeight = nullptr; }
    if (fDijetPtAveWeight) { delete fDijetPtAveWeight; fDijetPtAveWeight = nullptr; }
}

//________________
void DiJetAnalysis::init() {
    // Initialize analysis
    if ( fVerbose ) {
        std::cout << "\nDiJetAnalysis::init -- begin" << std::endl;
    }

    // Print analysis setup
    print();

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
              << "ptHat range                 : " << fPtHatRange[0] << "-" << fPtHatRange[1] << std::endl;
              if ( fRecoJetCut ) {
                  std::cout << "Reco jet cut parameters     : " << std::endl;
                  fRecoJetCut->report();
              }
              if ( fGenJetCut ) {
                  std::cout << "Gen jet cut parameters      : " << std::endl;
                  fGenJetCut->report();
              }
              if ( fDiJetCut ) {
                  std::cout << "Di-jet cut parameters       : " << std::endl;
                  fDiJetCut->report();
              }
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
void DiJetAnalysis::findMcWeight(const float& ptLead, const float& ptSublead) {
    fMcReweight = {1};
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

    fHM->hRecoJetCollectionSize->Fill( event->recoJetCollection()->size() );
    fHM->hGenVsRecoJetCollectionSize->Fill( event->recoJetCollection()->size(), event->genJetCollection()->size() );

    // Jet counter
    int recoJetCounter{0};
    // std::cout << "Reco jet collection size: " << event->recoJetCollection()->size() << std::endl;
    // Loop over reconstructed jets and store indices of good jets
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        recoJetCounter++;

        // Check if reconstructed jet passes selection criteria (*jet, isCM, isMC, requireMatching)
        if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, true, false) ) continue; 

        fRecoPtSortedJetIds.push_back( recoJetCounter-1 );
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

    // Sort indices based on the jet corrected pT (from high to low)
        std::sort( fRecoPtSortedJetIds.begin(), fRecoPtSortedJetIds.end(), [&](int i, int j) { 
            return event->recoJetCollection()->at(i)->ptJECCorr() > event->recoJetCollection()->at(j)->ptJECCorr(); 
    } );

    if ( fVerbose ) {
        // Print sorted indices and corresponding jet pT    
        for (const auto& id : fRecoPtSortedJetIds) {
            std::cout << Form("Sorted reco jet index: %d | pT: %5.1f eta: %3.2f\n", id, event->recoJetCollection()->at(id)->ptJECCorr(), etaLab(event->recoJetCollection()->at(id)->eta()));
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

        fHM->hGenJetCollectionSize->Fill( event->genJetCollection()->size() );

        int genJetCounter{0};
        // Loop over generated jets and store indices of good jets
        for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {
            genJetCounter++;

            // Check gen jet passes the selection criteria (*genJet, isCM)
            if ( fGenJetCut && !fGenJetCut->pass(*genJetIter, false) ) continue;

            fGenPtSortedJetIds.push_back( genJetCounter-1 );
        } // for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ )

        // Sort indices based on the jet pT (from high to low)
        std::sort( fGenPtSortedJetIds.begin(), fGenPtSortedJetIds.end(), [&](int i, int j) { 
            return event->genJetCollection()->at(i)->pt() > event->genJetCollection()->at(j)->pt(); 
        } );

        if ( fVerbose ) {
            // Print sorted indices and corresponding jet pT
            for (const auto& id : fGenPtSortedJetIds) {
                std::cout << Form("Sorted gen jet index: %d | pT: %5.1f eta: %3.2f\n", id, event->genJetCollection()->at(id)->pt(), etaLab(event->genJetCollection()->at(id)->eta()));
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

            // Check selection criteria (*recoJet, isCM, isMC, requireMatching)
            if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, true, true) ) continue;

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
void DiJetAnalysis::processDijets(const Event* event, const double &weight) {
    // Process and analyze MC dijets
    if ( fIsMc ) {
        // Process and analyze gen dijets
        processGenDijets(event, weight);
        processRefDijets(event, weight);
    }

    // Process and analyze reco dijets  
    processRecoDijets(event, weight);
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
    float etaGenLeadLab = etaLab( leadJet->eta() );
    float etaGenLeadCM = boostEta2CM( leadJet->eta() );

    GenJet* subLeadJet = event->genJetCollection()->at( fGenIdSubLead );
    float ptGenSubLead = subLeadJet->pt();
    float phiGenSubLead = subLeadJet->phi();
    float etaGenSubLeadLab = etaLab( subLeadJet->eta() );
    float etaGenSubLeadCM = boostEta2CM( subLeadJet->eta() );
    // if ( fVerbose ) {
    //     std::cout << Form("Gen lead jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", 
    //                       ptGenLead, etaGenLeadLab, etaGenLeadCM, phiGenLead);
    //     std::cout << Form("Gen sublead jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", 
    //                       ptGenSubLead, etaGenSubLeadLab, etaGenSubLeadCM, phiGenSubLead);
    // }

    // Create dijet
    fGenDijet->cleanParameters();
    fGenDijet->setLeadJetPt( ptGenLead );
    fGenDijet->setLeadJetEtaLab( etaGenLeadLab );
    fGenDijet->setLeadJetEtaCM( etaGenLeadCM );
    fGenDijet->setLeadJetPhi( phiGenLead );
    fGenDijet->setSubLeadJetPt( ptGenSubLead );
    fGenDijet->setSubLeadJetEtaLab( etaGenSubLeadLab );
    fGenDijet->setSubLeadJetEtaCM( etaGenSubLeadCM );
    fGenDijet->setSubLeadJetPhi( phiGenSubLead );

    // if ( fVerbose ) {
    //     std::cout << "Inclusive gen dijet:\n"; 
    //     dijetGen.print();
    // }

    float dijetGenPtAve = fGenDijet->ptAve();
    float dijetGenDphi = fGenDijet->dPhi();
    float dijetGenEtaLab = fGenDijet->etaLab();
    float dijetGenEtaCM = fGenDijet->etaCM();
    float dijetGenDetaCM = fGenDijet->dEtaCM();
    float dijetGenPhi = fGenDijet->phi();

    if ( fabs(ptGenLead - fGenDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(ptGenSubLead - fGenDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaGenLeadLab - fGenDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaGenSubLeadLab - fGenDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaGenLeadCM - fGenDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaGenSubLeadCM - fGenDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiGenLead - fGenDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiGenSubLead - fGenDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
        std::cout << "Error: pt of lead or sublead gen jet does not match the dijet pt within rounding" << std::endl;
        return;
    }

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
    fHM->hGenInclusiveDijetEtaDetaCMPt->Fill( dijetGenEtaLab, dijetGenDetaCM, dijetGenPtAve, 1. );
    fHM->hGenInclusiveDijetEtaDetaCMPtWeighted->Fill( dijetGenEtaLab, dijetGenDetaCM, dijetGenPtAve, weight );
    fHM->hGenInclusiveDijetXPb->Fill( x_Pb, 1. );
    fHM->hGenInclusiveDijetXPbWeighted->Fill( x_Pb, weight );
    fHM->hGenInclusiveDijetXp->Fill( x_p, 1. );
    fHM->hGenInclusiveDijetXpWeighted->Fill( x_p, weight );
    fHM->hGenInclusiveDijetXPbOverXp->Fill( xPbOverXp, 1. );
    fHM->hGenInclusiveDijetXPbOverXpWeighted->Fill( xPbOverXp, weight );
    fHM->hGenInclusiveDijetXPbOverXpEta->Fill( xPbOverXp, dijetGenEtaLab, 1. );
    fHM->hGenInclusiveDijetXPbOverXpEtaWeighted->Fill( xPbOverXp, dijetGenEtaLab, weight );

    double genDijetLeadSubLead[11] {dijetGenPtAve, dijetGenEtaLab, dijetGenDetaCM, dijetGenDphi, dijetGenPhi,
                                    ptGenLead, etaGenLeadLab, phiGenLead, 
                                    ptGenSubLead, etaGenSubLeadLab, phiGenSubLead };

    double genDijetLeadSubLeadCM[11] {dijetGenPtAve, dijetGenEtaCM, dijetGenDetaCM, dijetGenDphi, dijetGenPhi,
                                      ptGenLead, etaGenLeadCM, phiGenLead, 
                                      ptGenSubLead, etaGenSubLeadCM, phiGenSubLead };

    //
    // Lab frame
    //

    bool fIsGenDijetLabFound {false};
    if ( !fDiJetCut ) {
        fIsGenDijetLabFound = true; // No cut, so dijet is always found
    } 
    else {
        // Check if the dijet passes the cut (dijetGen, isCM)
        fIsGenDijetLabFound = fDiJetCut->pass(fGenDijet, false);
    }
    if ( fVerbose ) {
        std::cout << Form("Gen dijet in lab frame is %s\n", ((fIsGenDijetLabFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze gen dijets in lab frame
    if ( fIsGenDijetLabFound ) {

        // Flush the eta values to reflect the frame

        if ( fVerbose ) {
            std::cout << "Gen dijet parameters in the lab frame\n";
            std::cout << Form("--> Lead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", 
                              fGenDijet->leadJetPt(), fGenDijet->leadJetEtaLab(), fGenDijet->leadJetPhi());
            std::cout << Form("--> Sublead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", 
                              fGenDijet->subLeadJetPt(), fGenDijet->subLeadJetEtaLab(), fGenDijet->subLeadJetPhi());
            std::cout << Form("--> Dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f delta eta: %5.2f x_Pb: %5.2f x_p: %5.2f xPbOverXp: %5.2f\n", 
                              dijetGenPtAve, dijetGenEtaLab, dijetGenDphi, dijetGenDetaCM, x_Pb, x_p, xPbOverXp) << std::endl;
        }

        // Distributions for selected gen dijets in the lab frame
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

        // Fill the histogram with the dijet parameters [11 dimensions]
        fHM->hGenDijetLeadSubLead->Fill(genDijetLeadSubLead, weight);
        
        // Lead vs. SubLead
        fHM->hGenPtLeadPtSublead->Fill( ptGenLead, ptGenSubLead, weight );
        fHM->hGenEtaLeadEtaSublead->Fill( etaGenLeadLab, etaGenSubLeadLab, weight );
        // fHM->hGenPtLeadPtSubleadMcReweight->Fill( ptGenLead, ptGenSubLead, weight * fMcReweight );
        // fHM->hGenEtaLeadEtaSubleadMcReweight->Fill( etaGenLeadLab, etaGenSubLeadLab, weight * fMcReweight );

        // Dijet distributions
        fHM->hGenDijetEta->Fill(dijetGenEtaLab, weight * fMcReweight );
        fHM->hGenDijetPtEtaPhi->Fill(dijetGenPtAve, dijetGenEtaLab, dijetGenPhi, 1.);
        fHM->hGenDijetPtEtaPhiWeighted->Fill(dijetGenPtAve, dijetGenEtaLab, dijetGenPhi, weight * fMcReweight );
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

    bool fIsGenDijetCMFound{false};
    if ( !fDiJetCut ) {
        fIsGenDijetCMFound = true; // No cut, so dijet is always found
    } 
    else {
        // Check if the dijet passes the cut (dijetGen, isCM)
        fIsGenDijetCMFound = fDiJetCut->pass(fGenDijet, true);
    }
    if ( fVerbose ) {
        std::cout << Form("Gen dijet in CM frame is %s\n", ((fIsGenDijetCMFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze gen dijets in CM frame
    if ( fIsGenDijetCMFound ) {

        if ( fVerbose ) {
            std::cout << "Gen dijet parameters in the C.M. frame\n";
            std::cout << Form("--> Lead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", 
                              fGenDijet->leadJetPt(), fGenDijet->leadJetEtaCM(), fGenDijet->leadJetPhi());
            std::cout << Form("--> Sublead jet pt: %5.2f eta: %5.2f phi: %5.2f\n", 
                              fGenDijet->subLeadJetPt(), fGenDijet->subLeadJetEtaCM(), fGenDijet->subLeadJetPhi());
            std::cout << Form("--> Dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f delta eta: %5.2f x_Pb: %5.2f x_p: %5.2f xPbOverXp: %5.2f\n", 
                              dijetGenPtAve, dijetGenEtaCM, dijetGenDphi, dijetGenDetaCM, x_Pb, x_p, xPbOverXp);
        }

        fHM->hGenDijetLeadSubLeadCM->Fill(genDijetLeadSubLeadCM, weight * fMcReweight );

        // Lead and SubLead
        fHM->hGenPtLeadPtSubleadCM->Fill( ptGenLead, ptGenSubLead, weight * fMcReweight  );
        fHM->hGenEtaCMLeadEtaCMSublead->Fill( etaGenLeadCM, etaGenSubLeadCM, weight * fMcReweight  );
        

        fHM->hGenDijetEtaCM->Fill(dijetGenEtaCM, weight * fMcReweight );
        fHM->hGenDijetPtEtaPhiCM->Fill(dijetGenPtAve, dijetGenEtaCM, dijetGenPhi, 1.);
        fHM->hGenDijetPtEtaPhiCMWeighted->Fill(dijetGenPtAve, dijetGenEtaCM, dijetGenPhi, weight * fMcReweight );
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

    fGenDijet->cleanParameters();

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
    float etaRecoLeadLab = etaLab( recoLeadJet->eta() );
    float etaRecoLeadCM = boostEta2CM( recoLeadJet->eta() );
    float phiRecoLead = recoLeadJet->phi();

    // Subleading jet
    RecoJet* recoSubLeadJet = event->recoJetCollection()->at( fRecoIdSubLead );
    float ptRecoSubLead = recoSubLeadJet->ptJECCorr();
    float etaRecoSubLeadLab = etaLab( recoSubLeadJet->eta() );
    float etaRecoSubLeadCM = boostEta2CM( recoSubLeadJet->eta() );
    float phiRecoSubLead = recoSubLeadJet->phi();


    // if ( fVerbose ) {
    //     std::cout << Form("Reco lead jet id: %d pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", 
    //                       fRecoIdLead, ptRecoLead, etaRecoLeadLab, etaRecoLeadCM, phiRecoLead);
    //     std::cout << Form("Reco sublead jet id: %d pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n",
    //                       fRecoIdSubLead, ptRecoSubLead, etaRecoSubLeadLab, etaRecoSubLeadCM, phiRecoSubLead);
    // }

    // Create reco dijet
    fRecoDijet->cleanParameters();

    fRecoDijet->setLeadJetPt( ptRecoLead );
    fRecoDijet->setLeadJetEtaLab( etaRecoLeadLab );
    fRecoDijet->setLeadJetEtaCM( etaRecoLeadCM );
    fRecoDijet->setLeadJetPhi( phiRecoLead );

    fRecoDijet->setSubLeadJetPt( ptRecoSubLead );
    fRecoDijet->setSubLeadJetEtaLab( etaRecoSubLeadLab );
    fRecoDijet->setSubLeadJetEtaCM( etaRecoSubLeadCM );
    fRecoDijet->setSubLeadJetPhi( phiRecoSubLead );

    // if ( fVerbose ) {
    //     std::cout << "Inclusive reco dijet:\n"; 
    //     dijetReco.print();
    // }

    if ( fabs( ptRecoLead - fRecoDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs( ptRecoSubLead - fRecoDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs( etaRecoLeadLab - fRecoDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs( etaRecoSubLeadLab - fRecoDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs( etaRecoLeadCM - fRecoDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs( etaRecoSubLeadCM - fRecoDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs( phiRecoLead - fRecoDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
         fabs( phiRecoSubLead - fRecoDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
        std::cout << "Error: Reco dijet parameters do not match the leading and subleading jet parameters within rounding\n";
        return;
    }

    float dijetRecoPtAve = fRecoDijet->ptAve();
    float dijetRecoDphi = fRecoDijet->dPhi();
    float dijetRecoEtaLab = fRecoDijet->etaLab();
    float dijetRecoEtaCM = fRecoDijet->etaCM();
    float dijetRecoDetaCM = fRecoDijet->dEtaCM();
    float dijetRecoPhi = fRecoDijet->phi();

    // Reference leading jet
    GenJet* refLeadJet = {nullptr};
    float ptRefLead{0.};
    float etaRefLeadLab{0.};
    float etaRefLeadCM{0.};
    float phiRefLead{0.};
    // Reference subleading jet
    GenJet* refSubLeadJet = {nullptr};
    float ptRefSubLead{0.};
    float etaRefSubLeadLab{0.};
    float etaRefSubLeadCM{0.};
    float phiRefSubLead{0.};

    // Create ref dijet
    fRefDijet->cleanParameters();
    float dijetRefPtAve{0.};
    float dijetRefEtaLab{0.};
    float dijetRefEtaCM{0.};
    float dijetRefDphi{0.};
    float dijetRefPhi{0.};

    // For Monte Carlo check if reco jets have matched gen jets. Skip reco dijet if not
    if ( fIsMc ) {

        // // Leading and subleading jets must have matching gen jets (just a protection, easy to comment out)
        // if ( !recoLeadJet->hasMatching() || !recoSubLeadJet->hasMatching() ) {
        //     // if ( fVerbose ) {
        //         GenJet* genLeadJet = event->genJetCollection()->at( fGenIdLead );
        //         GenJet* genSubLeadJet = event->genJetCollection()->at( fGenIdSubLead );
        //         std::cerr << Form("[ERROR] Unmatched dijet idLead: %d idSubLead: %d Lead matched: %s SubLead matched: %s - Skip\n", 
        //                         fRecoIdLead, fRecoIdSubLead, 
        //                         (recoLeadJet->hasMatching() ? "[true]" : "[false]"), 
        //                         (recoSubLeadJet->hasMatching() ? "[true]" : "[false]"));
        //         std::cerr << Form("Reco lead pt: %5.2f etaLab: %5.2f SubLead pt: %5.2f etaLab: %5.2f <-> Gen lead pt: %5.2f etaLab: %5.2f SubLead pt: %5.2f etaLab: %5.2f\n",
        //                             recoLeadJet->ptJECCorr(), etaRecoLeadLab, recoSubLeadJet->ptJECCorr(), etaRecoSubLeadLab,
        //                             genLeadJet ? genLeadJet->pt() : -1, genLeadJet ? etaLab(genLeadJet->eta()) : -1,
        //                             genSubLeadJet ? genSubLeadJet->pt() : -1, genSubLeadJet ? etaLab(genSubLeadJet->eta()) : -1);
        //         std::cerr << std::endl;
        //             std::cerr << "Sorted reco jets (id, ptCorr, ptRaw, eta, matched):" << std::endl;
        //             for (const auto& id : fRecoPtSortedJetIds) {
        //                 RecoJet* jet = event->recoJetCollection()->at(id);
        //                 std::cerr << Form("-->  id: %d, ptCorr: %5.2f, ptRaw: %5.2f, eta: %3.2f matched: %s", jet->id(), jet->ptJECCorr(), jet->pt(), etaLab(jet->eta()), (jet->hasMatching() ? "[true]" : "[false]")) << std::endl;
        //             }
        //             std::cerr << "Sorted gen jets (id, pt, eta):" << std::endl;
        //             for (const auto& id : fGenPtSortedJetIds) {
        //                 GenJet* jet = event->genJetCollection()->at(id);
        //                 std::cerr << Form("-->  id: %d, pt: %5.2f, eta: %3.2f", jet->id(), jet->pt(), etaLab(jet->eta())) << std::endl;
        //             }
        //     // }
        //     return;
        // }

        if ( recoLeadJet->hasMatching() ) {
            // Matching gen jet for leading reco jet
            refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
            if ( !refLeadJet ) {
                std::cerr << "Error: Leading jet has no matching gen jet\n";
                return;
            }

            ptRefLead = refLeadJet->pt();
            etaRefLeadLab = etaLab( refLeadJet->eta() );
            etaRefLeadCM = boostEta2CM( refLeadJet->eta() );
            phiRefLead = refLeadJet->phi();
        } 
        
        if ( recoSubLeadJet->hasMatching() ) {
            // Matching gen jet for subleading reco jet
            refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );
            if ( !refSubLeadJet ) {
                std::cerr << "Error: Subleading jet has no matching gen jet\n";
                return;
            }

            ptRefSubLead = refSubLeadJet->pt();
            etaRefSubLeadLab = etaLab( refSubLeadJet->eta() );
            etaRefSubLeadCM = boostEta2CM( refSubLeadJet->eta() );
            phiRefSubLead = refSubLeadJet->phi();
        }

        // Make reference dijet if both jets have matched gen jets
        if (recoLeadJet->hasMatching() && recoSubLeadJet->hasMatching()) {
            fRefDijet->setLeadJetPt( ptRefLead );
            fRefDijet->setLeadJetEtaLab( etaRefLeadLab );
            fRefDijet->setLeadJetEtaCM( etaRefLeadCM );
            fRefDijet->setLeadJetPhi( phiRefLead );

            fRefDijet->setSubLeadJetPt( ptRefSubLead );
            fRefDijet->setSubLeadJetEtaLab( etaRefSubLeadLab );
            fRefDijet->setSubLeadJetEtaCM( etaRefSubLeadCM );
            fRefDijet->setSubLeadJetPhi( phiRefSubLead );

            dijetRefPtAve = fRefDijet->ptAve();
            dijetRefDphi = fRefDijet->dPhi();
            dijetRefEtaLab = fRefDijet->etaLab();
            dijetRefEtaCM = fRefDijet->etaCM();
            dijetRefPhi = fRefDijet->phi();


            if ( fabs( ptRefLead - fRefDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
                fabs( ptRefSubLead - fRefDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefLeadLab - fRefDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefSubLeadLab - fRefDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefLeadCM - fRefDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefSubLeadCM - fRefDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
                fabs( phiRefLead - fRefDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
                fabs( phiRefSubLead - fRefDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
                std::cerr << "Error: Ref dijet reference parameters do not match the leading and subleading jet parameters within rounding\n";
                return;
            }

            // if ( fVerbose ) {
            //     std::cout << "Matching gen dijet:\n"; 
            //     dijetRef.print();
            // }
        }
    } // if ( fIsMc )

    double recoDijetLeadSubLead[11] { dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDetaCM, dijetRecoDphi, dijetRecoPhi,
                                      ptRecoLead, etaRecoLeadLab, phiRecoLead,
                                      ptRecoSubLead, etaRecoSubLeadLab, phiRecoSubLead };

    double recoDijetLeadSubLeadCM[11] { dijetRecoPtAve, dijetRecoEtaCM, dijetRecoDetaCM, dijetRecoDphi, dijetRecoPhi,
                                        ptRecoLead, etaRecoLeadCM, phiRecoLead,
                                        ptRecoSubLead, etaRecoSubLeadCM, phiRecoSubLead };

    double reco2refDijetLeadSubLead[12] { dijetRecoPtAve, dijetRecoEtaLab, 
                                          ptRecoLead, etaRecoLeadLab, 
                                          ptRecoSubLead, etaRecoSubLeadLab, 
                                          dijetRefPtAve, dijetRefEtaLab,
                                          ptRefLead, etaRefLeadLab,
                                          ptRefSubLead, etaRefSubLeadLab };   

    double reco2refDijetLeadSubLeadCM[12] { dijetRecoPtAve, dijetRecoEtaCM, 
                                            ptRecoLead, etaRecoLeadCM, 
                                            ptRecoSubLead, etaRecoSubLeadCM, 
                                            dijetRefPtAve, dijetRefEtaCM,
                                            ptRefLead, etaRefLeadCM,
                                            ptRefSubLead, etaRefSubLeadCM };   
    //
    // Lab frame
    //

    bool fIsRecoDijetLabFound{false}; 
    if ( !fDiJetCut ) {
        fIsRecoDijetLabFound = true; // No cut, so dijet is always found
    } 
    else {
        // Check if the dijet passes the cut (dijetReco, isCM)
        fIsRecoDijetLabFound = fDiJetCut->pass(fRecoDijet, false);
    }
    if ( fVerbose ) {
        std::cout << Form("Reco dijet in the lab frame is %s\n", ((fIsRecoDijetLabFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze reco dijets in lab frame
    if ( fIsRecoDijetLabFound ) {

        if ( fVerbose ) {
            std::cout << "Reco dijet parameters in the lab frame: " << std::endl;
            std::cout << Form("--> Lead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRecoDijet->leadJetPt(), fRecoDijet->leadJetEtaLab(), fRecoDijet->leadJetPhi()) << std::endl;
            std::cout << Form("--> Sublead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRecoDijet->subLeadJetPt(), fRecoDijet->subLeadJetEtaLab(), fRecoDijet->subLeadJetPhi()) << std::endl;
            std::cout << Form("--> Dijet ptAve: %5.2f dijet eta: %5.2f dijet dphi: %5.2f", 
                              dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDphi) << std::endl;
        }

        // Fill histograms for the reco dijet in the lab frame
        fHM->hRecoDijetLeadSubLead->Fill(recoDijetLeadSubLead, weight);

        // Lead and SubLead
        fHM->hRecoPtLeadPtSublead->Fill( ptRecoLead, ptRecoSubLead, weight );
        fHM->hRecoEtaLeadEtaSublead->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight );
        fHM->hRecoPtLeadPtSubleadMcReweight->Fill( ptRecoLead, ptRecoSubLead, weight * fMcReweight );
        fHM->hRecoEtaLeadEtaSubleadMcReweight->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight * fMcReweight );

        fHM->hRecoDijetEta->Fill( dijetRecoEtaLab, weight * fMcReweight);
        fHM->hRecoDijetPtEta->Fill( dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight);
        fHM->hRecoDijetPtEtaPhi->Fill( dijetRecoPtAve, dijetRecoEtaLab, dijetRecoPhi, 1. );
        fHM->hRecoDijetPtEtaPhiWeighted->Fill( dijetRecoPtAve, dijetRecoEtaLab, dijetRecoPhi, weight * fMcReweight);
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
        if ( fIsMc && recoLeadJet->hasMatching() && recoSubLeadJet->hasMatching() ) {

            if ( fVerbose ) {
                std::cout << "Ref dijet parameters in the lab frame: " << std::endl;
                std::cout << Form("--> Lead pT: %5.2f eta: %5.2f phi: %5.2f", 
                                fRefDijet->leadJetPt(), fRefDijet->leadJetEtaLab(), fRefDijet->leadJetPhi()) << std::endl;
                std::cout << Form("--> Sublead pT: %5.2f eta: %5.2f phi: %5.2f", 
                                fRefDijet->subLeadJetPt(), fRefDijet->subLeadJetEtaLab(), fRefDijet->subLeadJetPhi()) << std::endl;
                std::cout << Form("--> Dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f", 
                                dijetRefPtAve, dijetRefEtaLab, dijetRefDphi) << std::endl;
            }

            fHM->hReco2RefDijetLeadSubLead->Fill( reco2refDijetLeadSubLead, weight * fMcReweight );

            // Simple Reco to Ref correlation
            double dijetUnfold[4] = { dijetRecoPtAve, dijetRecoEtaLab, dijetRefPtAve, dijetRefEtaLab };
            fHM->hRecoDijetPtEtaRefDijetPtEta->Fill(dijetUnfold, weight * fMcReweight);

            // Lead and SubLead
            fHM->hRefPtLeadPtSublead->Fill( ptRefLead, ptRefSubLead, weight );
            fHM->hRefEtaLeadEtaSublead->Fill( etaRefLeadLab, etaRefSubLeadLab, weight );
            fHM->hRefPtLeadPtSubleadMcReweight->Fill( ptRefLead, ptRefSubLead, weight * fMcReweight );
            fHM->hRefEtaLeadEtaSubleadMcReweight->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

            fHM->hRefDijetEta->Fill( dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEtaLab, dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Fill( dijetRecoEtaLab, dijetRefEtaLab, dijetRecoPtAve, weight * fMcReweight);
            fHM->hRefDijetPtEtaPhi->Fill( dijetRefPtAve, dijetRefEtaLab, dijetRefPhi, 1. );
            fHM->hRefDijetPtEtaPhiWeighted->Fill( dijetRefPtAve, dijetRefEtaLab, dijetRefPhi, weight * fMcReweight );
            (dijetRefEtaLab >= 0) ? fHM->hRefDijetPtEtaForward->Fill(dijetRefPtAve, dijetRefEtaLab) : fHM->hRefDijetPtEtaBackward->Fill(dijetRefPtAve, TMath::Abs(dijetRefEtaLab));
            (dijetRefEtaLab >= 0) ? fHM->hRefDijetPtEtaForwardWeighted->Fill(dijetRefPtAve, dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetPtEtaBackwardWeighted->Fill(dijetRefPtAve, TMath::Abs(dijetRefEtaLab), weight * fMcReweight);

            // Find exact dijet pT bins
            int refPtAveBin = findDijetPtAveBin( dijetRefPtAve );
            int refPtAveOldBin = findDijetPtAveOldBin( dijetRefPtAve );

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
        } // if ( fIsMc && recoLeadJet->hasMatching() && recoSubLeadJet->hasMatching() )
    } // if ( fIsRecoDijetLabFound )

    //
    // CM frame
    // 

    bool fIsRecoDijetCMFound{false};
    if ( !fDiJetCut ) {
        fIsRecoDijetCMFound = true; // No cut, so dijet is always found
    } 
    else {
        // Check if the dijet passes the cut (dijetReco, isCM)
        fIsRecoDijetCMFound = fDiJetCut->pass(fRecoDijet, true);
    }
    if ( fVerbose ) {
        std::cout << Form("\nReco dijet in CM frame is %s\n", ((fIsRecoDijetCMFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze reco dijets in CM frame
    if ( fIsRecoDijetCMFound ) {

        if ( fVerbose ) {
            std::cout << "Reco dijet parameters in the C.M. frame: " << std::endl;
            std::cout << Form("--> Lead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRecoDijet->leadJetPt(), fRecoDijet->leadJetEtaCM(), fRecoDijet->leadJetPhi()) << std::endl;
            std::cout << Form("--> Sublead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRecoDijet->subLeadJetPt(), fRecoDijet->subLeadJetEtaCM(), fRecoDijet->subLeadJetPhi()) << std::endl;
            std::cout << Form("--> Dijet ptAve: %5.2f dijet eta (CM): %5.2f dijet dphi: %5.2f", 
                              dijetRecoPtAve, dijetRecoEtaCM, dijetRecoDphi) << std::endl;
        }

        fHM->hRecoDijetLeadSubLeadCM->Fill(recoDijetLeadSubLeadCM, weight);

        fHM->hRecoPtLeadPtSubleadCM->Fill( ptRecoLead, ptRecoSubLead, weight );
        fHM->hRecoEtaCMLeadEtaCMSublead->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight );
        fHM->hRecoPtLeadPtSubleadCMMcReweight->Fill( ptRecoLead, ptRecoSubLead, weight * fMcReweight );
        fHM->hRecoEtaCMLeadEtaCMSubleadMcReweight->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight * fMcReweight );

        fHM->hRecoDijetEtaCM->Fill( dijetRecoEtaCM, weight * fMcReweight);
        fHM->hRecoDijetPtEtaCM->Fill( dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight);
        fHM->hRecoDijetPtEtaPhiCM->Fill( dijetRecoPtAve, dijetRecoEtaCM, dijetRecoPhi, 1. );
        fHM->hRecoDijetPtEtaPhiCMWeighted->Fill( dijetRecoPtAve, dijetRecoEtaCM, dijetRecoPhi, weight * fMcReweight);
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
        if ( fIsMc && recoLeadJet->hasMatching() && recoSubLeadJet->hasMatching() ) {

            if ( fVerbose ) {
                std::cout << "Ref dijet parameters in the C.M. frame: " << std::endl;
                std::cout << Form("--> Lead pT: %5.2f eta: %5.2f phi: %5.2f", 
                                  fRefDijet->leadJetPt(), fRefDijet->leadJetEtaCM(), fRefDijet->leadJetPhi()) << std::endl;
                std::cout << Form("--> Sublead pT: %5.2f eta: %5.2f phi: %5.2f", 
                                  fRefDijet->subLeadJetPt(), fRefDijet->subLeadJetEtaCM(), fRefDijet->subLeadJetPhi()) << std::endl;
                std::cout << Form("--> Dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f", 
                                  dijetRefPtAve, dijetRefEtaCM, dijetRefDphi) << std::endl;
            }

            fHM->hReco2RefDijetLeadSubLeadCM->Fill( reco2refDijetLeadSubLeadCM, weight * fMcReweight );
            // Simple Reco to Ref correlation
            double dijetUnfold[4] = { dijetRecoPtAve, dijetRecoEtaCM, dijetRefPtAve, dijetRefEtaCM };
            fHM->hRecoDijetPtEtaRefDijetPtEtaCM->Fill(dijetUnfold, weight * fMcReweight);

            // Lead and SubLead
            fHM->hRefPtLeadPtSubleadCM->Fill( ptRefLead, ptRefSubLead, weight );
            fHM->hRefEtaCMLeadEtaCMSublead->Fill( etaRefLeadCM, etaRefSubLeadCM, weight );
            fHM->hRefPtLeadPtSubleadCMMcReweight->Fill( ptRefLead, ptRefSubLead, weight * fMcReweight );
            fHM->hRefEtaCMLeadEtaCMSubleadMcReweight->Fill( etaRefLeadCM, etaRefSubLeadCM, weight * fMcReweight );
    

            fHM->hRefDijetEtaCM->Fill( dijetRefEtaCM, weight );
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPtAve, weight);
            fHM->hRefDijetPtEtaPhiCM->Fill( dijetRefPtAve, dijetRefEtaCM, dijetRefPhi, 1. );
            fHM->hRefDijetPtEtaPhiCMWeighted->Fill( dijetRefPtAve, dijetRefEtaCM, dijetRefPhi, weight * fMcReweight );

            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForward->Fill(dijetRefPtAve, dijetRefEtaCM, 1.) : fHM->hRefDijetPtEtaCMBackward->Fill(dijetRefPtAve, TMath::Abs(dijetRefEtaCM), 1.);
            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForwardWeighted->Fill(dijetRefPtAve, dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetPtEtaCMBackwardWeighted->Fill(dijetRefPtAve, TMath::Abs(dijetRefEtaCM), weight * fMcReweight);

            // Find exact dijet ptAve bin
            int refPtAveBin = findDijetPtAveBin( dijetRefPtAve );
            int refPtAveOldBin = findDijetPtAveOldBin( dijetRefPtAve );

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

    fRecoDijet->cleanParameters();
    fRefDijet->cleanParameters();

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
    float etaRefLeadLab = etaLab( refLeadJet->eta() );
    float etaRefLeadCM = boostEta2CM( refLeadJet->eta() );
    float phiRefLead = refLeadJet->phi();

    float ptRefSubLead = refSubLeadJet->pt();
    float etaRefSubLeadLab = etaLab( refSubLeadJet->eta() );
    float etaRefSubLeadCM = boostEta2CM( refSubLeadJet->eta() );
    float phiRefSubLead = refSubLeadJet->phi();

    // if (fVerbose) {
    //     std::cout << "Ref-selected inclusive leading and subleading jets:" << std::endl;
    //     std::cout << Form("Ref leading jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRefLead, etaRefLeadLab, etaRefLeadCM, phiRefLead);
    //     std::cout << Form("Ref subleading jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRefSubLead, etaRefSubLeadLab, etaRefSubLeadCM, phiRefSubLead);
    // }

    fRefDijet->cleanParameters();
    
    fRefDijet->setLeadJetPt( ptRefLead );
    fRefDijet->setLeadJetEtaLab( etaRefLeadLab );
    fRefDijet->setLeadJetEtaCM( etaRefLeadCM );
    fRefDijet->setLeadJetPhi( phiRefLead );

    fRefDijet->setSubLeadJetPt( ptRefSubLead );
    fRefDijet->setSubLeadJetEtaLab( etaRefSubLeadLab );
    fRefDijet->setSubLeadJetEtaCM( etaRefSubLeadCM );
    fRefDijet->setSubLeadJetPhi( phiRefSubLead );

    float dijetRefPtAve = fRefDijet->ptAve();
    float dijetRefEtaLab = fRefDijet->etaLab();
    float dijetRefEtaCM = fRefDijet->etaCM();
    float dijetRefDphi = fRefDijet->dPhi();
    float dijetRefPhi = fRefDijet->phi();

    if ( fabs(ptRefLead - fRefDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(ptRefSubLead - fRefDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefLeadLab - fRefDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefSubLeadLab - fRefDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefLeadCM - fRefDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefSubLeadCM - fRefDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRefLead - fRefDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRefSubLead - fRefDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
        std::cout << "Error: Ref dijet reference parameters do not match the leading and subleading jet parameters within rounding\n";
        return;
    }

    // if (fVerbose) {
    //     std::cout << "Ref-selected inclusive dijet parameters:" << std::endl;
    //     dijetRef.print();
    // }

    // Reconstructed jet partners
    float ptRecoLead = recoLeadJet->ptJECCorr();
    float etaRecoLeadLab = etaLab( recoLeadJet->eta() );
    float etaRecoLeadCM = boostEta2CM( recoLeadJet->eta() );
    float phiRecoLead = recoLeadJet->phi();
        
    float ptRecoSubLead = recoSubLeadJet->ptJECCorr();
    float etaRecoSubLeadLab = etaLab( recoSubLeadJet->eta() );
    float etaRecoSubLeadCM = boostEta2CM( recoSubLeadJet->eta() );
    float phiRecoSubLead = recoSubLeadJet->phi();

    // if (fVerbose) {
    //     std::cout << "Ref-selected inclusive reconstructed leading and subleading jets:" << std::endl;
    //     std::cout << Form("Reco leading jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRecoLead, etaRecoLeadLab, etaRecoLeadCM, phiRecoLead);
    //     std::cout << Form("Reco subleading jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRecoSubLead, etaRecoSubLeadLab, etaRecoSubLeadCM, phiRecoSubLead);
    // }

    fRecoDijet->cleanParameters();

    fRecoDijet->setLeadJetPt( ptRecoLead );
    fRecoDijet->setLeadJetEtaLab( etaRecoLeadLab );
    fRecoDijet->setLeadJetEtaCM( etaRecoLeadCM );
    fRecoDijet->setLeadJetPhi( phiRecoLead );

    fRecoDijet->setSubLeadJetPt( ptRecoSubLead );
    fRecoDijet->setSubLeadJetEtaLab( etaRecoSubLeadLab );
    fRecoDijet->setSubLeadJetEtaCM( etaRecoSubLeadCM );
    fRecoDijet->setSubLeadJetPhi( phiRecoSubLead );

    float dijetRecoPtAve = fRecoDijet->ptAve();
    float dijetRecoEtaLab = fRecoDijet->etaLab();
    float dijetRecoEtaCM = fRecoDijet->etaCM();
    float dijetRecoDphi = fRecoDijet->dPhi();
    // float dijetRecoPhi = fRecoDijet->phi();

    if ( fabs(ptRecoLead - fRecoDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(ptRecoSubLead - fRecoDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoLeadLab - fRecoDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoSubLeadLab - fRecoDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoLeadCM - fRecoDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoSubLeadCM - fRecoDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRecoLead - fRecoDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRecoSubLead - fRecoDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
        std::cout << "Error: Reco dijet reference parameters do not match the leading and subleading jet parameters within rounding\n";
        return;
    }

    // if (fVerbose) {
    //     std::cout << "Ref-selected inclusive reconstructed dijet parameters:" << std::endl;
    //     dijetReco.print();
    // }


    //
    // Lab frame
    //

    bool fIsRefSelDijetLabFound{false}; 
    if ( !fDiJetCut ) {
        fIsRefSelDijetLabFound = true; // No cut, so dijet is always found
    } 
    else {
        // Check if the dijet passes the cut (dijetRef, isLab)
        fIsRefSelDijetLabFound = fDiJetCut->pass(fRefDijet, false);
    }
    if ( fVerbose ) {
        std::cout << Form("Ref dijet in lab frame is %s\n", ((fIsRefSelDijetLabFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze ref-selected dijets in lab frame
    if ( fIsRefSelDijetLabFound ) {

        if ( fVerbose ) {
            std::cout << "Ref dijet parameters in the lab frame: " << std::endl;
            std::cout << Form("Ref lead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco lead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRefDijet->leadJetPt(), fRefDijet->leadJetEtaLab(), fRefDijet->leadJetPhi(), 
                              fRecoDijet->leadJetPt(), fRecoDijet->leadJetEtaLab(), fRecoDijet->leadJetPhi()) << std::endl;
            std::cout << Form("Ref sublead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco sublead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRefDijet->subLeadJetPt(), fRefDijet->subLeadJetEtaLab(), fRefDijet->subLeadJetPhi(), 
                              fRecoDijet->subLeadJetPt(), fRecoDijet->subLeadJetEtaLab(), fRecoDijet->subLeadJetPhi()) << std::endl;
            std::cout << Form("Ref dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f --> Reco dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f\n", 
                              dijetRefPtAve, dijetRefEtaLab, dijetRefDphi, 
                              dijetRecoPtAve, dijetRecoEtaLab, dijetRecoDphi) << std::endl;
        }

        fHM->hRefSelDijetEta->Fill(dijetRefEtaLab, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaPhi->Fill(dijetRefPtAve, dijetRefEtaLab, dijetRefPhi, 1.);
        fHM->hRefSelDijetPtEtaPhiWeighted->Fill(dijetRefPtAve, dijetRefEtaLab, dijetRefPhi, weight * fMcReweight );

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

    bool fIsRefSelDijetCMFound{false};
    if ( !fDiJetCut ) {
        fIsRefSelDijetCMFound = true; // No cut, so dijet is always found
    } 
    else {
        // Check if the dijet passes the cut (dijetRef, isCM)
        fIsRefSelDijetCMFound = fDiJetCut->pass(fRefDijet, true);
    }
    if ( fVerbose ) {
        std::cout << Form("Ref dijet in CM frame is %s\n", ((fIsRefSelDijetCMFound) ? "[good]" : "[bad]") ); 
    }

    // Analyze ref-selected dijets in CM frame
    if ( fIsRefSelDijetCMFound ) {

        if ( fVerbose ) {
            std::cout << "Ref dijet parameters in the C.M. frame: " << std::endl;
            std::cout << Form("Ref lead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco lead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRefDijet->leadJetPt(), fRefDijet->leadJetEtaCM(), fRefDijet->leadJetPhi(),
                              fRecoDijet->leadJetPt(), fRecoDijet->leadJetEtaCM(), fRecoDijet->leadJetPhi()) << std::endl;
            std::cout << Form("Ref sublead pt: %5.2f eta: %5.2f phi: %5.2f --> Reco sublead pt: %5.2f eta: %5.2f phi: %5.2f", 
                              fRefDijet->subLeadJetPt(), fRefDijet->subLeadJetEtaCM(), fRefDijet->subLeadJetPhi(),
                              fRecoDijet->subLeadJetPt(), fRecoDijet->subLeadJetEtaCM(), fRecoDijet->subLeadJetPhi()) << std::endl;
            std::cout << Form("Ref dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f --> Reco dijet ptAve: %5.2f eta: %5.2f dphi: %5.2f", 
                              dijetRefPtAve, dijetRefEtaCM, dijetRefDphi, 
                              dijetRecoPtAve, dijetRecoEtaCM, dijetRecoDphi) << std::endl;
        }

        fHM->hRefSelDijetEtaCM->Fill(dijetRefEtaCM, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaPhiCM->Fill(dijetRefPtAve, dijetRefEtaCM, dijetRefPhi, 1.);
        fHM->hRefSelDijetPtEtaPhiCMWeighted->Fill(dijetRefPtAve, dijetRefEtaCM, dijetRefPhi, weight * fMcReweight );

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

    fRefDijet->cleanParameters();
    fRecoDijet->cleanParameters();

    if ( fVerbose ) {
        std::cout << "DiJetAnalysis::processRefDijets -- end" << std::endl;
    }
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
        return;
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
    // float centrality = event->centrality();
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

    // Fill event histograms
    fHM->hVz->Fill( vz,  1. );
    fHM->hVzWeighted->Fill( vz, weight );

    fHM->hPtHat->Fill( ptHat, 1. );
    fHM->hPtHatWeighted->Fill( ptHat, weight );

    fHM->hHiBin->Fill( event->hiBin(), 1. );
    fHM->hHiBinWeighted->Fill( event->hiBin(), weight );

    // Process and analyze dijets
    processDijets(event, weight);

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
        std::cout << "\nDiJetAnalysis::processEvent -- end" << std::endl;
    }
}

//________________
void DiJetAnalysis::finish() {
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
