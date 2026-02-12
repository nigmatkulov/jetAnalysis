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
    fMixBufferSize{10}, fMixBufferOlga{}, fMixBuffer{} {

    fPtHatRange[0] = {0};
    fPtHatRange[1] = {100000000};
    for (int i=0; i<fJetPtBins; i++) {
        for (int j=0; j<fJetPtBins; j++) {
            fJetPtLeadPtSubleadReweightMatrix[i][j] = 1;
        }
    } // for (int i=0; i<fJetPtBins; i++)

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

    // pT Lead, pT SubLead weighting matrix
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
void DiJetAnalysis::addEventToMixBufferOlga(const double &ptAve, const DiJet& dijet) {
    // float binLow = 50.f;
    // float binHigh = 500.f;
    // float step = 50.f;
    // int bin = int( (ptAve - binLow) / step ) ;
    // int bin = findDijetPtAveBin(ptAve);
    // if ( bin < 0 || bin >= fMixBuffer.size() ) {
    // }
    // else {
    //     fMixBuffer[bin].push_back(dijet);
    // }
}

//________________
void DiJetAnalysis::addEventToMixBuffer(const double &vz, const DiJet& dijet) {

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
        // Indices of Lead and SubLead reco jets
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
            // Indices of Lead and SubLead gen jets
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
            // Indices of Lead and SubLead ref-selected reco jets
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

    // Loop over reconstructed jets
    for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

        float pt = (*recoJetIter)->ptJECCorr();
        float eta = etaLab( (*recoJetIter)->eta() );
        float etaCM = boostEta2CM( (*recoJetIter)->eta() );
        float phi = (*recoJetIter)->phi();
        float ptRaw = (*recoJetIter)->rawPt();

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
        if ( fabs( eta ) <= 2.4 ) { dummyIter = {0}; }
        else if ( fabs( eta ) <= 2.7 ) { dummyIter = {1}; }
        else if ( fabs( eta ) <= 3.0 ) { dummyIter = {2}; }
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

        // Check if jet passes the selection criteria (*recoJet, isCM, isMC, requireMatching)
        if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, false, false) ) continue;

        unsigned int runId = event->runId();

        if ( 40. < pt && pt < 90. ) {
            int idx = -1;
            if ( runId == 285480 ) { idx = 1; }
            else if ( runId == 285505 ) { idx = 2; }
            else if ( runId == 285517 ) { idx = 3; }
            else if ( runId == 285832 ) { idx = 4; }
            else if ( runId == 285993 ) { idx = 5; }

            // For all runs
            fHM->hRecoSubLeadJetEtaRun[0]->Fill( eta, weight );
            if ( idx > 0 ) fHM->hRecoSubLeadJetEtaRun[idx]->Fill( eta, weight );
            // Lead jet
            if ( (fRecoIdLead >= 0) && ((recoJetCounter-1) == fRecoIdLead) ) {
                fHM->hRecoLeadJetEtaRun[0]->Fill( eta, weight );
                if ( idx > 0 ) fHM->hRecoLeadJetEtaRun[idx]->Fill( eta, weight );
            }
            if ( (fRecoIdSubLead >= 0) && ((recoJetCounter-1) == fRecoIdSubLead) ) {
                fHM->hRecoSubLeadJetEtaRun[0]->Fill( eta, weight );
                if ( idx > 0 ) fHM->hRecoSubLeadJetEtaRun[idx]->Fill( eta, weight );
            }
        } // if ( 40. < pt && pt < 90. )

        //Inclusive (matched+unmatched) jets
        fHM->hRecoInclusiveAllJetPt->Fill(pt, weight);
        fHM->hRecoInclusiveAllJetEta->Fill(eta, weight);
        fHM->hRecoInclusiveAllJetEtaUnweighted->Fill(eta, 1.);
        fHM->hRecoInclusiveAllJetPtEta->Fill(eta, pt, weight);
        fHM->hRecoInclusiveAllJetPtEtaCM->Fill(etaCM, pt, weight );
        fHM->hRecoInclusiveAllJetPtRawEtaStdBins->Fill( eta, ptRaw, weight);
        fHM->hRecoInclusiveAllJetPtEtaStdBins->Fill( eta, pt, weight);
        fHM->hRecoInclusiveAllJetPtRawEta->Fill(eta, ptRaw, weight);
        fHM->hRecoInclusiveAllJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
        fHM->hRecoInclusiveAllJetPtRawEtaPtHat->Fill(eta, ptRaw, ptHat, weight);

        // Inclusive (matched+unmatched) Lead jet
        if ( (fRecoIdLead >= 0) && ((recoJetCounter - 1) == fRecoIdLead) ) {
            fHM->hRecoLeadAllJetPtEta->Fill( eta, pt, weight );
            fHM->hRecoLeadAllJetPtEtaCM->Fill( etaCM, pt, weight );
            fHM->hRecoLeadAllJetPtRawEtaStdBins->Fill( eta, ptRaw, weight);
            fHM->hRecoLeadAllJetPtEtaStdBins->Fill( eta, pt, weight);
            fHM->hRecoLeadAllJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );

        }

        // Inclusive (matched+unmatched) SubLead jet
        if ( (fRecoIdSubLead >= 0) && ((recoJetCounter - 1) == fRecoIdSubLead) ) {
            fHM->hRecoSubLeadAllJetPtEta->Fill( eta, pt, weight );
            fHM->hRecoSubLeadAllJetPtEtaCM->Fill( etaCM, pt, weight );
            fHM->hRecoSubLeadAllJetPtRawEtaStdBins->Fill( eta, ptRaw, weight);
            fHM->hRecoSubLeadAllJetPtEtaStdBins->Fill( eta, pt, weight);
            fHM->hRecoSubLeadAllJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
        }

        if ( pt > 30. ) {
            fHM->hRecoGoodInclusiveJetEtaLabFrame->Fill( eta, weight );
            fHM->hRecoGoodInclusiveJetEtaCMFrame->Fill( etaCM, weight );
        }

        // For the MC, check and study matching to gen reconstructed jets
        if ( fIsMc ) {

            // If reco jet has matching to gen jet
            if ( (*recoJetIter)->hasMatching() ) {

                GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
                if ( !matchedJet ) {
                    std::cerr << Form("Cannot retrieve gen jet with id: %d", (*recoJetIter)->genJetId() ) << std::endl;
                    continue;
                }
                
                float genPt = matchedJet->pt();
                float genEta = etaLab( matchedJet->eta() );
                float genEtaCM = boostEta2CM( matchedJet->eta() );
                float genPhi = matchedJet->phi();

                float JES = pt/genPt;

                double res[4]  = { JES, genPt, genEta, genPhi };
                double res1[4] = { JES, genPt, genEta, ptHat };
                double res2[4] = { JES, pt, eta, ptHat };
                double res3[4] = { JES, genPt, genEtaCM, ptHat };

                // Inclusive ref jets
                fHM->hRefInclusiveJetPt->Fill( genPt, weight );
                fHM->hRefInclusiveJetEta->Fill( genEta, weight );
                fHM->hRefInclusiveJetEtaUnweighted->Fill( genEta, 1. );
                fHM->hRefInclusiveJetPtEta->Fill( genEta, genPt, weight );
                fHM->hRefInclusiveJetPtEtaCM->Fill( genEtaCM, genPt, weight );
                // Must use reco jet pt to see the effect of JEC
                fHM->hRefInclusiveJetPtEtaPtHat->Fill( genEta, pt, ptHat, weight );

                // Jet energy scale
                fHM->hJESInclusiveJetPtEtaPhi->Fill(res, weight);

                // Inclusive matched reco jets
                fHM->hRecoInclusiveMatchedJetPt->Fill( pt, weight );
                fHM->hRecoInclusiveMatchedJetPtEta->Fill( eta, pt, weight );
                fHM->hRecoInclusiveMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );

                // fHM->hRecoInclusiveJetJECFactorVsPtEta->Fill( pt / ptRaw, genPt, genEta, weight );
                // fHM->hRecoInclusiveJetPtRawOverPtRefVsPtEta->Fill( ptRaw/genPt, genPt, genEta, weight );
                // fHM->hRecoInclusiveJetPtRawOverPtRefVsPtEtaStdBinning->Fill( ptRaw/genPt, genPt, genEta, weight );
                // fHM->hRecoInclusiveJetPtRawOverPtRefVsRecoPtEtaStdBinning->Fill( ptRaw/genPt, ptRaw, eta, weight );

                // Fill JES vs pt for |eta| < 1.4 (midrapidity)
                if ( fabs( genEta ) < 1.4 ) {
                    fHM->hInclusiveJetJESVsPtGen->Fill( genPt, JES, weight );
                }
                fHM->hInclusiveJetJESGenPtGenEtaPtHatWeighted->Fill( res1, weight );
                fHM->hInclusiveJetJESGenPtGenEtaCMPtHatWeighted->Fill( res3, weight );
                fHM->hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Fill( res2, weight );

                double correl[5] { pt, ptRaw, genPt, eta, genEta };
                fHM->hRecoInclusiveJetReco2Ref->Fill( correl, weight );

                // Lead jet
                if ( (fRecoIdLead >= 0) && ((recoJetCounter-1) == fRecoIdLead) ) {
                    fHM->hLeadJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                    fHM->hLeadJetJESGenPtEtaCMPtHatWeighted->Fill( res3, weight );
                    fHM->hRecoLeadJetReco2Ref->Fill( correl, weight );
                    // Lead matched reco jet
                    fHM->hRecoLeadMatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoLeadMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                    // Lead ref jet
                    fHM->hRefLeadJetPtEta->Fill( genEta, genPt, weight );
                    fHM->hRefLeadJetPtEtaCM->Fill( genEtaCM, genPt, weight );
                    fHM->hRefLeadJetPtEtaPtHat->Fill( genEta, pt, ptHat, weight );

                    if ( (*recoJetIter)->genJetId() == fGenIdLead ) {
                        fHM->hRefLeadUnswappedJetPtEta->Fill( genEta, genPt, weight );
                        fHM->hRefLeadUnswappedJetPtEtaPtHat->Fill( genEta, pt, ptHat, weight );
                    }

                }

                // SubLead jet
                if ( (fRecoIdSubLead >= 0) && ((recoJetCounter-1) == fRecoIdSubLead) ) {
                    fHM->hSubLeadJetJESGenPtEtaPtHatWeighted->Fill( res1, weight );
                    fHM->hSubLeadJetJESGenPtEtaCMPtHatWeighted->Fill( res3, weight );
                    fHM->hRecoSubLeadJetReco2Ref->Fill( correl, weight );
                    // SubLead matched reco jet
                    fHM->hRecoSubLeadMatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoSubLeadMatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                    // SubLead ref jet
                    fHM->hRefSubLeadJetPtEta->Fill( genEta, genPt, weight );
                    fHM->hRefSubLeadJetPtEtaCM->Fill( genEtaCM, genPt, weight );
                    fHM->hRefSubLeadJetPtEtaPtHat->Fill( genEta, pt, ptHat, weight );

                    if ( (*recoJetIter)->genJetId() == fGenIdSubLead ) {
                        fHM->hRefSubLeadUnswappedJetPtEta->Fill( genEta, genPt, weight );
                        fHM->hRefSubLeadUnswappedJetPtEtaPtHat->Fill( genEta, pt, ptHat, weight );
                    }
                }

            } // if ( (*recoJetIter)->hasMatching() )
            else {
                // Fill unmatched reco jets
                fHM->hRecoInclusiveUnmatchedJetPtEta->Fill(eta, pt, weight);
                fHM->hRecoInclusiveUnmatchedJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

                // Lead unmatched reco jet
                if ( (fRecoIdLead >= 0) && ((recoJetCounter-1) == fRecoIdLead) ) {
                    fHM->hRecoLeadUnmatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoLeadUnmatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );
                }

                // SubLead unmatched reco jet
                if ( (fRecoIdSubLead >= 0) && ((recoJetCounter-1) == fRecoIdSubLead) ) {
                    fHM->hRecoSubLeadUnmatchedJetPtEta->Fill( eta, pt, weight);
                    fHM->hRecoSubLeadUnmatchedJetPtEtaPtHat->Fill( eta, pt, ptHat, weight );

                }
            } // else
        } // if ( fIsMc )
    } // for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ )

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

        // Loop over generated jets and search for Lead and SubLead jets
        for ( genJetIter = event->genJetCollection()->begin(); genJetIter != event->genJetCollection()->end(); genJetIter++ ) {

            float pt = (*genJetIter)->pt();
            float eta = etaLab( (*genJetIter)->eta() );
            float etaCM = boostEta2CM( (*genJetIter)->eta() );
            float phi = (*genJetIter)->phi();
    
            if ( fVerbose ) {
                std::cout << Form("Gen jet #%d: pt: %5.2f eta: %5.2f phi: %5.2f\n", genJetCounter, pt, eta, phi);
            }

            genJetCounter++;

            // Selection criteria (*genJet, isCM)
            if ( fGenJetCut && !fGenJetCut->pass(*genJetIter, false) ) continue;

            // Inclusive gen jet
            fHM->hGenInclusiveJetPt->Fill(pt, weight);
            fHM->hGenInclusiveJetEta->Fill(eta, weight);
            fHM->hGenInclusiveJetEtaUnweighted->Fill(eta, 1.);
            fHM->hGenInclusiveJetPtEta->Fill(eta, pt, weight);
            fHM->hGenInclusiveJetPtEtaCM->Fill( etaCM, pt, weight);
            fHM->hGenInclusiveJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);

            // Lead gen jet
            if ( (fGenIdLead >= 0) && ((genJetCounter - 1) ==  fGenIdLead) ) {
                fHM->hGenLeadJetPtEta->Fill(eta, pt, weight);
                fHM->hGenLeadJetPtEtaCM->Fill( etaCM, pt, weight);
                fHM->hGenLeadJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
            }

            // SubLead gen jet
            if ( (fGenIdSubLead >= 0) && ((genJetCounter - 1) ==  fGenIdSubLead) ) {
                fHM->hGenSubLeadJetPtEta->Fill(eta, pt, weight);
                fHM->hGenSubLeadJetPtEtaCM->Fill( etaCM, pt, weight);
                fHM->hGenSubLeadJetPtEtaPtHat->Fill(eta, pt, ptHat, weight);
            }
    
            if ( pt > 30. ) {
                fHM->hGenGoodInclusiveJetEtaLabFrame->Fill( eta, weight);
                fHM->hGenGoodInclusiveJetEtaCMFrame->Fill( etaCM, weight );
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

        // Loop over reconstructed jets and search for Lead and SubLead gen-matched jets
        for ( recoJetIter = event->recoJetCollection()->begin(); recoJetIter != event->recoJetCollection()->end(); recoJetIter++ ) {

            refSelJetCounter++;

            // Check selection criteria (*recoJet, isCM, isMC, requireMatching)
            if ( fRecoJetCut && !fRecoJetCut->pass(*recoJetIter, false, true, true) ) continue; 

            // Retrieve matched gen jet
            GenJet *matchedJet = event->genJetCollection()->at( (*recoJetIter)->genJetId() );
            float genPt = matchedJet->pt();
            float genEta = etaLab( matchedJet->eta() );
            float genEtaCM = boostEta2CM( matchedJet->eta() );
            float genPhi = matchedJet->phi();

            if ( fVerbose ) {
                std::cout << Form("Ref jet #%d pt: %5.2f eta: %5.2f phi: %5.2f --> Reco jet #%d pt: %5.2f eta: %5.2f phi: %5.2f\n", 
                                  (*recoJetIter)->genJetId(), genPt, genEta, genPhi, refSelJetCounter-1, (*recoJetIter)->ptJECCorr(), (*recoJetIter)->eta(), (*recoJetIter)->phi());
            }

            // Inclusive ref-selected ref jet
            fHM->hRefSelInclusiveJetPt->Fill( genPt, weight * fMcReweight );
            fHM->hRefSelInclusiveJetEta->Fill( genEta, weight * fMcReweight );
            fHM->hRefSelInclusiveJetEtaUnweighted->Fill( genEta, 1. );
            fHM->hRefSelInclusiveJetPtEta->Fill(genEta, genPt, weight * fMcReweight);
            fHM->hRefSelInclusiveJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight * fMcReweight);

            // Lead ref-selected ref jet
            if ( (fRefSelRecoIdLead >= 0) && ((refSelJetCounter - 1) == fRefSelRecoIdLead) ) {
                fHM->hRefSelLeadJetPtEta->Fill(genEta, genPt, weight * fMcReweight);
                fHM->hRefSelLeadJetPtEtaPtHat->Fill(genEta, genPt, ptHat, weight * fMcReweight);
            }

            // SubLead ref-selected ref jet
            if ( (fRefSelRecoIdSubLead >= 0) && ((refSelJetCounter - 1) == fRefSelRecoIdSubLead) ) {
                fHM->hRefSelSubLeadJetPtEta->Fill(genEta, genPt, weight * fMcReweight);
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

        fHM->hRecoLeadJetPtOverPtHatVsLeadJetPt->Fill( ptRecoLead, ptRecoLead/ptHat, 1. );
        fHM->hRecoLeadJetPtOverPtHatVsLeadJetPtWeighted->Fill( ptRecoLead, ptRecoLead/ptHat, weight );
        
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

        fHM->hGenLeadJetPtOverPtHatVsLeadJetPt->Fill( ptGenLead, ptGenLead/ptHat, 1. );
        fHM->hGenLeadJetPtOverPtHatVsLeadJetPtWeighted->Fill( ptGenLead, ptGenLead/ptHat, weight );

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

    // if ( ptGenSubLead > 30.f && fabs(dijetGenDphi) > TMath::PiOver2() ) {
    //     double genDijetLeadSublead[10] { dijetGenPtAve, dijetGenEtaLab, dijetGenEtaCM, dijetGenPhi, 
    //                                      ptGenLead, etaGenLeadLab, etaGenLeadCM, 
    //                                      ptGenSubLead, etaGenSubLeadLab, etaGenSubLeadCM };

    //     fHM->hGenDijetInfo->Fill(genDijetLeadSublead, weight * fMcReweight);
    // }

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

    // Fill forward and backward distributions for gen dijets in the CM frame
    if ( ptGenLead > 50. && ptGenSubLead > 40. && fabs(dijetGenDphi) > (TMath::TwoPi() / 3) ) {
        for (int iEtaCut{0}; iEtaCut<6; iEtaCut++) {
            const float etaCut = 1.4 + iEtaCut * 0.1;
            if ( fabs(etaGenLeadCM) < etaCut && fabs(etaGenSubLeadCM) < etaCut ) {
                (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaForwardArr[iEtaCut]->Fill(dijetGenPtAve, dijetGenEtaCM, weight * fMcReweight) :
                                       fHM->hGenDijetPtEtaBackwardArr[iEtaCut]->Fill(dijetGenPtAve, fabs(dijetGenEtaCM), weight * fMcReweight);
            }
        }
    } // if ( ptGenLead > 50. && ptGenSubLead > 40. && fabs(dijetGenDphi) > (TMath::TwoPi() / 3) )

    //
    // Lab frame
    //

    fIsGenDijetLabFound = {false};
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

        fHM->hGenDijetEta->Fill(dijetGenEtaLab, weight * fMcReweight );
        fHM->hGenDijetPtEta->Fill(dijetGenPtAve, dijetGenEtaLab, 1.);
        fHM->hGenDijetPtEtaWeighted->Fill(dijetGenPtAve, dijetGenEtaLab, weight * fMcReweight );
        fHM->hGenDijetPtEtaCMInLab->Fill(dijetGenPtAve, dijetGenEtaCM, weight * fMcReweight);
        (dijetGenEtaLab >= 0) ? fHM->hGenDijetPtEtaForward->Fill(dijetGenPtAve, dijetGenEtaLab) : fHM->hGenDijetPtEtaBackward->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaLab));
        (dijetGenEtaLab >= 0) ? fHM->hGenDijetPtEtaForwardWeighted->Fill(dijetGenPtAve, dijetGenEtaLab, weight * fMcReweight) : fHM->hGenDijetPtEtaBackwardWeighted->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaLab), weight * fMcReweight);
        (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaForwardCMInLab->Fill(dijetGenPtAve, dijetGenEtaCM, weight * fMcReweight) : fHM->hGenDijetPtEtaBackwardCMInLab->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaCM), weight * fMcReweight);

        fHM->hGenDijetPtAveLeadPtSubLeadPt->Fill( dijetGenPtAve, ptGenLead, ptGenSubLead, weight );
        fHM->hGenDijetPtAveLeadEtaSubLeadEta->Fill( dijetGenPtAve, etaGenLeadLab, etaGenSubLeadLab, weight );
        fHM->hGenDijetEtaLeadEtaSubLeadEta->Fill( dijetGenEtaLab,  etaGenLeadLab, etaGenSubLeadLab, weight );
    } // if ( goodDijetLab )

    //
    // CM frame
    //

    fIsGenDijetCMFound = {false};
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

        fHM->hGenEtaCMLeadEtaCMSublead->Fill( etaGenLeadCM, etaGenSubLeadCM, weight );

        fHM->hGenDijetEtaCM->Fill(dijetGenEtaCM, weight * fMcReweight );
        fHM->hGenDijetPtEtaCM->Fill(dijetGenPtAve, dijetGenEtaCM, 1.);
        fHM->hGenDijetPtEtaCMWeighted->Fill(dijetGenPtAve, dijetGenEtaCM, weight * fMcReweight );
        fHM->hGenDijetPtEtaLabInCM->Fill(dijetGenPtAve, dijetGenEtaLab, weight * fMcReweight);
        (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForward->Fill(dijetGenPtAve, dijetGenEtaCM) : fHM->hGenDijetPtEtaCMBackward->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaCM));
        (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaCMForwardWeighted->Fill(dijetGenPtAve, dijetGenEtaCM, weight * fMcReweight) : fHM->hGenDijetPtEtaCMBackwardWeighted->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaCM), weight * fMcReweight);
        (dijetGenEtaCM >= 0) ? fHM->hGenDijetPtEtaForwardLabInCM->Fill(dijetGenPtAve, dijetGenEtaLab, weight * fMcReweight) : fHM->hGenDijetPtEtaBackwardLabInCM->Fill(dijetGenPtAve, TMath::Abs(dijetGenEtaLab), weight * fMcReweight);

        fHM->hGenDijetPtAveLeadPtSubLeadPtCM->Fill( dijetGenPtAve, ptGenLead, ptGenSubLead, weight );
        fHM->hGenDijetPtAveLeadEtaSubLeadEtaCM->Fill( dijetGenPtAve, etaGenLeadCM, etaGenSubLeadCM, weight );
        fHM->hGenDijetEtaLeadEtaSubLeadEtaCM->Fill( dijetGenEtaCM,  etaGenLeadCM, etaGenSubLeadCM, weight );
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
        std::cout << Form("Event weight: %f", weight) << std::endl;
    }

    if ( weight <= 0. ) {
        std::cerr << "Error: weight is not positive: " << weight << std::endl;
        return;
    }

    fMcReweight = {1.};
    fIsRecoDijetLabFound = {false};
    fIsRecoDijetCMFound = {false};

    // ptHat value
    // float ptHat = event->ptHat();

    if ( fRecoIdLead < 0 || fRecoIdSubLead < 0 ) {
        if ( fVerbose ) {
            std::cout << "Reco dijet not found" << std::endl;
            std::cout << "DiJetAnalysis::processRecoDijets -- end" << std::endl;
        }
        return;
    }

    unsigned int runId = event->runId();

    // Lead jet
    RecoJet* recoLeadJet = event->recoJetCollection()->at( fRecoIdLead );
    float ptRawRecoLead = recoLeadJet->rawPt();
    float ptRecoLead = recoLeadJet->ptJECCorr();
    float etaRecoLeadLab = etaLab( recoLeadJet->eta() );
    float etaRecoLeadCM = boostEta2CM( recoLeadJet->eta() );
    float phiRecoLead = recoLeadJet->phi();

    // SubLead jet
    RecoJet* recoSubLeadJet = event->recoJetCollection()->at( fRecoIdSubLead );
    float ptRawRecoSubLead = recoSubLeadJet->rawPt();
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
        std::cout << "Error: Reco dijet parameters do not match the Lead and SubLead jet parameters within rounding\n";
        return;
    }

    float dijetRecoPtAveRaw = (ptRawRecoLead + ptRawRecoSubLead) / 2.f;
    float dijetRecoPtAve = fRecoDijet->ptAve();
    float dijetRecoDphi = fRecoDijet->dPhi();
    float dijetRecoEtaLab = fRecoDijet->etaLab();
    float dijetRecoEtaCM = fRecoDijet->etaCM();
    // float dijetRecoDetaCM = fRecoDijet->dEtaCM();
    // float dijetRecoPhi = fRecoDijet->phi();

    // This will allow to perform an entire analysis offline
    // if ( ptRecoSubLead > 30.f && fabs(dijetRecoDphi) > TMath::PiOver2() ) {
    //     double dijetRecoInfo[10] { dijetRecoPtAve, dijetRecoEtaLab, dijetRecoEtaCM, dijetRecoDphi,
    //                                ptRecoLead, etaRecoLeadLab, etaRecoLeadCM,
    //                                ptRecoSubLead, etaRecoSubLeadLab, etaRecoSubLeadCM };
    //     fHM->hRecoDijetInfo->Fill(dijetRecoInfo, weight * fMcReweight);
    // }

    // Reference Lead jet
    GenJet* refLeadJet = {nullptr};
    float ptRefLead{0.};
    float etaRefLeadLab{0.};
    float etaRefLeadCM{0.};
    float phiRefLead{0.};
    // Reference SubLead jet
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
    // float dijetRefPhi{0.};

    // For Monte Carlo check if reco jets have matched gen jets. Skip reco dijet if not
    if ( fIsMc ) {

        // // Lead and SubLead jets must have matching gen jets (just a protection, easy to comment out)
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
            // Matching gen jet for Lead reco jet
            refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
            if ( !refLeadJet ) {
                std::cerr << "Error: Lead jet has no matching gen jet\n";
                return;
            }

            ptRefLead = refLeadJet->pt();
            etaRefLeadLab = etaLab( refLeadJet->eta() );
            etaRefLeadCM = boostEta2CM( refLeadJet->eta() );
            phiRefLead = refLeadJet->phi();
        } 
        
        if ( recoSubLeadJet->hasMatching() ) {
            // Matching gen jet for SubLead reco jet
            refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );
            if ( !refSubLeadJet ) {
                std::cerr << "Error: SubLead jet has no matching gen jet\n";
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
            // dijetRefPhi = fRefDijet->phi();

            if ( fabs( ptRefLead - fRefDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
                fabs( ptRefSubLead - fRefDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefLeadLab - fRefDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefSubLeadLab - fRefDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefLeadCM - fRefDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
                fabs( etaRefSubLeadCM - fRefDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
                fabs( phiRefLead - fRefDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
                fabs( phiRefSubLead - fRefDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
                std::cerr << "Error: Ref dijet reference parameters do not match the Lead and SubLead jet parameters within rounding\n";
                return;
            }

            // if ( fVerbose ) {
            //     std::cout << "Matching gen dijet:\n"; 
            //     dijetRef.print();
            // }
        }
    } // if ( fIsMc )

    if (ptRecoLead > 50. && ptRecoSubLead > 40. && 
        fabs(etaRecoSubLeadLab)<1.2 && 
        fabs(dijetRecoDphi) > TMath::TwoPi()/3. &&
        50. < dijetRecoPtAve && dijetRecoPtAve < 90.) {

        float xj = ptRecoSubLead / ptRecoLead;
        if (xj > 1.) { xj = 1.07; }
        if ( -2.4 < etaRecoLeadLab && etaRecoLeadLab < -1.8) {
            fHM->hRecoDijetXj[0]->Fill( xj, weight * fMcReweight );
        }
        else if ( -1.2 < etaRecoLeadLab && etaRecoLeadLab < 1.2) {
            fHM->hRecoDijetXj[1]->Fill( xj, weight * fMcReweight );
        }
        else if ( 1.8 < etaRecoLeadLab && etaRecoLeadLab < 2.4) {
            fHM->hRecoDijetXj[2]->Fill( xj, weight * fMcReweight );
        }
    } // if (ptRecoLead > 50. && ptRecoSubLead > 40. && fabs(etaRecoSubLeadLab)<1.2 && fabs(dijetRecoDphi) > TMath::TwoPi()/3. && 50. < dijetRecoPtAve && dijetRecoPtAve < 90.)

    if (ptRecoLead > 50. && ptRecoSubLead > 40. && 
        fabs(etaRecoSubLeadCM)<1.2 && 
        fabs(dijetRecoDphi) > TMath::TwoPi()/3. &&
        50. < dijetRecoPtAve && dijetRecoPtAve < 90.) {

        float xj = ptRecoSubLead / ptRecoLead;
        if (xj > 1.) { xj = 1.07; }
        if ( -1.9 < etaRecoLeadCM && etaRecoLeadCM < -1.6) {
            fHM->hRecoDijetXjCM[0]->Fill( xj, weight * fMcReweight );
        }
        else if ( -1.2 < etaRecoLeadCM && etaRecoLeadCM < 1.2) {
            fHM->hRecoDijetXjCM[1]->Fill( xj, weight * fMcReweight );
        }
        else if ( 1.6 < etaRecoLeadCM && etaRecoLeadCM < 1.9) {
            fHM->hRecoDijetXjCM[2]->Fill( xj, weight * fMcReweight );
        }
    } // if (ptRecoLead > 50. && ptRecoSubLead > 40. && fabs(etaRecoSubLeadLab)<1.2 && fabs(dijetRecoDphi) > TMath::TwoPi()/3. && 50. < dijetRecoPtAve && dijetRecoPtAve < 90.)

    // Fill forward and backward distributions for reco dijets in the CM frame (with corrected pT)
    if ( ptRecoLead > 50. && ptRecoSubLead > 40. && fabs(dijetRecoDphi) > (TMath::TwoPi() / 3) ) {
        for (int iEtaCut{0}; iEtaCut < 6; iEtaCut++) {
            const float etaCut = 1.4 + iEtaCut * 0.1;
            if ( fabs(etaRecoLeadCM) < etaCut && fabs(etaRecoSubLeadCM) < etaCut ) {
                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaForwardArr[iEtaCut]->Fill(dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight) : 
                                        fHM->hRecoDijetPtEtaBackwardArr[iEtaCut]->Fill(dijetRecoPtAve, fabs(dijetRecoEtaCM), weight * fMcReweight);
            }
        } // for (int iEtaCut{0}; iEtaCut < 6; iEtaCut++)
    } // if ( ptRecoLead > 50. && ptRecoSubLead > 40. && fabs(dijetRecoDphi) > (TMath::TwoPi() / 3) )


    // Fill forward and backward distributions for reco dijets in the CM frame (with uncorrected pT)
    if ( ptRawRecoLead > 50. && ptRawRecoSubLead > 40. && fabs(dijetRecoDphi) > (TMath::TwoPi() / 3) ) {
        for (int iEtaCut{0}; iEtaCut < 6; iEtaCut++) {
            const float etaCut = 1.4 + iEtaCut * 0.1;

            if ( fabs(etaRecoLeadCM) < etaCut && fabs(etaRecoSubLeadCM) < etaCut ) {
                (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtRawEtaForwardArr[iEtaCut]->Fill(dijetRecoPtAveRaw, dijetRecoEtaCM, weight * fMcReweight) : 
                                        fHM->hRecoDijetPtRawEtaBackwardArr[iEtaCut]->Fill(dijetRecoPtAveRaw, fabs(dijetRecoEtaCM), weight * fMcReweight);
            }
        } // for (int iEtaCut{0}; iEtaCut < 6; iEtaCut++)
    } // if ( ptRawRecoLead > 50. && ptRawRecoSubLead > 40. && fabs(dijetRecoDphi) > (TMath::TwoPi() / 3) )


    //
    // Lab frame
    //

    fIsRecoDijetLabFound = {false}; 
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

        // Correlation between Lead and SubLead
        fHM->hRecoPtLeadPtSublead->Fill( ptRecoLead, ptRecoSubLead, weight );
        fHM->hRecoEtaLeadEtaSublead->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight );
        fHM->hRecoPtLeadPtSubleadMcReweight->Fill( ptRecoLead, ptRecoSubLead, weight * fMcReweight );
        fHM->hRecoEtaLeadEtaSubleadMcReweight->Fill( etaRecoLeadLab, etaRecoSubLeadLab, weight * fMcReweight );

        fHM->hRecoDijetEta->Fill( dijetRecoEtaLab, weight * fMcReweight);
        fHM->hRecoDijetPtEta->Fill( dijetRecoPtAve, dijetRecoEtaLab, 1. );
        fHM->hRecoDijetPtEtaWeighted->Fill( dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight);
        fHM->hRecoDijetPtEtaCMInLab->Fill( dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight);

        (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetPtEtaForward->Fill(dijetRecoPtAve, dijetRecoEtaLab, 1.) : fHM->hRecoDijetPtEtaBackward->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaLab), 1.);
        (dijetRecoEtaLab >= 0) ? fHM->hRecoDijetPtEtaForwardWeighted->Fill(dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight) : fHM->hRecoDijetPtEtaBackwardWeighted->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaLab), weight * fMcReweight);
        (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaForwardCMInLab->Fill(dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetPtEtaBackwardCMInLab->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);

        fHM->hRecoDijetPtAveLeadPtSubLeadPt->Fill( dijetRecoPtAve, ptRecoLead, ptRecoSubLead, weight );
        fHM->hRecoDijetPtAveLeadEtaSubLeadEta->Fill( dijetRecoPtAve, etaRecoLeadLab, etaRecoSubLeadLab, weight );
        fHM->hRecoDijetEtaLeadEtaSubLeadEta->Fill( dijetRecoEtaLab,  etaRecoLeadLab, etaRecoSubLeadLab, weight );

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

            fHM->hRecoDijetPtEtaMatched->Fill( dijetRecoPtAve, dijetRecoEtaLab, weight );

            fHM->hRefPtLeadPtSublead->Fill( ptRefLead, ptRefSubLead, weight );
            fHM->hRefEtaLeadEtaSublead->Fill( etaRefLeadLab, etaRefSubLeadLab, weight );
            fHM->hRefPtLeadPtSubleadMcReweight->Fill( ptRefLead, ptRefSubLead, weight * fMcReweight );
            fHM->hRefEtaLeadEtaSubleadMcReweight->Fill( etaRefLeadLab, etaRefSubLeadLab, weight * fMcReweight );

            double dijetRecoUnfold[12] = { dijetRecoPtAve, dijetRecoEtaLab,
                                           ptRecoLead, etaRecoLeadLab,
                                           ptRecoSubLead, etaRecoSubLeadLab,
                                           dijetRefPtAve, dijetRefEtaLab,
                                           ptRefLead, etaRefLeadLab,
                                           ptRefSubLead, etaRefSubLeadLab };

            double dijetUnfold[4] = { dijetRecoPtAve, dijetRecoEtaLab, dijetRefPtAve, dijetRefEtaLab };

            fHM->hRecoDijetPtEtaRefDijetPtEta->Fill(dijetUnfold, 1.);
            fHM->hRecoDijetPtEtaRefDijetPtEtaWeighted->Fill(dijetUnfold, weight * fMcReweight);

            // fHM->hReco2RefFull->Fill(dijetRecoUnfold, weight * fMcReweight );
            fHM->hRefDijetEta->Fill( dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetEtaVsRecoDijetEta->Fill( dijetRecoEtaLab, dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPt->Fill( dijetRecoEtaLab, dijetRefEtaLab, dijetRecoPtAve, 1.);
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtWeighted->Fill( dijetRecoEtaLab, dijetRefEtaLab, dijetRecoPtAve, weight * fMcReweight );
            // Dijet ptAve must be from RECO!!!
            fHM->hRefDijetPtEta->Fill( dijetRecoPtAve, dijetRefEtaLab, 1. );
            fHM->hRefDijetPtEtaWeighted->Fill( dijetRecoPtAve, dijetRefEtaLab, weight * fMcReweight );
            fHM->hRefDijetPtEtaCMInLab->Fill( dijetRecoPtAve, dijetRefEtaCM, weight * fMcReweight );

            (dijetRefEtaLab >= 0) ? fHM->hRefDijetPtEtaForward->Fill(dijetRecoPtAve, dijetRefEtaLab) : fHM->hRefDijetPtEtaBackward->Fill(dijetRecoPtAve, TMath::Abs(dijetRefEtaLab));
            (dijetRefEtaLab >= 0) ? fHM->hRefDijetPtEtaForwardWeighted->Fill(dijetRecoPtAve, dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetPtEtaBackwardWeighted->Fill(dijetRecoPtAve, TMath::Abs(dijetRefEtaLab), weight * fMcReweight);
            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaForwardCMInLab->Fill(dijetRecoPtAve, dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetPtEtaBackwardCMInLab->Fill(dijetRecoPtAve, TMath::Abs(dijetRefEtaCM), weight * fMcReweight);
        } // if ( fIsMc && recoLeadJet->hasMatching() && recoSubLeadJet->hasMatching() )
    } // if ( fIsRecoDijetLabFound )

    //
    // CM frame
    // 

    fIsRecoDijetCMFound = {false};
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

        fHM->hRecoEtaCMLeadEtaCMSublead->Fill( etaRecoLeadCM, etaRecoSubLeadCM, weight );
        fHM->hRecoDijetEtaCM->Fill( dijetRecoEtaCM, weight * fMcReweight);
        fHM->hRecoDijetPtEtaCM->Fill( dijetRecoPtAve, dijetRecoEtaCM, 1. );
        fHM->hRecoDijetPtEtaCMWeighted->Fill( dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight);
        fHM->hRecoDijetPtEtaLabInCM->Fill( dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight);

        (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaCMForward->Fill(dijetRecoPtAve, dijetRecoEtaCM, 1.) : fHM->hRecoDijetPtEtaCMBackward->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaCM), 1.);
        (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaCMForwardWeighted->Fill(dijetRecoPtAve, dijetRecoEtaCM, weight * fMcReweight) : fHM->hRecoDijetPtEtaCMBackwardWeighted->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaCM), weight * fMcReweight);
        (dijetRecoEtaCM >= 0) ? fHM->hRecoDijetPtEtaForwardLabInCM->Fill(dijetRecoPtAve, dijetRecoEtaLab, weight * fMcReweight) : fHM->hRecoDijetPtEtaBackwardLabInCM->Fill(dijetRecoPtAve, TMath::Abs(dijetRecoEtaLab), weight * fMcReweight);

        fHM->hRecoDijetPtAveLeadPtSubLeadPtCM->Fill( dijetRecoPtAve, ptRecoLead, ptRecoSubLead, weight );
        fHM->hRecoDijetPtAveLeadEtaSubLeadEtaCM->Fill( dijetRecoPtAve, etaRecoLeadCM, etaRecoSubLeadCM, weight );
        fHM->hRecoDijetEtaLeadEtaSubLeadEtaCM->Fill( dijetRecoEtaCM,  etaRecoLeadCM, etaRecoSubLeadCM, weight );

        fHM->hRecoDijetLeadPtEta->Fill( ptRecoLead, etaRecoLeadLab, weight );
        fHM->hRecoDijetLeadPtEtaStdBins->Fill( ptRecoLead, etaRecoLeadLab, weight );
        fHM->hRecoDijetSubLeadPtEta->Fill( ptRecoSubLead, etaRecoSubLeadLab, weight );
        fHM->hRecoDijetSubLeadPtEtaStdBins->Fill( ptRecoSubLead, etaRecoSubLeadLab, weight );
        
        if ( 50. < dijetRecoPtAve && dijetRecoPtAve < 90. ) {
            int idx = -1;
            if ( runId == 285480 ) { idx = 1; }
            else if ( runId == 285505 ) { idx = 2; }
            else if ( runId == 285517 ) { idx = 3; }
            else if ( runId == 285832 ) { idx = 4; }
            else if ( runId == 285993 ) { idx = 5; }

            // For all runs
            fHM->hRecoDijetEtaCMRun[0]->Fill( dijetRecoEtaCM, weight );
            if ( idx > 0 ) fHM->hRecoDijetEtaCMRun[idx]->Fill( dijetRecoEtaCM, weight );
        } // if ( 50. < dijetRecoPtAve && dijetRecoPtAve < 90. )

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

            fHM->hRefEtaCMLeadEtaCMSublead->Fill( etaRefLeadCM, etaRefSubLeadCM, weight );

            fHM->hRecoDijetPtEtaCMMatched->Fill( dijetRecoPtAve, dijetRecoEtaCM, weight );

            fHM->hRefDijetEtaCM->Fill( dijetRefEtaCM, weight );
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCM->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPtAve, 1.);
            fHM->hRefDijetEtaVsRecoDijetEtaVsRecoDijetPtCMWeighted->Fill( dijetRecoEtaCM, dijetRefEtaCM, dijetRecoPtAve, weight * fMcReweight );
            // Dijet ptAve must be from RECO!!!
            fHM->hRefDijetPtEtaCM->Fill( dijetRecoPtAve, dijetRefEtaCM, 1. );
            fHM->hRefDijetPtEtaCMWeighted->Fill( dijetRecoPtAve, dijetRefEtaCM, weight * fMcReweight );
            fHM->hRefDijetPtEtaLabInCM->Fill( dijetRecoPtAve, dijetRefEtaLab, weight * fMcReweight );

            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForward->Fill(dijetRecoPtAve, dijetRefEtaCM, 1.) : fHM->hRefDijetPtEtaCMBackward->Fill(dijetRecoPtAve, TMath::Abs(dijetRefEtaCM), 1.);
            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaCMForwardWeighted->Fill(dijetRecoPtAve, dijetRefEtaCM, weight * fMcReweight) : fHM->hRefDijetPtEtaCMBackwardWeighted->Fill(dijetRecoPtAve, TMath::Abs(dijetRefEtaCM), weight * fMcReweight);
            (dijetRefEtaCM >= 0) ? fHM->hRefDijetPtEtaForwardLabInCM->Fill(dijetRecoPtAve, dijetRefEtaLab, weight * fMcReweight) : fHM->hRefDijetPtEtaBackwardLabInCM->Fill(dijetRecoPtAve, TMath::Abs(dijetRefEtaLab), weight * fMcReweight);
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

    // ptHat value
    // float ptHat = event->ptHat();

    if ( fRefSelRecoIdLead < 0 || fRefSelRecoIdSubLead < 0 ) {
        if ( fVerbose ) {
            std::cout << "Ref-selected reco dijet not found\n";
            std::cout << "DiJetAnalysis::processRefDijets -- end" << std::endl;
        }
        return;
    }

    // Retrieve Lead and SubLead jets
    RecoJet *recoLeadJet = event->recoJetCollection()->at( fRefSelRecoIdLead );
    RecoJet *recoSubLeadJet = event->recoJetCollection()->at( fRefSelRecoIdSubLead );
    GenJet *refLeadJet = event->genJetCollection()->at( recoLeadJet->genJetId() );
    GenJet *refSubLeadJet = event->genJetCollection()->at( recoSubLeadJet->genJetId() );

    // Retrieve kinematic information for reference Lead and SubLead jets, dijets
    float ptRefLead = refLeadJet->pt();
    float etaRefLeadLab = etaLab( refLeadJet->eta() );
    float etaRefLeadCM = boostEta2CM( refLeadJet->eta() );
    float phiRefLead = refLeadJet->phi();

    float ptRefSubLead = refSubLeadJet->pt();
    float etaRefSubLeadLab = etaLab( refSubLeadJet->eta() );
    float etaRefSubLeadCM = boostEta2CM( refSubLeadJet->eta() );
    float phiRefSubLead = refSubLeadJet->phi();

    // if (fVerbose) {
    //     std::cout << "Ref-selected inclusive Lead and SubLead jets:" << std::endl;
    //     std::cout << Form("Ref Lead jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRefLead, etaRefLeadLab, etaRefLeadCM, phiRefLead);
    //     std::cout << Form("Ref SubLead jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRefSubLead, etaRefSubLeadLab, etaRefSubLeadCM, phiRefSubLead);
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
    // float dijetRefDphi = fRefDijet->dPhi();
    float dijetRefPhi = fRefDijet->phi();

    if ( fabs(ptRefLead - fRefDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(ptRefSubLead - fRefDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefLeadLab - fRefDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefSubLeadLab - fRefDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefLeadCM - fRefDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRefSubLeadCM - fRefDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRefLead - fRefDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRefSubLead - fRefDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
        std::cout << "Error: Ref dijet reference parameters do not match the Lead and SubLead jet parameters within rounding\n";
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
    //     std::cout << "Ref-selected inclusive reconstructed Lead and SubLead jets:" << std::endl;
    //     std::cout << Form("Reco Lead jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRecoLead, etaRecoLeadLab, etaRecoLeadCM, phiRecoLead);
    //     std::cout << Form("Reco SubLead jet pt: %5.2f etaLab: %5.2f etaCM: %5.2f phi: %5.2f\n", ptRecoSubLead, etaRecoSubLeadLab, etaRecoSubLeadCM, phiRecoSubLead);
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
    // float dijetRecoDphi = fRecoDijet->dPhi();
    float dijetRecoPhi = fRecoDijet->phi();

    if ( fabs(ptRecoLead - fRecoDijet->leadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(ptRecoSubLead - fRecoDijet->subLeadJetPt()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoLeadLab - fRecoDijet->leadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoSubLeadLab - fRecoDijet->subLeadJetEtaLab()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoLeadCM - fRecoDijet->leadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(etaRecoSubLeadCM - fRecoDijet->subLeadJetEtaCM()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRecoLead - fRecoDijet->leadJetPhi()) > std::numeric_limits<float>::epsilon() ||
         fabs(phiRecoSubLead - fRecoDijet->subLeadJetPhi()) > std::numeric_limits<float>::epsilon() ) {
        std::cout << "Error: Reco dijet reference parameters do not match the Lead and SubLead jet parameters within rounding\n";
        return;
    }

    // if (fVerbose) {
    //     std::cout << "Ref-selected inclusive reconstructed dijet parameters:" << std::endl;
    //     dijetReco.print();
    // }


    //
    // Lab frame
    //

    fIsRefSelDijetLabFound = {false};
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
            std::cout << Form("Ref dijet ptAve: %5.2f eta: %5.2f phi: %5.2f --> Reco dijet ptAve: %5.2f eta: %5.2f phi: %5.2f\n", 
                              dijetRefPtAve, dijetRefEtaLab, dijetRefPhi, 
                              dijetRecoPtAve, dijetRecoEtaLab, dijetRecoPhi) << std::endl;
        }

        // Dijet reco vs ref for unfolding
        double dijetRecoUnfold[12] = { dijetRecoPtAve, dijetRecoEtaLab,
                                       ptRecoLead, etaRecoLeadLab,
                                       ptRecoSubLead, etaRecoSubLeadLab,
                                       dijetRefPtAve, dijetRefEtaLab,
                                       ptRefLead, etaRefLeadLab,
                                       ptRefSubLead, etaRefSubLeadLab };    

        // fHM->hRefSel2RecoFull->Fill(dijetRecoUnfold, weight * fMcReweight );
        fHM->hRefSelDijetEta->Fill(dijetRefEtaLab, weight * fMcReweight );
        fHM->hRefSelDijetPtEta->Fill(dijetRefPtAve, dijetRefEtaLab, 1.);
        fHM->hRefSelDijetPtEtaWeighted->Fill(dijetRefPtAve, dijetRefEtaLab, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaCMInLab->Fill(dijetRefPtAve, dijetRefEtaCM, weight * fMcReweight );

    } // if ( fIsRefSelDijetLabFound )

    //
    // CM frame
    //

    fIsRefSelDijetCMFound = {false};
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
            std::cout << Form("Ref dijet ptAve: %5.2f eta: %5.2f phi: %5.2f --> Reco dijet ptAve: %5.2f eta: %5.2f phi: %5.2f", 
                              dijetRefPtAve, dijetRefEtaCM, dijetRefPhi, 
                              dijetRecoPtAve, dijetRecoEtaCM, dijetRecoPhi) << std::endl;
        }

        fHM->hRefSelDijetEtaCM->Fill(dijetRefEtaCM, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaCM->Fill(dijetRefPtAve, dijetRefEtaCM, 1.);
        fHM->hRefSelDijetPtEtaCMWeighted->Fill(dijetRefPtAve, dijetRefEtaCM, weight * fMcReweight );
        fHM->hRefSelDijetPtEtaLabInCM->Fill(dijetRefPtAve, dijetRefEtaLab, weight * fMcReweight );
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

    if ( fVerbose ) {
        std::cout << "Event weight: " << weight << std::endl;
    }

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

    //
    // Loop over jet collections (reco, gen, ref)
    //
    processInclusiveJets(event, weight);

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
