#pragma once
#include "Constants.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <memory>
namespace jetanalysis {
struct UnfoldingPairHistograms {
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMMatchedGenSelected[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMMatchedGenFound[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMSelectedRecoFound[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMSelectedRecoNotFound[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMatchedRecoSelected[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMSelectedGenNotFound[nEtaCuts];
};

// Extended truth is the same-event global gen-leading pair, irrespective of
// gen cuts. The category distinguishes outside-fiducial (1) and fiducial (2).
struct ExtendedUnfoldingHistograms {
    std::unique_ptr<THnSparseD> response[nEtaCuts];
    std::unique_ptr<THnSparseD> truth[nEtaCuts];
    std::unique_ptr<THnSparseD> matchedTruth[nEtaCuts];
    std::unique_ptr<THnSparseD> miss[nEtaCuts];
    std::unique_ptr<TH2D> measured[nEtaCuts];
    std::unique_ptr<TH2D> matchedMeasured[nEtaCuts];
    std::unique_ptr<TH2D> fake[nEtaCuts];
};

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

    std::unique_ptr<TH2D> hGenDijetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hGenDijetPtEtaLab;
    std::unique_ptr<TH2D> hGenDijetPtEtaCM;

    std::unique_ptr<TH2D> hGenDijetPtEtaLabArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<THnSparseD> hGenDijetPtEtaCMXPbArr[nEtaCuts];
    std::unique_ptr<THnSparseD> hGenDijetPtEtaCMXpArr[nEtaCuts];
    std::unique_ptr<THnSparseD> hGenDijetPtEtaCMXPbOverXpArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMControlledSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMStudiedSampleArr[nEtaCuts];
#ifdef DO_MIXED_PRIOR
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMixedPriorArr[nEtaCuts];
#endif

    // Gen dijet smeared histograms
    std::unique_ptr<TH2D> hGenDijetDefPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefExtraPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefDoublePtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetDefNinetyPtEtaCMArr[nEtaCuts];

    ExtendedUnfoldingHistograms unfoldingExtended;
    ExtendedUnfoldingHistograms unfoldingExtendedJerDefExtra;
    ExtendedUnfoldingHistograms unfoldingExtendedControlled;
    ExtendedUnfoldingHistograms unfoldingExtendedStudied;

    UnfoldingPairHistograms unfoldingPairs;
    UnfoldingPairHistograms unfoldingPairsJerDefExtra;
    UnfoldingPairHistograms unfoldingPairsControlled;
    UnfoldingPairHistograms unfoldingPairsStudied;

    // Unfolding histograms (response matrices and missing histograms)
    std::unique_ptr<THnSparseF> hGenDijetPtEtaCMVsRecoPtEtaCMArr[nEtaCuts];
    std::unique_ptr<THnSparseF> hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr[nEtaCuts];
    std::unique_ptr<THnSparseF> hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr[nEtaCuts];
    // std::unique_ptr<THnSparseF> hRefDijetPtEtaCMVsRecoPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMissArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMissControlledSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMissStudiedSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMFakeArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMFakeControlledSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMFakeStudiedSampleArr[nEtaCuts];
    std::unique_ptr<THnSparseF> hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hGenDijetPtEtaCMMissJerDefExtraArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMFakeJerDefExtraArr[nEtaCuts];
    std::unique_ptr<TH1D> hUnfoldingPairClassificationArr[nEtaCuts];
    std::unique_ptr<TH1D> hUnfoldingPairClassificationJerDefExtraArr[nEtaCuts];

    //
    // Ref-selected level histograms
    //
    std::unique_ptr<TH2D> hRefSelDijetPtEtaCMArr[nEtaCuts];

    //
    // Reco level histograms
    //

    // Reco jet histograms
    std::unique_ptr<TH2D> hRecoInclusiveJetPtEtaLabUnflippedUnweighted;
    std::unique_ptr<TH2D> hRecoInclusiveJetRawPtEtaLabUnflipped;
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
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMRunArr[nDiagnosticRunHistograms];
    std::unique_ptr<TH2D> hRecoDijetPtEtaLabArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMControlledSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMStudiedSampleArr[nEtaCuts];
    
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefNoSF[nEtaCuts];             // JER no SF variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefDoubleSF[nEtaCuts];         // JER default variation with double SF for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefNinetySF[nEtaCuts];         // JER default variation with 90% SF for reco dijet histograms

    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerUp[nEtaCuts];         // JER up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDown[nEtaCuts];       // JER down variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDef[nEtaCuts];        // JER default variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDownExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerUpExtra[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefUnfold[nEtaCuts];        // JER default variation for reco dijet histograms for unfolding
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMJerDefUnfoldControlledSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMJerDefUnfoldStudiedSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJerDefExtraUnfold[nEtaCuts];

    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJeuUp[nEtaCuts];         // JEU up variation for reco dijet histograms
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMArrJeuDown[nEtaCuts];       // JEU down variation for reco dijet histograms

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
    std::unique_ptr<TH2D> hRefSelInclusiveJetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefSelInclusiveJetPtEtaLab;
    std::unique_ptr<TH2D> hRefSelInclusiveJetPtEtaCM;
    std::unique_ptr<TH2D> hRefSelInclusiveJetTrkMaxPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefSelInclusiveJetTrkMaxPtEtaLab;
    std::unique_ptr<TH2D> hRefSelInclusiveJetTrkMaxPtEtaCM;
    std::unique_ptr<TH2D> hRefSelInclusiveJetNoSelPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefSelInclusiveJetNoSelPtEtaLab;
    std::unique_ptr<TH2D> hRefSelInclusiveJetNoSelPtEtaCM;
    std::unique_ptr<TH2D> hRefInclusiveJetPtEtaStdBins;

    // Ref dijet histograms

    std::unique_ptr<TH2D> hRefDijetPtEtaLabUnflipped;
    std::unique_ptr<TH2D> hRefDijetPtEtaLab;
    std::unique_ptr<TH2D> hRefDijetPtEtaCM;

    std::unique_ptr<TH2D> hRefDijetPtEtaLabArr[nEtaCuts];
    std::unique_ptr<TH2D> hRecoDijetPtEtaLabPointingResArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaLabPointingResArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaLabPointingResRefPtAveArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaCMArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaCMControlledSampleArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaCMStudiedSampleArr[nEtaCuts];

    std::unique_ptr<TH2D> hRefDijetPtEtaCMArrJerUp[nEtaCuts];         // JER up variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaCMArrJerDown[nEtaCuts];       // JER down variation for ref dijet histograms
    std::unique_ptr<TH2D> hRefDijetPtEtaCMArrJerDef[nEtaCuts];        // JER default variation for ref dijet histograms

    // MC-only pointing diagnostic. Reco and legacy Ref use reconstructed pTave and
    // require a matched dijet pair. Reco applies the acceptance to reconstructed
    // jet eta, while Ref applies it to matched-reference jet eta.
    // RefPtAve retains that Ref population with matched-reference pTave on x.
    std::unique_ptr<TH2D> hRecoDijetPtEtaCMPointingResArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaCMPointingResArr[nEtaCuts];
    std::unique_ptr<TH2D> hRefDijetPtEtaCMPointingResRefPtAveArr[nEtaCuts];

};


void createHistograms(Histograms&, const bool&);
}
