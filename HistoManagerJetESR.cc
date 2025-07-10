/**
 * @file HistoManagerJetESR.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Histograms for JES and JER studies
 * @version 1.01
 * @date 2025-01-06
 * 
 * @copyright Copyright (c) 2025
 * 
 */

// Jet analysis headers
#include "HistoManagerJetESR.h"

// ROOT headers
#include "TObject.h"
#include "TIterator.h"
#include "TString.h"
#include "TClass.h"
#include "TKey.h"
#include "TROOT.h"
#include "TSystem.h"

//________________
HistoManagerJetESR::HistoManagerJetESR() :
    BaseHistoManager(),

    //
    // Event histograms
    //
    hVz{nullptr}, 
    hVzCentWeighted{nullptr},
    hVzPtHatWeighted{nullptr},
    hVzWeighted{nullptr},
    hPtHat{nullptr},
    hPtHatWeighted{nullptr},
    hHiBin{nullptr},
    hHiBinWeighted{nullptr},
    hPtHatWeight{nullptr},
    hPtHatWeightWeighted{nullptr},
    hCentrality{nullptr},
    hCentralityWeighted{nullptr},
    hVzPtHatCent{nullptr},
    hVzPtHatCentWeighted{nullptr},

    //
    // Gen jet histograms
    //

    hGenJetCollectionSize{nullptr},    
    hGenVsRecoJetCollectionSize{nullptr},

    hGenLeadingJetPtOverPtHatVsLeadingJetPt{nullptr},
    hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted{nullptr},
    hGenDijetPtOverPtHatVsDijetPt{nullptr},
    hGenDijetPtOverPtHatVsDijetPtWeighted{nullptr},
    hGenDijetPtAveOverPtHatVsDijetPtAve{nullptr},
    hGenDijetPtAveOverPtHatVsDijetPtAveWeighted{nullptr},

    hGenInclusiveJetPt{nullptr},
    hGenInclusiveJetEta{nullptr},
    hGenInclusiveJetEtaUnweighted{nullptr},
    hGenInclusiveJetPtEta{nullptr},
    hGenInclusiveJetPtEtaPtHat{nullptr},
    hGenLeadJetPtEta{nullptr},
    hGenLeadJetPtEtaPtHat{nullptr},
    hGenSubLeadJetPtEta{nullptr},
    hGenSubLeadJetPtEtaPtHat{nullptr},
    hGenGoodInclusiveJetEtaLabFrame{nullptr},
    hGenGoodInclusiveJetEtaCMFrame{nullptr},

    //
    // Reco jet histograms
    //

    hRecoJetCollectionSize{nullptr},

    hRecoLeadingJetPtOverPtHatVsLeadingJetPt{nullptr},
    hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted{nullptr},
    hRecoDijetPtOverPtHatVsDijetPt{nullptr},
    hRecoDijetPtOverPtHatVsDijetPtWeighted{nullptr},
    hRecoDijetPtAveOverPtHatVsDijetPtAve{nullptr},
    hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted{nullptr},

    hRecoInclusiveJetNHF{nullptr},
    hRecoInclusiveJetNEmF{nullptr},
    hRecoInclusiveJetNumOfConst{nullptr},
    hRecoInclusiveJetMUF{nullptr},
    hRecoInclusiveJetCHF{nullptr},
    hRecoInclusiveJetChargedMult{nullptr},
    hRecoInclusiveJetCEmF{nullptr},
    hRecoInclusiveJetNumOfNeutPart{nullptr},

    hRecoInclusiveAllJetPtRawEta{nullptr},

    hRecoInclusiveAllJetPt{nullptr},
    hRecoInclusiveAllJetEta{nullptr},
    hRecoInclusiveAllJetEtaUnweighted{nullptr},
    hRecoInclusiveAllJetPtEta{nullptr},
    hRecoInclusiveAllJetPtEtaPtHat{nullptr},
    hRecoLeadAllJetPtEta{nullptr},
    hRecoLeadAllJetPtEtaPtHat{nullptr},
    hRecoSubLeadAllJetPtEta{nullptr},
    hRecoSubLeadAllJetPtEtaPtHat{nullptr},

    hRecoGoodInclusiveJetEtaLabFrame{nullptr},
    hRecoGoodInclusiveJetEtaCMFrame{nullptr},

    //
    // Ref jet histograms
    //

    hRefInclusiveJetPt{nullptr},
    hRefInclusiveJetEta{nullptr},
    hRefInclusiveJetEtaUnweighted{nullptr},
    hRefInclusiveJetPtEta{nullptr},
    hRefInclusiveJetPtEtaPtHat{nullptr},
    hRefLeadJetPtEta{nullptr},
    hRefLeadJetPtEtaPtHat{nullptr},
    hRefLeadUnswappedJetPtEta{nullptr},
    hRefLeadUnswappedJetPtEtaPtHat{nullptr},
    hRefSubLeadJetPtEta{nullptr},
    hRefSubLeadJetPtEtaPtHat{nullptr},
    hRefSubLeadUnswappedJetPtEta{nullptr},
    hRefSubLeadUnswappedJetPtEtaPtHat{nullptr},

    hRecoInclusiveMatchedJetPt{nullptr},
    hRecoInclusiveMatchedJetPtEta{nullptr},
    hRecoInclusiveMatchedJetPtEtaPtHat{nullptr},

    hRecoInclusiveUnmatchedJetPt{nullptr},
    hRecoInclusiveUnmatchedJetPtEta{nullptr},
    hRecoInclusiveUnmatchedJetPtEtaPtHat{nullptr},

    hRecoLeadMatchedJetPtEta{nullptr},
    hRecoLeadMatchedJetPtEtaPtHat{nullptr},
    hRecoLeadUnmatchedJetPtEta{nullptr},
    hRecoLeadUnmatchedJetPtEtaPtHat{nullptr},

    hRecoSubLeadMatchedJetPtEta{nullptr},
    hRecoSubLeadMatchedJetPtEtaPtHat{nullptr},
    hRecoSubLeadUnmatchedJetPtEta{nullptr},
    hRecoSubLeadUnmatchedJetPtEtaPtHat{nullptr},

    hInclusiveJetJESVsPtGen{nullptr},
    hInclusiveJetJESGenPtGenEtaPtHatWeighted{nullptr},
    hInclusiveJetJESRecoPtRecoEtaPtHatWeighted{nullptr},
    hLeadingJetJESGenPtEtaPtHatWeighted{nullptr},
    hSubleadingJetJESGenPtEtaPtHatWeighted{nullptr},

    hInclusiveReco2RefJetPtEtaPhiPtHat{nullptr},
    hLeadReco2RefJetPtEtaPhiPtHat{nullptr},
    hSubLeadReco2RefJetPtEtaPhiPtHat{nullptr},

    //
    // Ref-selected jet histograms
    //

    hRefSelInclusiveJetPt{nullptr},
    hRefSelInclusiveJetEta{nullptr},
    hRefSelInclusiveJetEtaUnweighted{nullptr},
    hRefSelInclusiveJetPtEta{nullptr},
    hRefSelInclusiveJetPtEtaPtHat{nullptr},
    hRefSelLeadJetPtEta{nullptr},
    hRefSelLeadJetPtEtaPtHat{nullptr},
    hRefSelSubLeadJetPtEta{nullptr},
    hRefSelSubLeadJetPtEtaPtHat{nullptr},

    //
    // Private members
    //

    fPtHatBins{100}, fPtHatRange{15., 1015.},
    fCentBins{101}, fCentRange{-0.5, 100.5},
    fPtBins{150}, fPtRange{5., 1505.}, 
    fEtaBins{52}, fEtaRange{-5.2, 5.2},
    fPhiBins{32}, fPhiRange{-TMath::Pi(), TMath::Pi()},
    fJESBins{100}, fJESRange{0., 2.},
    fVerbose{false}, fIsMc{false} {    
    // Empty constructor
}

//________________
HistoManagerJetESR::~HistoManagerJetESR() {
    // Destructor

    //
    // Event histograms
    //
    if (hVz) { delete hVz; hVz = nullptr; }
    if (hVzCentWeighted) { delete hVzCentWeighted; hVzCentWeighted = nullptr; }
    if (hVzPtHatWeighted) { delete hVzPtHatWeighted; hVzPtHatWeighted = nullptr; }
    if (hVzWeighted) { delete hVzWeighted; hVzWeighted = nullptr; }
    if (hPtHat) { delete hPtHat; hPtHat = nullptr; }
    if (hPtHatWeighted) { delete hPtHatWeighted; hPtHatWeighted = nullptr; }
    if (hHiBin) { delete hHiBin; hHiBin = nullptr; }
    if (hHiBinWeighted) { delete hHiBinWeighted; hHiBinWeighted = nullptr; }
    if (hPtHatWeight) { delete hPtHatWeight; hPtHatWeight = nullptr; }
    if (hPtHatWeightWeighted) { delete hPtHatWeightWeighted; hPtHatWeightWeighted = nullptr; }
    if (hCentrality) { delete hCentrality; hCentrality = nullptr; }
    if (hCentralityWeighted) { delete hCentralityWeighted; hCentralityWeighted = nullptr; }
    if (hVzPtHatCent) { delete hVzPtHatCent; hVzPtHatCent = nullptr; }
    if (hVzPtHatCentWeighted) { delete hVzPtHatCentWeighted; hVzPtHatCentWeighted = nullptr; }

    //
    // Gen jet histograms
    //
    if ( fIsMc ) {
        if (hGenJetCollectionSize) { delete hGenJetCollectionSize; hGenJetCollectionSize = nullptr; }
        if (hGenVsRecoJetCollectionSize) { delete hGenVsRecoJetCollectionSize; hGenVsRecoJetCollectionSize = nullptr; }

        if (hGenLeadingJetPtOverPtHatVsLeadingJetPt) { delete hGenLeadingJetPtOverPtHatVsLeadingJetPt; hGenLeadingJetPtOverPtHatVsLeadingJetPt = nullptr; }
        if (hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted) { delete hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted; hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted = nullptr; }
        if (hGenDijetPtOverPtHatVsDijetPt) { delete hGenDijetPtOverPtHatVsDijetPt; hGenDijetPtOverPtHatVsDijetPt = nullptr; }
        if (hGenDijetPtOverPtHatVsDijetPtWeighted) { delete hGenDijetPtOverPtHatVsDijetPtWeighted; hGenDijetPtOverPtHatVsDijetPtWeighted = nullptr; }
        if (hGenDijetPtAveOverPtHatVsDijetPtAve) { delete hGenDijetPtAveOverPtHatVsDijetPtAve; hGenDijetPtAveOverPtHatVsDijetPtAve = nullptr; }
        if (hGenDijetPtAveOverPtHatVsDijetPtAveWeighted) { delete hGenDijetPtAveOverPtHatVsDijetPtAveWeighted; hGenDijetPtAveOverPtHatVsDijetPtAveWeighted = nullptr; }

        if (hGenInclusiveJetPt) { delete hGenInclusiveJetPt; hGenInclusiveJetPt = nullptr; }
        if (hGenInclusiveJetEta) { delete hGenInclusiveJetEta; hGenInclusiveJetEta = nullptr; }
        if (hGenInclusiveJetEtaUnweighted) { delete hGenInclusiveJetEtaUnweighted; hGenInclusiveJetEtaUnweighted = nullptr; }
        if (hGenInclusiveJetPtEta) { delete hGenInclusiveJetPtEta; hGenInclusiveJetPtEta = nullptr; }
        if (hGenInclusiveJetPtEtaPtHat) { delete hGenInclusiveJetPtEtaPtHat; hGenInclusiveJetPtEtaPtHat = nullptr; }
        if (hGenLeadJetPtEta) { delete hGenLeadJetPtEta; hGenLeadJetPtEta = nullptr; }
        if (hGenLeadJetPtEtaPtHat) { delete hGenLeadJetPtEtaPtHat; hGenLeadJetPtEtaPtHat = nullptr; }
        if (hGenSubLeadJetPtEta) { delete hGenSubLeadJetPtEta; hGenSubLeadJetPtEta = nullptr; }
        if (hGenSubLeadJetPtEtaPtHat) { delete hGenSubLeadJetPtEtaPtHat; hGenSubLeadJetPtEtaPtHat = nullptr; }
        if (hGenGoodInclusiveJetEtaLabFrame) { delete hGenGoodInclusiveJetEtaLabFrame; hGenGoodInclusiveJetEtaLabFrame = nullptr; }
        if (hGenGoodInclusiveJetEtaCMFrame) { delete hGenGoodInclusiveJetEtaCMFrame; hGenGoodInclusiveJetEtaCMFrame = nullptr; }
    }

    //
    // Reco jet histograms
    //

    if (hRecoJetCollectionSize) { delete hRecoJetCollectionSize; hRecoJetCollectionSize = nullptr; }

    for (int i = 0; i < 4; ++i) {
        if (hRecoInclusiveJetNHF[i]) { delete hRecoInclusiveJetNHF[i]; hRecoInclusiveJetNHF[i] = nullptr; }
        if (hRecoInclusiveJetNEmF[i]) { delete hRecoInclusiveJetNEmF[i]; hRecoInclusiveJetNEmF[i] = nullptr; }
        if (hRecoInclusiveJetNumOfConst[i]) { delete hRecoInclusiveJetNumOfConst[i]; hRecoInclusiveJetNumOfConst[i] = nullptr; }
        if (hRecoInclusiveJetMUF[i]) { delete hRecoInclusiveJetMUF[i]; hRecoInclusiveJetMUF[i] = nullptr; }
        if (hRecoInclusiveJetCHF[i]) { delete hRecoInclusiveJetCHF[i]; hRecoInclusiveJetCHF[i] = nullptr; }
        if (hRecoInclusiveJetChargedMult[i]) { delete hRecoInclusiveJetChargedMult[i]; hRecoInclusiveJetChargedMult[i] = nullptr; }
        if (hRecoInclusiveJetCEmF[i]) { delete hRecoInclusiveJetCEmF[i]; hRecoInclusiveJetCEmF[i] = nullptr; }
        if (hRecoInclusiveJetNumOfNeutPart[i]) { delete hRecoInclusiveJetNumOfNeutPart[i]; hRecoInclusiveJetNumOfNeutPart[i] = nullptr; }
    }

    if (hRecoInclusiveAllJetPtRawEta) { delete hRecoInclusiveAllJetPtRawEta; hRecoInclusiveAllJetPtRawEta = nullptr; }

    if (hRecoInclusiveAllJetPt) { delete hRecoInclusiveAllJetPt; hRecoInclusiveAllJetPt = nullptr; }
    if (hRecoInclusiveAllJetEta) { delete hRecoInclusiveAllJetEta; hRecoInclusiveAllJetEta = nullptr; }
    if (hRecoInclusiveAllJetEtaUnweighted) { delete hRecoInclusiveAllJetEtaUnweighted; hRecoInclusiveAllJetEtaUnweighted = nullptr; }
    if (hRecoInclusiveAllJetPtEta) { delete hRecoInclusiveAllJetPtEta; hRecoInclusiveAllJetPtEta = nullptr; }
    if (hRecoInclusiveAllJetPtEtaPtHat) { delete hRecoInclusiveAllJetPtEtaPtHat; hRecoInclusiveAllJetPtEtaPtHat = nullptr; }
    if (hRecoLeadAllJetPtEta) { delete hRecoLeadAllJetPtEta; hRecoLeadAllJetPtEta = nullptr; }
    if (hRecoLeadAllJetPtEtaPtHat) { delete hRecoLeadAllJetPtEtaPtHat; hRecoLeadAllJetPtEtaPtHat = nullptr; }
    if (hRecoSubLeadAllJetPtEta) { delete hRecoSubLeadAllJetPtEta; hRecoSubLeadAllJetPtEta = nullptr; }
    if (hRecoSubLeadAllJetPtEtaPtHat) { delete hRecoSubLeadAllJetPtEtaPtHat; hRecoSubLeadAllJetPtEtaPtHat = nullptr; }

    if (hRecoGoodInclusiveJetEtaLabFrame) { delete hRecoGoodInclusiveJetEtaLabFrame; hRecoGoodInclusiveJetEtaLabFrame = nullptr; }
    if (hRecoGoodInclusiveJetEtaCMFrame) { delete hRecoGoodInclusiveJetEtaCMFrame; hRecoGoodInclusiveJetEtaCMFrame = nullptr; }

    //
    // Ref jet histograms
    //

    if ( fIsMc ) {
        if (hRefInclusiveJetPt) { delete hRefInclusiveJetPt; hRefInclusiveJetPt = nullptr; }
        if (hRefInclusiveJetEta) { delete hRefInclusiveJetEta; hRefInclusiveJetEta = nullptr; }
        if (hRefInclusiveJetEtaUnweighted) { delete hRefInclusiveJetEtaUnweighted; hRefInclusiveJetEtaUnweighted = nullptr; }
        if (hRefInclusiveJetPtEta) { delete hRefInclusiveJetPtEta; hRefInclusiveJetPtEta = nullptr; }
        if (hRefInclusiveJetPtEtaPtHat) { delete hRefInclusiveJetPtEtaPtHat; hRefInclusiveJetPtEtaPtHat = nullptr; }
        if (hRefLeadJetPtEta) { delete hRefLeadJetPtEta; hRefLeadJetPtEta = nullptr; }
        if (hRefLeadJetPtEtaPtHat) { delete hRefLeadJetPtEtaPtHat; hRefLeadJetPtEtaPtHat = nullptr; }
        if (hRefLeadUnswappedJetPtEta) { delete hRefLeadUnswappedJetPtEta; hRefLeadUnswappedJetPtEta = nullptr; }
        if (hRefLeadUnswappedJetPtEtaPtHat) { delete hRefLeadUnswappedJetPtEtaPtHat; hRefLeadUnswappedJetPtEtaPtHat = nullptr; }
        if (hRefSubLeadJetPtEta) { delete hRefSubLeadJetPtEta; hRefSubLeadJetPtEta = nullptr; }
        if (hRefSubLeadJetPtEtaPtHat) { delete hRefSubLeadJetPtEtaPtHat; hRefSubLeadJetPtEtaPtHat = nullptr; }
        if (hRefSubLeadUnswappedJetPtEta) { delete hRefSubLeadUnswappedJetPtEta; hRefSubLeadUnswappedJetPtEta = nullptr; }
        if (hRefSubLeadUnswappedJetPtEtaPtHat) { delete hRefSubLeadUnswappedJetPtEtaPtHat; hRefSubLeadUnswappedJetPtEtaPtHat = nullptr; }

        if (hRecoLeadingJetPtOverPtHatVsLeadingJetPt) { delete hRecoLeadingJetPtOverPtHatVsLeadingJetPt; hRecoLeadingJetPtOverPtHatVsLeadingJetPt = nullptr; };
        if (hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted) { delete hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted; hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted = nullptr; };
        if (hRecoDijetPtOverPtHatVsDijetPt) { delete hRecoDijetPtOverPtHatVsDijetPt; hRecoDijetPtOverPtHatVsDijetPt = nullptr; };
        if (hRecoDijetPtOverPtHatVsDijetPtWeighted) { delete hRecoDijetPtOverPtHatVsDijetPtWeighted; hRecoDijetPtOverPtHatVsDijetPtWeighted = nullptr; };
        if (hRecoDijetPtAveOverPtHatVsDijetPtAve) { delete hRecoDijetPtAveOverPtHatVsDijetPtAve; hRecoDijetPtAveOverPtHatVsDijetPtAve = nullptr; };
        if (hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted) { delete hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted; hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted = nullptr; };

        if (hRecoInclusiveMatchedJetPt) { delete hRecoInclusiveMatchedJetPt; hRecoInclusiveMatchedJetPt = nullptr; }
        if (hRecoInclusiveMatchedJetPtEta) { delete hRecoInclusiveMatchedJetPtEta; hRecoInclusiveMatchedJetPtEta = nullptr; }
        if (hRecoInclusiveMatchedJetPtEtaPtHat) { delete hRecoInclusiveMatchedJetPtEtaPtHat; hRecoInclusiveMatchedJetPtEtaPtHat = nullptr; }

        if (hRecoInclusiveUnmatchedJetPt) { delete hRecoInclusiveUnmatchedJetPt; hRecoInclusiveUnmatchedJetPt = nullptr; }
        if (hRecoInclusiveUnmatchedJetPtEta) { delete hRecoInclusiveUnmatchedJetPtEta; hRecoInclusiveUnmatchedJetPtEta = nullptr; }
        if (hRecoInclusiveUnmatchedJetPtEtaPtHat) { delete hRecoInclusiveUnmatchedJetPtEtaPtHat; hRecoInclusiveUnmatchedJetPtEtaPtHat = nullptr; }

        if (hRecoLeadMatchedJetPtEta) { delete hRecoLeadMatchedJetPtEta; hRecoLeadMatchedJetPtEta = nullptr; }
        if (hRecoLeadMatchedJetPtEtaPtHat) { delete hRecoLeadMatchedJetPtEtaPtHat; hRecoLeadMatchedJetPtEtaPtHat = nullptr; }
        if (hRecoLeadUnmatchedJetPtEta) { delete hRecoLeadUnmatchedJetPtEta; hRecoLeadUnmatchedJetPtEta = nullptr; }
        if (hRecoLeadUnmatchedJetPtEtaPtHat) { delete hRecoLeadUnmatchedJetPtEtaPtHat; hRecoLeadUnmatchedJetPtEtaPtHat = nullptr; }

        if (hRecoSubLeadMatchedJetPtEta) { delete hRecoSubLeadMatchedJetPtEta; hRecoSubLeadMatchedJetPtEta = nullptr; }
        if (hRecoSubLeadMatchedJetPtEtaPtHat) { delete hRecoSubLeadMatchedJetPtEtaPtHat; hRecoSubLeadMatchedJetPtEtaPtHat = nullptr; }
        if (hRecoSubLeadUnmatchedJetPtEta) { delete hRecoSubLeadUnmatchedJetPtEta; hRecoSubLeadUnmatchedJetPtEta = nullptr; }
        if (hRecoSubLeadUnmatchedJetPtEtaPtHat) { delete hRecoSubLeadUnmatchedJetPtEtaPtHat; hRecoSubLeadUnmatchedJetPtEtaPtHat = nullptr; }

        if (hInclusiveJetJESVsPtGen) { delete hInclusiveJetJESVsPtGen; hInclusiveJetJESVsPtGen = nullptr; }
        if (hInclusiveJetJESGenPtGenEtaPtHatWeighted) { delete hInclusiveJetJESGenPtGenEtaPtHatWeighted; hInclusiveJetJESGenPtGenEtaPtHatWeighted = nullptr; }
        if (hInclusiveJetJESRecoPtRecoEtaPtHatWeighted) { delete hInclusiveJetJESRecoPtRecoEtaPtHatWeighted; hInclusiveJetJESRecoPtRecoEtaPtHatWeighted = nullptr; }
        if (hLeadingJetJESGenPtEtaPtHatWeighted) { delete hLeadingJetJESGenPtEtaPtHatWeighted; hLeadingJetJESGenPtEtaPtHatWeighted = nullptr; }
        if (hSubleadingJetJESGenPtEtaPtHatWeighted) { delete hSubleadingJetJESGenPtEtaPtHatWeighted; hSubleadingJetJESGenPtEtaPtHatWeighted = nullptr; }

        if (hInclusiveReco2RefJetPtEtaPhiPtHat) { delete hInclusiveReco2RefJetPtEtaPhiPtHat; hInclusiveReco2RefJetPtEtaPhiPtHat = nullptr; }
        if (hLeadReco2RefJetPtEtaPhiPtHat) { delete hLeadReco2RefJetPtEtaPhiPtHat; hLeadReco2RefJetPtEtaPhiPtHat = nullptr; }
        if (hSubLeadReco2RefJetPtEtaPhiPtHat) { delete hSubLeadReco2RefJetPtEtaPhiPtHat; hSubLeadReco2RefJetPtEtaPhiPtHat = nullptr; }
    }

    //
    // Ref-selected jet histograms
    //

    if ( fIsMc ) {
        if (hRefSelInclusiveJetPt) { delete hRefSelInclusiveJetPt; hRefSelInclusiveJetPt = nullptr; }
        if (hRefSelInclusiveJetEta) { delete hRefSelInclusiveJetEta; hRefSelInclusiveJetEta = nullptr; }
        if (hRefSelInclusiveJetEtaUnweighted) { delete hRefSelInclusiveJetEtaUnweighted; hRefSelInclusiveJetEtaUnweighted = nullptr; }
        if (hRefSelInclusiveJetPtEta) { delete hRefSelInclusiveJetPtEta; hRefSelInclusiveJetPtEta = nullptr; }
        if (hRefSelInclusiveJetPtEtaPtHat) { delete hRefSelInclusiveJetPtEtaPtHat; hRefSelInclusiveJetPtEtaPtHat = nullptr; }
        if (hRefSelLeadJetPtEta) { delete hRefSelLeadJetPtEta; hRefSelLeadJetPtEta = nullptr; }
        if (hRefSelLeadJetPtEtaPtHat) { delete hRefSelLeadJetPtEtaPtHat; hRefSelLeadJetPtEtaPtHat = nullptr; }
        if (hRefSelSubLeadJetPtEta) { delete hRefSelSubLeadJetPtEta; hRefSelSubLeadJetPtEta = nullptr; }
        if (hRefSelSubLeadJetPtEtaPtHat) { delete hRefSelSubLeadJetPtEtaPtHat; hRefSelSubLeadJetPtEtaPtHat = nullptr; }
    } // if ( fIsMc )

}

//________________
void HistoManagerJetESR::init() {

    int    vzBins{320};
    double vzRange[2] = {-31., 31.};
    // int    multBins{1800};
    // double multRange[2] {-0.5, 1799.5};
    int    hiBinBins{203};
    double hiBinRange[2] = {-1.5, 201.5};
    // int    weightBins{110};
    // double weightRange[2] = {-0.05, 1.05};
    int    ptHatBins{100};
    double ptHatRange[2] = {0., 1000.};
    int    fMultBins{32}; 
    double fMultRange[2] = {-0.5, 31.5};
    int    ptHatWeightBins{100};
    double ptHatWeightRange[2] = {0., 10.};
    int    fFracBins{100};
    double fFracRange[2] = {0., 1.};

    int    prescale = 1;

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

    // JES and JER bins [4 dimensions]
    // 0 - JES, 1 - pt, 2 - eta, 3 - ptHat
    int    bins4D_jet_JESPtEtaPtHat[4] { fJESBins, fPtBins, fEtaBins, fPtHatBins };
    double xmin4D_jet_JESPtEtaPtHat[4] { fJESRange[0], fPtRange[0], fEtaRange[0], fPtHatRange[0] };
    double xmax4D_jet_JESPtEtaPtHat[4] { fJESRange[1], fPtRange[1], fEtaRange[1], fPtHatRange[1] };

    // Reco 2 Ref correlations for inclusive jets [7 dimensions]
    // 0 - reco ptCorr, 1- reco eta, 2 - reco phi,
    // 3 - ref pt, 4 - ref eta, 5 - ref phi, 6 - ptHat
    int    bins7D_reco2ref_jet_PtEtaPhiPtHat[7] { fPtBins, fEtaBins, fPhiBins, fPtBins, fEtaBins, fPhiBins, fPtHatBins };
    double xmin7D_reco2ref_jet_PtEtaPhiPtHat[7] { fPtRange[0], fEtaRange[0], fPhiRange[0], fPtRange[0], fEtaRange[0], fPhiRange[0], fPtHatRange[0] };
    double xmax7D_reco2ref_jet_PtEtaPhiPtHat[7] { fPtRange[1], fEtaRange[1], fPhiRange[1], fPtRange[1], fEtaRange[1], fPhiRange[1], fPtHatRange[1] };
    

    //
    // Event histograms
    //

    hVz = new TH1D("hVz","Vertex z position;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVz->Sumw2();
    hVzPtHatWeighted = new TH1D("hVzPtHatWeighted","Vertex z position with #hat{p_{T}} weight;vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzPtHatWeighted->Sumw2();
    hVzCentWeighted = new TH1D("hVzCentWeighted", "Vertex z position (centrality weighted);vz (cm);Entries", vzBins, vzRange[0], vzRange[1]);
    hVzCentWeighted->Sumw2();
    hVzWeighted = new TH1D("hVzWeighted","Vertex z position;vz (cm);Entries", 400, -50., 50.);
    hVzWeighted->Sumw2();
    hPtHat = new TH1D("hPtHat","#hat{p_{T}};#hat{p_{T}} (GeV);Entries", 
                      ptHatBins, ptHatRange[0], ptHatRange[1]);
    hPtHat->Sumw2();
    hPtHatWeighted = new TH1D("hPtHatWeighted","#hat{p_{T}} with #hat{p_{T}} weight;#hat{p_{T}} (GeV);Entries", 
                              ptHatBins, ptHatRange[0], ptHatRange[1]);
    hPtHatWeighted->Sumw2();
    hHiBin = new TH1D("hHiBin","HiBin a.k.a. centrality;HiBin;Entries", 
                      hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBin->Sumw2();
    hHiBinWeighted = new TH1D("hHiBinWeighted","HiBin a.k.a. centrality with #hat{p_{T}} weight;HiBin;Entries", 
                              hiBinBins, hiBinRange[0], hiBinRange[1]);
    hHiBinWeighted->Sumw2();
    hPtHatWeight = new TH1D("hPtHatWeight","#hat{p_{T}} weight;#hat{p_{T}}^{weight};Entries", ptHatWeightBins, ptHatWeightRange[0], ptHatWeightRange[1]);
    hPtHatWeight->Sumw2();
    hPtHatWeightWeighted = new TH1D("hPtHatWeightWeighted","#hat{p_{T}} weight (#hat{p_{T}} weighted);#hat{p_{T}}^{weight};Entries", ptHatWeightBins, ptHatWeightRange[0], ptHatWeightRange[1]);
    hPtHatWeightWeighted->Sumw2();

    hCentrality = new TH1D("hCentrality","Collision centrality;Centrality (%);Entries", fCentBins, fCentRange[0], fCentRange[1]);
    hCentrality->Sumw2();
    hCentralityWeighted = new TH1D("hCentralityWeighted","Collision centrality (weighted);Centrality (%);Entries", fCentBins, fCentRange[0], fCentRange[1]);
    hCentralityWeighted->Sumw2();

    hVzPtHatCent = new TH3D("hVzPtHatCent","Vertex and centrality;vz (cm);#hat{p_{T}} (GeV);centrality (%)",
                            vzBins, vzRange[0], vzRange[1], fPtHatBins, fPtHatRange[0], fPtHatRange[1], fCentBins, fCentRange[0], fCentRange[1]);
    hVzPtHatCent->Sumw2();
    hVzPtHatCentWeighted = new TH3D("hVzPtHatCentWeighted","Vertex, #hat{p_{T}} and centrality (weighted);vz (cm);#hat{p_{T}} (GeV);centrality (%)",
                                    vzBins, vzRange[0], vzRange[1], fPtHatBins, fPtHatRange[0], fPtHatRange[1], fCentBins, fCentRange[0], fCentRange[1]);
    hVzPtHatCentWeighted->Sumw2();


    //
    // Reco jet histograms
    //

    // Inclusive jet histograms
    hRecoJetCollectionSize = new TH1D("hRecoJetCollectionSize","Reco jet collection size;Number of jets;Entries", 100, -0.5, 99.5);
    hRecoJetCollectionSize->Sumw2();

    // For 4 eta ranges: <=2.4, <=2.7, <=3, >3
    for (int i{0}; i<4; i++) {
        float low{0}, hi{0};
        if      ( i == 0 ) { low = {-2.4f}; hi = {2.4f}; }
        else if ( i == 1 ) { low = {2.4f}; hi = {2.7f}; }
        else if ( i == 2 ) { low = {2.7f}; hi = {3.0f}; }
        else if ( i == 3 ) { low = {3.f}; hi = {100.0f}; }
        hRecoInclusiveJetNHF[i] = new TH1D(Form("hRecoInclusiveJetNHF_%d",i), Form("Neutral hadron fraction for %3.1f< #eta<=%3.1f;Neutral hadron fraction;1/N dN/dNHF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetNHF[i]->Sumw2();
        hRecoInclusiveJetNEmF[i] = new TH1D(Form("hRecoInclusiveJetNEmF_%d",i), Form("Neutral EM fraction for %3.1f< #eta<=%3.1f;Neutral EM fraction;1/N dN/dNEF", low, hi), 
                            fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetNEmF[i]->Sumw2();
        hRecoInclusiveJetNumOfConst[i] = new TH1D(Form("hRecoInclusiveJetNumOfConst_%d",i), Form("Number of constituents for %3.1f< #eta<=%3.1f;Number of constituents;1/N dN/dNconst", low, hi), 
                                  fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetNumOfConst[i]->Sumw2();
        hRecoInclusiveJetMUF[i] = new TH1D(Form("hRecoInclusiveJetMUF_%d",i), Form("Muon fraction for %3.1f< #eta<=%3.1f;Muon fraction;1/N dN/dMUF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetMUF[i]->Sumw2();
        hRecoInclusiveJetCHF[i] = new TH1D(Form("hRecoInclusiveJetCHF_%d",i), Form("Charged hadron fraction for %3.1f< #eta<=%3.1f;Charged hadron fraction;1/N dN/dCHF", low, hi), 
                           fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetCHF[i]->Sumw2();
        hRecoInclusiveJetChargedMult[i] = new TH1D(Form("hRecoInclusiveJetChargedMult_%d",i), Form("Charged particle multiplicity for %3.1f< #eta<=%3.1f;Charged multiplicity;1/N dN/dChMult", low, hi), 
                                   fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetChargedMult[i]->Sumw2();
        hRecoInclusiveJetCEmF[i] = new TH1D(Form("hRecoInclusiveJetCEmF_%d",i), Form("Charged EM fraction for %3.1f< #eta<=%3.1f;Charged EM fraction;1/N dN/dChEmF", low, hi), 
                            fFracBins, fFracRange[0], fFracRange[1]);
        hRecoInclusiveJetCEmF[i]->Sumw2();
        hRecoInclusiveJetNumOfNeutPart[i] = new TH1D(Form("hRecoInclusiveJetNumOfNeutPart_%d",i), Form("Number of neutral particles for %3.1f< #eta<=%3.1f;Number of neutral particles;1/N dN/dNumOfNeutrals", low, hi), 
                                     fMultBins, fMultRange[0], fMultRange[1]);
        hRecoInclusiveJetNumOfNeutPart[i]->Sumw2();
    } // for (int i{0}; i<4; i++)



    hRecoLeadAllJetPtEta = new TH2D("hRecoLeadAllJetPtEta","Leading jet all p_{T} vs #eta;#eta;p_{T} (GeV)", 
                                        fEtaBins, fEtaRange[0], fEtaRange[1], 
                                        fPtBins, fPtRange[0], fPtRange[1]);
    hRecoLeadAllJetPtEta->Sumw2();

    hRecoLeadAllJetPtEtaPtHat = new TH3D("hRecoLeadAllJetPtEtaPtHat","Leading jet (matched+unmatched) p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                         prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                         fPtBins, fPtRange[0], fPtRange[1],
                                         fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hRecoLeadAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoLeadAllJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
    hRecoLeadAllJetPtEtaPtHat->Sumw2();

    hRecoInclusiveAllJetPtRawEta = new TH2D("hRecoInclusiveAllJetPtRawEta","Reco jet raw p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtRawEta->Sumw2();

    hRecoSubLeadAllJetPtEta = new TH2D("hRecoSubLeadAllJetPtEta","Subleading jet all p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1], 
                                            fPtBins, fPtRange[0], fPtRange[1]);
    hRecoSubLeadAllJetPtEta->Sumw2();

    hRecoSubLeadAllJetPtEtaPtHat = new TH3D("hRecoSubLeadAllJetPtEtaPtHat","Subleading jet (matched+unmatched) p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hRecoSubLeadAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoSubLeadAllJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
    hRecoSubLeadAllJetPtEtaPtHat->Sumw2();


    hRecoInclusiveAllJetPt = new TH1D("hRecoInclusiveAllJetPt","Reco jet p_{T} (matched+unmatched);Reco p_{T} (GeV);Entries",
                                   fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPt->Sumw2();
    hRecoInclusiveAllJetEta = new TH1D("hRecoInclusiveAllJetEta","Reco jet #eta (matched+unmatched);#eta;Entries",
                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoInclusiveAllJetEta->Sumw2();
    hRecoInclusiveAllJetEtaUnweighted = new TH1D("hRecoInclusiveAllJetEtaUnweighted","Reco jet #eta (matched+unmatched) unweighted;#eta;Entries",
                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoInclusiveAllJetEtaUnweighted->Sumw2();
    hRecoInclusiveAllJetPtEta = new TH2D("hRecoInclusiveAllJetPtEta","Reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                      fEtaBins, fEtaRange[0], fEtaRange[1],
                                      fPtBins, fPtRange[0], fPtRange[1]);
    hRecoInclusiveAllJetPtEta->Sumw2();
    hRecoInclusiveAllJetPtEtaPtHat = new TH3D("hRecoInclusiveAllJetPtEtaPtHat","Reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                           prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                           fPtBins, fPtRange[0], fPtRange[1],
                                           fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
    hRecoInclusiveAllJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
    hRecoInclusiveAllJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
    hRecoInclusiveAllJetPtEtaPtHat->Sumw2();


    hRecoGoodInclusiveJetEtaLabFrame = new TH1D("hRecoGoodInclusiveJetEtaLabFrame","Reco good jet #eta in lab frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaLabFrame->Sumw2();
    hRecoGoodInclusiveJetEtaCMFrame = new TH1D("hRecoGoodInclusiveJetEtaCMFrame","Reco good jet #eta in CM frame;#eta;Entries",
                                                fEtaBins, fEtaRange[0], fEtaRange[1]);
    hRecoGoodInclusiveJetEtaCMFrame->Sumw2();


    //
    // Monte Carlo histograms
    //

    if ( fIsMc ) {

        //
        // Reco (matched and unmatched) histograms
        //

        hRecoLeadingJetPtOverPtHatVsLeadingJetPt = new TH2D("hRecoLeadingJetPtOverPtHatVsLeadingJetPt", "Reco leading jet p_{T}/#hat{p}_{T} vs reco leading jet p_{T}^{reco};Reco leading jet p_{T} (GeV);Reco leading jet p_{T}/#hat{p}_{T}",
                                                            fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Sumw2();
        hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted = new TH2D("hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted", "Reco leading jet p_{T}/#hat{p}_{T} vs reco leading jet p_{T}^{reco} weighted;Reco leading jet p_{T} (GeV);Reco leading jet p_{T}/#hat{p}_{T}",
                                                                    fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Sumw2();
        hRecoDijetPtOverPtHatVsDijetPt = new TH2D("hRecoDijetPtOverPtHatVsDijetPt", "Reco dijet p_{T}/#hat{p}_{T} vs dijet p_{T};Reco dijet p_{T} (GeV);Reco dijet p_{T}/#hat{p}_{T}",
                                                  fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtOverPtHatVsDijetPt->Sumw2();
        hRecoDijetPtOverPtHatVsDijetPtWeighted = new TH2D("hRecoDijetPtOverPtHatVsDijetPtWeighted", "Reco dijet p_{T}/#hat{p}_{T} vs reco dijet p_{T} weighted;Reco dijet p_{T}^{reco} (GeV);Reco dijet p_{T}^{reco}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtOverPtHatVsDijetPtWeighted->Sumw2();
        hRecoDijetPtAveOverPtHatVsDijetPtAve = new TH2D("hRecoDijetPtAveOverPtHatVsDijetPtAve", "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs reco dijet p_{T}^{ave};Reco dijet p_{T}^{ave} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtAveOverPtHatVsDijetPtAve->Sumw2();
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted = new TH2D("hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted", "Reco dijet p_{T}^{ave}/#hat{p}_{T} vs reco dijet p_{T}^{ave} weighted;Reco dijet p_{T}^{ave} (GeV);Reco dijet p_{T}^{ave}/#hat{p}_{T}",
                                                                fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Sumw2();


        hRecoInclusiveMatchedJetPt = new TH1D("hRecoInclusiveMatchedJetPt","Matched reco jet p_{T} (matched);Reco p_{T} (GeV);Entries",
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveMatchedJetPt->Sumw2();
        hRecoInclusiveMatchedJetPtEta = new TH2D("hRecoInclusiveMatchedJetPtEta","Matched reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveMatchedJetPtEta->Sumw2();
        hRecoInclusiveMatchedJetPtEtaPtHat = new TH3D("hRecoInclusiveMatchedJetPtEtaPtHat","Matched reco jet p_{T} vs #eta vs #hat{p_{T}};#eta;p_{T} (GeV);#hat{p_{T}} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRecoInclusiveMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveMatchedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRecoInclusiveMatchedJetPtEtaPtHat->Sumw2();
        
        hRecoInclusiveUnmatchedJetPt = new TH1D("hRecoInclusiveUnmatchedJetPt","Unmatched reco jet p_{T} (unmatched);Reco p_{T} (GeV);Entries",
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveUnmatchedJetPt->Sumw2();
        hRecoInclusiveUnmatchedJetPtEta = new TH2D("hRecoInclusiveUnmatchedJetPtEta","Unmatched reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoInclusiveUnmatchedJetPtEta->Sumw2();
        hRecoInclusiveUnmatchedJetPtEtaPtHat = new TH3D("hRecoInclusiveUnmatchedJetPtEtaPtHat","Unmatched reco jet p_{T} vs #eta vs #hat{p}_{T}};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRecoInclusiveUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoInclusiveUnmatchedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRecoInclusiveUnmatchedJetPtEtaPtHat->Sumw2();

        hRecoLeadMatchedJetPtEta = new TH2D("hRecoLeadMatchedJetPtEta","Leading matched reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoLeadMatchedJetPtEta->Sumw2();
        hRecoLeadMatchedJetPtEtaPtHat = new TH3D("hRecoLeadMatchedJetPtEtaPtHat","Leading matched reco jet p_{T} vs #eta vs #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRecoLeadMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoLeadMatchedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRecoLeadMatchedJetPtEtaPtHat->Sumw2();

        hRecoLeadUnmatchedJetPtEta = new TH2D("hRecoLeadUnmatchedJetPtEta","Leading unmatched reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoLeadUnmatchedJetPtEta->Sumw2();
        hRecoLeadUnmatchedJetPtEtaPtHat = new TH3D("hRecoLeadUnmatchedJetPtEtaPtHat","Leading unmatched reco jet p_{T} vs #eta vs #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRecoLeadUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoLeadUnmatchedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRecoLeadUnmatchedJetPtEtaPtHat->Sumw2();

        hRecoSubLeadMatchedJetPtEta = new TH2D("hRecoSubLeadMatchedJetPtEta","Subleading matched reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoSubLeadMatchedJetPtEta->Sumw2();
        hRecoSubLeadMatchedJetPtEtaPtHat = new TH3D("hRecoSubLeadMatchedJetPtEtaPtHat","Subleading matched reco jet p_{T} vs #eta vs #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRecoSubLeadMatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoSubLeadMatchedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRecoSubLeadMatchedJetPtEtaPtHat->Sumw2();
        hRecoSubLeadUnmatchedJetPtEta = new TH2D("hRecoSubLeadUnmatchedJetPtEta","Subleading unmatched reco jet p_{T} vs #eta;#eta;p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRecoSubLeadUnmatchedJetPtEta->Sumw2();
        hRecoSubLeadUnmatchedJetPtEtaPtHat = new TH3D("hRecoSubLeadUnmatchedJetPtEtaPtHat","Subleading unmatched reco jet p_{T} vs #eta vs #hat{p}_{T};#eta;p_{T} (GeV);#hat{p}_{T} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRecoSubLeadUnmatchedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRecoSubLeadUnmatchedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRecoSubLeadUnmatchedJetPtEtaPtHat->Sumw2();

        // JES and JER for inclusive, leading, and subleading jets

        // pt corr/pt gen vs pt gen at midrapidity -1.4 < eta < 1.4 
        hInclusiveJetJESVsPtGen = new TH2D("hInclusiveJetJESVsPtGen","Inclusive jet JES vs p_{T}^{Gen};p_{T}^{Gen} (GeV);p_{T}^{corr}/p_{T}^{Gen}",
                                                        fPtBins, fPtRange[0], fPtRange[1], 
                                                        fJESBins, fJESRange[0], fJESRange[1]);
        hInclusiveJetJESVsPtGen->Sumw2();
        // pt corr/pt, gen pt, gen eta, ptHat
        hInclusiveJetJESGenPtGenEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESGenPtGenEtaPtHatWeighted","Inclusive jet JES vs p_{T}^{Gen}, #eta, #hat{p}_{T};p_{T}^{Gen} (GeV);#eta;#hat{p}_{T} (GeV)",
                                                        4, bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        // pt corr/pt, reco pt, reco eta, ptHat
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted = new THnSparseD("hInclusiveJetJESRecoPtRecoEtaPtHatWeighted","Inclusive jet JES vs p_{T}^{Reco}, #eta, #hat{p}_{T};p_{T}^{Reco} (GeV);#eta;#hat{p}_{T} (GeV)",
                                                        4, bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        // pt corr/pt gen vs pt gen, gen eta, ptHat
        hLeadingJetJESGenPtEtaPtHatWeighted = new THnSparseD("hLeadingJetJESGenPtEtaPtHatWeighted","Leading jet JES vs p_{T}^{Gen}, #eta, #hat{p}_{T};p_{T}^{Gen} (GeV);#eta;#hat{p}_{T} (GeV)",
                                                        4, bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);
        // pt corr/pt gen vs pt gen, gen eta, ptHat
        hSubleadingJetJESGenPtEtaPtHatWeighted = new THnSparseD("hSubleadingJetJESGenPtEtaPtHatWeighted","Subleading jet JES vs p_{T}^{Gen}, #eta, #hat{p}_{T};p_{T}^{Gen} (GeV);#eta;#hat{p}_{T} (GeV)",
                                                        4, bins4D_jet_JESPtEtaPtHat, xmin4D_jet_JESPtEtaPtHat, xmax4D_jet_JESPtEtaPtHat);

        // Reco 2 Ref correlations for inclusive jets [7 dimensions]
        // 0 - reco ptCorr, 1- reco eta, 2 - reco phi,
        // 3 - ref pt, 4 - ref eta, 5 - ref phi, 6 - ptHat
        hInclusiveReco2RefJetPtEtaPhiPtHat = new THnSparseD("hInclusiveReco2RefJetPtEtaPhiPtHat","Inclusive reco 2 ref jet p_{T} vs #eta vs #phi vs #hat{p}_{T};Reco p_{T}^{corr} (GeV);Reco #eta;Reco #phi (rad);Ref p_{T} (GeV);Ref #eta;Ref #phi (rad);#hat{p}_{T} (GeV)",
                                                        7, bins7D_reco2ref_jet_PtEtaPhiPtHat, xmin7D_reco2ref_jet_PtEtaPhiPtHat, xmax7D_reco2ref_jet_PtEtaPhiPtHat);
        hLeadReco2RefJetPtEtaPhiPtHat = new THnSparseD("hLeadReco2RefJetPtEtaPhiPtHat","Leading jet reco 2 ref jet p_{T} vs #eta vs #phi vs #hat{p}_{T};Reco p_{T}^{corr} (GeV);Reco #eta;Reco #phi (rad);Ref p_{T} (GeV);Ref #eta;Ref #phi (rad);#hat{p}_{T} (GeV)",
                                                        7, bins7D_reco2ref_jet_PtEtaPhiPtHat, xmin7D_reco2ref_jet_PtEtaPhiPtHat, xmax7D_reco2ref_jet_PtEtaPhiPtHat);
        hSubLeadReco2RefJetPtEtaPhiPtHat = new THnSparseD("hSubLeadReco2RefJetPtEtaPhiPtHat","Subleading jet reco 2 ref jet p_{T} vs #eta vs #phi vs #hat{p}_{T};Reco p_{T}^{corr} (GeV);Reco #eta;Reco #phi (rad);Ref p_{T} (GeV);Ref #eta;Ref #phi (rad);#hat{p}_{T} (GeV)",
                                                        7, bins7D_reco2ref_jet_PtEtaPhiPtHat, xmin7D_reco2ref_jet_PtEtaPhiPtHat, xmax7D_reco2ref_jet_PtEtaPhiPtHat);

        //
        // Gen inclusive jets
        //

        hGenJetCollectionSize = new TH1D("hGenJetCollectionSize","Gen jet collection size;Gen jet collection size;Entries", 100, -0.5, 99.5);
        hGenJetCollectionSize->Sumw2();

        // Histograms to check event overweighing
        hGenVsRecoJetCollectionSize = new TH2D("hGenVsRecoJetCollectionSize","Reco vs Gen jet collection size;Reco jet collection size;Gen jet collection size", 100, -0.5, 99.5, 100, -0.5, 99.5);
        hGenVsRecoJetCollectionSize->Sumw2();

        hGenLeadingJetPtOverPtHatVsLeadingJetPt = new TH2D("hGenLeadingJetPtOverPtHatVsLeadingJetPt", "Leading jet p_{T}^{Gen}/#hat{p}_{T} vs leading jet p_{T}^{Gen};p_{T}^{Gen} (GeV);p_{T}^{Gen}/#hat{p}_{T}",
                                                            fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenLeadingJetPtOverPtHatVsLeadingJetPt->Sumw2();
        hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted = new TH2D("hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted", "Leading jet p_{T}^{Gen}/#hat{p}_{T} vs leading jet p_{T}^{Gen} weighted;p_{T}^{Gen} (GeV);p_{T}^{Gen}/#hat{p}_{T}",
                                                                    fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Sumw2();
        hGenDijetPtOverPtHatVsDijetPt = new TH2D("hGenDijetPtOverPtHatVsDijetPt", "Gen dijet p_{T}/#hat{p}_{T} vs gen dijet p_{T};Gen dijet p_{T} (GeV);Gen dijet p_{T}/#hat{p}_{T}",
                                                  fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtOverPtHatVsDijetPt->Sumw2();
        hGenDijetPtOverPtHatVsDijetPtWeighted = new TH2D("hGenDijetPtOverPtHatVsDijetPtWeighted", "Gen dijet p_{T}/#hat{p}_{T} vs gen dijet p_{T} weighted;Gen dijet p_{T} (GeV);Gen dijet p_{T}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtOverPtHatVsDijetPtWeighted->Sumw2();
        hGenDijetPtAveOverPtHatVsDijetPtAve = new TH2D("hGenDijetPtAveOverPtHatVsDijetPtAve", "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs gen dijet p_{T}^{ave};Gen dijet p_{T}^{ave} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}",
                                                          fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtAveOverPtHatVsDijetPtAve->Sumw2();
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted = new TH2D("hGenDijetPtAveOverPtHatVsDijetPtAveWeighted", "Gen dijet p_{T}^{ave}/#hat{p}_{T} vs gen dijet p_{T}^{ave} weighted;Gen dijet p_{T}^{ave} (GeV);Gen dijet p_{T}^{ave}/#hat{p}_{T}",
                                                                fPtBins, fPtRange[0], fPtRange[1], 350, 0., 3.5);
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Sumw2();   
        
        // Single-jet distributions (inclusive, leading, subleading)
        hGenInclusiveJetPt = new TH1D("hGenInclusiveJetPt","Inclusive gen jet;Gen p_{T}^{inclusive} (GeV)",
                                    fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPt->Sumw2();
        hGenInclusiveJetEta = new TH1D("hGenInclusiveJetEta","Inclusive gen jet;Gen #eta^{inclusive}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenInclusiveJetEta->Sumw2();
        hGenInclusiveJetEtaUnweighted = new TH1D("hGenInclusiveJetEtaUnweighted","Inclusive gen jet unweighted;Gen #eta^{inclusive}",
                                        fEtaBins, fEtaRange[0], fEtaRange[1] );
        hGenInclusiveJetEtaUnweighted->Sumw2();
        hGenInclusiveJetPtEta = new TH2D("hGenInclusiveJetPtEta","Gen inclusive jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenInclusiveJetPtEta->Sumw2();
        hGenInclusiveJetPtEtaPtHat = new TH3D("hGenInclusiveJetPtEtaPtHat","Gen inclusive jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        hGenInclusiveJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hGenInclusiveJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hGenInclusiveJetPtEtaPtHat->Sumw2();
        hGenLeadJetPtEta = new TH2D("hGenLeadJetPtEta","Gen leading jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenLeadJetPtEta->Sumw2();
        hGenLeadJetPtEtaPtHat = new TH3D("hGenLeadJetPtEtaPtHat","Gen leading jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        hGenLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hGenLeadJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hGenLeadJetPtEtaPtHat->Sumw2();
        hGenSubLeadJetPtEta = new TH2D("hGenSubLeadJetPtEta","Gen subleading jet acceptance;Gen #eta;Gen p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1] );
        hGenSubLeadJetPtEta->Sumw2();
        hGenSubLeadJetPtEtaPtHat = new TH3D("hGenSubLeadJetPtEtaPtHat","Gen subleading jet acceptance vs pT hat;Gen #eta;Gen p_{T} (GeV);p_{T}^{hat} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1] );
        hGenSubLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hGenSubLeadJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hGenSubLeadJetPtEtaPtHat->Sumw2();

        hGenGoodInclusiveJetEtaLabFrame = new TH1D("hGenGoodInclusiveJetEtaLabFrame","Gen good inclusive jet #eta in lab frame;#eta^{Inclusive}",
                                                   fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenGoodInclusiveJetEtaLabFrame->Sumw2();
        hGenGoodInclusiveJetEtaCMFrame = new TH1D("hGenGoodInclusiveJetEtaCMFrame","Gen good inclusive jet #eta in CM frame;#eta^{Inclusive}_{CM}",
                                                  fEtaBins, fEtaRange[0], fEtaRange[1]);
        hGenGoodInclusiveJetEtaCMFrame->Sumw2();


        //
        // Ref jet histograms
        //

        // Single-jet distributions (inclusive, leading, subleading)
        hRefInclusiveJetPt = new TH1D("hRefInclusiveJetPt","Ref jet p_{T};Ref p_{T} (GeV);Entries",
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPt->Sumw2();
        hRefInclusiveJetEta = new TH1D("hRefInclusiveJetEta","Ref jet #eta;Ref #eta;Entries",
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefInclusiveJetEta->Sumw2();
        hRefInclusiveJetEtaUnweighted = new TH1D("hRefInclusiveJetEtaUnweighted","Ref jet #eta (unweighted);Ref #eta;Entries",
                                        fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefInclusiveJetEtaUnweighted->Sumw2();
        hRefInclusiveJetPtEta = new TH2D("hRefInclusiveJetPtEta","Ref jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefInclusiveJetPtEta->Sumw2();
        hRefInclusiveJetPtEtaPtHat = new TH3D("hRefInclusiveJetPtEtaPtHat","Ref jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefInclusiveJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);        
        hRefInclusiveJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefInclusiveJetPtEtaPtHat->Sumw2();
        hRefLeadJetPtEta = new TH2D("hRefLeadJetPtEta","Ref leading jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefLeadJetPtEta->Sumw2();
        hRefLeadJetPtEtaPtHat = new TH3D("hRefLeadJetPtEtaPtHat","Ref leading jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefLeadJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefLeadJetPtEtaPtHat->Sumw2();
        hRefLeadUnswappedJetPtEta = new TH2D("hRefLeadUnswappedJetPtEta","Ref leading unswapped jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefLeadUnswappedJetPtEta->Sumw2();
        hRefLeadUnswappedJetPtEtaPtHat = new TH3D("hRefLeadUnswappedJetPtEtaPtHat","Ref leading unswapped jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefLeadUnswappedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefLeadUnswappedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefLeadUnswappedJetPtEtaPtHat->Sumw2();
        hRefSubLeadJetPtEta = new TH2D("hRefSubLeadJetPtEta","Ref subleading jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSubLeadJetPtEta->Sumw2();
        hRefSubLeadJetPtEtaPtHat = new TH3D("hRefSubLeadJetPtEtaPtHat","Ref subleading jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefSubLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefSubLeadJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefSubLeadJetPtEtaPtHat->Sumw2();
        hRefSubLeadUnswappedJetPtEta = new TH2D("hRefSubLeadJetUnswappedPtEta","Ref subleading unswapped jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSubLeadUnswappedJetPtEta->Sumw2();
        hRefSubLeadUnswappedJetPtEtaPtHat = new TH3D("hRefSubLeadJetUnswappedPtEtaPtHat","Ref subleading unswapped jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefSubLeadUnswappedJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefSubLeadUnswappedJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefSubLeadUnswappedJetPtEtaPtHat->Sumw2();

        //
        // Ref-selected jet histograms
        //

        hRefSelInclusiveJetPt = new TH1D("hRefSelInclusiveJetPt","Ref selected inclusive jet p_{T};Ref p_{T} (GeV);Entries",
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelInclusiveJetPt->Sumw2();                                         
        hRefSelInclusiveJetEta = new TH1D("hRefSelInclusiveJetEta","Ref selected inclusive jet #eta;Ref #eta;Entries",
                                            fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefSelInclusiveJetEta->Sumw2();
        hRefSelInclusiveJetEtaUnweighted = new TH1D("hRefSelInclusiveJetEtaUnweighted","Ref selected inclusive jet #eta (unweighted);Ref #eta;Entries",
                                            fEtaBins, fEtaRange[0], fEtaRange[1]);
        hRefSelInclusiveJetEtaUnweighted->Sumw2();
        hRefSelInclusiveJetPtEta = new TH2D("hRefSelInclusiveJetPtEta","Ref selected inclusive jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                            fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelInclusiveJetPtEta->Sumw2();
        hRefSelInclusiveJetPtEtaPtHat = new TH3D("hRefSelInclusiveJetPtEtaPtHat","Ref selected inclusive jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                            prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                            fPtBins, fPtRange[0], fPtRange[1],
                                            fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefSelInclusiveJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefSelInclusiveJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefSelInclusiveJetPtEtaPtHat->Sumw2();
        hRefSelLeadJetPtEta = new TH2D("hRefSelLeadJetPtEta","Ref selected leading jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelLeadJetPtEta->Sumw2();
        hRefSelLeadJetPtEtaPtHat = new TH3D("hRefSelLeadJetPtEtaPtHat","Ref selected leading jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefSelLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefSelLeadJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefSelLeadJetPtEtaPtHat->Sumw2();
        hRefSelSubLeadJetPtEta = new TH2D("hRefSelSubLeadJetPtEta","Ref selected subleading jet p_{T} vs #eta;Ref #eta;Ref p_{T} (GeV)",
                                        fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1]);
        hRefSelSubLeadJetPtEta->Sumw2();
        hRefSelSubLeadJetPtEtaPtHat = new TH3D("hRefSelSubLeadJetPtEtaPtHat","Ref selected subleading jet p_{T} vs #eta vs #hat{p}_{T};Ref #eta;Ref p_{T} (GeV);#hat{p}_{T} (GeV)",
                                        prescale * fEtaBins, fEtaRange[0], fEtaRange[1],
                                        fPtBins, fPtRange[0], fPtRange[1],
                                        fPtHatBins, fPtHatRange[0], fPtHatRange[1]);
        hRefSelSubLeadJetPtEtaPtHat->GetXaxis()->Set(jetEtaL2L3StdBins, jetEtaL2L3StdVals);
        hRefSelSubLeadJetPtEtaPtHat->SetBinsLength(-1); // Recalculate fNcells
        hRefSelSubLeadJetPtEtaPtHat->Sumw2();

    } // if ( fIsMc )
}

//________________
void HistoManagerJetESR::writeOutput() {

    //
    // Event information
    //
    hVz->Write();
    hVzCentWeighted->Write();
    hVzPtHatWeighted->Write();
    hVzWeighted->Write();
    hPtHat->Write();
    hPtHatWeighted->Write();
    hHiBin->Write();
    hHiBinWeighted->Write();
    hPtHatWeight->Write();
    hPtHatWeightWeighted->Write();
    hCentrality->Write();
    hCentralityWeighted->Write();
    hVzPtHatCent->Write();
    hVzPtHatCentWeighted->Write();

    //
    // Reco jets
    //

    hRecoJetCollectionSize->Write();

    for (int i=0; i<4; i++) {
        hRecoInclusiveJetNHF[i]->Write();
        hRecoInclusiveJetNEmF[i]->Write();
        hRecoInclusiveJetNumOfConst[i]->Write();
        hRecoInclusiveJetMUF[i]->Write();
        hRecoInclusiveJetCHF[i]->Write();
        hRecoInclusiveJetChargedMult[i]->Write();
        hRecoInclusiveJetCEmF[i]->Write();
        hRecoInclusiveJetNumOfNeutPart[i]->Write();
    }

    hRecoInclusiveAllJetPtRawEta->Write();

    hRecoInclusiveAllJetPt->Write();
    hRecoInclusiveAllJetEta->Write();
    hRecoInclusiveAllJetEtaUnweighted->Write();
    hRecoInclusiveAllJetPtEta->Write();
    hRecoInclusiveAllJetPtEtaPtHat->Write();
    hRecoLeadAllJetPtEta->Write();
    hRecoLeadAllJetPtEtaPtHat->Write(); 
    hRecoSubLeadAllJetPtEta->Write();
    hRecoSubLeadAllJetPtEtaPtHat->Write();
    hRecoGoodInclusiveJetEtaLabFrame->Write();
    hRecoGoodInclusiveJetEtaCMFrame->Write();

    //
    // Monte Carlo information
    //

    if ( fIsMc ) {

        //
        // Reco (matched+unmatched) jets
        //

        hRecoLeadingJetPtOverPtHatVsLeadingJetPt->Write();
        hRecoLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Write();
        hRecoDijetPtOverPtHatVsDijetPt->Write();
        hRecoDijetPtOverPtHatVsDijetPtWeighted->Write();
        hRecoDijetPtAveOverPtHatVsDijetPtAve->Write();
        hRecoDijetPtAveOverPtHatVsDijetPtAveWeighted->Write();

        hRecoInclusiveMatchedJetPt->Write();
        hRecoInclusiveMatchedJetPtEta->Write();
        hRecoInclusiveMatchedJetPtEtaPtHat->Write();
        hRecoInclusiveUnmatchedJetPt->Write();
        hRecoInclusiveUnmatchedJetPtEta->Write();
        hRecoInclusiveUnmatchedJetPtEtaPtHat->Write();
        hRecoLeadMatchedJetPtEta->Write();
        hRecoLeadMatchedJetPtEtaPtHat->Write();
        hRecoLeadUnmatchedJetPtEta->Write();
        hRecoLeadUnmatchedJetPtEtaPtHat->Write();
        hRecoSubLeadMatchedJetPtEta->Write();
        hRecoSubLeadMatchedJetPtEtaPtHat->Write();
        hRecoSubLeadUnmatchedJetPtEta->Write();
        hRecoSubLeadUnmatchedJetPtEtaPtHat->Write();

        hInclusiveJetJESVsPtGen->Write();
        hInclusiveJetJESGenPtGenEtaPtHatWeighted->Write();
        hInclusiveJetJESRecoPtRecoEtaPtHatWeighted->Write();
        hLeadingJetJESGenPtEtaPtHatWeighted->Write();
        hSubleadingJetJESGenPtEtaPtHatWeighted->Write();
        hInclusiveReco2RefJetPtEtaPhiPtHat->Write();
        hLeadReco2RefJetPtEtaPhiPtHat->Write();
        hSubLeadReco2RefJetPtEtaPhiPtHat->Write();

        //
        // Gen jets
        //

        hGenJetCollectionSize->Write();
        hGenVsRecoJetCollectionSize->Write();

        hGenLeadingJetPtOverPtHatVsLeadingJetPt->Write();
        hGenLeadingJetPtOverPtHatVsLeadingJetPtWeighted->Write();
        hGenDijetPtOverPtHatVsDijetPt->Write();
        hGenDijetPtOverPtHatVsDijetPtWeighted->Write();
        hGenDijetPtAveOverPtHatVsDijetPtAve->Write();
        hGenDijetPtAveOverPtHatVsDijetPtAveWeighted->Write();

        hGenInclusiveJetPt->Write();
        hGenInclusiveJetEta->Write();
        hGenInclusiveJetEtaUnweighted->Write();
        hGenInclusiveJetPtEta->Write();
        hGenInclusiveJetPtEtaPtHat->Write();
        hGenLeadJetPtEta->Write();
        hGenLeadJetPtEtaPtHat->Write();
        hGenSubLeadJetPtEta->Write();
        hGenSubLeadJetPtEtaPtHat->Write();
        hGenGoodInclusiveJetEtaLabFrame->Write();
        hGenGoodInclusiveJetEtaCMFrame->Write();

        //
        // Ref jet histograms
        //

        hRefInclusiveJetPt->Write();
        hRefInclusiveJetEta->Write();
        hRefInclusiveJetEtaUnweighted->Write();
        hRefInclusiveJetPtEta->Write();
        hRefInclusiveJetPtEtaPtHat->Write();
        hRefLeadJetPtEta->Write();
        hRefLeadJetPtEtaPtHat->Write();
        hRefLeadUnswappedJetPtEta->Write();
        hRefLeadUnswappedJetPtEtaPtHat->Write();
        hRefSubLeadJetPtEta->Write();
        hRefSubLeadJetPtEtaPtHat->Write();
        hRefSubLeadUnswappedJetPtEta->Write();
        hRefSubLeadUnswappedJetPtEtaPtHat->Write();

        //
        // Ref-selected jet histograms
        //

        hRefSelInclusiveJetPt->Write();
        hRefSelInclusiveJetEta->Write();
        hRefSelInclusiveJetEtaUnweighted->Write();
        hRefSelInclusiveJetPtEta->Write();
        hRefSelInclusiveJetPtEtaPtHat->Write();
        hRefSelLeadJetPtEta->Write();
        hRefSelLeadJetPtEtaPtHat->Write();
        hRefSelSubLeadJetPtEta->Write();
        hRefSelSubLeadJetPtEtaPtHat->Write();
    } // if ( fIsMc )
}
