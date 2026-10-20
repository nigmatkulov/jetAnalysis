#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TNamed.h"
#include "TString.h"
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
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>
#include <cstdint>
#include <cstring>


#include "Histograms.h"
#include "jetanalysis/Output.h"
namespace jetanalysis {
enum class JERSmearingMode : int {
    Hybrid = 0,
    Stochastic = 1,
    TruthToPseudoReco = 2
};
void writeOutput(const std::string &output, Histograms &hs, const bool isMc,
                 const int mode, const bool useEtaDependentJER, const OutputSelection& selection) {

    TString oFileName(output);
    const auto jerSmearingMode=static_cast<JERSmearingMode>(mode);
    std::cout << "Writing output to file: " << oFileName.Data();
    // Output file
    int compressionSetting = 208; // LZMA compression
    // TFile *fOut = TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting);
    auto fOut = std::unique_ptr<TFile>( TFile::Open( oFileName.Data(), "RECREATE", "", compressionSetting) );
    if (!fOut || fOut->IsZombie()) throw std::runtime_error("Cannot create ROOT output: "+output);
    if (isMc) {
        const char *jerModeName = jerSmearingMode == JERSmearingMode::Hybrid
            ? "hybrid"
            : (jerSmearingMode == JERSmearingMode::Stochastic
                ? "stochastic"
                : "truth-to-pseudo-reco");
        TNamed("unfoldingExtendedContract",
               "v1: same-event global leading gen/reco pairs; no geometric matching; "
               "truth category 1=outside fiducial,2=fiducial; "
               "response axes=genPt,genEta,recoPt,recoEta,truthCategory; "
               "truth/miss axes=genPt,genEta,truthCategory; preserve coordinate flow bins").Write();
        TNamed("jerSmearingMode", jerModeName).Write();
        TNamed("useEtaDependentJER", useEtaDependentJER ? "true" : "false").Write();
    }

    //
    // Event level histograms
    //
    selection.write(hs.hVzRaw.get());
    selection.write(hs.hPtHatUnweighted.get());
    selection.write(hs.hPtHat.get());
    selection.write(hs.hVzUnweighted.get());
    selection.write(hs.hVz.get());

    //
    // Gen level histograms
    //
    if (isMc) {
        // Event histograms
        selection.write(hs.hGenDijetPtAveOverPtHatVsPtHat.get());
        selection.write(hs.hGenDijetPtAveOverPtHatVsPtHatPass.get());
        selection.write(hs.hGenLeadJetPtOverPtHatVsPtHat.get());
        selection.write(hs.hGenLeadJetPtOverPtHatVsPtHatPass.get());
        selection.write(hs.hRecoDijetPtAveOverPtHatVsPtHat.get());
        selection.write(hs.hRecoDijetPtAveOverPtHatVsPtHatPass.get());
        selection.write(hs.hRecoLeadJetPtOverPtHatVsPtHat.get());
        selection.write(hs.hRecoLeadJetPtOverPtHatVsPtHatPass.get());

        // Gen jets
        selection.write(hs.hGenInclusiveJetPtEtaLabUnflippedUnweighted.get());
        selection.write(hs.hGenInclusiveJetPtEtaLabUnflipped.get());
        selection.write(hs.hGenInclusiveJetPtEtaLab.get());
        selection.write(hs.hGenInclusiveJetPtEtaCM.get());
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            selection.write(hs.hGenInclusiveJetEtaCMShiftedUnweighted[iShift].get());
            selection.write(hs.hGenInclusiveJetEtaCMShifted[iShift].get());
        }
        selection.write(hs.hGenInclusiveJetPtEtaStdBins.get());
        selection.write(hs.hGenInclusiveJetPtPtHat.get());

        // Gen dijets
        selection.write(hs.hGenDijetPtAve.get());
        for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
            selection.write(hs.hGenDijetEtaCMShiftedUnweighted[iShift].get());
            selection.write(hs.hGenDijetEtaCMShifted[iShift].get());
        }

        selection.write(hs.hGenDijetPtEtaLabUnflipped.get());
        selection.write(hs.hGenDijetPtEtaLab.get());
        selection.write(hs.hGenDijetPtEtaCM.get());

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
            selection.write(hs.hGenDijetPtEtaLabArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMXPbArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMXpArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMXPbOverXpArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMControlledSampleArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMStudiedSampleArr[iCut].get());
#ifdef DO_MIXED_PRIOR
            selection.write(hs.hGenDijetPtEtaCMMixedPriorArr[iCut].get());
#endif

            // Smearing studies
            selection.write(hs.hGenDijetDefPtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetDefExtraPtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetDefDoublePtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetDefNinetyPtEtaCMArr[iCut].get());

            for (auto *extended : {&hs.unfoldingExtended, &hs.unfoldingExtendedJerDefExtra,
                                   &hs.unfoldingExtendedControlled, &hs.unfoldingExtendedStudied}) {
                extended->response[iCut]->Write();
                extended->truth[iCut]->Write();
                extended->matchedTruth[iCut]->Write();
                extended->miss[iCut]->Write();
                extended->measured[iCut]->Write();
                extended->matchedMeasured[iCut]->Write();
                extended->fake[iCut]->Write();
            }

            for (auto *pairs : {&hs.unfoldingPairs, &hs.unfoldingPairsJerDefExtra,
                                &hs.unfoldingPairsControlled, &hs.unfoldingPairsStudied}) {
                pairs->hRecoDijetPtEtaCMMatchedGenSelected[iCut]->Write();
                pairs->hRecoDijetPtEtaCMMatchedGenFound[iCut]->Write();
                pairs->hGenDijetPtEtaCMSelectedRecoFound[iCut]->Write();
                pairs->hGenDijetPtEtaCMSelectedRecoNotFound[iCut]->Write();
                pairs->hGenDijetPtEtaCMMatchedRecoSelected[iCut]->Write();
                pairs->hRecoDijetPtEtaCMSelectedGenNotFound[iCut]->Write();
            }

            // Response matrix 
            selection.write(hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMVsRecoPtEtaCMControlledSampleArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMVsRecoPtEtaCMStudiedSampleArr[iCut].get());
            // selection.write(hs.hRefDijetPtEtaCMVsRecoPtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMMissArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMMissControlledSampleArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMMissStudiedSampleArr[iCut].get());
            selection.write(hs.hRecoDijetPtEtaCMFakeArr[iCut].get());
            selection.write(hs.hRecoDijetPtEtaCMFakeControlledSampleArr[iCut].get());
            selection.write(hs.hRecoDijetPtEtaCMFakeStudiedSampleArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMVsRecoJerDefExtraPtEtaCMArr[iCut].get());
            selection.write(hs.hGenDijetPtEtaCMMissJerDefExtraArr[iCut].get());
            selection.write(hs.hRecoDijetPtEtaCMFakeJerDefExtraArr[iCut].get());
            selection.write(hs.hUnfoldingPairClassificationArr[iCut].get());
            selection.write(hs.hUnfoldingPairClassificationJerDefExtraArr[iCut].get());
        }

        selection.write(hs.hGenInclusiveJetJESDefPtEta.get());
        selection.write(hs.hGenInclusiveJetJESDefExtraPtEta.get());
        selection.write(hs.hGenInclusiveJetJESDefDoublePtEta.get());
        selection.write(hs.hGenInclusiveJetJESDefNinetyPtEta.get());

        selection.write(hs.hGenInclusiveJetDefPtEtaFlipped.get());
        selection.write(hs.hGenInclusiveJetDefExtraPtEtaFlipped.get());
        selection.write(hs.hGenInclusiveJetDefDoublePtEtaFlipped.get());
        selection.write(hs.hGenInclusiveJetDefNinetyPtEtaFlipped.get());

        //
        // Ref level histograms
        //
        
        // Ref jets
        selection.write(hs.hRefInclusiveJetNoSelPtEtaLabUnflipped.get());
        selection.write(hs.hRefInclusiveJetNoSelPtEtaLab.get());
        selection.write(hs.hRefInclusiveJetNoSelPtEtaCM.get());
        selection.write(hs.hRefInclusiveJetTrkMaxPtEtaLabUnflipped.get());
        selection.write(hs.hRefInclusiveJetTrkMaxPtEtaLab.get());
        selection.write(hs.hRefInclusiveJetTrkMaxPtEtaCM.get());
        selection.write(hs.hRefInclusiveJetPtEtaLabUnflipped.get());
        selection.write(hs.hRefInclusiveJetPtEtaLab.get());
        selection.write(hs.hRefInclusiveJetPtEtaCM.get());
        selection.write(hs.hRefSelInclusiveJetNoSelPtEtaLabUnflipped.get());
        selection.write(hs.hRefSelInclusiveJetNoSelPtEtaLab.get());
        selection.write(hs.hRefSelInclusiveJetNoSelPtEtaCM.get());
        selection.write(hs.hRefSelInclusiveJetTrkMaxPtEtaLabUnflipped.get());
        selection.write(hs.hRefSelInclusiveJetTrkMaxPtEtaLab.get());
        selection.write(hs.hRefSelInclusiveJetTrkMaxPtEtaCM.get());
        selection.write(hs.hRefSelInclusiveJetPtEtaLabUnflipped.get());
        selection.write(hs.hRefSelInclusiveJetPtEtaLab.get());
        selection.write(hs.hRefSelInclusiveJetPtEtaCM.get());

        selection.write(hs.hRefInclusiveJetPtEtaStdBins.get());

        // Ref dijets
        selection.write(hs.hRefDijetPtEtaLabUnflipped.get());
        selection.write(hs.hRefDijetPtEtaLab.get());
        selection.write(hs.hRefDijetPtEtaCM.get());

        for (int iCut{0}; iCut < nEtaCuts; ++iCut) {

            //
            // Ref level histograms
            //
            selection.write(hs.hRefDijetPtEtaLabArr[iCut].get());
            selection.write(hs.hRecoDijetPtEtaLabPointingResArr[iCut].get());
            selection.write(hs.hRefDijetPtEtaLabPointingResArr[iCut].get());
            selection.write(hs.hRefDijetPtEtaLabPointingResRefPtAveArr[iCut].get());
            selection.write(hs.hRefDijetPtEtaCMArr[iCut].get());
            selection.write(hs.hRefDijetPtEtaCMControlledSampleArr[iCut].get());
            selection.write(hs.hRefDijetPtEtaCMStudiedSampleArr[iCut].get());

            //
            // Ref-selected level histograms
            //
            selection.write(hs.hRefSelDijetPtEtaCMArr[iCut].get());
        }

    } // if (isMc)

    //
    // Reco level histograms
    //
    
    // Reco jets
    selection.write(hs.hRecoInclusiveJetNoSelPtEtaLabUnflipped.get());
    selection.write(hs.hRecoInclusiveJetNoSelPtEtaLab.get());
    selection.write(hs.hRecoInclusiveJetNoSelPtEtaCM.get());
    selection.write(hs.hRecoInclusiveJetTrkMaxPtEtaLabUnflipped.get());
    selection.write(hs.hRecoInclusiveJetTrkMaxPtEtaLab.get());
    selection.write(hs.hRecoInclusiveJetTrkMaxPtEtaCM.get());
    selection.write(hs.hRecoInclusiveJetPtEtaLabUnflippedUnweighted.get());
    selection.write(hs.hRecoInclusiveJetPtEtaLabUnflipped.get());
    selection.write(hs.hRecoInclusiveJetPtEtaLab.get());
    selection.write(hs.hRecoInclusiveJetPtEtaCM.get());

    // for (int iShift = 0; iShift < nEtaShifts; ++iShift) {
    //     selection.write(hs.hRecoInclusiveJetEtaCMShiftedUnweighted[iShift].get());
    //     selection.write(hs.hRecoInclusiveJetEtaCMShifted[iShift].get());
    // }
    selection.write(hs.hRecoInclusiveJetPtEtaStdBins.get());
    if (isMc) {
        selection.write(hs.hRecoInclusiveJetJESPtEta.get());
        selection.write(hs.hRecoInclusiveJetPtEtaLabMatched.get());
        selection.write(hs.hRecoInclusiveJetPtEtaLabUnmatched.get());
        selection.write(hs.hRecoInclusiveJetJESDefNoSFPtEta.get());
        selection.write(hs.hRecoInclusiveJetJESDefDoublePtEta.get());
        selection.write(hs.hRecoInclusiveJetJESDefNinetyPtEta.get());
        selection.write(hs.hRecoInclusiveJetJESDefPtEta.get());
        selection.write(hs.hRecoInclusiveJetJESUpPtEta.get());
        selection.write(hs.hRecoInclusiveJetJESDownPtEta.get());
        selection.write(hs.hRecoInclusiveJetPtPtHat.get());
    }
    
    // Reco dijets
    selection.write(hs.hRecoDijetPtEtaLabUnflipped.get());
    selection.write(hs.hRecoDijetPtEtaLab.get());
    selection.write(hs.hRecoDijetPtEtaCM.get());
    selection.write(hs.hRecoInclusiveJetRawPtEtaLabUnflipped.get());
    if (!isMc) {
        for (int iRun{0}; iRun < nDiagnosticRunHistograms; ++iRun) {
            selection.write(hs.hRecoDijetPtEtaCMRunArr[iRun].get());
        }
    }

    for (int iCut{0}; iCut < nEtaCuts; ++iCut) {
        selection.write(hs.hRecoDijetPtEtaLabArr[iCut].get());
        selection.write(hs.hRecoDijetPtEtaCMArr[iCut].get());
        if (isMc) {
            selection.write(hs.hRecoDijetPtEtaCMControlledSampleArr[iCut].get());
            selection.write(hs.hRecoDijetPtEtaCMStudiedSampleArr[iCut].get());
            // Reco
            if (hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDefNoSF[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDefDoubleSF[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDefNinetySF[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerUp[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerUp[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDown[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDown[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDef[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDef[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDownExtra[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDownExtra[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDefExtra[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDefExtra[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerUpExtra[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerUpExtra[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDefUnfold[iCut].get());
            if (hs.hRecoDijetPtEtaCMJerDefUnfoldControlledSampleArr[iCut]) selection.write(hs.hRecoDijetPtEtaCMJerDefUnfoldControlledSampleArr[iCut].get());
            if (hs.hRecoDijetPtEtaCMJerDefUnfoldStudiedSampleArr[iCut]) selection.write(hs.hRecoDijetPtEtaCMJerDefUnfoldStudiedSampleArr[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJerDefExtraUnfold[iCut].get());
            // Ref
            if (hs.hRefDijetPtEtaCMArrJerUp[iCut]) selection.write(hs.hRefDijetPtEtaCMArrJerUp[iCut].get());
            if (hs.hRefDijetPtEtaCMArrJerDown[iCut]) selection.write(hs.hRefDijetPtEtaCMArrJerDown[iCut].get());
            if (hs.hRefDijetPtEtaCMArrJerDef[iCut]) selection.write(hs.hRefDijetPtEtaCMArrJerDef[iCut].get());
            if (hs.hRecoDijetPtEtaCMPointingResArr[iCut]) selection.write(hs.hRecoDijetPtEtaCMPointingResArr[iCut].get());
            if (hs.hRefDijetPtEtaCMPointingResArr[iCut]) selection.write(hs.hRefDijetPtEtaCMPointingResArr[iCut].get());
            if (hs.hRefDijetPtEtaCMPointingResRefPtAveArr[iCut]) selection.write(hs.hRefDijetPtEtaCMPointingResRefPtAveArr[iCut].get());
        }

        if (!isMc) {
            if (hs.hRecoDijetPtEtaCMArrJeuUp[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJeuUp[iCut].get());
            if (hs.hRecoDijetPtEtaCMArrJeuDown[iCut]) selection.write(hs.hRecoDijetPtEtaCMArrJeuDown[iCut].get());
        }
    }

    fOut->Close();
    std::cout << GREEN << "\t[DONE]" << RESET << std::endl;
}


}
