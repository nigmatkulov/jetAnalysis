// ROOT headers
#include "TF1.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TROOT.h"

// C++ headers
#include <iostream>

// Jet analysis headers
#include "JetAnalysis.h"

ClassImp(JetAnalysis)

//________________
JetAnalysis::JetAnalysis() : BaseAnalysis(), fDebug(kFALSE), fHM(nullptr) {
    /* Empty */
}

//________________
JetAnalysis::~JetAnalysis() {
    if (fHM) {delete fHM; fHM = nullptr; }
}

//________________
void JetAnalysis::init() {
    // Initialize analysis
    //std::cout << "JetAnalysis::init" << std::endl;
}

//________________
void JetAnalysis::processEvent(const Event* event) {
    // Perform the analysis
    //std::cout << "JetAnalysis::processEvent" << std::endl;

    if (fHM) {
        float ptHatW = event->ptHatWeight();
        fHM->hVz->Fill( event->vz() );
        fHM->hVzWeighted->Fill(event->vz(), ptHatW);
        fHM->hHiBin->Fill(event->hiBin());
        fHM->hHiBinWieghted->Fill(event->hiBin(), ptHatW);
        fHM->hPtHat->Fill(event->ptHat());
        fHM->hPtHatWeighted->Fill(event->ptHat(), ptHatW);
        fHM->hPtHatWeight->Fill(ptHatW);
        Double_t centrality = event->centrality();
        fHM->hCentrality->Fill(centrality);

        fHM->hNRecoJets->Fill( event->pfJetCollection()->size() );
        PartFlowJetIterator pfJetIter;
        for ( pfJetIter = event->pfJetCollection()->begin();
              pfJetIter != event->pfJetCollection()->end();
              pfJetIter++ ) {

            // For speed up purpose here
            if ( !(*pfJetIter)->hasMatching() ) continue;

            // Reco jets
            Double_t recoPtCorr = (*pfJetIter)->recoJetPtJECCorr();
            Double_t recoEta = (*pfJetIter)->recoJetEta();
            Double_t recoPhi = (*pfJetIter)->recoJetPhi();
            fHM->hRecoJetPtRaw->Fill( (*pfJetIter)->recoJetPt() );
            fHM->hRecoJetPt->Fill( recoPtCorr );
            fHM->hRecoJetPtWeighted->Fill( recoPtCorr, ptHatW );
            fHM->hRecoJetPtCorrVsPtRaw->Fill( (*pfJetIter)->recoJetPt(), recoPtCorr );
            fHM->hRecoJetEta->Fill( recoEta );
            fHM->hRecoJetPhi->Fill( (*pfJetIter)->recoJetPhi() );

            Double_t recoJetRaw[4]{(Double_t)(*pfJetIter)->recoJetPt(), recoEta, recoPhi, centrality};
            Double_t recoJetCorr[4]{recoPtCorr, recoEta, recoPhi, centrality};
            fHM->hRecoJet->Fill( recoJetRaw );
            fHM->hRecoJetCorr->Fill( recoJetCorr );
            fHM->hRecoJetCorrWeighted->Fill( recoJetCorr, ptHatW );

            // Ref jets
            Double_t genPt = (*pfJetIter)->refJetPt();
            Double_t genEta = (*pfJetIter)->refJetEta();
            Double_t genPhi = (*pfJetIter)->refJetPhi();
            fHM->hRefJetPt->Fill( genPt );
            fHM->hRefJetPtWeighted->Fill( genPt, ptHatW );
            fHM->hRefJetEta->Fill( genEta );
            fHM->hRefJetPhi->Fill( genPhi );

            Double_t genJet[4]{genPt, genEta, genPhi, centrality};
            fHM->hRefJet->Fill( genJet );
            fHM->hRefJetWeighted->Fill( genJet, ptHatW );

            Double_t JESRaw = (*pfJetIter)->recoJetPt() / genPt;
            Double_t JES = recoPtCorr / genPt;
            Double_t flavorFromB = (*pfJetIter)->refFlavorForB();

            Double_t jesRaw[4]{JESRaw, genPt, flavorFromB, centrality};
            Double_t jesCorr[4]{JES, genPt, flavorFromB, centrality};
            fHM->hJESRaw->Fill( jesRaw );
            fHM->hJESRawWeighted->Fill( jesRaw, ptHatW );
            fHM->hJESReco->Fill( jesCorr );
            fHM->hJESRecoWeighted->Fill( jesCorr, ptHatW );

        } // for ( pfJetIter = event->pfJetCollection()->begin();
    } // if (fHM)
}

//________________
void JetAnalysis::finish() {
    // Save data and close files
    std::cout << "JetAnalysis::finish" << std::endl;
}

//________________
void JetAnalysis::report() {
    // Force to report everyone
}

//________________
TList* JetAnalysis::getOutputList() {
    TList *outputList = new TList();

    // Add list of settings for cuts

    return outputList;
}