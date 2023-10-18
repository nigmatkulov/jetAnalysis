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
        fHM->hHiBinWeighted->Fill(event->hiBin(), ptHatW);
        fHM->hPtHat->Fill(event->ptHat());
        fHM->hPtHatWeighted->Fill(event->ptHat(), ptHatW);
        fHM->hPtHatWeight->Fill(ptHatW);
        Double_t centrality = event->centrality();
        fHM->hCentrality->Fill(centrality);

        Double_t ptHat = event->ptHat();

        std::cout << "HiBin: " << event->hiBin() << " centrality: " << centrality << std::endl;

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
            Double_t recoJetCorrPtHat[4];
            if ( centrality < 10 ) {
                recoJetCorrPtHat[0] = recoPtCorr;
                recoJetCorrPtHat[1] = recoEta;
                recoJetCorrPtHat[2] = recoPhi; 
                recoJetCorrPtHat[3] = ptHat;
            }
            fHM->hRecoJet->Fill( recoJetRaw );
            fHM->hRecoJetCorr->Fill( recoJetCorr );
            fHM->hRecoJetCorrWeighted->Fill( recoJetCorr, ptHatW );
            if ( centrality < 10 ) {
                fHM->hRecoJetCorrEtaPhiPtHat->Fill( recoJetCorrPtHat);
                fHM->hRecoJetCorrEtaPhiPtHatWeighted->Fill( recoJetCorrPtHat, ptHatW );
            }

            // Ref jets
            Double_t genPt = (*pfJetIter)->refJetPt();
            Double_t genEta = (*pfJetIter)->refJetEta();
            Double_t genPhi = (*pfJetIter)->refJetPhi();
            fHM->hRefJetPt->Fill( genPt );
            fHM->hRefJetPtWeighted->Fill( genPt, ptHatW );
            fHM->hRefJetEta->Fill( genEta );
            fHM->hRefJetPhi->Fill( genPhi );

            fHM->hRecoJetPtCorrVsGenPt->Fill(genPt, recoPtCorr);

            Double_t genJet[4]{genPt, genEta, genPhi, centrality};
            Double_t genJetPtHat[4];
            if ( centrality < 10 ) {
                genJetPtHat[0] = genPt;
                genJetPtHat[1] = genEta;
                genJetPtHat[2] = genPhi; 
                genJetPtHat[3] = ptHat;
            }
            fHM->hRefJet->Fill( genJet );
            fHM->hRefJetWeighted->Fill( genJet, ptHatW );
            if ( centrality < 10 ) {
                fHM->hRefJetEtaPhiPtHat->Fill(genJetPtHat);
                fHM->hRefJetEtaPhiPtHatWeighted->Fill(genJetPtHat, ptHatW);
            }

            Double_t JESRaw = (*pfJetIter)->recoJetPt() / genPt;
            Double_t JES = recoPtCorr / genPt;
            Double_t JERRaw = ((*pfJetIter)->recoJetPt() - genPt) / genPt;
            Double_t JER = (recoPtCorr - genPt) / genPt;
            Double_t flavorFromB = (*pfJetIter)->refFlavorForB();

            Double_t jesRaw[4]{JESRaw, genPt, flavorFromB, centrality};
            Double_t jesCorr[4]{JES, genPt, flavorFromB, centrality};
            Double_t jesCorrPtHat[4];
            if ( centrality < 10 ) {
                jesCorrPtHat[0] = JES;
                jesCorrPtHat[1] = genPt;
                jesCorrPtHat[2] = genEta;
                jesCorrPtHat[3] = ptHat;
            }
            fHM->hJESRaw->Fill( jesRaw );
            fHM->hJESRawWeighted->Fill( jesRaw, ptHatW );
            fHM->hJESReco->Fill( jesCorr );
            fHM->hJESRecoWeighted->Fill( jesCorr, ptHatW );
            if ( centrality < 10 ) {
                fHM->hJESRecoEtaPhiPtHat->Fill(jesCorrPtHat);
                fHM->hJESRecoEtaPhiPtHatWeighted->Fill(jesCorrPtHat, ptHatW);
            }

            Double_t jerRaw[4]{JERRaw, genPt, flavorFromB, centrality};
            Double_t jerCorr[4]{JER, genPt, flavorFromB, centrality};
            fHM->hJERRaw->Fill(jerRaw);
            fHM->hJERRawWeighted->Fill(jerRaw, ptHatW);
            fHM->hJERReco->Fill(jerCorr);
            fHM->hJERRecoWeighted->Fill(jerCorr, ptHatW);

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