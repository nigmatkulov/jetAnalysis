// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>
void test_unfolding_pair_populations() {
    TH1::AddDirectory(false);
    Histograms hs;
    createHistograms(hs, true);
    auto reco = [](float pt, float eta, float phi, bool ref = true) {
        RecoJet j{};
        j.recoPt = j.refPtSmeared = pt;
        j.recoEta = eta; j.recoPhi = phi; j.refPt = ref ? pt : -1.f;
        return j;
    };
    auto check = [](bool ok, const char *why) { if (!ok) throw std::runtime_error(why); };
    auto components = [](UnfoldingPairHistograms &p, int i) {
        return std::vector<TH2D*>{p.hGenDijetPtEtaCMMatchedRecoSelected[i].get(),
            p.hGenDijetPtEtaCMSelectedRecoFound[i].get(), p.hGenDijetPtEtaCMSelectedRecoNotFound[i].get(),
            p.hRecoDijetPtEtaCMMatchedGenSelected[i].get(), p.hRecoDijetPtEtaCMMatchedGenFound[i].get(),
            p.hRecoDijetPtEtaCMSelectedGenNotFound[i].get()};
    };
    int cases = 0;
    for (bool pp : {false, true}) for (bool pb : {false, true})
    for (bool controlled : {false, true}) for (float weight : {2.5f, -0.75f})
    for (int scenario = 0; scenario < 11; ++scenario) {
        std::vector<GenJet> gen{{100.f,.2f,0.f,-1.f},{60.f,-.2f,3.14f,-1.f}};
        std::vector<RecoJet> measured{reco(100,.2,0),reco(60,-.2,3.14)};
        // Expected nonzero component indices, gen side and reco side (-1=unselected).
        int g = 0, r = 3;
        if (scenario == 1) { measured[1].refPtSmeared = 30; g = 1; r = -1; } // threshold fail
        if (scenario == 2) { measured.pop_back(); g = 2; r = -1; }
        if (scenario == 3) { gen[1].pt = 30; g = -1; r = 4; }
        if (scenario == 4) { gen.pop_back(); g = -1; r = 5; }
        if (scenario == 5) { // true counterparts exist, but third reco jet displaces one
            measured.push_back(reco(80,.8,3.14)); g = 1; r = 5;
        }
        if (scenario == 6) { // selected reco pair corresponds to lower-ranked truth jets
            gen.push_back({80.f,.8f,3.14f,-1.f}); g = 2; r = 4;
        }
        if (scenario == 7) { measured[0].refPt = -1; g = 2; r = 5; }
        if (scenario == 8) { measured[0].refPtSmeared = 60; measured[1].refPtSmeared = 100; }
        if (scenario == 9) { gen.clear(); measured.clear(); g = r = -1; }
        if (scenario == 10) { // one good reco jet cannot supply two truth counterparts
            measured[0].recoEta = .2; measured[1].recoEta = 1.;
            gen[1].eta = .2; gen[1].phi = .1; g = -1; r = 5;
        }
        for (int i = 0; i < nEtaCuts; ++i) {
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[i]->Reset();
            hs.hGenDijetPtEtaCMMissArr[i]->Reset(); hs.hRecoDijetPtEtaCMFakeArr[i]->Reset();
            for (auto *p : {&hs.unfoldingPairs, &hs.unfoldingPairsControlled, &hs.unfoldingPairsStudied})
                for (auto *h : components(*p,i)) h->Reset();
        }
        isFillControlledSample = controlled;
        prepareUnfolding(pp,pb,true,weight,gen,measured,
            hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr,hs.hGenDijetPtEtaCMMissArr,
            hs.hRecoDijetPtEtaCMFakeArr,hs.hUnfoldingPairClassificationArr,
            nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,-1,
            &hs.unfoldingPairs,&hs.unfoldingPairsControlled,&hs.unfoldingPairsStudied);
        for (int i = 0; i < nEtaCuts; ++i) {
            auto h = components(hs.unfoldingPairs,i);
            auto selected = components(controlled ? hs.unfoldingPairsControlled : hs.unfoldingPairsStudied,i);
            auto other = components(controlled ? hs.unfoldingPairsStudied : hs.unfoldingPairsControlled,i);
            for (int k = 0; k < 6; ++k) {
                const double expected = k == g || k == r ? weight : 0.;
                check(h[k]->Integral(0,-1,0,-1) == expected,"population membership/weight");
                check(selected[k]->Integral(0,-1,0,-1) == expected,"split membership");
                check(other[k]->GetEntries() == 0,"split leakage");
                check(h[k]->GetEntries() == (expected != 0 ? 1 : 0),"population exclusivity");
                for (int b = 0; b < h[k]->GetNcells(); ++b) {
                    check(std::abs(h[k]->GetBinError(b)*h[k]->GetBinError(b) -
                          (h[k]->GetBinContent(b) != 0 ? weight*weight : 0.)) < 1e-10,"Sumw2");
                }
            }
            for (int b = 0; b < h[0]->GetNcells(); ++b) {
                check(h[1]->GetBinContent(b)+h[2]->GetBinContent(b)==hs.hGenDijetPtEtaCMMissArr[i]->GetBinContent(b),"miss partition");
                check(h[4]->GetBinContent(b)+h[5]->GetBinContent(b)==hs.hRecoDijetPtEtaCMFakeArr[i]->GetBinContent(b),"fake partition");
            }
            check(h[0]->GetEntries()==hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr[i]->GetEntries(),"matched response membership");
        }
        ++cases;
    }
    // Verify production wiring and first-replica-only fills in both families.
    JERSmearingHelper jer;
    for (auto mode : {JERSmearingMode::Hybrid, JERSmearingMode::Stochastic,
                      JERSmearingMode::TruthToPseudoReco}) for (bool etaDependent : {false,true}) {
        std::vector<GenJet> gen{{100.f,.2f,0.f,-1.f},{60.f,-.2f,3.14f,-1.f}};
        std::vector<RecoJet> measured{reco(100,.2,0),reco(60,-.2,3.14)};
        double before[2][nEtaCuts][3];
        for (int f = 0; f < 2; ++f) for (int i = 0; i < nEtaCuts; ++i) {
            auto h = components(f ? hs.unfoldingPairsJerDefExtra : hs.unfoldingPairs,i);
            before[f][i][0] = h[0]->GetEntries();
            before[f][i][1] = h[1]->GetEntries()+h[2]->GetEntries();
            before[f][i][2] = h[4]->GetEntries()+h[5]->GetEntries();
        }
        processRecoDijets(true,true,true,1.f,measured,gen,hs,&jer,mode,etaDependent);
        for (int f = 0; f < 2; ++f) for (int i = 0; i < nEtaCuts; ++i) {
            auto h = components(f ? hs.unfoldingPairsJerDefExtra : hs.unfoldingPairs,i);
            check(h[0]->GetEntries()+h[1]->GetEntries()+h[2]->GetEntries() -
                  before[f][i][0]-before[f][i][1] == 1.,"caller truth diagnostic replica count");
            check(h[3]->GetEntries()+h[4]->GetEntries()+h[5]->GetEntries() -
                  before[f][i][0]-before[f][i][2] == 1.,"caller reco diagnostic replica count");
        }
    }
    TString output = "/tmp/jetAnalysis-unfolding-pair-populations.root";
    writeOutput(output,hs,true,JERSmearingMode::Hybrid,false);
    TFile file(output,"READ");
    for (auto *p : {&hs.unfoldingPairs,&hs.unfoldingPairsJerDefExtra,&hs.unfoldingPairsControlled,&hs.unfoldingPairsStudied})
        for (int i = 0; i < nEtaCuts; ++i) for (auto *h : components(*p,i))
            check(file.Get(h->GetName()) != nullptr,"missing persisted diagnostic key");
    std::cout << "unfolding populations: " << cases << " cases, all cuts, signed weights, partitions and persisted keys passed\n";
}
