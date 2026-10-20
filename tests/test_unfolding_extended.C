// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>
void test_unfolding_extended() {
    TH1::AddDirectory(false);
    Histograms hs; createHistograms(hs,true);
    auto check = [](bool ok,const char *why) { if (!ok) throw std::runtime_error(why); };
    auto makeReco = [](float pt,float eta,float phi) {
        RecoJet j{}; j.recoPt=j.refPtSmeared=pt; j.recoEta=eta; j.recoPhi=phi; j.refPt=pt; return j;
    };
    auto reset = [&]() {
        for (auto *e : {&hs.unfoldingExtended,&hs.unfoldingExtendedJerDefExtra,
                        &hs.unfoldingExtendedControlled,&hs.unfoldingExtendedStudied})
            for (int i=0;i<nEtaCuts;++i) {
                for (auto *h : {e->truth[i].get(),e->matchedTruth[i].get(),e->miss[i].get(),e->response[i].get()}) {
                    h->Reset(); h->Sumw2();
                }
                for (auto *h : {e->measured[i].get(),e->matchedMeasured[i].get(),e->fake[i].get()}) h->Reset();
            }
    };
    auto verify = [&](ExtendedUnfoldingHistograms &e,int i,bool truth,bool fid,bool reco,float w) {
        const bool matched=truth&&reco;
        check(e.truth[i]->GetEntries()==(truth?1:0),"extended truth membership");
        check(e.response[i]->GetEntries()==(matched?1:0),"extended response membership");
        check(e.miss[i]->GetEntries()==(truth&&!reco?1:0),"extended misses");
        check(e.measured[i]->GetEntries()==(reco?1:0),"extended measured");
        check(e.fake[i]->GetEntries()==(reco&&!truth?1:0),"extended fakes");
        int dims[3]={0,1,4};
        auto responseTruth=std::unique_ptr<THnSparse>(e.response[i]->Projection(3,dims,"E"));
        for (Long64_t b=0;b<e.truth[i]->GetNbins();++b) {
            int c[3]; const double v=e.truth[i]->GetBinContent(b,c);
            check(c[2]==(fid?2:1),"truth fiducial category");
            check(v==w && e.truth[i]->GetBinError2(b)==w*w,"truth weight/Sumw2");
            auto value = [&](THnSparse *h) { auto k=h->GetBin(c,false); return k<0?0.:h->GetBinContent(k); };
            auto error2 = [&](THnSparse *h) { auto k=h->GetBin(c,false); return k<0?0.:h->GetBinError2(k); };
            check(v==value(e.matchedTruth[i].get())+value(e.miss[i].get()),"truth accounting");
            check(value(responseTruth.get())==value(e.matchedTruth[i].get()),"response truth marginal");
            check(error2(responseTruth.get())==error2(e.matchedTruth[i].get()),"response truth Sumw2");
            check(e.truth[i]->GetBinError2(b)==error2(e.matchedTruth[i].get())+error2(e.miss[i].get()),"miss Sumw2 partition");
        }
        auto responseReco=std::unique_ptr<TH2D>(e.response[i]->Projection(3,2,"E"));
        for (int b=0;b<e.measured[i]->GetNcells();++b) {
            check(e.measured[i]->GetBinContent(b)==e.matchedMeasured[i]->GetBinContent(b)+e.fake[i]->GetBinContent(b),"measured accounting");
            check(responseReco->GetBinContent(b)==e.matchedMeasured[i]->GetBinContent(b),"response reco marginal");
            check(std::abs(e.measured[i]->GetBinError(b)*e.measured[i]->GetBinError(b)-
                e.matchedMeasured[i]->GetBinError(b)*e.matchedMeasured[i]->GetBinError(b)-
                e.fake[i]->GetBinError(b)*e.fake[i]->GetBinError(b))<1e-10,"fake Sumw2 partition");
        }
    };
    int cases=0;
    for (bool pp : {false,true}) for (bool pb : {false,true})
    for (bool controlled : {false,true}) for (float w : {2.5f,-.75f})
    for (int scenario=0;scenario<18;++scenario) {
        reset(); isFillControlledSample=controlled;
        std::vector<GenJet> gen{{100.f,.2f,0.f,-1.f},{60.f,-.2f,3.14f,-1.f}};
        std::vector<RecoJet> reco{makeReco(100,.2,0),makeReco(60,-.2,3.14)};
        if(scenario==1) gen[1].pt=30; // nonfiducial threshold, reco selected
        if(scenario==2) gen[1].phi=.1; // nonfiducial angular cut, reco selected
        if(scenario==3) gen[0].eta=1.8; // category varies across acceptance cuts
        if(scenario==4) reco.clear();
        if(scenario==5) reco.resize(1);
        if(scenario==6) gen.clear();
        if(scenario==7) gen.resize(1);
        if(scenario==8) { reco.push_back(makeReco(80,.8,3.14)); } // rank displacement
        if(scenario==9) reco[0].refPt=reco[1].refPt=-1; // no reference gating in extended
        if(scenario==10) { gen[0].eta=gen[1].eta=6.; } // coordinate overflow, nonfiducial
        if(scenario==11) { gen[0].pt=1500;gen[1].pt=1200; } // coordinate overflow, fiducial
        if(scenario==12) { gen[1].pt=30; reco.clear(); } // nonfiducial inefficiency
        if(scenario==13) { gen[0].eta=gen[1].eta=6.; reco.clear(); } // overflow miss
        if(scenario==14) { gen.clear(); reco.clear(); } // neither population exists
        if(scenario==15) reco[1].refPtSmeared=30;
        if(scenario==16) reco[1].recoPhi=.1;
        if(scenario==17) reco[0].recoEta=1.8;
        prepareUnfolding(pp,pb,true,w,gen,reco,hs.hGenDijetPtEtaCMVsRecoPtEtaCMArr,
            hs.hGenDijetPtEtaCMMissArr,hs.hRecoDijetPtEtaCMFakeArr,hs.hUnfoldingPairClassificationArr,
            nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,-1,nullptr,nullptr,nullptr,
            &hs.unfoldingExtended,&hs.unfoldingExtendedControlled,&hs.unfoldingExtendedStudied);
        for(int i=0;i<nEtaCuts;++i) {
            bool fid=gen.size()>=2 && analysis_contract::passesDijetSelection(gen[0].pt,gen[1].pt,
                etaCM(gen[0].eta,analysis_contract::kNominalEtaShift,pb,true,pp),
                etaCM(gen[1].eta,analysis_contract::kNominalEtaShift,pb,true,pp),
                deltaPhi(gen[0].phi,gen[1].phi),etaCuts[i]);
            const bool selected=reco.size()>=2 && analysis_contract::passesDijetSelection(
                reco[0].refPtSmeared,reco[1].refPtSmeared,
                etaCM(reco[0].recoEta,analysis_contract::kNominalEtaShift,pb,true,pp),
                etaCM(reco[1].recoEta,analysis_contract::kNominalEtaShift,pb,true,pp),
                deltaPhi(reco[0].recoPhi,reco[1].recoPhi),etaCuts[i]);
            verify(hs.unfoldingExtended,i,gen.size()>=2,fid,selected,w);
            verify(controlled?hs.unfoldingExtendedControlled:hs.unfoldingExtendedStudied,i,gen.size()>=2,fid,selected,w);
            verify(controlled?hs.unfoldingExtendedStudied:hs.unfoldingExtendedControlled,i,false,false,false,w);
        }
        ++cases;
    }
    JERSmearingHelper jer;
    for(auto mode : {JERSmearingMode::Hybrid,JERSmearingMode::Stochastic,JERSmearingMode::TruthToPseudoReco})
    for(bool eta : {false,true}) for(int nr : {0,1,2}) {
        reset(); isFillControlledSample=true;
        std::vector<GenJet> gen{{100.f,.2f,0.f,-1.f},{60.f,-.2f,3.14f,-1.f}};
        std::vector<RecoJet> reco{makeReco(100,.2,0),makeReco(60,-.2,3.14)};reco.resize(nr);
        processRecoDijets(true,true,true,1.f,reco,gen,hs,&jer,mode,eta);
        for(int i=0;i<nEtaCuts;++i) for(auto *e : {&hs.unfoldingExtended,&hs.unfoldingExtendedJerDefExtra,&hs.unfoldingExtendedControlled})
            verify(*e,i,true,true,nr==2,1.f);
    }
    TString path="/tmp/jetAnalysis-unfolding-extended.root";
    writeOutput(path,hs,true,JERSmearingMode::Hybrid,false);
    TFile file(path,"READ");
    check(file.Get("unfoldingExtendedContract")!=nullptr,"missing extended metadata");
    int persisted=0;
    for(auto *key : *file.GetListOfKeys())
        if(std::string(key->GetName()).find("hUnfoldingExtended")==0) ++persisted;
    check(persisted==224,"missing extended output keys");
    std::cout<<"extended unfolding: "<<cases<<" cases, flows, category marginals, Sumw2, caller modes and 224 keys passed\n";
}
