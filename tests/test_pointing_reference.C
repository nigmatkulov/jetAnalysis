// Load JetCorrector.cc, JetUncertainty.cc and processForestSimple.C first.
#include <stdexcept>
void test_pointing_reference(const char *path = "/tmp/pointing-reference.root") {
 Histograms hs; createHistograms(hs,true);
 std::vector<GenJet> gen; JERSmearingHelper jer;

 // Reco selection and matching govern both references, with reference eta acceptance.
 for (int scenario=0;scenario<6;++scenario) {
   RecoJet a{}, b{};
   a.recoPt=80; b.recoPt=60; a.refPt=20; b.refPt=10;
   a.recoEta=b.recoEta=2.8f; a.refEta=b.refEta=.25f;
   a.recoPhi=0; b.recoPhi=3.14f;
   if (scenario==1) b.refPt=-1;  // unmatched
   if (scenario==2) b.recoPt=39; // reco threshold
   if (scenario==3) b.recoPhi=1; // reco delta phi
   if (scenario==4) a.refEta=3.5f; // reference acceptance
   if (scenario==5) a.recoPt=50, b.recoPt=40; // inclusive thresholds
   const double expected=(scenario==0 || scenario==5) ? -0.75 : 0.;
   std::vector<RecoJet> reco{a,b};
   processRecoDijets(true,true,true,-.75f,reco,gen,hs,&jer,JERSmearingMode::Hybrid,false);
   for (int i=0;i<nEtaCuts;++i) {
     for (auto pair : {std::make_pair(hs.hRefDijetPtEtaCMPointingResArr[i].get(),
                                    hs.hRefDijetPtEtaCMPointingResRefPtAveArr[i].get()),
                       std::make_pair(hs.hRefDijetPtEtaLabPointingResArr[i].get(),
                                    hs.hRefDijetPtEtaLabPointingResRefPtAveArr[i].get())}) {
       auto old=pair.first; auto alt=pair.second;
       if (old->Integral(0,-1,0,-1)!=expected || alt->Integral(0,-1,0,-1)!=expected)
         throw std::runtime_error("Pointing selection or signed weight changed");
       if (expected && (alt->GetBinContent(alt->FindBin(15.,.25))!=expected ||
                       old->GetBinContent(old->FindBin(scenario==5 ? 45. : 70.,.25))!=expected))
         throw std::runtime_error("Pointing fill coordinates changed");
       if (old->GetNbinsX()!=100 || old->GetNbinsY()!=72 ||
           alt->GetNbinsX()!=100 || alt->GetNbinsY()!=72)
         throw std::runtime_error("Pointing binning changed");
       old->Reset(); alt->Reset();
     }
   }
 }
 for (bool pb : {false,true}) for (int p=55;p<1000;p+=10)
 for (int e=0;e<58;++e) {
   float eta=-2.85f+e*.1f;
   RecoJet a{},b{};
   a.recoPt=p+10; b.recoPt=p-10;
   a.refPt=a.refPtSmeared=.84f*a.recoPt;
   b.refPt=b.refPtSmeared=.84f*b.recoPt;
   a.recoEta=eta+.05f; b.recoEta=eta-.05f;
   a.refEta=.9f*eta+.07f; b.refEta=.9f*eta-.07f;
   a.recoPhi=a.refPhi=0; b.recoPhi=b.refPhi=3.14f;
   std::vector<RecoJet> reco{a,b};
   processRecoDijets(false,pb,true,1.f + .02f*e + .0001f*e*e,reco,gen,hs,&jer,JERSmearingMode::Hybrid,false);
 }

 // The same events, eta bins and weights must survive a change of pT coordinate.
 for (int i=0;i<nEtaCuts;++i) {
   for (auto pair : {std::make_pair(hs.hRefDijetPtEtaCMPointingResArr[i].get(),
                                  hs.hRefDijetPtEtaCMPointingResRefPtAveArr[i].get()),
                     std::make_pair(hs.hRefDijetPtEtaLabPointingResArr[i].get(),
                                  hs.hRefDijetPtEtaLabPointingResRefPtAveArr[i].get())}) {
     auto old=pair.first; auto alt=pair.second;
     if (old->GetEntries()!=alt->GetEntries() || std::abs(old->Integral(0,-1,0,-1)-alt->Integral(0,-1,0,-1))>1e-7)
       throw std::runtime_error("Reference populations or weights differ");
     if (!(alt->GetMean(1)<old->GetMean(1)))
       throw std::runtime_error("Reference pT coordinate did not change");
     auto a=std::unique_ptr<TH1D>(old->ProjectionY("oldEta",0,-1));
     auto b=std::unique_ptr<TH1D>(alt->ProjectionY("refPtEta",0,-1));
     for (int k=0;k<a->GetNcells();++k)
       if (std::abs(a->GetBinContent(k)-b->GetBinContent(k))>1e-8 || std::abs(a->GetBinError(k)-b->GetBinError(k))>1e-8)
         throw std::runtime_error("Reference eta or Sumw2 differs");
     std::cout << old->GetName() << " entries=" << old->GetEntries()
               << " weight=" << old->Integral(0,-1,0,-1)
               << " mean pT=" << old->GetMean(1) << " vs " << alt->GetMean(1) << std::endl;
   }
 }
 TString output(path); writeOutput(output,hs,true,JERSmearingMode::Hybrid,false);
}
