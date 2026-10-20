#define main reference_main
#include "processForestSimple.C"
#undef main
#include "jetanalysis/CMSForestReader.h"
#include <cstring>
void equal(bool ok,const char* name) { if(!ok) throw std::runtime_error(std::string("Reader differs: ")+name); }
template<class T> bool bits(T a,T b) {return std::memcmp(&a,&b,sizeof(T))==0;}
int main(int argc,char** argv) {
    if(argc!=4) return 2;
    try {
        jetanalysis::ReaderConfig c; c.input=argv[1];c.mcType=std::stoi(argv[2]);c.pp=std::string(argv[3])=="pp";
        const bool mc=c.mcType!=0;
        TChain hlt("hltanalysis/HltTree"),events("hiEvtAnalyzer/HiTree"),skim("skimanalysis/HltTree"),jets("ak4PFJetAnalyzer/t");
        setupInput(c.input.c_str(),hlt,events,skim,jets,mc,c.pp);
        jetanalysis::CMSForestReader reader(c);
        equal(reader.entries()==hlt.GetEntries(),"entries");
        jetanalysis::Event event;
        std::uint64_t reco=0,gen=0;
        for(Long64_t i=0;i<hlt.GetEntries();++i) {
            equal(hlt.GetEntry(i)>0,"reference GetEntry");equal(reader.readNext(event),"readNext");
            if(true) equal(event.info.vz==vz,"vz");
            if(true) equal(event.info.run==run,"run");
            if(true) equal(event.info.lumi==lumi,"lumi");
            if(true) equal(event.info.evt==evt,"evt");
            if(mc) equal(event.info.ptHat==ptHat,"ptHat");
            if(mc && c.pp) equal(event.info.inputEventWeight==inputEventWeight,"inputEventWeight");
            if(!c.pp) equal(event.info.pBeamScrapingFilter==pBeamScrapingFilter,"pBeamScrapingFilter");
            if(!c.pp) equal(event.info.pPAprimaryVertexFilter==pPAprimaryVertexFilter,"pPAprimaryVertexFilter");
            if(!c.pp) equal(event.info.HBHENoiseFilterResultRun2Loose==HBHENoiseFilterResultRun2Loose,"HBHENoiseFilterResultRun2Loose");
            if(!c.pp) equal(event.info.phfCoincFilter==phfCoincFilter,"phfCoincFilter");
            if(!c.pp) equal(event.info.pVertexFilterCutdz1p0==pVertexFilterCutdz1p0,"pVertexFilterCutdz1p0");
            if(!c.pp) equal(event.info.pVertexFilterCutGplus==pVertexFilterCutGplus,"pVertexFilterCutGplus");
            if(!c.pp) equal(event.info.pVertexFilterCutVtx1==pVertexFilterCutVtx1,"pVertexFilterCutVtx1");
            if(!mc) equal(event.info.HLT_PAAK4PFJet60_Eta5p1_v4==HLT_PAAK4PFJet60_Eta5p1_v4,"HLT_PAAK4PFJet60_Eta5p1_v4");
            if(!mc) equal(event.info.HLT_PAAK4PFJet80_Eta5p1_v3==HLT_PAAK4PFJet80_Eta5p1_v3,"HLT_PAAK4PFJet80_Eta5p1_v3");
            if(!mc) equal(event.info.HLT_PAAK4PFJet100_Eta5p1_v3==HLT_PAAK4PFJet100_Eta5p1_v3,"HLT_PAAK4PFJet100_Eta5p1_v3");
            equal(bits(event.info.eventRho,eventRho),"rho");
            equal(event.recoJets.size()==static_cast<std::size_t>(nRecoJets),"reco count");
            reco+=event.recoJets.size();
            for(int j=0;j<nRecoJets;++j) {
                equal(bits(event.recoJets[j].recoPtRaw,recoJetPtRaw[j]),"rawpt");
                equal(bits(event.recoJets[j].recoTrackMaxPt,recoJetTrackMaxPt[j]),"trackMax");
                equal(bits(event.recoJets[j].recoEta,recoJetEta[j]),"jteta");
                equal(bits(event.recoJets[j].recoPhi,recoJetPhi[j]),"jtphi");
                equal(bits(event.recoJets[j].recoPfNHF,recoJetPfNHF[j]),"jtPfNHF");
                equal(bits(event.recoJets[j].recoPfNEF,recoJetPfNEF[j]),"jtPfNEF");
                equal(bits(event.recoJets[j].recoPfCHF,recoJetPfCHF[j]),"jtPfCHF");
                equal(bits(event.recoJets[j].recoPfMUF,recoJetPfMUF[j]),"jtPfMUF");
                equal(bits(event.recoJets[j].recoPfCEF,recoJetPfCEF[j]),"jtPfCEF");
                equal(bits(event.recoJets[j].recoPfCHM,recoJetPfCHM[j]),"jtPfCHM");
                equal(bits(event.recoJets[j].recoPfCEM,recoJetPfCEM[j]),"jtPfCEM");
                equal(bits(event.recoJets[j].recoPfNHM,recoJetPfNHM[j]),"jtPfNHM");
                equal(bits(event.recoJets[j].recoPfNEM,recoJetPfNEM[j]),"jtPfNEM");
                equal(bits(event.recoJets[j].recoPfMUM,recoJetPfMUM[j]),"jtPfMUM");
                if(mc) equal(bits(event.recoJets[j].refPt,refJetPt[j]),"refpt");
                if(mc) equal(bits(event.recoJets[j].refEta,refJetEta[j]),"refeta");
                if(mc) equal(bits(event.recoJets[j].refPhi,refJetPhi[j]),"refphi");
            }
            if(mc) {
                equal(event.genJets.size()==static_cast<std::size_t>(nGenJets),"gen count");gen+=event.genJets.size();
                for(int j=0;j<nGenJets;++j) {
                    equal(bits(event.genJets[j].pt,genJetPt[j]),"genpt");
                    equal(bits(event.genJets[j].eta,genJetEta[j]),"geneta");
                    equal(bits(event.genJets[j].phi,genJetPhi[j]),"genphi");
                }
            }
        }
        equal(!reader.readNext(event),"end of input");
        std::cout<<"EXACT READER MATCH: "<<reader.entries()<<" events, "<<reco<<" reco jets, "<<gen<<" gen jets\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
