#include "jetanalysis/CMSForestReader.h"
#include "jetanalysis/RootBranchContract.h"
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TClass.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <functional>
#include <algorithm>
#include <optional>
namespace jetanalysis {
namespace {
template<class T> using Value=std::unique_ptr<TTreeReaderValue<T>>;
template<class T> using Array=std::unique_ptr<TTreeReaderArray<T>>;
template<class T> Value<T> scalar(TTreeReader& reader,const char* name) {
    requireLeaf<T>(*reader.GetTree(),name);
    return std::make_unique<TTreeReaderValue<T>>(reader,name);
}
template<class T> Array<T> array(TTreeReader& reader,const char* name,const char* count) {
    if (count) requireLeaf<T>(*reader.GetTree(),name,count);
    else {
        auto* b=reader.GetTree()->GetBranch(name);
        TClass* cl=nullptr; EDataType type=kOther_t;
        if (!b || b->GetExpectedType(cl,type)!=0 || cl!=TClass::GetClass(typeid(std::vector<T>)))
            throw std::runtime_error(std::string("Expected exact vector branch: ")+name);
    }
    return std::make_unique<TTreeReaderArray<T>>(reader,name);
}
template<class T> T get(const Value<T>& value) {
    if (!value || value->GetSetupStatus()<0) throw std::runtime_error("Invalid scalar branch reader");
    return **value;
}
template<class T> void size(const Array<T>& value,std::size_t n) {
    if (!value || value->GetSetupStatus()<0 || value->GetSize()!=n)
        throw std::runtime_error("Array size or type differs from validated schema");
}
struct FileState {
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTreeReader> events, skim, hlt, jets, tracks, particles;
    ReaderConfig config;
    std::string rhoName;
    Long64_t entries=0;
    Value<Float_t> v_vz;
    Value<UInt_t> v_run;
    Value<UInt_t> v_lumi;
    Value<ULong64_t> v_evt;
    Value<Float_t> v_pthat;
    Value<Float_t> v_weight;
    Value<Int_t> v_pBeamScrapingFilter;
    Value<Int_t> v_pPAprimaryVertexFilter;
    Value<Int_t> v_HBHENoiseFilterResultRun2Loose;
    Value<Int_t> v_phfCoincFilter;
    Value<Int_t> v_pVertexFilterCutdz1p0;
    Value<Int_t> v_pVertexFilterCutGplus;
    Value<Int_t> v_pVertexFilterCutVtx1;
    Value<Int_t> v_HLT_PAAK4PFJet60_Eta5p1_v4;
    Value<Int_t> v_HLT_PAAK4PFJet80_Eta5p1_v3;
    Value<Int_t> v_HLT_PAAK4PFJet100_Eta5p1_v3;
    Value<Float_t> rho;
    Value<Int_t> nref, ngen, nTrk;
    Array<Float_t> a_rawpt;
    Array<Float_t> a_trackMax;
    Array<Float_t> a_jteta;
    Array<Float_t> a_jtphi;
    Array<Float_t> a_jtPfNHF;
    Array<Float_t> a_jtPfNEF;
    Array<Float_t> a_jtPfCHF;
    Array<Float_t> a_jtPfMUF;
    Array<Float_t> a_jtPfCEF;
    Array<Int_t> a_jtPfCHM;
    Array<Int_t> a_jtPfCEM;
    Array<Int_t> a_jtPfNHM;
    Array<Int_t> a_jtPfNEM;
    Array<Int_t> a_jtPfMUM;
    Array<Float_t> a_refpt;
    Array<Float_t> a_refeta;
    Array<Float_t> a_refphi;
    Array<Float_t> a_genpt;
    Array<Float_t> a_geneta;
    Array<Float_t> a_genphi;
    Array<Float_t> a_trkPt;
    Array<Float_t> a_trkEta;
    Array<Float_t> a_trkPhi;
    Array<Float_t> a_pt;
    Array<Float_t> a_eta;
    Array<Float_t> a_phi;
    Array<Int_t> a_trkCharge;
    Array<Int_t> a_chg;
    Array<Int_t> a_pdg;
    Array<Int_t> a_sube;
    Array<Bool_t> a_highPurity;
    FileState(const std::string& path,const ReaderConfig& c): file(TFile::Open(path.c_str(),"READ")), config(c) {
        if (!file || file->IsZombie() || file->GetNkeys()==0) throw std::runtime_error("Unusable ROOT file: "+path);
        auto open=[&](const std::string& name) {
            auto* t=file->Get<TTree>(name.c_str());
            if (!t) throw std::runtime_error("Missing required TTree: "+name+" in "+path);
            if (events && t->GetEntries()!=entries) throw std::runtime_error("Unaligned trees in "+path);
            return std::make_unique<TTreeReader>(t);
        };
        events=open("hiEvtAnalyzer/HiTree"); entries=events->GetTree()->GetEntries();
        // Retain the reference forest's required tree inventory by default.
        hlt=open("hltanalysis/HltTree"); skim=open("skimanalysis/HltTree");
        const bool mc=config.mcType!=0;
        if (true) v_vz=scalar<Float_t>(*events,"vz");
        if (true) v_run=scalar<UInt_t>(*events,"run");
        if (true) v_lumi=scalar<UInt_t>(*events,"lumi");
        if (true) v_evt=scalar<ULong64_t>(*events,"evt");
        if (mc) v_pthat=scalar<Float_t>(*events,"pthat");
        if (mc && config.pp) v_weight=scalar<Float_t>(*events,"weight");
        if (!config.pp) v_pBeamScrapingFilter=scalar<Int_t>(*skim,"pBeamScrapingFilter");
        if (!config.pp) v_pPAprimaryVertexFilter=scalar<Int_t>(*skim,"pPAprimaryVertexFilter");
        if (!config.pp) v_HBHENoiseFilterResultRun2Loose=scalar<Int_t>(*skim,"HBHENoiseFilterResultRun2Loose");
        if (!config.pp) v_phfCoincFilter=scalar<Int_t>(*skim,"phfCoincFilter");
        if (!config.pp) v_pVertexFilterCutdz1p0=scalar<Int_t>(*skim,"pVertexFilterCutdz1p0");
        if (!config.pp) v_pVertexFilterCutGplus=scalar<Int_t>(*skim,"pVertexFilterCutGplus");
        if (!config.pp) v_pVertexFilterCutVtx1=scalar<Int_t>(*skim,"pVertexFilterCutVtx1");
        if (!mc) v_HLT_PAAK4PFJet60_Eta5p1_v4=scalar<Int_t>(*hlt,"HLT_PAAK4PFJet60_Eta5p1_v4");
        if (!mc) v_HLT_PAAK4PFJet80_Eta5p1_v3=scalar<Int_t>(*hlt,"HLT_PAAK4PFJet80_Eta5p1_v3");
        if (!mc) v_HLT_PAAK4PFJet100_Eta5p1_v3=scalar<Int_t>(*hlt,"HLT_PAAK4PFJet100_Eta5p1_v3");
        if (events->GetTree()->GetBranch("rho")) {rhoName="rho";rho=scalar<Float_t>(*events,"rho");}
        else if (events->GetTree()->GetBranch("fixedGridRhoFastjetAll")) {rhoName="fixedGridRhoFastjetAll";rho=scalar<Float_t>(*events,"fixedGridRhoFastjetAll");}
        else std::cerr<<"[WARNING] No event-rho branch found; deterministic JER hash uses rho=0.\n";
        if (config.readRecoJets || (mc && config.readGenJets)) jets=open(config.jetTree);
        if (config.readRecoJets) {
            nref=scalar<Int_t>(*jets,"nref");
            a_rawpt=array<Float_t>(*jets,"rawpt","nref");
            a_trackMax=array<Float_t>(*jets,"trackMax","nref");
            a_jteta=array<Float_t>(*jets,"jteta","nref");
            a_jtphi=array<Float_t>(*jets,"jtphi","nref");
            a_jtPfNHF=array<Float_t>(*jets,"jtPfNHF","nref");
            a_jtPfNEF=array<Float_t>(*jets,"jtPfNEF","nref");
            a_jtPfCHF=array<Float_t>(*jets,"jtPfCHF","nref");
            a_jtPfMUF=array<Float_t>(*jets,"jtPfMUF","nref");
            a_jtPfCEF=array<Float_t>(*jets,"jtPfCEF","nref");
            a_jtPfCHM=array<Int_t>(*jets,"jtPfCHM","nref");
            a_jtPfCEM=array<Int_t>(*jets,"jtPfCEM","nref");
            a_jtPfNHM=array<Int_t>(*jets,"jtPfNHM","nref");
            a_jtPfNEM=array<Int_t>(*jets,"jtPfNEM","nref");
            a_jtPfMUM=array<Int_t>(*jets,"jtPfMUM","nref");
            if (mc) a_refpt=array<Float_t>(*jets,"refpt","nref");
            if (mc) a_refeta=array<Float_t>(*jets,"refeta","nref");
            if (mc) a_refphi=array<Float_t>(*jets,"refphi","nref");
        }
        if (mc && config.readGenJets) {
            ngen=scalar<Int_t>(*jets,"ngen");
            a_genpt=array<Float_t>(*jets,"genpt","ngen");
            a_geneta=array<Float_t>(*jets,"geneta","ngen");
            a_genphi=array<Float_t>(*jets,"genphi","ngen");
        }
        if (config.readTracks) {
            tracks=open(config.trackTree); nTrk=scalar<Int_t>(*tracks,"nTrk");
            a_trkPt=array<Float_t>(*tracks,"trkPt","nTrk");
            a_trkEta=array<Float_t>(*tracks,"trkEta","nTrk");
            a_trkPhi=array<Float_t>(*tracks,"trkPhi","nTrk");
            a_trkCharge=array<Int_t>(*tracks,"trkCharge","nTrk");
            a_highPurity=array<Bool_t>(*tracks,"highPurity","nTrk");
        }
        if (config.readParticles) {
            if (!mc) throw std::runtime_error("Generator particles require MC input");
            particles=open(config.particleTree);
            a_pt=array<Float_t>(*particles,"pt",nullptr);
            a_eta=array<Float_t>(*particles,"eta",nullptr);
            a_phi=array<Float_t>(*particles,"phi",nullptr);
            a_chg=array<Int_t>(*particles,"chg",nullptr);
            a_pdg=array<Int_t>(*particles,"pdg",nullptr);
            a_sube=array<Int_t>(*particles,"sube",nullptr);
        }
    }
    void read(Long64_t index,Event& event) {
        for (auto* r : {events.get(),skim.get(),hlt.get(),jets.get(),tracks.get(),particles.get()})
            if (r && r->SetEntry(index)!=TTreeReader::kEntryValid) throw std::runtime_error("Failed to read entry "+std::to_string(index));
        event.info=EventInfo{};
        const bool mc=config.mcType!=0;
        if (true) event.info.vz=get(v_vz);
        if (true) event.info.run=get(v_run);
        if (true) event.info.lumi=get(v_lumi);
        if (true) event.info.evt=get(v_evt);
        if (mc) event.info.ptHat=get(v_pthat);
        if (mc && config.pp) event.info.inputEventWeight=get(v_weight);
        if (!config.pp) event.info.pBeamScrapingFilter=get(v_pBeamScrapingFilter);
        if (!config.pp) event.info.pPAprimaryVertexFilter=get(v_pPAprimaryVertexFilter);
        if (!config.pp) event.info.HBHENoiseFilterResultRun2Loose=get(v_HBHENoiseFilterResultRun2Loose);
        if (!config.pp) event.info.phfCoincFilter=get(v_phfCoincFilter);
        if (!config.pp) event.info.pVertexFilterCutdz1p0=get(v_pVertexFilterCutdz1p0);
        if (!config.pp) event.info.pVertexFilterCutGplus=get(v_pVertexFilterCutGplus);
        if (!config.pp) event.info.pVertexFilterCutVtx1=get(v_pVertexFilterCutVtx1);
        if (!mc) event.info.HLT_PAAK4PFJet60_Eta5p1_v4=get(v_HLT_PAAK4PFJet60_Eta5p1_v4);
        if (!mc) event.info.HLT_PAAK4PFJet80_Eta5p1_v3=get(v_HLT_PAAK4PFJet80_Eta5p1_v3);
        if (!mc) event.info.HLT_PAAK4PFJet100_Eta5p1_v3=get(v_HLT_PAAK4PFJet100_Eta5p1_v3);
        event.info.eventRho=rho ? get(rho) : 0.f;
        event.sample=static_cast<SampleType>(config.mcType);
        event.hasRecoJets=config.readRecoJets; event.hasGenJets=mc && config.readGenJets;
        event.hasTracks=config.readTracks; event.hasParticles=config.readParticles;
        event.recoJets.clear(); event.genJets.clear(); event.tracks.clear(); event.particles.clear();
        if (event.hasRecoJets) {
            const Int_t n=get(nref);
            if (n<0 || n>1500) throw std::runtime_error("nref exceeds reference capacity 1500");
            size(a_rawpt,n);
            size(a_trackMax,n);
            size(a_jteta,n);
            size(a_jtphi,n);
            size(a_jtPfNHF,n);
            size(a_jtPfNEF,n);
            size(a_jtPfCHF,n);
            size(a_jtPfMUF,n);
            size(a_jtPfCEF,n);
            size(a_jtPfCHM,n);
            size(a_jtPfCEM,n);
            size(a_jtPfNHM,n);
            size(a_jtPfNEM,n);
            size(a_jtPfMUM,n);
            if (mc) size(a_refpt,n);
            if (mc) size(a_refeta,n);
            if (mc) size(a_refphi,n);
            event.recoJets.resize(n);
            for (Int_t i=0;i<n;++i) {
                auto& jet=event.recoJets[i]; jet=RecoJet{};
                jet.recoPtJeuUp=jet.recoPtJeuDown=jet.refPtSmeared=-1.f;
                jet.refPt=-1.f;
                jet.recoPtRaw=(*a_rawpt)[i];
                jet.recoTrackMaxPt=(*a_trackMax)[i];
                jet.recoEta=(*a_jteta)[i];
                jet.recoPhi=(*a_jtphi)[i];
                jet.recoPfNHF=(*a_jtPfNHF)[i];
                jet.recoPfNEF=(*a_jtPfNEF)[i];
                jet.recoPfCHF=(*a_jtPfCHF)[i];
                jet.recoPfMUF=(*a_jtPfMUF)[i];
                jet.recoPfCEF=(*a_jtPfCEF)[i];
                jet.recoPfCHM=(*a_jtPfCHM)[i];
                jet.recoPfCEM=(*a_jtPfCEM)[i];
                jet.recoPfNHM=(*a_jtPfNHM)[i];
                jet.recoPfNEM=(*a_jtPfNEM)[i];
                jet.recoPfMUM=(*a_jtPfMUM)[i];
                if (mc) jet.refPt=(*a_refpt)[i];
                if (mc) jet.refEta=(*a_refeta)[i];
                if (mc) jet.refPhi=(*a_refphi)[i];
            }
        }
        if (event.hasGenJets) {
            const Int_t n=get(ngen);
            if (n<0 || n>1500) throw std::runtime_error("ngen exceeds reference capacity 1500");
            size(a_genpt,n); size(a_geneta,n); size(a_genphi,n);
            event.genJets.resize(n);
            for(Int_t i=0;i<n;++i) event.genJets[i]={(*a_genpt)[i],(*a_geneta)[i],(*a_genphi)[i],-1.f};
        }
        if (event.hasTracks) {
            const Int_t n=get(nTrk); if(n<0) throw std::runtime_error("Negative track count");
            size(a_trkPt,n); size(a_trkEta,n); size(a_trkPhi,n); size(a_trkCharge,n); size(a_highPurity,n);
            event.tracks.resize(n);
            for(Int_t i=0;i<n;++i) event.tracks[i]={(*a_trkPt)[i],(*a_trkEta)[i],(*a_trkPhi)[i],(*a_trkCharge)[i],bool((*a_highPurity)[i])};
        }
        if(event.hasParticles) {
            const auto n=a_pt->GetSize();
            size(a_pt,n); size(a_eta,n); size(a_phi,n); size(a_chg,n); size(a_pdg,n); size(a_sube,n);
            event.particles.resize(n);
            for(std::size_t i=0;i<n;++i) event.particles[i]={(*a_pt)[i],(*a_eta)[i],(*a_phi)[i],(*a_chg)[i],(*a_pdg)[i],(*a_sube)[i]};
        }
    }
};
} // namespace
struct CMSForestReader::Impl {
    ReaderConfig config;
    std::vector<std::string> paths;
    std::unique_ptr<FileState> state;
    std::int64_t total=0, local=0;
    std::size_t next=0;
    explicit Impl(const ReaderConfig& c):config(c) {
        std::vector<std::string> candidates;
        if (c.mcType<0 || c.mcType>2) throw std::runtime_error("Invalid MC type");
        if (c.input.size()>=5 && c.input.substr(c.input.size()-5)==".root") candidates.push_back(c.input);
        else {
            std::ifstream stream(c.input); if(!stream) throw std::runtime_error("Cannot open input list: "+c.input);
            std::string line,path;
            while(std::getline(stream,line)) {
                std::istringstream fields(line); if(!(fields>>path) || path[0]=='#') continue;
                if(path.find("Forest")==std::string::npos || path.find(".root")==std::string::npos)
                    throw std::runtime_error("Invalid HiForest list entry: "+line);
                candidates.push_back(path);
            }
        }
        std::optional<std::string> rhoSchema;
        for(const auto& path:candidates) {
            // Match the reference's skip policy for unreadable/incomplete forests;
            // schema mismatches remain fatal rather than silently losing files.
            std::unique_ptr<TFile> probe(TFile::Open(path.c_str(),"READ"));
            bool usable=probe && !probe->IsZombie() && probe->GetNkeys()>0;
            if(usable) for(const auto& name: std::vector<std::string>{"hiEvtAnalyzer/HiTree","hltanalysis/HltTree","skimanalysis/HltTree"})
                usable=usable && probe->Get<TTree>(name.c_str());
            if(usable && (config.readRecoJets || (config.mcType && config.readGenJets)))
                usable=probe->Get<TTree>(config.jetTree.c_str());
            if(!usable) {std::cerr<<"SKIPPING INPUT FILE: unreadable or incomplete forest: "<<path<<'\n';continue;}
            probe.reset();
            FileState checked(path,config);
            if (rhoSchema && *rhoSchema!=checked.rhoName)
                throw std::runtime_error("Event-rho branch layout changes between files; refusing ambiguous legacy friend-chain behavior");
            rhoSchema=checked.rhoName;
            total+=checked.entries; paths.push_back(path);
        }
        std::cout<<"Input file summary: accepted="<<paths.size()<<" skipped="<<candidates.size()-paths.size()<<'\n';
        if(paths.empty() || total<=0) throw std::runtime_error("No usable input ROOT entries remain");
    }
    bool read(Event& event) {
        while(!state || local==state->entries) {
            if(next==paths.size()) return false;
            state=std::make_unique<FileState>(paths[next++],config); local=0;
        }
        state->read(local++,event); return true;
    }
};
CMSForestReader::CMSForestReader(const ReaderConfig& c):impl_(std::make_unique<Impl>(c)) {}
CMSForestReader::~CMSForestReader()=default;
std::int64_t CMSForestReader::entries() const {return impl_->total;}
bool CMSForestReader::readNext(Event& event) {return impl_->read(event);}
} // namespace jetanalysis
