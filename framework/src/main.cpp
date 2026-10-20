#include "jetanalysis/CMSForestReader.h"
#include "jetanalysis/DiJetAnalysis.h"
#include "jetanalysis/AnalysisManager.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <regex>
namespace {
void usage() {
    std::cout << "Usage: jetAnalysis INPUT OUTPUT MC_TYPE DIRECTION PTHAT TRIGGER JET_ID VERTEX [pPb|pp] [JER_MODE] [ETA_JER]\n"
        "Options: --jet-collection AK4PF|AK4CSPF|AK5PF|... --jet-tree PATH\n"
        " --jec FILE (repeatable) --jeu FILE --aux-directory DIR\n"
        " --read-tracks --read-particles --track-tree PATH --particle-tree PATH\n"
        " --no-reco-jets --no-gen-jets --reader-only --profile-stages --exclude-output PREFIX (repeatable)\n";
}
int number(const std::string& s) {
    std::size_t used=0; int value=std::stoi(s,&used);
    if(used!=s.size()) throw std::runtime_error("Invalid integer: "+s);
    return value;
}
}
int main(int argc,char** argv) {
    if(argc==2 && std::string(argv[1])=="--help") {usage();return 0;}
    try {
        if(argc<9) {usage();return 2;}
        jetanalysis::AnalysisConfig c;
        c.reader.input=argv[1]; c.output=argv[2]; c.reader.mcType=number(argv[3]);
        c.isPbGoingDir=number(argv[4]); c.ptHatSample=number(argv[5]); c.triggerId=number(argv[6]);
        c.recoJetSelMethod=number(argv[7]); c.vertexFilterSelection=number(argv[8]);
        int i=9;
        if(i<argc && std::string(argv[i]).rfind("--",0)!=0) {
            const std::string system=argv[i++];
            if(system!="pp" && system!="pPb") throw std::runtime_error("Expected pPb or pp");
            c.reader.pp=system=="pp";
        }
        if(i<argc && std::string(argv[i]).rfind("--",0)!=0) c.jerSmearingMode=number(argv[i++]);
        if(i<argc && std::string(argv[i]).rfind("--",0)!=0) c.useEtaDependentJER=number(argv[i++]);
        bool readerOnly=false;
        for(;i<argc;++i) {
            const std::string flag=argv[i];
            auto value=[&]() {if(++i>=argc) throw std::runtime_error("Missing value for "+flag);return std::string(argv[i]);};
            if(flag=="--jet-tree") c.reader.jetTree=value();
            else if(flag=="--jet-collection") {
                std::string n=value();
                std::smatch match;
                if(!std::regex_match(n,match,std::regex("AK([0-9]+)(CS)?PF")))
                    throw std::runtime_error("Expected AK4PF, AK4CSPF, AK5PF, etc.; use --jet-tree for an explicit tree path");
                c.reader.jetTree="ak"+std::string(match[2].matched ? "Cs" : "")+match[1].str()+"PFJetAnalyzer/t";
            }
            else if(flag=="--jec") c.jecFiles.push_back(value());
            else if(flag=="--jeu") c.jeuFile=value();
            else if(flag=="--aux-directory") c.auxDirectory=value();
            else if(flag=="--exclude-output") {auto p=value();if(p.empty()) throw std::runtime_error("Empty output prefix");c.excludeOutputPrefixes.push_back(p);}
            else if(flag=="--track-tree") c.reader.trackTree=value();
            else if(flag=="--particle-tree") c.reader.particleTree=value();
            else if(flag=="--read-tracks") c.reader.readTracks=true;
            else if(flag=="--read-particles") c.reader.readParticles=true;
            else if(flag=="--no-reco-jets") c.reader.readRecoJets=false;
            else if(flag=="--no-gen-jets") c.reader.readGenJets=false;
            else if(flag=="--profile-stages") c.profileStages=true;
            else if(flag=="--reader-only") readerOnly=true;
            else throw std::runtime_error("Unknown option: "+flag);
        }
        if(c.reader.mcType<0 || c.reader.mcType>2 || c.isPbGoingDir<0 || c.isPbGoingDir>1 ||
           c.triggerId<0 || c.triggerId>3 || c.recoJetSelMethod<0 || c.recoJetSelMethod>3 ||
           c.vertexFilterSelection<0 || c.vertexFilterSelection>2 || c.jerSmearingMode<0 || c.jerSmearingMode>2 ||
           c.useEtaDependentJER<0 || c.useEtaDependentJER>1) throw std::runtime_error("Analysis argument outside supported range");
        if(!readerOnly && c.reader.jetTree!="ak4PFJetAnalyzer/t")
            throw std::runtime_error("DiJetAnalysis is restricted to AK4PF until collection-specific JEC/JER are approved; use --reader-only");
        if(!readerOnly && (!c.reader.readRecoJets || (c.reader.mcType && !c.reader.readGenJets)))
            throw std::runtime_error("DiJetAnalysis requires reco jets and, for MC, gen jets; use --reader-only for other collection combinations");
        if(c.reader.input.empty()) {
            const std::string gen=c.reader.mcType==2 ? "pythia" : "embedding";
            const std::string dir=c.isPbGoingDir ? "Pbgoing" : "pgoing";
            const std::string tag=c.reader.mcType==2 ? "unembedded" : "embedded";
            c.reader.input="~/cernbox/ana/pPb8160/"+gen+"/"+dir+"/forest/HiForestSkim_pPb_MC_pthat"+
                std::to_string(c.ptHatSample)+"_"+dir+"_"+tag+".root";
        }
        jetanalysis::CMSForestReader reader(c.reader);
        if(readerOnly) {
            jetanalysis::Event e; std::uint64_t events=0,reco=0,gen=0,tracks=0,particles=0;
            while(reader.readNext(e)) {++events;reco+=e.recoJets.size();gen+=e.genJets.size();tracks+=e.tracks.size();particles+=e.particles.size();}
            std::cout<<"Read "<<events<<" events, reco="<<reco<<" gen="<<gen<<" tracks="<<tracks<<" particles="<<particles<<'\n';
        } else {
            jetanalysis::AnalysisManager manager(reader,c.profileStages);
            manager.addTask(std::make_unique<jetanalysis::DiJetAnalysis>(std::move(c)));
            manager.run();
        }
        return 0;
    } catch(const std::exception& error) {std::cerr<<"jetAnalysis: "<<error.what()<<'\n';return 1;}
}
