#include <jetanalysis/AnalysisManager.h>
#include <jetanalysis/CMSForestReader.h>
#include <jetanalysis/DiJetAnalysis.h>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    if(argc!=3 && argc!=4) {std::cerr<<"Usage: dijet_task INPUT.root OUTPUT.root [AUX_DIRECTORY]\n";return 2;}
    try {
        jetanalysis::AnalysisConfig config;
        config.reader.input=argv[1];config.output=argv[2];
        if(argc==4)config.auxDirectory=argv[3];
        config.reader.mcType=1;config.isPbGoingDir=0;config.ptHatSample=540;
        jetanalysis::CMSForestReader reader(config.reader);
        jetanalysis::AnalysisManager manager(reader);
        manager.addTask(std::make_unique<jetanalysis::DiJetAnalysis>(config));
        manager.run();
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
