#include "jetanalysis/AnalysisManager.h"
#include "jetanalysis/CMSForestReader.h"
#include "jetanalysis/DiJetAnalysis.h"
#include <stdexcept>
#include <iostream>
class InputObserver final : public jetanalysis::AnalysisTask {
    std::int64_t count_=0,expected_=0;
public:
    void initialize(const jetanalysis::RunContext& c) override {expected_=c.entries;}
    void process(const jetanalysis::Event& e) override {
        const auto i=e.info.evt-((std::uint64_t{1}<<63)+17);
        if(e.info.vz!=(i==1 ? 16.f : -2.f) || !e.hasTracks || e.tracks.size()!=2 || e.tracks[0].pt!=1.5f)
            throw std::runtime_error("An earlier task changed the shared input event");
        ++count_;
    }
    void finish() override {
        if(count_!=expected_)throw std::runtime_error("Observer missed events rejected by DiJetAnalysis");
    }
};
int main(int argc,char** argv) {
    if(argc!=4)return 2;
    try {
        jetanalysis::AnalysisConfig c;
        c.reader.input=argv[1];c.reader.mcType=1;c.reader.readTracks=true;
        c.isPbGoingDir=1;c.ptHatSample=540;c.vertexFilterSelection=1;c.jerSmearingMode=1;c.useEtaDependentJER=1;
        jetanalysis::CMSForestReader reader(c.reader);
        jetanalysis::AnalysisManager manager(reader);
        c.output=argv[2];manager.addTask(std::make_unique<jetanalysis::DiJetAnalysis>(c));
        c.output=argv[3];manager.addTask(std::make_unique<jetanalysis::DiJetAnalysis>(c));
        manager.addTask(std::make_unique<InputObserver>());
        manager.run();
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
