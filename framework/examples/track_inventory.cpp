// A compilable nonphysics example of an external task using only the core library.
#include <jetanalysis/AnalysisManager.h>
#include <jetanalysis/CMSForestReader.h>
#include <iostream>
#include <stdexcept>
class TrackInventory final : public jetanalysis::AnalysisTask {
    std::uint64_t events_=0,tracks_=0;
public:
    void initialize(const jetanalysis::RunContext&) override {}
    void process(const jetanalysis::Event& event) override {
        if(!event.hasTracks)throw std::runtime_error("TrackInventory requires tracks");
        ++events_;tracks_+=event.tracks.size();
    }
    void finish() override {std::cout<<events_<<" events, "<<tracks_<<" tracks\n";}
};
int main(int argc,char** argv) {
    if(argc!=2) {std::cerr<<"Usage: track_inventory INPUT.root\n";return 2;}
    try {
        jetanalysis::ReaderConfig config;
        config.input=argv[1];config.mcType=0;config.readRecoJets=false;config.readGenJets=false;config.readTracks=true;
        jetanalysis::CMSForestReader reader(config);
        jetanalysis::AnalysisManager manager(reader);
        manager.addTask(std::make_unique<TrackInventory>());
        manager.run();
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
