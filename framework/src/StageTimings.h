#pragma once
#include <array>
#include <chrono>
#include <iostream>
namespace jetanalysis {
// Opt-in instrumentation: disabled runs do not read the clock per event.
struct StageTimings {
    enum Stage {Read,RecoLoad,GenLoad,GenJets,GenDijets,RefDijets,RecoJets,RecoDijets,Output,Count};
    bool enabled=false;
    std::array<double,Count> seconds{};
    template<class F> void measure(Stage stage,F operation) {
        if(!enabled) {operation();return;}
        const auto start=std::chrono::steady_clock::now();
        operation();
        seconds[stage]+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    }
    void report() const {
        if(!enabled)return;
        const char* names[]={"read_and_convert","reco_load_and_jec","gen_load_sort","gen_jets","gen_dijets","reference_dijets","reco_jets","reco_dijets_and_responses","output"};
        for(int i=1;i<Count;++i) std::cout<<"STAGE "<<names[i]<<" "<<seconds[i]<<" seconds\n";
    }
};
}
