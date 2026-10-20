#include "jetanalysis/AnalysisManager.h"
#include <chrono>
#include <iostream>
#include <stdexcept>
namespace jetanalysis {
namespace {
RunSummary execute(EventReader& reader,const std::vector<AnalysisTask*>& tasks,bool profile) {
    if(tasks.empty()) throw std::invalid_argument("AnalysisManager requires at least one task");
    const RunContext context{reader.entries()};
    if(context.entries<0) throw std::runtime_error("Negative reader entry count");
    for(auto* task:tasks) task->initialize(context);
    RunSummary result;
    Event event;
    for(std::int64_t i=0;i<context.entries;++i) {
        bool loaded;
        if(profile) {
            const auto start=std::chrono::steady_clock::now();
            loaded=reader.readNext(event);
            result.readSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        } else loaded=reader.readNext(event);
        if(!loaded) throw std::runtime_error("Reader ended before its declared entry count");
        for(auto* task:tasks) task->process(event);
        ++result.processedEvents;
    }
    // Failed runs do not invoke normal finish/output callbacks.
    for(auto* task:tasks) task->finish();
    if(profile) std::cout<<"STAGE read_and_convert "<<result.readSeconds<<" seconds\n";
    return result;
}
}
AnalysisManager::AnalysisManager(EventReader& reader,bool profile):reader_(reader),profileReads_(profile) {}
void AnalysisManager::addTask(std::unique_ptr<AnalysisTask> task) {
    if(started_) throw std::logic_error("Cannot add tasks after a run starts");
    if(!task) throw std::invalid_argument("Null analysis task");
    tasks_.push_back(std::move(task));
}
RunSummary AnalysisManager::run() {
    if(started_) throw std::logic_error("AnalysisManager is single-use");
    started_=true;
    std::vector<AnalysisTask*> tasks;
    for(auto& task:tasks_) tasks.push_back(task.get());
    return execute(reader_,tasks,profileReads_);
}
RunSummary runTask(EventReader& reader,AnalysisTask& task,bool profile) {
    return execute(reader,{&task},profile);
}
} // namespace jetanalysis
