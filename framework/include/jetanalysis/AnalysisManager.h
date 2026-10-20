#pragma once
#include "AnalysisTask.h"
#include "EventReader.h"
#include <memory>
#include <vector>
namespace jetanalysis {
struct RunSummary {
    std::int64_t processedEvents=0;
    double readSeconds=0; // Measured only when profiling is enabled.
};
class AnalysisManager {
public:
    // Reader is borrowed and must outlive run(); tasks are exclusively owned.
    explicit AnalysisManager(EventReader& reader, bool profileReads=false);
    void addTask(std::unique_ptr<AnalysisTask> task);
    RunSummary run();
private:
    EventReader& reader_;
    bool profileReads_, started_=false;
    std::vector<std::unique_ptr<AnalysisTask>> tasks_;
};
// Single borrowed-task convenience, using the same lifecycle and event loop.
RunSummary runTask(EventReader&, AnalysisTask&, bool profileReads=false);
} // namespace jetanalysis
