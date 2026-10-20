#pragma once
#include "AnalysisTask.h"
#include "Config.h"
#include "EventReader.h"
#include <memory>
namespace jetanalysis {
class DiJetAnalysis final : public AnalysisTask {
public:
    explicit DiJetAnalysis(AnalysisConfig config);
    ~DiJetAnalysis() override;
    void initialize(const RunContext&) override;
    void process(const Event&) override;
    void finish() override;
    void run(EventReader& reader); // Compatible single-task application shortcut.
private:
    AnalysisConfig config_;
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace jetanalysis
