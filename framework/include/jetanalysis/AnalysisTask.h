#pragma once
#include "Event.h"
#include <cstdint>
namespace jetanalysis {
struct RunContext {
    std::int64_t entries=0; // Complete invocation count, before task-specific cuts.
};
// Tasks receive the same immutable input event. Each owns its corrections,
// selections, derived collections, output, and any event-mixing history.
class AnalysisTask {
public:
    virtual ~AnalysisTask() = default;
    virtual void initialize(const RunContext&) = 0;
    virtual void process(const Event&) = 0;
    virtual void finish() = 0;
};
} // namespace jetanalysis
