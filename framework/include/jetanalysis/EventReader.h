#pragma once
#include "Event.h"
#include <cstdint>
namespace jetanalysis {
class EventReader {
public:
    virtual ~EventReader() = default;
    virtual std::int64_t entries() const = 0;
    virtual bool readNext(Event&) = 0;
};
} // namespace jetanalysis
