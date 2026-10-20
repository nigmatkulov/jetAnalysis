#pragma once
#include "Config.h"
#include "EventReader.h"
#include <memory>
namespace jetanalysis {
class CMSForestReader final : public EventReader {
public:
    explicit CMSForestReader(const ReaderConfig&);
    ~CMSForestReader() override;
    std::int64_t entries() const override;
    bool readNext(Event&) override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace jetanalysis
