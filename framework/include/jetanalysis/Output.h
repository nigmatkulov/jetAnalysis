#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <TObject.h>
namespace jetanalysis {
struct Histograms;
struct OutputSelection {
    std::vector<std::string> excludedPrefixes;
    void write(TObject* object) const {
        const std::string name=object->GetName();
        for (const auto& prefix : excludedPrefixes)
            if (name.compare(0, prefix.size(), prefix)==0) return;
        if (object->Write() <= 0) throw std::runtime_error("Failed to write ROOT object: "+name);
    }
};
void writeOutput(const std::string&, Histograms&, bool, int, bool, const OutputSelection&);
}
