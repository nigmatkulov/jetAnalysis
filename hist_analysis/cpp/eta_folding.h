#ifndef JETANALYSIS_ETA_FOLDING_H
#define JETANALYSIS_ETA_FOLDING_H

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TString.h>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace eta_folding {
// Match python/forward_backward.py: mirror complete bins, never interpolate.
inline std::unique_ptr<TH1> fold(const TH1 &full, bool backward, const char *name) {
    const auto *two = dynamic_cast<const TH2 *>(&full);
    if (full.GetDimension() != 1 && full.GetDimension() != 2)
        throw std::invalid_argument("Eta folding requires TH1 or TH2");
    const TAxis &eta = two ? *two->GetYaxis() : *full.GetXaxis();
    const int n = eta.GetNbins(), half = n / 2;
    if (n % 2 || std::abs(eta.GetBinLowEdge(half+1)) > 1.e-9)
        throw std::invalid_argument("Eta folding requires zero as a bin edge");
    for (int i = 1; i <= n+1; ++i)
        if (std::abs(eta.GetBinLowEdge(i) + eta.GetBinLowEdge(n+2-i)) > 1.e-9)
            throw std::invalid_argument("Eta folding requires symmetric bins");
    std::vector<double> edges(half+1);
    for (int i = 1; i <= half; ++i) edges[i] = eta.GetBinUpEdge(half+i);
    std::unique_ptr<TH1> result;
    if (two) {
        std::vector<double> x(two->GetNbinsX()+1);
        for (int i = 0; i <= two->GetNbinsX(); ++i) x[i] = two->GetXaxis()->GetBinLowEdge(i+1);
        result = std::make_unique<TH2D>(name, full.GetTitle(), two->GetNbinsX(), x.data(), half, edges.data());
    } else result = std::make_unique<TH1D>(name, full.GetTitle(), half, edges.data());
    result->SetDirectory(nullptr);
    result->Sumw2();
    result->GetXaxis()->SetTitle(two ? full.GetXaxis()->GetTitle() : "|#eta_{CM}|");
    result->GetYaxis()->SetTitle(two ? "|#eta_{CM}|" : full.GetYaxis()->GetTitle());
    for (int bin = 1; bin <= half+1; ++bin) {
        const int source = backward ? (bin <= half ? half-bin+1 : 0)
                                    : (bin <= half ? half+bin : n+1);
        const int lastX = two ? two->GetNbinsX()+1 : 0;
        for (int x = 0; x <= lastX; ++x) {
            const int from = two ? full.GetBin(x, source) : source;
            const int to = two ? result->GetBin(x, bin) : bin;
            result->SetBinContent(to, full.GetBinContent(from));
            result->SetBinError(to, full.GetBinError(from));
        }
    }
    return result;
}

// Legacy plotting macros keep their local hemisphere variables, but always
// derive them from full inputs, including when opening older ROOT files.
inline TH1 *get(TFile *file, TString key) {
    const bool backward = key.Contains("Backward");
    if (!backward && !key.Contains("Forward")) return dynamic_cast<TH1 *>(file->Get(key));
    TString fullKey = key;
    fullKey.ReplaceAll("PtEtaForward", "PtEtaCM");
    fullKey.ReplaceAll("PtEtaBackward", "PtEtaCM");
    fullKey.ReplaceAll("EtaCMForwardShifted", "EtaCMShifted");
    fullKey.ReplaceAll("EtaCMBackwardShifted", "EtaCMShifted");
    fullKey.ReplaceAll("CMArr_", "CM_");
    auto *full = dynamic_cast<TH1 *>(file->Get(fullKey));
    if (!full) return nullptr;
    auto result = fold(*full, backward, (key + "_from_full").Data());
    result->SetDirectory(file);
    return result.release();
}
} // namespace eta_folding
#endif
