#ifndef JETANALYSIS_HIST_ANALYSIS_REGULARIZATION_SCAN_HPP
#define JETANALYSIS_HIST_ANALYSIS_REGULARIZATION_SCAN_HPP

#include <cstddef>
#include <vector>

namespace jetanalysis::regularization {

struct ScanInput {
  std::size_t nMeasuredBins{};
  std::size_t nTruthBins{};
  std::size_t nToys{};
  std::size_t maxIterations{};

  // Row-major P(measured i | truth j), indexed as i*nTruthBins + j.
  // Every truth column must have a positive sum. Columns are normalized by the
  // kernel, allowing callers to pass either probabilities or migration counts.
  std::vector<double> response;

  // RooUnfold Bayesian training prior and external factorized efficiency.
  std::vector<double> trainingPrior;
  std::vector<double> truthEfficiency;

  // Row-major toy contents, indexed as toy*nMeasuredBins + measured bin.
  std::vector<double> measuredToys;
};

struct ScanMoments {
  std::size_t nTruthBins{};
  std::size_t nToys{};
  std::size_t completedIterations{};

  // Row-major iteration/truth-bin arrays. These are sufficient to reconstruct
  // the ensemble mean and unbiased variance without retaining unfolded toys.
  std::vector<double> sums;
  std::vector<double> sumsSquared;
};

// Advance every toy through all Bayesian iterations once. Independent toy
// chunks run in worker threads over plain numeric buffers; no ROOT object is
// shared between workers. workerCount=0 selects hardware_concurrency().
ScanMoments scanBayesMoments(const ScanInput& input,
                             std::size_t workerCount = 0);

}  // namespace jetanalysis::regularization

extern "C" {

// Stable C ABI used by the Python ctypes bridge. Returns zero on success and
// copies a diagnostic into errorBuffer on failure.
int jetanalysis_regularization_scan_moments(
    std::size_t nMeasuredBins, std::size_t nTruthBins, std::size_t nToys,
    std::size_t maxIterations, const double* response,
    const double* trainingPrior, const double* truthEfficiency,
    const double* measuredToys, std::size_t workerCount, double* sums,
    double* sumsSquared, char* errorBuffer, std::size_t errorBufferSize);

}

#endif
