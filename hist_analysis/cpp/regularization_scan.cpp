#include "regularization_scan.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

namespace jetanalysis::regularization {
namespace {

struct WorkerMoments {
  std::vector<long double> sums;
  std::vector<long double> sumsSquared;
};

void requireFiniteNonnegative(const std::vector<double>& values,
                              const char* label) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index]) || values[index] < 0.0) {
      throw std::invalid_argument(std::string(label) + " contains an invalid value at " +
                                  std::to_string(index));
    }
  }
}

void validate(const ScanInput& input) {
  if (input.nMeasuredBins == 0 || input.nTruthBins == 0 || input.nToys < 2 ||
      input.maxIterations == 0) {
    throw std::invalid_argument(
        "positive bin/iteration counts and at least two toys are required");
  }
  if (input.response.size() != input.nMeasuredBins * input.nTruthBins) {
    throw std::invalid_argument("response dimensions do not match ScanInput");
  }
  if (input.trainingPrior.size() != input.nTruthBins ||
      input.truthEfficiency.size() != input.nTruthBins) {
    throw std::invalid_argument("truth-vector dimensions do not match ScanInput");
  }
  if (input.measuredToys.size() != input.nToys * input.nMeasuredBins) {
    throw std::invalid_argument("toy dimensions do not match ScanInput");
  }
  requireFiniteNonnegative(input.response, "response");
  requireFiniteNonnegative(input.trainingPrior, "trainingPrior");
  requireFiniteNonnegative(input.measuredToys, "measuredToys");
  double priorTotal = 0.0;
  for (std::size_t truth = 0; truth < input.nTruthBins; ++truth) {
    priorTotal += input.trainingPrior[truth];
    const double efficiency = input.truthEfficiency[truth];
    if (!std::isfinite(efficiency) || efficiency < 0.0) {
      throw std::invalid_argument("truthEfficiency must be finite and nonnegative");
    }
    double columnTotal = 0.0;
    for (std::size_t measured = 0; measured < input.nMeasuredBins; ++measured) {
      columnTotal += input.response[measured * input.nTruthBins + truth];
    }
    if (!(columnTotal > 0.0) &&
        (input.trainingPrior[truth] != 0.0 || efficiency != 0.0)) {
      throw std::invalid_argument(
          "an empty response truth column requires zero prior and efficiency");
    }
  }
  if (!(priorTotal > 0.0)) {
    throw std::invalid_argument("trainingPrior must have positive total content");
  }
}

std::vector<double> normalizedResponse(const ScanInput& input) {
  std::vector<double> result(input.response);
  for (std::size_t truth = 0; truth < input.nTruthBins; ++truth) {
    double total = 0.0;
    for (std::size_t measured = 0; measured < input.nMeasuredBins; ++measured) {
      total += result[measured * input.nTruthBins + truth];
    }
    if (total > 0.0) {
      for (std::size_t measured = 0; measured < input.nMeasuredBins; ++measured) {
        result[measured * input.nTruthBins + truth] /= total;
      }
    }
  }
  return result;
}

void normalize(std::vector<double>& values) {
  double total = 0.0;
  for (const double value : values) total += value;
  if (!(total > 0.0) || !std::isfinite(total)) {
    throw std::runtime_error("Bayesian prior acquired a non-positive total");
  }
  for (double& value : values) value /= total;
}

void scanToyRange(const ScanInput& input, const std::vector<double>& response,
                  const std::size_t firstToy, const std::size_t lastToy,
                  WorkerMoments& output) {
  const std::size_t nTruth = input.nTruthBins;
  const std::size_t nMeasured = input.nMeasuredBins;
  std::vector<double> prior(input.trainingPrior);
  std::vector<double> denominator(nMeasured);
  std::vector<double> unfolded(nTruth);

  for (std::size_t toy = firstToy; toy < lastToy; ++toy) {
    prior = input.trainingPrior;
    normalize(prior);
    const double* measuredToy = input.measuredToys.data() + toy * nMeasured;

    for (std::size_t iteration = 0; iteration < input.maxIterations; ++iteration) {
      std::fill(denominator.begin(), denominator.end(), 0.0);
      for (std::size_t measured = 0; measured < nMeasured; ++measured) {
        const double* row = response.data() + measured * nTruth;
        for (std::size_t truth = 0; truth < nTruth; ++truth) {
          denominator[measured] += row[truth] * prior[truth];
        }
      }

      std::fill(unfolded.begin(), unfolded.end(), 0.0);
      for (std::size_t measured = 0; measured < nMeasured; ++measured) {
        if (!(denominator[measured] > 0.0)) continue;
        const double factor = measuredToy[measured] / denominator[measured];
        const double* row = response.data() + measured * nTruth;
        for (std::size_t truth = 0; truth < nTruth; ++truth) {
          unfolded[truth] += row[truth] * prior[truth] * factor;
        }
      }

      const std::size_t offset = iteration * nTruth;
      for (std::size_t truth = 0; truth < nTruth; ++truth) {
        const double efficiency = input.truthEfficiency[truth];
        if (efficiency == 0.0) {
          if (unfolded[truth] != 0.0) {
            throw std::runtime_error(
                "zero-efficiency truth bin acquired non-zero unfolded content");
          }
          continue;
        }
        const double inclusive = unfolded[truth] / efficiency;
        if (!std::isfinite(inclusive)) {
          throw std::runtime_error("unfolded toy contains a non-finite value");
        }
        output.sums[offset + truth] += inclusive;
        output.sumsSquared[offset + truth] += inclusive * inclusive;
      }
      prior = unfolded;
      normalize(prior);
    }
  }
}

}  // namespace

ScanMoments scanBayesMoments(const ScanInput& input, std::size_t workerCount) {
  validate(input);
  const std::vector<double> response = normalizedResponse(input);
  if (workerCount == 0) {
    workerCount = std::max<std::size_t>(1, std::thread::hardware_concurrency());
  }
  workerCount = std::min(workerCount, input.nToys);
  const std::size_t momentSize = input.maxIterations * input.nTruthBins;
  std::vector<WorkerMoments> workers(workerCount);
  for (WorkerMoments& worker : workers) {
    worker.sums.assign(momentSize, 0.0L);
    worker.sumsSquared.assign(momentSize, 0.0L);
  }

  std::vector<std::thread> threads;
  std::vector<std::exception_ptr> failures(workerCount);
  threads.reserve(workerCount);
  for (std::size_t worker = 0; worker < workerCount; ++worker) {
    const std::size_t first = input.nToys * worker / workerCount;
    const std::size_t last = input.nToys * (worker + 1) / workerCount;
    threads.emplace_back([&, worker, first, last] {
      try {
        scanToyRange(input, response, first, last, workers[worker]);
      } catch (...) {
        failures[worker] = std::current_exception();
      }
    });
  }
  for (std::thread& thread : threads) thread.join();
  for (const std::exception_ptr& failure : failures) {
    if (failure) std::rethrow_exception(failure);
  }

  ScanMoments result;
  result.nTruthBins = input.nTruthBins;
  result.nToys = input.nToys;
  result.completedIterations = input.maxIterations;
  result.sums.assign(momentSize, 0.0);
  result.sumsSquared.assign(momentSize, 0.0);
  // Merge in a fixed worker order for reproducibility at a fixed worker count.
  for (std::size_t index = 0; index < momentSize; ++index) {
    long double sum = 0.0L;
    long double sumSquared = 0.0L;
    for (const WorkerMoments& worker : workers) {
      sum += worker.sums[index];
      sumSquared += worker.sumsSquared[index];
    }
    result.sums[index] = static_cast<double>(sum);
    result.sumsSquared[index] = static_cast<double>(sumSquared);
  }
  return result;
}

}  // namespace jetanalysis::regularization

extern "C" int jetanalysis_regularization_scan_moments(
    const std::size_t nMeasuredBins, const std::size_t nTruthBins,
    const std::size_t nToys, const std::size_t maxIterations,
    const double* response, const double* trainingPrior,
    const double* truthEfficiency, const double* measuredToys,
    const std::size_t workerCount, double* sums, double* sumsSquared,
    char* errorBuffer, const std::size_t errorBufferSize) {
  auto reportError = [errorBuffer, errorBufferSize](const char* message) {
    if (errorBuffer != nullptr && errorBufferSize > 0) {
      std::snprintf(errorBuffer, errorBufferSize, "%s", message);
    }
  };
  try {
    if (response == nullptr || trainingPrior == nullptr ||
        truthEfficiency == nullptr || measuredToys == nullptr || sums == nullptr ||
        sumsSquared == nullptr) {
      throw std::invalid_argument("null array passed to compiled scan");
    }
    jetanalysis::regularization::ScanInput input;
    input.nMeasuredBins = nMeasuredBins;
    input.nTruthBins = nTruthBins;
    input.nToys = nToys;
    input.maxIterations = maxIterations;
    input.response.assign(response, response + nMeasuredBins * nTruthBins);
    input.trainingPrior.assign(trainingPrior, trainingPrior + nTruthBins);
    input.truthEfficiency.assign(truthEfficiency, truthEfficiency + nTruthBins);
    input.measuredToys.assign(measuredToys, measuredToys + nToys * nMeasuredBins);
    const auto result =
        jetanalysis::regularization::scanBayesMoments(input, workerCount);
    std::copy(result.sums.begin(), result.sums.end(), sums);
    std::copy(result.sumsSquared.begin(), result.sumsSquared.end(), sumsSquared);
    reportError("");
    return 0;
  } catch (const std::exception& error) {
    reportError(error.what());
    return 1;
  } catch (...) {
    reportError("unknown compiled regularization scan failure");
    return 2;
  }
}
