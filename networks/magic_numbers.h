#ifndef MAGIC_NUMBERS_H_
#define MAGIC_NUMBERS_H_

// Magic numbers

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <future>
#include <numbers>
#include <string_view>
#include <tuple>
#include <vector>

//! This is very ugly and should be changed.
constexpr std::string_view kDataDirectory = "Uni/Exponential_Networks_cpp";

// Numerical constants.
const double pi = std::numbers::pi;
constexpr std::complex<double> J = std::complex<double>(0, 1);
const std::complex<double> kDummyRoot = std::complex<double>(1, 0.3);
const std::complex<double> kZeta3 = std::exp(2.0 * pi * J / 3.0);
const double kNumOffset = 1e-10;

// Precision constants.
const double kZerosPrecisions = 1e-20;
const double kOdeAbsError = 1e-13;
const double kOdeRelError = 1e-13;
const double kFiberCompTolerance = 1e-8;

// Iteration constants.
const uint32_t kInitialSteps = 1000;
const uint32_t kZerosMaxIterations = 300;
const uint32_t kLineStepsPerUnit = 30;

// Map constants.
const uint32_t kMapResolutionReal = 1000;
const uint32_t kMapResolutionImag = 1000;
const std::array<double, 2> kMapRangeReal = {-5.0, 5.0};
const std::array<double, 2> kMapRangeImag = {-4.0, 4.0};

//D0-D4 constants.
const double kD0Mass = 4*pi*pi;
const double kD4angle = -0.01;
const double kD4Cutoff = 201.9;

// ODE integration constants.
const double kIntegratePeriod = 5.0;
const double kInitialStepSize = 1e-8;
const double kCutoff = 200.0;
const double kLowCutoffX = 1e-17;
const double kHighCutoffX = 1000;
const double kCutoffY = 10000;
const uint32_t kMaxSteps = 100000;

// Datapoint index constants.
constexpr uint32_t kIndexLowerBound = 0;
constexpr uint32_t kIndexUpperBound = 1;

#endif  // MAGIC_NUMBERS_H_
