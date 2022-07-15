#ifndef MAGIC_NUMBERS_H_
#define MAGIC_NUMBERS_H_

// Magic numbers

#include <cstdint>
#include <string_view>
#include <complex>
#include <numbers> 
#include <vector>
#include <cmath>
#include <array>
#include <tuple>
#include <future>

using namespace std::complex_literals;

//! This is very ugly and should be changed.
constexpr std::string_view kDataDirectory = "Uni/Exponential_Networks";

// Numerical constants.
const double pi = std::numbers::pi;
constexpr std::complex<double> J = std::complex<double>(0,1);
const std::complex<double> kDummyRoot = std::complex<double>(1,0.3);
const std::complex<double> kZeta3 = std::exp(2.0 *pi*J / 3.0);
const double kNumOffset = 1e-10;

// Precision constants.
const double kZerosPrecisions = 1e-15;
const double kOdeAbsError = 1e-13;
const double kOdeRelError = 1e-13;
const double kFiberCompTolerance = 1e-8;

// Iteration constants.
const uint32_t kInitialSteps = 500;
const uint32_t kZerosMaxIterations = 200;
const uint32_t kLineStepsPerUnit = 20;

// Map constants.
const uint32_t kMapResolutionReal = 3000;
const uint32_t kMapResolutionImag = 3000;
const std::array<double,2> kMapRangeReal = {-5.0, 5.0};
const std::array<double,2> kMapRangeImag = {-4.0, 4.0};


// ODE integration constants.
const double kIntegratePeriod = 10.0;
const double kInitialStepSize = 1e-6;
const double kCutoff = 200.0;
const uint32_t kMaxSteps = 100000;

// Datapoint index constants.
constexpr uint32_t kIndexLowerBound = 0;
constexpr uint32_t kIndexUpperBound = 1;

#endif  // MAGIC_NUMBERS_H_
