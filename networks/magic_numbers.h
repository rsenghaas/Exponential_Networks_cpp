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
using namespace std::numbers;

//! This is very ugly and should be changed.
constexpr std::string_view kDataDirectory = "Uni/Exponential_Networks";

// Numerical constants.
constexpr std::complex<double> J = std::complex<double>(0,1);
const std::complex<double> kDummyRoot = std::complex<double>(1,0.3);
const std::complex<double> kZeta3 = std::exp(2.0 *pi*J / 3.0);

// Precision constants.
const double kZerosPrecisions = 1e-15;
const double kOdeAbsError = 1e-13;
const double kOdeRelError = 1e-13;

// Iteration constants.
const uint32_t kInitialEulerSteps = 100;
const uint32_t kZerosMaxIterations = 200;

// ODE integration constants.
const double kIntegratePeriod = 10.0;
const double kInitialStepSize = 1e-5;
const double kCutoff = 40.0;

// Datapoint index constants.
constexpr uint32_t kIndexX = 0; // Base coordinate.
constexpr uint32_t kIndexY = 1; // Fiber coordinate.

constexpr uint32_t kIndexY1 = 1; // Fiber coordinate 1 for ode state.
constexpr uint32_t kIndexY2 = 2; // Fiber coordinate 2 for ode state.



#endif  // MAGIC_NUMBERS_H_
