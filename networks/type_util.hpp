#ifndef TYPE_UTIL_HPP_
#define TYPE_UTIL_HPP_

#include <array>
#include <complex>
#include <future>
#include <tuple>
#include <vector>

#include "magic_numbers.h"
#include "spdlog/spdlog.h"

// ! Type definition should be in an extra header.
// Type definitions.
using cplx = std::complex<double>;

auto complex_to_string(const cplx &z) -> std::string;

constexpr uint32_t kIndexX = 0;  // Base coordinate.
constexpr uint32_t kIndexY = 1;  // Fiber coordinate.

constexpr uint32_t kIndexY1 = 1;  // Fiber coordinate 1 for ode state.
constexpr uint32_t kIndexY2 = 2;  // Fiber coordinate 2 for ode state.

using state_type = std::array<cplx, 3>;
using future_type = std::tuple<std::future<std::vector<state_type>>,
                               std::future<std::vector<double>>>;

// Time indices
constexpr uint32_t kIndexStartTime = 0;
constexpr uint32_t kIndexEndTime = 1;

struct path_point {
  uint32_t id;
  uint32_t pp_vec_index;
  int32_t coordinate_real;
  int32_t coordinate_imag;
  std::array<uint32_t, 2> t;
};

// Path indices
constexpr uint32_t kIndexFirstPath = 0;
constexpr uint32_t kIndexSecondPath = 1;

constexpr uint32_t kIndexCoordReal = 0;
constexpr uint32_t kIndexCoordImag = 1;

struct intersection {
  std::array<uint32_t, 2> ids;
  std::vector<std::array<int32_t, 2>> coordinates;
  std::array<std::array<uint32_t, 2>, 2> times;
};

auto print_intersection(const intersection &intersect) -> void;
auto print_path_point(const path_point &pp) -> void;
auto print_state_type(const state_type &v) -> void;

inline auto operator+(const state_type &v1, const state_type &v2)
    -> state_type {
  return state_type{v1.at(kIndexX) + v2.at(kIndexX),
                    v1.at(kIndexY1) + v2.at(kIndexY1),
                    v1.at(kIndexY2) + v2.at(kIndexY2)};
}
inline auto operator*(const double &lambda, const state_type &v) -> state_type {
  return state_type{lambda * v.at(kIndexX), lambda * v.at(kIndexY1),
                    lambda * v.at(kIndexY2)};
}

#endif  // TYPE_UTIL_HPP_
