#ifndef TYPE_UTIL_HPP_
#define TYPE_UTIL_HPP_

#include <complex>
#include <vector>
#include <array>
#include <future>
#include <tuple>

#include "magic_numbers.h"

// ! Type definition should be in an extra header.
// Type definitions.
using cplx = std::complex<double>;
typedef std::array<cplx, 3> state_type;
typedef std::tuple<std::future<std::vector<state_type>>, std::future<std::vector<double>>> future_type; 

inline state_type operator+(const state_type &v1, const state_type &v2){return state_type{v1.at(kIndexX) + v2.at(kIndexX), v1.at(kIndexY1) + v2.at(kIndexY1), v1.at(kIndexY2) + v2.at(kIndexY2)};}
inline state_type operator*(const double& lambda, const state_type& v){return state_type{lambda * v.at(kIndexX), lambda * v.at(kIndexY1), lambda * v.at(kIndexY2)};}

#endif  // TYPE_UTIL_HPP_
