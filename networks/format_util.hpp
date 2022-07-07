#ifndef FORMAT_UTIL_HPP_
#define FORMAT_UTIL_HPP_

#include <fmt/core.h>
#include <complex>
#include <cmath>
#include <string>
#include <array>

#include "magic_numbers.h"
#include "type_util.hpp"
auto complex_to_string(const cplx &z) -> std::string;

#endif  // FORMAT_UTIL_HPP_