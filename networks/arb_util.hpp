#ifndef ARB_UTIL_HPP_
#define ARB_UTIL_HPP_

#include <acb.h>
#include <arb.h>

#include <string>

#include "type_util.hpp"

auto arb_to_double(arb_t x) -> double;
auto acb_to_cplx(acb_t z) -> cplx;

#endif  // ARB_UTIL_HPP_
