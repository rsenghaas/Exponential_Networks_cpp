#ifndef ARB_UTIL_HPP_
#define ARB_UTIL_HPP_

#include <flint/acb.h>
#include <flint/arb.h>

#include <string>

#include "type_util.hpp"

auto arb_to_double(const arb_t x) -> double;
auto acb_to_cplx(const acb_t z) -> cplx;
auto acb_abs_diff_small(const acb_t a,
                        const acb_t b,
                        double eps) -> bool;
bool acb_abs_too_small_or_large(const acb_t a,
                                double bound);


#endif  // ARB_UTIL_HPP_
