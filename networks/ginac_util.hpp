#ifndef GINAC_UTIL_HPP_
#define GINAC_UTIL_HPP_

#include <ginac/ginac.h>
#include <spdlog/spdlog.h>

#include <array>
#include <cmath>
#include <complex>
#include <string>

#include "magic_numbers.h"
#include "type_util.hpp"

auto numeric_to_complex(const GiNaC::numeric &z) -> cplx;
auto complex_to_ex(const cplx &z) -> GiNaC::ex;
auto eval_ex_to_complex(const GiNaC::ex &f, const GiNaC::symbol &x,
                        const cplx &z) -> cplx;

auto sylvester_matrix(const GiNaC::ex &f, const GiNaC::ex &g,
                      const GiNaC::symbol &x) -> GiNaC::matrix;
auto discriminant(const GiNaC::ex &f, const GiNaC::symbol &x) -> GiNaC::ex;

auto roots(const GiNaC::ex &f, const GiNaC::symbol &x, bool initial_schur) -> std::vector<cplx>;
auto schur_roots(const GiNaC::ex &f, const GiNaC::symbol &x) -> std::vector<cplx>;
auto refine_root_newton(const GiNaC::ex& f, const GiNaC::symbol& x, cplx z0) -> cplx;

#endif  // GINAC_UTIL_HPP_
