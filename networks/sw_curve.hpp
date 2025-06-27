#ifndef SW_CURVE_HPP_
#define SW_CURVE_HPP_

#include <ginac/ginac.h>

#include <array>
#include <boost/range/combine.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <mutex>
#include <numbers>
#include <sstream>
#include <vector>

#include "arb_util.hpp"
#include "ginac_util.hpp"
#include "magic_numbers.h"
#include "type_util.hpp"

class SW_curve {
 public:
  SW_curve(GiNaC::ex (*func)(const GiNaC::symbol &, const GiNaC::symbol &),
           std::string diff_mode)
      : x_("x"),
        y_("y"),
        H_(func(x_, y_)),
        mode(std::move(diff_mode)),
        dH_dx_(),
        dH_dy_(),
        d2H_dy2_() {
    compute_derivatives();
  }

  auto eval_H(const cplx &x, const cplx &y) -> cplx;
  auto eval_dH_dx(const cplx &x, const cplx &y) -> cplx;
  auto eval_dH_dy(const cplx &x, const cplx &y) -> cplx;
  auto eval_d2H_dy2(const cplx &x, const cplx &y) -> cplx;

  auto sw_differential(const state_type &v, state_type &dv) -> void;
  auto elliptic_differential(const state_type &v, state_type &dv) -> void;

  auto get_branch_points() -> std::vector<cplx>;
  auto get_ramification_points() -> std::vector<std::array<cplx, 2>>;
  auto get_fiber(const cplx &x) -> std::vector<cplx>;
  auto get_branched_sheet(const cplx &x) -> cplx;
  auto match_fiber(state_type &v) -> void;
  auto save_branch_points(std::vector<cplx> branch_points) -> void;
  auto refine_root(const cplx& y, const cplx &x) -> cplx;
  auto get_kappa(const cplx &x, const cplx &y, uint32_t n) -> cplx;

  std::mutex sw_mutex;
  std::string mode;

 private:
  GiNaC::symbol x_, y_;
  GiNaC::ex H_;
  GiNaC::ex dH_dx_;
  GiNaC::ex dH_dy_;
  GiNaC::ex d2H_dy2_;

  std::vector<GiNaC::ex> branch_pts_;
  std::vector<GiNaC::ex> sing_pts;

  auto compute_derivatives() -> void;
};

#endif  // SW_CURVE_HPP_
