#include "sw_curve.hpp"

#include <acb.h>
#include <acb_elliptic.h>
#include <arb.h>

#include <fstream>

auto SW_curve::compute_derivatives() -> void {
  dH_dx_ = H_.diff(x_);
  dH_dy_ = H_.diff(y_);
  d2H_dy2_ = H_.diff(y_, 2);
  std::stringstream H_ss, dH_dx_ss, dH_dy_ss, d2H_dy2_ss;
  H_ss << H_;
  dH_dx_ss << dH_dx_;
  dH_dy_ss << dH_dy_;
  d2H_dy2_ss << d2H_dy2_;
  spdlog::info("The Network is initialized with H = {}.", H_ss.str());
  spdlog::info(
      "The derivatives are calculated as:\n * dH/dx = {}\n * dH/dy = {}\n * "
      "d²H/dy² = {}.",
      dH_dx_ss.str(), dH_dy_ss.str(), d2H_dy2_ss.str());
}

auto SW_curve::eval_H(const cplx &x, const cplx &y) -> cplx {
  GiNaC::ex x_ex = complex_to_ex(x);
  GiNaC::ex y_ex = complex_to_ex(y);
  return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(
      H_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::eval_dH_dx(const cplx &x, const cplx &y) -> cplx {
  GiNaC::ex x_ex = complex_to_ex(x);
  GiNaC::ex y_ex = complex_to_ex(y);
  return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(
      dH_dx_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::eval_dH_dy(const cplx &x, const cplx &y) -> cplx {
  GiNaC::ex x_ex = complex_to_ex(x);
  GiNaC::ex y_ex = complex_to_ex(y);
  return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(
      dH_dy_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::eval_d2H_dy2(const cplx &x, const cplx &y) -> cplx {
  GiNaC::ex x_ex = complex_to_ex(x);
  GiNaC::ex y_ex = complex_to_ex(y);
  return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(
      d2H_dy2_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::get_branch_points() -> std::vector<cplx> {
  GiNaC::ex disc = discriminant(H_, y_).expand();
  std::cout << disc.is_polynomial(x_) << std::endl;
  GiNaC::ex red = disc;
  GiNaC::ex x_zero = complex_to_ex(0);
  while (GiNaC::ex_to<GiNaC::numeric>(
             red.subs(GiNaC::lst{x_ == x_zero}).evalf()) == 0) {
    red = (red * GiNaC::pow(x_, -1)).expand();
  }
  std::stringstream disc_ss;
  disc_ss << red;
  spdlog::info("The reduced discriminant is {}.", disc_ss.str());
  std::vector<cplx> branch_points = roots(red, x_);
  save_branch_points(branch_points);
  return branch_points;
}

auto SW_curve::get_ramification_points() -> std::vector<std::array<cplx, 2>> {
  std::vector<cplx> branch_points = get_branch_points();
  std::vector<std::array<cplx, 2>> ramification_points;
  for (auto &b : branch_points) {
    std::array<cplx, 2> r;
    r.at(kIndexX) = b;
    r.at(kIndexY) = get_branched_sheet(b);
    ramification_points.push_back(r);
  }
  return ramification_points;
}

// At the moment this only does exponential networks
// It might to return to expressions and build the differential outside the
// class.
auto SW_curve::sw_differential(const state_type &v, state_type &dv) -> void {
  dv.at(kIndexX) = v.at(kIndexX) / (v.at(kIndexY2) - v.at(kIndexY1));
  dv.at(kIndexY1) = -eval_dH_dx(v.at(kIndexX), std::exp(v.at(kIndexY1))) /
                    eval_dH_dy(v.at(kIndexX), std::exp(v.at(kIndexY1))) *
                    dv.at(kIndexX) / std::exp(v.at(kIndexY1));
  dv.at(kIndexY2) = -eval_dH_dx(v.at(kIndexX), std::exp(v.at(kIndexY2))) /
                    eval_dH_dy(v.at(kIndexX), std::exp(v.at(kIndexY2))) *
                    dv.at(kIndexX) / std::exp(v.at(kIndexY2));
}

auto SW_curve::elliptic_differential(const state_type &v, state_type &dv)
    -> void {
  acb_t p_y1, p_y2, tau, y1, y2;
  acb_init(y1);
  acb_init(y2);
  acb_init(tau);
  acb_init(p_y1);
  acb_init(p_y2);
  acb_set_d_d(y1, v.at(kIndexY1).real(), v.at(kIndexY1).imag());
  acb_set_d_d(y2, v.at(kIndexY2).real(), v.at(kIndexY2).imag());
  acb_onei(tau);
  acb_elliptic_p(p_y1, y1, tau, 50);
  acb_elliptic_p(p_y2, y2, tau, 50);

  cplx p_y1_d = acb_to_cplx(p_y1);
  cplx p_y2_d = acb_to_cplx(p_y2);

  dv.at(kIndexX) = v.at(kIndexX) / (v.at(kIndexY2) - v.at(kIndexY1));
  dv.at(kIndexY1) = -eval_dH_dx(v.at(kIndexX), p_y1_d) /
                    eval_dH_dy(v.at(kIndexX), p_y1_d) * dv.at(kIndexX) / p_y1_d;
  dv.at(kIndexY2) = -eval_dH_dx(v.at(kIndexX), p_y2_d) /
                    eval_dH_dy(v.at(kIndexX), p_y2_d) * dv.at(kIndexX) / p_y2_d;
}

auto SW_curve::get_fiber(const cplx &x) -> std::vector<cplx> {
  GiNaC::ex fiber_poly = H_.subs({x_ == complex_to_ex(x)});
  std::vector<cplx> fiber = roots(fiber_poly, y_);
  return fiber;
}

auto SW_curve::match_fiber(state_type &v) -> void {
  cplx x = v.at(kIndexX);
  std::vector<cplx> fiber = get_fiber(x);
  cplx nearest_fiber_y1 = fiber.at(0);
  cplx nearest_fiber_y2 = fiber.at(0);
  for (auto &f : fiber) {
    if (std::abs(std::exp(v.at(kIndexY1)) - f) <
        std::abs(std::exp(v.at(kIndexY1)) - std::exp(nearest_fiber_y1))) {
      nearest_fiber_y1 = std::log(f);
    }
    if (std::abs(std::exp(v.at(kIndexY2)) - f) <
        std::abs(std::exp(v.at(kIndexY2)) - std::exp(nearest_fiber_y2))) {
      nearest_fiber_y2 = std::log(f);
    }
  }
  cplx dv_y1 = v.at(kIndexY1) - nearest_fiber_y1;
  cplx dv_y2 = v.at(kIndexY2) - nearest_fiber_y2;
  auto k1 = static_cast<int32_t>(std::round(((dv_y1) / (2 * pi * J)).real()));
  auto k2 = static_cast<int32_t>(std::round(((dv_y2) / (2 * pi * J)).real()));
  v.at(kIndexY1) = nearest_fiber_y1 + 2 * pi * J * static_cast<double>(k1);
  v.at(kIndexY2) = nearest_fiber_y2 + 2 * pi * J * static_cast<double>(k2);
}

auto SW_curve::get_branched_sheet(const cplx &x) -> cplx {
  std::vector<cplx> fiber = get_fiber(x);
  double diff = std::numeric_limits<double>::max();
  cplx sheet;
  for (auto it1 = std::rbegin(fiber); it1 != std::rend(fiber); it1++) {
    for (auto it2 = it1 + 1; it2 != std::rend(fiber); it2++) {
      double new_diff = std::abs(*it1 - *it2);
      if (new_diff < diff) {
        diff = new_diff;
        sheet = (*it1 + *it2) / 2.0;
      }
    }
  }
  return sheet;
}

auto SW_curve::save_branch_points(std::vector<cplx> branch_points) -> void {
  std::fstream data_file;
  std::string sep = ",";
  std::string filename = "data/map_data/branch_points.csv";
  data_file.open(filename, std::ios::out);
  if (!data_file) {
    spdlog::debug("{} could not be created.", filename);
  } else {
    for (auto &pt : branch_points) {
      std::string output_line = complex_to_string(pt) + "\n";
      data_file << output_line;
    }
    spdlog::debug("Branch points saved to {}.", filename);
    data_file.close();
  }
}