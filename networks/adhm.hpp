#ifndef ADHM_BPS_HPP_
#define ADHM_BPS_HPP_

#include <ginac/ginac.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

#include "maps.hpp"
#include "path.hpp"
#include "sw_curve.hpp"
#include "type_util.hpp"

auto H_c3(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;

class ADHM {
 public:
  explicit ADHM(double theta) : theta_(theta) {
    curve_ = std::make_shared<SW_curve>(H_c3);
    start_paths();
  }

  auto BPS_state() -> void;

 private:
  double theta_;

  // SW_curve data.
  std::shared_ptr<SW_curve> curve_;
  std::vector<std::array<cplx, 2>> ramification_points_;

  // Path stuff.
  auto start_paths() -> void;
  auto evolve_path(std::vector<Path>::iterator path_it, double cutoff) -> void;
  auto add_new_path(state_type start_point, uint32_t path_id) -> void;

  std::vector<Path> new_paths_;

  // Map.
  Map map_;
  auto draw_map(std::vector<Path>::iterator path_it) -> void;

  // Intersection stuff.
  auto handle_self_intersections(
      std::vector<std::vector<path_point>> intersection_candidates,
      std::vector<intersection> &new_intersections,
      std::vector<path_point>::iterator current_pp_it,
      std::vector<Path>::iterator current_path_it) -> bool;
  auto handle_new_intersections() -> void;
};

#endif  // ADHM_BPS_HPP_
