#ifndef NETS_HPP_
#define NETS_HPP_

#include <fmt/core.h>
#include <ginac/ginac.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include "magic_numbers.h"
#include "maps.hpp"
#include "ode_integrator.hpp"
#include "path.hpp"
#include "spdlog/spdlog.h"
#include "sw_curve.hpp"
#include "type_util.hpp"

struct computed_intersection {
  std::array<state_type, 2> states;
  std::array<uint32_t, 2> ids;
  std::array<uint32_t, 2> times;
};

class Network {
 public:
  Network(GiNaC::ex (*func)(const GiNaC::symbol &, const GiNaC::symbol &),
          double theta)
      : theta_(theta) {
    curve_ = std::make_shared<SW_curve>(func, "");
    ramification_points_ = curve_->get_ramification_points();
    start_path();
  }

  // IO methods.
  auto print_ramification_points() -> void;
  auto evolution_step() -> void;
  auto draw_map() -> void;
  auto save_intersections() -> void;

  auto engineer_bound_states() -> void;

 private:
  // Network parameter.
  double theta_;
  uint32_t next_id_{0};

  // SW_curve data.
  std::shared_ptr<SW_curve> curve_;
  std::vector<std::array<cplx, 2>> ramification_points_;

  // Path stuff.
  auto start_path() -> void;
  auto determine_sign(const state_type &r, state_type &v) -> void;

  auto move_to_evolved(std::vector<Path>::iterator) -> void;
  auto add_new_paths() -> void;
  std::vector<Path> new_paths_;
  std::vector<Path> evolved_paths_;

  // Map stuff.
  Map map_;

  // Intersection stuff.
  std::vector<intersection> new_intersections_;
  std::vector<computed_intersection> computed_intersections_;
  auto handle_new_intersections(
      std::vector<std::vector<path_point>> intersection_candidates,
      std::vector<path_point>::iterator current_pp_it,
      std::vector<Path>::iterator current_path_it) -> bool;
  auto compute_intersection_points() -> void;
};

#endif  // NETS_HPP_
