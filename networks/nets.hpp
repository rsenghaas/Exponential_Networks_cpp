#ifndef NETS_HPP_
#define NETS_HPP_

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
    start_paths();
  }

  // IO methods.
  auto print_ramification_points() -> void;
  auto evolution_step() -> void;
  auto draw_map() -> void;
  auto save_intersections() -> void;

 protected:
  // Network parameter.
  double theta_;
  uint32_t next_id_{0};

  // SW_curve data.
  std::shared_ptr<SW_curve> curve_;
  std::vector<std::array<cplx, 2>> ramification_points_;

  // Path stuff.
  auto start_paths() -> void;
  auto evolve_path(std::vector<Path>::iterator path_it, double cutoff) -> void;
  auto add_new_path(state_type start_point) -> void;
  auto save_data(uint32_t id) -> void;
  auto get_iterator_by_id(std::vector<Path>& path_vec, uint32_t id)
    -> std::vector<Path>::iterator;
  // auto determine_sign(const state_type &r, state_type &v) -> void;

  auto move_to_evolved(std::vector<Path>::iterator) -> void;
  auto add_new_paths() -> void;
  std::vector<Path> new_paths_;
  std::vector<Path> evolved_paths_;

  // Map stuff.
  auto draw_map(std::vector<Path>::iterator path_it) -> void;
  auto add_to_map(std::vector<Path>::iterator path_it) -> void;

  // Intersection stuff.
  auto self_intersection_handler(uint32_t id, bool truncate, int32_t n,
                                 uint32_t intersection_number, bool shift,
                                 bool swap) -> void;

  auto self_intersections(std::vector<Path>::iterator path_it)
      -> std::vector<intersection>;

  auto two_path_intersection_handler(uint32_t id_A, uint32_t id_B,
                                     bool truncate_A, bool truncate_B,
                                     int32_t n, uint32_t intersection_number,
                                     bool shift, bool swap) -> void;
  auto two_path_intersections(std::vector<Path>::iterator path_A_it,
                              std::vector<Path>::iterator path_B_it)
      -> std::vector<intersection>;

  auto compute_intersection_points(intersection &inter,
                                   std::vector<Path>::iterator path_A_it,
                                   std::vector<Path>::iterator path_B_it,
                                   int32_t n, state_type &new_state) -> bool;
  
  
  // *********Old Stuff*********** //
  Map map_;
  std::vector<intersection> new_intersections_;
  std::vector<computed_intersection> computed_intersections_;
  auto handle_new_intersections(
      std::vector<std::vector<path_point>> intersection_candidates,
      std::vector<path_point>::iterator current_pp_it,
      std::vector<Path>::iterator current_path_it) -> bool;
  auto compute_intersection_points_old() -> void;
  
  // Path cutoff stuff.
  uint32_t path_2_endtime_;
  uint32_t path_2_end_partner_;
};

#endif  // NETS_HPP_
