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
  explicit ADHM(double theta) : theta_(theta), path_2_endtime_(0), path_2_end_partner_(2) {
    curve_ = std::make_shared<SW_curve>(H_c3, "");
    start_paths();
  }

  auto BPS_state(std::vector<uint32_t> pattern_vec) -> void;
  auto backwards(std::vector<uint32_t> pattern_vec) -> void;
  auto custom_BPS() -> void;

 private:
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

  std::vector<Path> new_paths_;

  // Map.
  Map map_;
  auto draw_map(std::vector<Path>::iterator path_it) -> void;
  auto add_to_map(std::vector<Path>::iterator path_it) -> void;

  // Intersection stuff.
  auto handle_self_intersections(
      std::vector<std::vector<path_point>> intersection_candidates,
      std::vector<intersection> &new_intersections,
      std::vector<path_point>::iterator current_pp_it,
      std::vector<Path>::iterator current_path_it) -> bool;
  auto handle_new_intersections() -> void;
  // WARN: This probably should go somewhere else, since this is generally
  // associated to path and doesn't use any member variables.
  //
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

  auto shift_loop(uint32_t path_index, int32_t n, bool swap) -> void;
  auto trivial_loop(uint32_t path_index) -> void;

  auto compute_intersection_points(intersection &inter,
                                   std::vector<Path>::iterator path_A_it,
                                   std::vector<Path>::iterator path_B_it,
                                   int32_t n, state_type &new_state) -> bool;

  // Path cutoff stuff.
  uint32_t path_2_endtime_;
  uint32_t path_2_end_partner_;
};

#endif  // ADHM_BPS_HPP_
