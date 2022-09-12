#include "adhm.hpp"

#include <fmt/core.h>
#include <ginac/ginac.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include "Eigen/Dense"
#include "maps.hpp"

auto determine_sign(SW_curve& curve, const state_type& r, state_type& v)
    -> void {
  state_type dv;
  curve.sw_differential(v, dv);
  cplx dx = v.at(kIndexX) - r.at(kIndexX);
  cplx next_dx = dv.at(kIndexX);
  if (next_dx.real() * dx.real() + next_dx.imag() * dx.imag() < 0) {
    cplx temp = v.at(kIndexY1);
    v.at(kIndexY1) = v.at(kIndexY2);
    v.at(kIndexY2) = temp;
  }
}

auto H_c3(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return -x + GiNaC::pow(y, 2) + y;
}

auto ADHM::start_paths() -> void {
  ramification_points_ = curve_->get_ramification_points();
  for (auto& r : ramification_points_) {
    cplx b = r.at(kIndexX);
    cplx y = r.at(kIndexY);

    state_type start_state;
    start_state.at(kIndexX) = b;
    start_state.at(kIndexY1) = std::log(y);
    start_state.at(kIndexY2) = std::log(y);

    // Compute the expansion around a branch point.
    // rsenghaas TODO: I should move this to the SW_curve class and compute it
    // generally (kappa = 1 / n! d^n H / dy^n * (1 / d H / dx))
    cplx kappa =
        -1.0 / 2 * curve_->eval_d2H_dy2(b, y) / curve_->eval_dH_dx(b, y);

    cplx dx;
    cplx dy;
    for (uint32_t k = 0; k < 3; k++) {
      std::vector<state_type> v;
      std::vector<double> masses;
      state_type next_state;

      v.push_back(start_state);
      masses.push_back(0);

      dx = std::pow(std::pow(3.0 / 4 * std::pow(kappa, 1.0 / 2) * y * b *
                                 std::exp(J * theta_) * kInitialStepSize,
                             2.0),
                    1.0 / 3) *
           std::pow(kZeta3,
                    k);  // ! This is only valid for exponential networks.
      dy = 3.0 / 4.0 * y * b * std::exp(J * theta_) * kInitialStepSize / dx;

      next_state.at(kIndexX) = b + dx;
      next_state.at(kIndexY1) = start_state.at(kIndexY1) + dy / y;
      next_state.at(kIndexY2) = start_state.at(kIndexY2) - dy / y;
      curve_->match_fiber(v.back());

      determine_sign(*curve_, start_state, next_state);
      v.push_back(next_state);

      spdlog::debug(
          "Numerical Check: \n |H(x, y_1)| = {}\n |H(x, y_2)| = {}.",
          complex_to_string(std::abs(curve_->eval_H(
              v.back().at(kIndexX), std::exp(v.back().at(kIndexY1))))),
          complex_to_string(std::abs(curve_->eval_H(
              v.back().at(kIndexX), std::exp(v.back().at(kIndexY2))))));

      cplx dlog_y1 = next_state.at(kIndexY1) - start_state.at(kIndexY1);
      cplx dlog_y2 = next_state.at(kIndexY2) - start_state.at(kIndexY2);
      state_type dv{dx, dlog_y1, dlog_y2};
      masses.push_back(compute_dm(next_state, dv));
      // First step is done.

      // Doing some runge kutta steps before going into the boost::odeint
      // integration.
      for (uint32_t i = 0; i < kInitialSteps; i++) {
        ODE_runge_kutta_step(curve_, v, masses, kInitialStepSize, theta_);
      }

      Path path(v, masses, new_paths_.size());
      new_paths_.push_back(std::move(path));
      next_id_++;

      spdlog::debug("Path {} appended.", new_paths_.back().path_id_);
      spdlog::debug(
          "Numerical Check: \n |H(x, y_1)| = {}\n |H(x, y_2)| = {}.",
          complex_to_string(std::abs(curve_->eval_H(
              v.back().at(kIndexX), std::exp(v.back().at(kIndexY1))))),
          complex_to_string(std::abs(curve_->eval_H(
              v.back().at(kIndexX), std::exp(v.back().at(kIndexY2))))));
    }
    spdlog::debug("Paths are started.");
  }
}

auto get_iterator_by_id(std::vector<Path>& path_vec, uint32_t id)
    -> std::vector<Path>::iterator {
  auto ret_it = path_vec.begin();
  while (ret_it != path_vec.end()) {
    if (ret_it->path_id_ == id) {
      break;
    }
    ret_it++;
  }
  return ret_it;
}

const state_type cutoffPoint = {
    -1.8897512312374315e-09 - 1.8882623082801205e-10 * J,
    -20.08185327082071 + 3.2411832963911094 * J,
    -1.8897512104510875e-09 + 3.1415926534009646 * J};

auto ADHM::backwards(std::vector<uint32_t> pattern_vec) -> void {
  auto path_it = get_iterator_by_id(new_paths_, 0);
  std::vector<state_type> v_cutoff = {cutoffPoint};
  std::vector<double> masses_cutoff = {0};
  Path cutoff_path(v_cutoff, masses_cutoff, 0);
  *path_it = cutoff_path;
  evolve_path(path_it, 4 * kD4Cutoff);
  state_type end_pt = path_it->get_endpoint();
  spdlog::debug("Cutoff point");
  // curve_->match_fiber(end_pt);
  print_state_type(end_pt);

  // TODO: This should eventually be an input parameter.
  if (pattern_vec.front() != 1) {
    spdlog::error("Illegal pattern");
    return;
  }

  for (auto vec_it = pattern_vec.begin(); vec_it + 1 != pattern_vec.end();
       vec_it++) {
    spdlog::debug("{}", *vec_it);
    if (*vec_it != *(vec_it + 1) && *vec_it + 1 != *(vec_it + 1) &&
        *vec_it - 1 != *(vec_it + 1)) {
      spdlog::error("Illegal pattern");
      return;
    }
  }

  // std::reverse(pattern_vec.begin(), pattern_vec.end());

  spdlog::debug("{}", pattern_vec.front());

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kD4Cutoff);

  if (pattern_vec.size() == 1) {
    save_data(0);
    save_data(1);
    return;
  }
  std::reverse(pattern_vec.begin(), pattern_vec.end());

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, kCutoff);

  spdlog::debug("New paths size {}", new_paths_.size());

  if (pattern_vec.at(0) == pattern_vec.at(1)) {
    self_intersection_handler(1, true, -1, 0, false, false);

    path_it = get_iterator_by_id(new_paths_, 1);
    save_data(1);

    path_it = get_iterator_by_id(new_paths_, 3);
    evolve_path(path_it, kCutoff);

    two_path_intersection_handler(2, 3, false, true, 0, 0, false, false);
  } else {
    spdlog::debug("Shift loop");
    int32_t n = pattern_vec.at(1);
    bool swap = true;
    if (pattern_vec.at(0) < pattern_vec.at(1)) {
      n = pattern_vec.at(0);
      swap = false;
    }

    self_intersection_handler(1, true, -n, 0, true, swap);

    path_it = get_iterator_by_id(new_paths_, 1);
    save_data(1);

    path_it = get_iterator_by_id(new_paths_, 2);
    evolve_path(path_it, kCutoff);

    path_it = get_iterator_by_id(new_paths_, 3);
    evolve_path(path_it, kCutoff);

    two_path_intersection_handler(2, 3, false, true, 0, 0, true, false);
  }
  save_data(2);
  save_data(3);

  uint32_t current_index = 4;
  for (auto it = pattern_vec.begin() + 1; it != pattern_vec.end() - 1; it++) {
    spdlog::debug("Drawing loops.");
    if (*it == *(it + 1)) {
      trivial_loop(current_index);

    } else if (*it > *(it + 1)) {
      shift_loop(current_index, -*(it + 1), true);
    } else {
      shift_loop(current_index, -*(it), false);
    }
    current_index += 2;
    spdlog::debug("Finished loop.");
  }

  double M_total;
  for (auto& pth : new_paths_) {
    M_total += pth.get_endmass();
  }
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  save_data(current_index);

  two_path_intersection_handler(0, 2, true, false, 1, 0, false, false);
  save_data(0);
  current_index = new_paths_.size() - 1;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, current_index, false, true, 0, 0, false,
                                false);
  save_data(current_index);
  current_index++;

  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  save_data(current_index);
  current_index++;

  for (auto it = pattern_vec.begin() + 1; it != pattern_vec.end() - 1; it++) {
    if (*it == *(it + 1)) {
      auto path_it = get_iterator_by_id(new_paths_, current_index);
      evolve_path(path_it, kCutoff);

      two_path_intersection_handler(1, current_index, false, true, 1, 0, false,
                                    false);
      save_data(current_index);

      path_it = get_iterator_by_id(new_paths_, current_index);
      evolve_path(path_it, kCutoff);

      two_path_intersection_handler(2, current_index, false, true, 0, 0, false,
                                    false);
      save_data(current_index);
    } else if (*it > *(it + 1)) {
      shift_loop(current_index, -*(it + 1), true);
    } else {
      shift_loop(current_index, -*(it), false);
    }
    current_index += 2;
  }
  path_it = get_iterator_by_id(new_paths_, 2);
  path_it->truncate(0, path_2_endtime_ + 1);
  save_data(2);
}

auto ADHM::custom_BPS() -> void {
  auto path_it = get_iterator_by_id(new_paths_, 0);
  std::vector<state_type> v_cutoff = {cutoffPoint};
  std::vector<double> masses_cutoff = {0};
  Path cutoff_path(v_cutoff, masses_cutoff, 0);
  *path_it = cutoff_path;
  evolve_path(path_it, 4 * kD4Cutoff);
  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);
  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(0, 2, true, false, 3, 0, true, true);
  next_id_--;
  new_paths_.pop_back();
  save_data(0);

  self_intersection_handler(1, false, -1, 0, false, false);
  self_intersection_handler(1, true, -1, 0, true, true);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, kCutoff);

  path_it = get_iterator_by_id(new_paths_, 4);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(2, 3, true, true, 0, 0, true, false);
  two_path_intersection_handler(2, 4, false, true, 0, 0, true, false);

  save_data(2);
  save_data(3);

  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 5, false, true, -1, 0, true, false);
  save_data(5);

  path_it = get_iterator_by_id(new_paths_, 7);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(4, 7, true, true, -1, 0, false, false);
  save_data(4);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 8, false, true, 0, 0, true, false);
  save_data(8);

  path_it = get_iterator_by_id(new_paths_, 9);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 9, false, true, -2, 0, true, true);
  save_data(9);

  path_it = get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 10, false, true, 0, 0, true, false);
  save_data(10);

  path_it = get_iterator_by_id(new_paths_, 11);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 11, false, true, -1, 0, true, true);
  save_data(11);

  path_it = get_iterator_by_id(new_paths_, 12);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 12, false, true, 0, 0, true, false);
  save_data(12);

  /* self_intersection_handler(1, true, -1, 0, false, false);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, kCutoff);

  path_it = get_iterator_by_id(new_paths_, 4);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(2,3,false,true,0,0,true,false);
  two_path_intersection_handler(2,4,true,true,0,0,true,false);

  save_data(2);
  save_data(3);
  save_data(4);

  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,6,false,true,-1,0,false,false);
  save_data(6);

  path_it = get_iterator_by_id(new_paths_, 7);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2,7,false,true,0,0,true,false);
  save_data(7);

  path_it =get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,8,false,true,-1,0,true,true);
  save_data(8);

  path_it =get_iterator_by_id(new_paths_, 9);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2,9,false,true,-1,0,true,false);
  save_data(9); */

  /*path_it =get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,10,false,true,-1,0,true,false);
  save_data(10);*/
}

auto ADHM::BPS_state(std::vector<uint32_t> pattern_vec) -> void {
  auto path_it = get_iterator_by_id(new_paths_, 0);
  std::vector<state_type> v_cutoff = {cutoffPoint};
  std::vector<double> masses_cutoff = {0};
  Path cutoff_path(v_cutoff, masses_cutoff, 0);
  *path_it = cutoff_path;
  evolve_path(path_it, 4 * kD4Cutoff);
  state_type end_pt = path_it->get_endpoint();
  spdlog::debug("Cutoff point");
  // curve_->match_fiber(end_pt);
  print_state_type(end_pt);

  // TODO: This should eventually be an input parameter.
  if (pattern_vec.front() != 1) {
    spdlog::error("Illegal pattern");
    return;
  }

  for (auto vec_it = pattern_vec.begin(); vec_it + 1 != pattern_vec.end();
       vec_it++) {
    spdlog::debug("{}", *vec_it);
    if (*vec_it != *(vec_it + 1) && *vec_it + 1 != *(vec_it + 1) &&
        *vec_it - 1 != *(vec_it + 1)) {
      spdlog::error("Illegal pattern");
      return;
    }
  }

  std::reverse(pattern_vec.begin(), pattern_vec.end());

  spdlog::debug("{}", pattern_vec.front());

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kD4Cutoff);

  save_data(1);
  if (pattern_vec.size() == 1) {
    save_data(1);
    return;
  }

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(0, 2, true, false, 0, 0, false, false);
  next_id_--;
  save_data(0);
  new_paths_.pop_back();
  spdlog::debug("New paths size {}", new_paths_.size());

  /* self_intersection_handler(1, true, -1, 0, true, true);
  path_it = get_iterator_by_id(new_paths_, 1);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, kD4Cutoff);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 4*kD4Cutoff);

  two_path_intersection_handler(2, 3, true, true, 0, 0, true, false);
  save_data(2);
  save_data(3); */

  if (pattern_vec.at(0) == pattern_vec.at(1)) {
    self_intersection_handler(1, true, -1, 0, false, false);

    path_it = get_iterator_by_id(new_paths_, 1);
    save_data(1);

    path_it = get_iterator_by_id(new_paths_, 3);
    evolve_path(path_it, kCutoff);

    two_path_intersection_handler(2, 3, false, true, 0, 0, false, false);
  } else {
    spdlog::debug("Shift loop");
    int32_t n = pattern_vec.at(1);
    bool swap = true;
    if (pattern_vec.at(0) < pattern_vec.at(1)) {
      n = pattern_vec.at(0);
      swap = false;
    }

    self_intersection_handler(1, true, -n, 0, true, swap);

    path_it = get_iterator_by_id(new_paths_, 1);
    save_data(1);

    path_it = get_iterator_by_id(new_paths_, 2);
    evolve_path(path_it, kCutoff);

    path_it = get_iterator_by_id(new_paths_, 3);
    evolve_path(path_it, kCutoff);

    two_path_intersection_handler(2, 3, false, true, 0, 0, true, false);
  }
  save_data(2);
  save_data(3);

  uint32_t current_index = 4;
  for (auto it = pattern_vec.begin() + 1; it != pattern_vec.end() - 1; it++) {
    spdlog::debug("Drawing loops.");
    if (*it == *(it + 1)) {
      trivial_loop(current_index);
    } else if (*it > *(it + 1)) {
      shift_loop(current_index, -*(it + 1), true);
    } else {
      shift_loop(current_index, -*(it), false);
    }
    current_index += 2;
    spdlog::debug("Finished loop.");
  }

  path_it = get_iterator_by_id(new_paths_, 2);
  path_it->truncate(0, path_2_endtime_ + 1);
  save_data(2);
  double M_total;
  for (auto& pth : new_paths_) {
    M_total += pth.get_endmass();
  }
  spdlog::debug("Mass: {}", M_total);

  // evolve_path(path_it, 4*kD4Cutoff);
  // save_data(current_index);

  /* path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, 4*kD4Cutoff);
  // save_data(5);

  // two_path_intersection_handler(2, 5, false, true, 1, 0, true);
  two_path_intersection_handler(2, 5, false, true, 0, 0, false);
  // two_path_intersection_handler(2, 3, false, false, 1);
  // two_path_intersection_handler(2, 3, true, true, 2);
  save_data(5);
  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, 4*kD4Cutoff);

  save_data(6);*/

  /* path_it = get_iterator_by_id(new_paths_, 7);
  evolve_path(path_it, 4*kD4Cutoff);
  two_path_intersection_handler(2, 7, false, true, 0, 0, false);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, 4*kD4Cutoff);
  save_data(8); */

  // two_path_intersection_handler(2, 4, true, true, 0);
  // save_data(4);

  // path_it = get_iterator_by_id(new_paths_, 5);
  // evolve_path(path_it, kD4Cutoff);
  // save_data(5);

  /* path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, 2* kD4Cutoff);
  path_it-> save_data();

  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, 3*kD4Cutoff);
  path_it-> save_data(); */

  // two_path_intersection_handler(1, 4, false, true, 0);
  // two_path_intersection_handler(1, 5, false, true, 0);
  // two_path_intersection_handler(1, 6, false, true, 0);

  // two_path_intersection_handler(2, 3, true, true, 0);
  // two_path_intersection_handler(2, 3, true, true, -1);

  // path_it = get_iterator_by_id(new_paths_, 2);
  // path_it-> save_data();

  // path_it = get_iterator_by_id(new_paths_, 3);
  // path_it-> save_data();

  // path_it = get_iterator_by_id(new_paths_, 4);
  // evolve_path(path_it, 3*kD4Cutoff);
  // two_path_intersection_handler(1, 4, false, true, -1);

  // path_it = get_iterator_by_id(new_paths_, 4);
  // path_it->save_data();

  // path_it = get_iterator_by_id(new_paths_, 5);
  // evolve_path(path_it, kD4Cutoff);
  // two_path_intersection_handler(2, 5, false, true, 0);

  // path_it = get_iterator_by_id(new_paths_, 5);
  // path_it->save_data();

  /* path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kD4Cutoff);

  path_it = get_iterator_by_id(new_paths_, 6);
  path_it->save_data(); */

  new_paths_.clear();
}

auto ADHM::shift_loop(uint32_t path_index, int32_t n, bool swap) -> void {
  spdlog::debug("Shift loop {}.", n);
  auto path_it = get_iterator_by_id(new_paths_, path_index);
  evolve_path(path_it, 4 * kCutoff);

  two_path_intersection_handler(1, path_index, false, true, n, 0, true, swap);
  save_data(path_index);

  path_index = new_paths_.size() - 1;
  path_it = get_iterator_by_id(new_paths_, path_index);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(2, path_index, false, true, 0, 0, true, false);
  save_data(path_index);
}

auto ADHM::trivial_loop(uint32_t path_index) -> void {
  spdlog::debug("Fluid loop.");
  auto path_it = get_iterator_by_id(new_paths_, path_index);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(1, path_index, false, true, -1, 0, false,
                                false);
  save_data(path_index);

  path_index = new_paths_.size() - 1;
  path_it = get_iterator_by_id(new_paths_, path_index);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(2, path_index, false, true, 0, 0, false, false);
  save_data(path_index);
}

auto ADHM::save_data(uint32_t id) -> void {
  auto path_it = get_iterator_by_id(new_paths_, id);
  path_it->save_data();
}

auto ADHM::self_intersection_handler(uint32_t id, bool truncate, int32_t n,
                                     uint32_t intersection_number, bool shift,
                                     bool swap) -> void {
  auto path_it = get_iterator_by_id(new_paths_, id);
  std::vector<intersection> intersections = self_intersections(path_it);
  spdlog::debug("{} self intersections", intersections.size());
  auto inter_it = intersections.begin();
  inter_it += intersection_number;
  state_type next_state;
  state_type pt_A =
      path_it->get_point(inter_it->times.at(0).at(kIndexStartTime));
  print_state_type(pt_A);
  int32_t state_A_k = get_log_sheet(pt_A);
  spdlog::debug("Path_A has k {}", state_A_k);
  state_type pt_B =
      path_it->get_point(inter_it->times.at(1).at(kIndexStartTime));
  print_state_type(pt_B);
  int32_t state_B_k = get_log_sheet(pt_B);
  spdlog::debug("Path_B has k {}", state_B_k);
  if (compute_intersection_points(*inter_it, path_it, path_it, n, next_state)) {
    print_state_type(next_state);
    if (truncate) {
      path_it->truncate(0,
                        inter_it->times.at(kIndexSecondPath).at(kIndexEndTime));
      // BUG: This doesn't take care of the y's, one should probably make this
      // better via match_fiber!
      path_it->add_single_point(next_state);
    }
    if (shift) {
      state_type shift_state = path_it->get_point(
          inter_it->times.at(kIndexFirstPath).at(kIndexStartTime));
      if (swap) {
        cplx temp = shift_state.at(kIndexY1);
        shift_state.at(kIndexY1) = shift_state.at(kIndexY2);
        shift_state.at(kIndexY2) = temp;
      }
      shift_state.at(kIndexX) = next_state.at(kIndexX);
      curve_->match_fiber(shift_state);
      shift_state.at(kIndexY1) = std::log(std::exp(shift_state.at(kIndexY1))) +
                                 2 * pi * J * static_cast<double>(n);
      add_new_path(shift_state);
    } else {
      add_new_path(next_state);
    }
  }
}

auto ADHM::two_path_intersection_handler(uint32_t id_A, uint32_t id_B,
                                         bool truncate_A, bool truncate_B,
                                         int32_t n,
                                         uint32_t intersection_number,
                                         bool shift, bool swap) -> void {
  auto path_A_it = get_iterator_by_id(new_paths_, id_A);
  auto path_B_it = get_iterator_by_id(new_paths_, id_B);
  if (id_A > id_B) {
    path_A_it = get_iterator_by_id(new_paths_, id_B);
    path_B_it = get_iterator_by_id(new_paths_, id_A);
  }
  std::vector<intersection> intersections =
      two_path_intersections(path_A_it, path_B_it);
  auto inter_it = intersections.begin();
  if (inter_it->times.at(kIndexSecondPath).at(kIndexStartTime) == 0) {
    inter_it++;
  }
  inter_it += intersection_number;
  print_intersection(*inter_it);
  state_type pt_A =
      path_A_it->get_point(inter_it->times.at(0).at(kIndexStartTime));
  int32_t state_A_k = get_log_sheet(pt_A);
  spdlog::debug("Path_A has k {}", state_A_k);
  state_type pt_B =
      path_B_it->get_point(inter_it->times.at(1).at(kIndexStartTime));
  int32_t state_B_k = get_log_sheet(pt_B);
  spdlog::debug("Path_B has k {}", state_B_k);

  state_type next_state;
  if (compute_intersection_points(*inter_it, path_A_it, path_B_it, n,
                                  next_state)) {
    spdlog::debug("Two path intersection.");
    if (id_A == 2 && inter_it->times.at(kIndexFirstPath).at(kIndexEndTime) >
                         path_2_endtime_) {
      path_2_endtime_ = inter_it->times.at(kIndexFirstPath).at(kIndexEndTime);
    }
    if (id_B == 2 && inter_it->times.at(kIndexSecondPath).at(kIndexEndTime) >
                         path_2_endtime_) {
      path_2_endtime_ = inter_it->times.at(kIndexSecondPath).at(kIndexEndTime);
    }

    if (truncate_A) {
      path_A_it->truncate(
          0, inter_it->times.at(kIndexFirstPath).at(kIndexEndTime));
      path_A_it->add_single_point(next_state);
    }
    if (truncate_B) {
      path_B_it->truncate(
          0, inter_it->times.at(kIndexSecondPath).at(kIndexEndTime));
      path_B_it->add_single_point(next_state);
    }
    if (shift) {
      auto shift_path_it = get_iterator_by_id(new_paths_, id_A);
      state_type shift_state = shift_path_it->get_point(
          inter_it->times.at(kIndexFirstPath).at(kIndexStartTime));
      if (swap) {
        cplx temp = shift_state.at(kIndexY1);
        shift_state.at(kIndexY1) = shift_state.at(kIndexY2);
        shift_state.at(kIndexY2) = temp;
      }
      shift_state.at(kIndexX) = next_state.at(kIndexX);
      curve_->match_fiber(shift_state);
      shift_state.at(kIndexY1) = std::log(std::exp(shift_state.at(kIndexY1))) +
                                 2 * pi * J * static_cast<double>(n);
      print_state_type(shift_state);
      add_new_path(shift_state);
    } else {
      add_new_path(next_state);
    }
  }
}

auto ADHM::add_new_path(state_type start_point) -> void {
  std::vector<state_type> v;
  v.push_back(start_point);

  std::vector<double> masses;
  masses.push_back(0);

  Path path(v, masses, next_id_);
  new_paths_.push_back(std::move(path));
  next_id_++;
}

auto ADHM::evolve_path(std::vector<Path>::iterator path_it, double cutoff)
    -> void {
  future_type future = path_it->integrate(curve_, theta_, cutoff);
  path_it->update(std::ref(future));
}

/* auto ADHM::draw_map(std::vector<Path>::iterator path_it) -> void {
  Map temp_map;
  std::vector<intersection> new_intersections;
  path_it->compute_map_points();
  for (auto pp_it = path_it->pp_vec.begin();
       pp_it != std::prev(path_it->pp_vec.end()); ++pp_it) {
    handle_self_intersections(temp_map.draw_line(*pp_it, *std::next(pp_it)),
                              new_intersections, pp_it, path_it);
  }
  for (auto& inter : new_intersections) {
    print_intersection(inter);
  }
} */

auto compare_states(state_type v1, state_type v2) -> bool {
  if ((std::abs(std::exp(v1.at(kIndexY1) - v1.at(kIndexY2)) - 1.0) <
       kFiberCompTolerance) &&
      (std::abs(std::exp(v2.at(kIndexY1) - v2.at(kIndexY2)) - 1.0) <
       kFiberCompTolerance)) {
    return true;
  }
  if (std::abs((v1.at(kIndexY1) - v1.at(kIndexY2)) -
               (v2.at(kIndexY2) - v2.at(kIndexY1))) < kFiberCompTolerance) {
    return true;
  }
  if (std::abs((v1.at(kIndexY1) - v1.at(kIndexY2)) -
               (v2.at(kIndexY1) - v2.at(kIndexY2))) < kFiberCompTolerance) {
    return true;
  }
  return false;
}

/* auto neighbour_pixel(std::array<int32_t, 2> coord_arr1,
                     std::array<int32_t, 2> coord_arr2) -> bool {
  return !(std::max(std::abs(coord_arr1.at(kIndexCoordReal) -
                             coord_arr2.at(kIndexCoordReal)),
                    std::abs(coord_arr1.at(kIndexCoordImag) -
                             coord_arr2.at(kIndexCoordImag))) > 1);
} */

/* auto ADHM::handle_self_intersections(
    std::vector<std::vector<path_point>> intersection_candidates,
    std::vector<intersection>& new_intersections,
    std::vector<path_point>::iterator current_pp_it,
    std::vector<Path>::iterator current_path_it) -> bool {
  for (auto& ic : intersection_candidates) {
    for (auto& pp : ic) {
      bool candidate_inserted = false;
      std::array<int32_t, 2> coord_arr{pp.coordinate_real, pp.coordinate_imag};
      for (auto& prev_ic : new_intersections) {
        if (!neighbour_pixel(prev_ic.coordinates.back(), coord_arr)) {
          continue;
        }
        prev_ic.coordinates.push_back(coord_arr);
        prev_ic.times.at(kIndexSecondPath).at(kIndexEndTime) =
            current_pp_it->t.at(kIndexEndTime);
        if (prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) ==
                pp.t.at(kIndexEndTime) + 1 ||
            prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) ==
                pp.t.at(kIndexEndTime)) {
          prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) =
              pp.t.at(kIndexStartTime);
        } else {
          prev_ic.times.at(kIndexFirstPath).at(kIndexEndTime) =
              pp.t.at(kIndexEndTime);
        }
        candidate_inserted = true;
        break;
      }
      if (!candidate_inserted) {
        state_type v1 = current_path_it->get_point(pp.t.at(kIndexStartTime));
        state_type v2 =
            current_path_it->get_point(current_pp_it->t.at(kIndexStartTime));
        v1.at(kIndexX) = path_point_to_complex(pp);
        v2.at(kIndexX) = path_point_to_complex(pp);
        print_state_type(v1);
        print_state_type(v2);
        curve_->match_fiber(v1);
        curve_->match_fiber(v2);
        if (compare_states(v1, v2)) {
          return false;
        }
        intersection new_candidate;
        new_candidate.ids = std::array<uint32_t, 2>{pp.id, current_pp_it->id};
        new_candidate.coordinates.push_back(coord_arr);
        new_candidate.times.at(kIndexFirstPath) = pp.t;
        new_candidate.times.at(kIndexSecondPath) = current_pp_it->t;
        new_intersections.push_back(new_candidate);
      }
    }
  }
} */

auto ADHM::self_intersections(std::vector<Path>::iterator path_it)
    -> std::vector<intersection> {
  SinglePathMap path_map = SinglePathMap(path_it);
  path_map.get_self_intersections();
  return path_map.intersections;
}

auto ADHM::two_path_intersections(std::vector<Path>::iterator path_A_it,
                                  std::vector<Path>::iterator path_B_it)
    -> std::vector<intersection> {
  SinglePathMap path_A_map = SinglePathMap(path_A_it);
  SinglePathMap path_B_map = SinglePathMap(path_B_it);
  Map two_path_map;
  if (path_B_it->path_id_ == 3) {
    for (auto& pp : path_B_map.pp_vec) {
      if (pp.t.at(kIndexStartTime) < 35) {
        continue;
      }
      // print_path_point(pp);
      if (pp.t.at(kIndexStartTime) > 38) {
        break;
      }
    }
  }
  two_path_map.add_path(path_A_map.pp_vec);
  return two_path_map.get_intersections(path_B_map.pp_vec);
}

auto ADHM::add_to_map(std::vector<Path>::iterator path_it) -> void {}

auto ADHM::handle_new_intersections() -> void {}

auto line_intersection(std::array<cplx, 2> A, std::array<cplx, 2> B, cplx& z)
    -> bool {
  // spdlog::debug("Compute line intersection.");
  Eigen::Matrix2d M;
  M(0, 0) = A.at(1).real() - A.at(0).real();
  M(1, 0) = A.at(1).imag() - A.at(0).imag();
  M(0, 1) = -(B.at(1).real() - B.at(0).real());
  M(1, 1) = -(B.at(1).imag() - B.at(0).imag());
  Eigen::Vector2d A0;
  Eigen::Vector2d B0;
  A0 << A.at(0).real(), A.at(0).imag();
  B0 << B.at(0).real(), B.at(0).imag();
  if (M.determinant() != 0) {
    Eigen::Matrix2d M_inv = M.inverse();
    Eigen::Vector2d t = M_inv * (B0 - A0);
    if ((t(0) >= 0 && t(0) <= 1) && (t(1) >= 0 && t(1) <= 1)) {
      z = (A0(0) + t(0) * M(0, 0)) + J * (A0(1) + t(0) * M(1, 0));
      // spdlog::debug("Here is an intersection: {}!", complex_to_string(z));
      return true;
    }
  }
  return false;
}

auto intersect_states(state_type state_A, state_type state_B,
                      state_type& new_state) -> bool {
  if (std::abs(std::exp(state_A.at(kIndexY1)) -
               std::exp(state_B.at(kIndexY2))) < kFiberCompTolerance) {
    new_state.at(kIndexX) = state_A.at(kIndexX);
    new_state.at(kIndexY1) = std::log(std::exp(state_B.at(kIndexY1)));
    new_state.at(kIndexY2) = std::log(std::exp(state_A.at(kIndexY2)));
    return true;
  }
  if (std::abs(std::exp(state_A.at(kIndexY2)) -
               std::exp(state_B.at(kIndexY1))) < kFiberCompTolerance) {
    new_state.at(kIndexX) = state_A.at(kIndexX);
    new_state.at(kIndexY1) = std::log(std::exp(state_A.at(kIndexY1)));
    new_state.at(kIndexY2) = std::log(std::exp(state_B.at(kIndexY2)));
    return true;
  }
  if (std::abs(std::exp(state_A.at(kIndexY1)) -
               std::exp(state_B.at(kIndexY1))) < kFiberCompTolerance &&
      std::abs(std::exp(state_A.at(kIndexY2)) -
               std::exp(state_B.at(kIndexY2))) < kFiberCompTolerance) {
    new_state.at(kIndexX) = state_A.at(kIndexX);
    new_state.at(kIndexY1) = std::log(std::exp(state_A.at(kIndexY1)));
    new_state.at(kIndexY2) = std::log(std::exp(state_B.at(kIndexY1)));
    return true;
  }
  print_state_type(state_A);
  print_state_type(state_B);
  return false;
}

auto ADHM::compute_intersection_points(intersection& inter,
                                       std::vector<Path>::iterator path_A_it,
                                       std::vector<Path>::iterator path_B_it,
                                       int32_t n, state_type& new_state)
    -> bool {
  cplx z;
  // uint32_t A_id{inter.ids.at(kIndexFirstPath)};
  // uint32_t B_id{inter.ids.at(kIndexSecondPath)};

  uint32_t A_start_time{inter.times.at(kIndexFirstPath).at(kIndexStartTime)};
  uint32_t A_end_time{inter.times.at(kIndexFirstPath).at(kIndexEndTime)};

  uint32_t B_start_time{inter.times.at(kIndexSecondPath).at(kIndexStartTime)};
  uint32_t B_end_time{inter.times.at(kIndexSecondPath).at(kIndexEndTime)};

  std::array<cplx, 2> pt_A;
  std::array<cplx, 2> pt_B;

  for (uint32_t t_A = A_start_time - 10; t_A < A_end_time + 10; t_A++) {
    if (t_A < 0) {
      continue;
    }
    if (t_A > path_A_it->get_length()) {
      return false;
    }

    for (uint32_t t_B = B_start_time - 10; t_B < B_end_time + 10; t_B++) {
      if (t_B < 0) {
        continue;
      }
      if (t_B > path_B_it->get_length()) {
        return false;
      }

      pt_A = {path_A_it->get_point(t_A).at(kIndexX),
              path_A_it->get_point(t_A + 1).at(kIndexX)};
      pt_B = {path_B_it->get_point(t_B).at(kIndexX),
              path_B_it->get_point(t_B + 1).at(kIndexX)};

      if (line_intersection(pt_A, pt_B, z)) {
        spdlog::debug("Line intersection found at {}", complex_to_string(z));
        state_type state_A = path_A_it->get_point(t_A);
        state_type state_B = path_B_it->get_point(t_B);
        inter.times.at(kIndexFirstPath).at(kIndexStartTime) = t_A;
        inter.times.at(kIndexFirstPath).at(kIndexEndTime) = t_A + 1;
        inter.times.at(kIndexSecondPath).at(kIndexStartTime) = t_B;
        inter.times.at(kIndexSecondPath).at(kIndexEndTime) = t_B + 1;
        state_A.at(kIndexX) = z;
        state_B.at(kIndexX) = z;
        curve_->match_fiber(state_A);
        curve_->match_fiber(state_B);
        if (intersect_states(state_A, state_B, new_state)) {
          new_state.at(kIndexY1) += 2 * pi * J * static_cast<double>(n);
          spdlog::debug("Return state!");
          return true;
        }
      }
    }
  }
  return false;
}
