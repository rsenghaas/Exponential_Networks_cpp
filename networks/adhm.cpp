#include "adhm.hpp"

#include <fmt/core.h>
#include <ginac/ginac.h>
#include <spdlog/spdlog.h>

#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

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

auto ADHM::BPS_state() -> void {
  auto path_it = new_paths_.begin();
  evolve_path(path_it, kD4Cutoff);
  // path_it->save_data();

  // state_type end_point = path_it->get_endpoint();
  // print_state_type(end_point);
  // std::swap(end_point.at(kIndexY1), end_point.at(kIndexY2));
  // print_state_type(end_point);
  // end_point.at(kIndexY1);

  path_it++;
  evolve_path(path_it, kD4Cutoff);
  path_it->save_data();

  path_it++;
  evolve_path(path_it, kD4Cutoff);
  path_it->save_data();

  new_paths_.clear();
  // add_new_path(end_point, 3);
}

auto ADHM::add_new_path(state_type start_point, uint32_t path_id) -> void {
  std::vector<state_type> v;
  v.push_back(start_point);

  std::vector<double> masses;
  masses.push_back(0);

  Path path(v, masses, path_id);
  new_paths_.push_back(std::move(path));
}

auto ADHM::evolve_path(std::vector<Path>::iterator path_it, double cutoff)
    -> void {
  future_type future;
  future = path_it->integrate(curve_, theta_, cutoff);
  path_it->update(std::ref(future));
}

auto ADHM::draw_map(std::vector<Path>::iterator path_it) -> void {
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
}

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

auto neighbour_pixel(std::array<int32_t, 2> coord_arr1,
                     std::array<int32_t, 2> coord_arr2) -> bool {
  return !(std::max(std::abs(coord_arr1.at(kIndexCoordReal) -
                             coord_arr2.at(kIndexCoordReal)),
                    std::abs(coord_arr1.at(kIndexCoordImag) -
                             coord_arr2.at(kIndexCoordImag))) > 1);
}

auto ADHM::handle_self_intersections(
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
}

auto ADHM::handle_new_intersections() -> void {}
