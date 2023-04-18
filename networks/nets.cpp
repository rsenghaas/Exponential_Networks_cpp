#include "nets.hpp"

#include <fstream>
#include <sstream>

#include "Eigen/Dense"
#include "magic_numbers.h"

// * Log has range (-i*pi,i*pi]

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

auto Network::get_iterator_by_id(std::vector<Path>& path_vec, uint32_t id)
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

auto Network::start_paths() -> void {
  for (auto& r : ramification_points_) {
    if(std::abs(r.at(kIndexX)) < 0.00001) {
      continue;
    }
    cplx b = r.at(kIndexX);
    spdlog::debug("Branch point at x = {}.", complex_to_string(b));
    cplx y = r.at(kIndexY);
    
    state_type start_state;
    start_state.at(kIndexX) = b;
    start_state.at(kIndexY1) = std::log(y);
    start_state.at(kIndexY2) = std::log(y);

    // Compute the expansion around a branch point.
    // TODO: I should move this to the SW_curve class and compute it generally
    // (kappa = 1 / n! d^n H / dy^n * (1 / d H / dx))
    cplx kappa =
        -1.0 / 2 * curve_->eval_d2H_dy2(b, y) / curve_->eval_dH_dx(b, y);
    // spdlog::debug("kappa = {}", complex_to_string(kappa));

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

auto Network::evolve_path(std::vector<Path>::iterator path_it, double cutoff) -> void {
  future_type future = path_it->integrate(curve_, theta_, cutoff);
  path_it->update(std::ref(future));
}

auto Network::add_new_path(state_type start_point) -> void {
  std::vector<state_type> v;
  v.push_back(start_point);

  std::vector<double> masses;
  masses.push_back(0);

  Path path(v, masses, next_id_);
  new_paths_.push_back(std::move(path));
  next_id_++;
}

auto Network::save_data(uint32_t id) -> void {
  auto path_it = get_iterator_by_id(new_paths_, id);
  path_it->save_data();
}

auto Network::self_intersection_handler(uint32_t id, bool truncate, int32_t n,
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
      auto endpoint = path_it->get_endpoint();
      endpoint.at(kIndexX) = next_state.at(kIndexX);
      curve_->match_fiber(endpoint);
      path_it->add_single_point(endpoint);
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

auto Network::self_intersections(std::vector<Path>::iterator path_it)
    -> std::vector<intersection> {
  SinglePathMap path_map = SinglePathMap(path_it);
  path_map.get_self_intersections();
  return path_map.intersections;
}

auto Network::two_path_intersections(std::vector<Path>::iterator path_A_it,
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

auto Network::two_path_intersection_handler(uint32_t id_A, uint32_t id_B,
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
      path_2_end_partner_ = id_B;
      path_2_endtime_ = inter_it->times.at(kIndexFirstPath).at(kIndexEndTime);
    }
    if (id_B == 2 && inter_it->times.at(kIndexSecondPath).at(kIndexEndTime) >
                         path_2_endtime_) {
      path_2_end_partner_ = id_A;
      path_2_endtime_ = inter_it->times.at(kIndexSecondPath).at(kIndexEndTime);
    }

    if (truncate_A) {
      path_A_it->truncate(
          0, inter_it->times.at(kIndexFirstPath).at(kIndexEndTime));
      auto endpoint = path_A_it->get_endpoint();
      endpoint.at(kIndexX) = next_state.at(kIndexX);
      curve_->match_fiber(endpoint);
      path_A_it->add_single_point(endpoint);
    }
    if (truncate_B) {
      path_B_it->truncate(
          0, inter_it->times.at(kIndexSecondPath).at(kIndexEndTime));
      auto endpoint = path_B_it->get_endpoint();
      endpoint.at(kIndexX) = next_state.at(kIndexX);
      curve_->match_fiber(endpoint);
      path_B_it->add_single_point(endpoint);
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

auto Network::compute_intersection_points(intersection& inter,
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

// One shoudl probably restrict that to groups of 4 at a time.
/* auto Network::evolution_step() -> void {
  std::vector<future_type> futures;
  for (auto& new_path : new_paths_) {
    futures.emplace_back(new_path.integrate(curve_, theta_, kD4Cutoff));
  }
  auto new_path_it = new_paths_.begin();
  for (auto& f : futures) {
    if (new_path_it == new_paths_.end()) {
      spdlog::debug("Iterator in future loop overshoots!");
    }
    new_path_it->update(std::ref(f));
    new_path_it->compute_map_points();
    for (auto pp_it = new_path_it->pp_vec.begin();
         pp_it != std::prev(new_path_it->pp_vec.end()); ++pp_it) {
      if (!handle_new_intersections(map_.draw_line(*pp_it, *std::next(pp_it)),
                                    pp_it, new_path_it)) {
        spdlog::debug("Path {} crashed into other path", new_path_it->path_id_);
        break;
      }
    }

    new_path_it->save_data();
    evolved_paths_.push_back(std::move(*new_path_it));

    new_path_it++;
  }
  // map_.print_map_data();
  compute_intersection_points_old();
  // new_paths_.clear();
  // add_new_paths();
} */

// WARN: This is somewhat misplaced.
// WARN: This is defined twice (also in maps.cpp).

auto net_neighbour_pixel(std::array<int32_t, 2> coord_arr1,
                     std::array<int32_t, 2> coord_arr2) -> bool {
  return !(std::max(std::abs(coord_arr1.at(kIndexCoordReal) -
                             coord_arr2.at(kIndexCoordReal)),
                    std::abs(coord_arr1.at(kIndexCoordImag) -
                             coord_arr2.at(kIndexCoordImag))) > 1);
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

// This will either become a bool or and uint32_t to return some flags.
auto Network::handle_new_intersections(
    std::vector<std::vector<path_point>> intersection_candidates,
    std::vector<path_point>::iterator current_pp_it,
    std::vector<Path>::iterator current_path_it) -> bool {
  for (auto& ic : intersection_candidates) {
    spdlog::debug("Intersection candidate at ({},{}) where x = {}.",
                  ic.front().coordinate_real, ic.front().coordinate_imag,
                  complex_to_string(path_point_to_complex(ic.at(0))));
    for (auto& pp : ic) {
      // spdlog::debug("Found intersection of {} and {} found", pp.id,
      // current_pp_it->id);
      bool candidate_inserted = false;
      std::array<int32_t, 2> coord_arr{pp.coordinate_real, pp.coordinate_imag};
      for (auto& prev_ic : new_intersections_) {
        if (pp.id != prev_ic.ids.at(kIndexFirstPath) ||
            current_pp_it->id != prev_ic.ids.at(kIndexSecondPath)) {
          // spdlog::debug("Id check failed!");
          // print_intersection(prev_ic);
          // print_path_point(pp);
          // print_path_point(*current_pp_it);
          continue;
        }
        // spdlog::debug("Id check passed.");
        if (!net_neighbour_pixel(prev_ic.coordinates.back(), coord_arr)) {
          // spdlog::debug("Neighbour check failed.");
          // print_intersection(prev_ic);
          // print_path_point(pp);
          // print_path_point(*current_pp_it);++ program in python
          continue;
        }
        // spdlog::debug("Neighbour check passed.");
        if ((prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) !=
                 pp.t.at(kIndexEndTime) + 1 &&
             prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) !=
                 pp.t.at(kIndexEndTime) &&
             prev_ic.times.at(kIndexFirstPath).at(kIndexEndTime) !=
                 pp.t.at(kIndexStartTime) - 1 &&
             prev_ic.times.at(kIndexFirstPath).at(kIndexEndTime) !=
                 pp.t.at(kIndexStartTime)) ||
            prev_ic.times.at(kIndexSecondPath).at(kIndexEndTime) !=
                current_pp_it->t.at(kIndexStartTime) - 1) {
          // spdlog::debug("Times check failed.");
          // print_intersection(prev_ic);
          // print_path_point(pp);
          // print_path_point(*current_pp_it);
          // continue;
        }
        // spdlog::debug("Time check passed.");
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
      // The Path id's are always ordered by size.
      if (!candidate_inserted) {
        state_type v1;
        spdlog::debug("Start initialize path points.");
        if (evolved_paths_.size() == pp.id) {
          v1 = current_path_it->get_point(pp.t.at(kIndexStartTime));
        } else {
          v1 = evolved_paths_.at(pp.id).get_point(pp.t.at(kIndexStartTime));
        }
        state_type v2 =
            current_path_it->get_point(current_pp_it->t.at(kIndexStartTime));
        v1.at(kIndexX) = path_point_to_complex(pp);
        v2.at(kIndexX) = path_point_to_complex(pp);
        spdlog::debug("Matching fibers.");
        print_state_type(v1);
        print_state_type(v2);
        curve_->match_fiber(v1);
        curve_->match_fiber(v2);
        spdlog::debug("After matching");
        print_state_type(v1);
        print_state_type(v2);
        spdlog::debug("Compare paths {} and {}.", pp.id, current_pp_it->id);
        if (compare_states(v1, v2)) {
          spdlog::debug("Overlapping paths.");
          print_path_point(pp);
          print_path_point(*current_pp_it);
          // return false;
          return true;
        }
        spdlog::debug("Candidate created.");
        intersection new_candidate;
        new_candidate.ids = std::array<uint32_t, 2>{pp.id, current_pp_it->id};
        new_candidate.coordinates.push_back(coord_arr);
        new_candidate.times.at(kIndexFirstPath) = pp.t;
        new_candidate.times.at(kIndexSecondPath) = current_pp_it->t;
        new_intersections_.push_back(new_candidate);
      }
    }
  }
  return true;
}

auto Network::add_new_paths() -> void {
  uint32_t id = evolved_paths_.size();
  for (auto& cpi : computed_intersections_) {
    state_type state_A = cpi.states.at(0);
    state_type state_B = cpi.states.at(1);
    state_type new_state;
    if (std::abs(std::exp(state_A.at(kIndexY1)) -
                 std::exp(state_B.at(kIndexY2))) < kFiberCompTolerance) {
      int32_t k = static_cast<int32_t>(std::round(
          ((state_A.at(kIndexY1) - state_B.at(kIndexY2)) / (2 * pi * J))
              .real()));
      new_state.at(kIndexX) = state_A.at(kIndexX);
      new_state.at(kIndexY1) =
          state_B.at(kIndexY1) + 2 * pi * J * static_cast<double>(k);
      new_state.at(kIndexY2) = state_A.at(kIndexY2);

    } else if (std::abs(std::exp(state_B.at(kIndexY1)) -
                        std::exp(state_A.at(kIndexY2))) < kFiberCompTolerance) {
      int32_t k = static_cast<int32_t>(std::round(
          ((state_B.at(kIndexY1) - state_A.at(kIndexY2)) / (2 * pi * J))
              .real()));
      new_state.at(kIndexX) = state_A.at(kIndexX);
      new_state.at(kIndexY1) =
          state_A.at(kIndexY2) + 2 * pi * J * static_cast<double>(k);
      new_state.at(kIndexY2) = state_B.at(kIndexY1);
    }
    //! This looks kind of ugly..., direction is not canonical at all!
    else if (std::abs(std::exp(state_A.at(kIndexY1)) -
                      std::exp(state_B.at(kIndexY1))) < kFiberCompTolerance) {
      int32_t k = static_cast<int32_t>(std::round(
          ((state_A.at(kIndexY1) - state_B.at(kIndexY1)) / (2 * pi * J))
              .real()));
      if (k != 0) {
        new_state.at(kIndexX) = state_A.at(kIndexX);
        new_state.at(kIndexY1) =
            state_A.at(kIndexY2) + 2 * pi * J * static_cast<double>(k);
        new_state.at(kIndexY2) = state_A.at(kIndexY2);
      } else {
        continue;
      }
    } else {
      continue;
    }
    std::vector<state_type> v{new_state};
    std::vector<double> m{0};
    Path path(v, m, id);
    new_paths_.push_back(std::move(path));
    id++;
  }
  computed_intersections_.clear();
}

auto save_string_to_file(std::string filename, std::string s) {
  std::fstream data_file;
  data_file.open(filename, std::ios::out);
  if (!data_file) {
    spdlog::debug("{} could not be created.", filename);
  } else {
    data_file << s;
  }
  spdlog::debug("Data saved to {}.", filename);
  data_file.close();
}

auto Network::compute_intersection_points_old() -> void {
  spdlog::debug("Computing the precise intersections!");
  std::string s = "";
  uint32_t found_intersections = 0;
  uint32_t unsuccess_id = 0;
  for (auto& ic : new_intersections_) {
    // print_intersection(ic);
    cplx z;
    uint32_t id_A = ic.ids.at(kIndexFirstPath);
    uint32_t id_B = ic.ids.at(kIndexSecondPath);

    uint32_t A_times0 = ic.times.at(kIndexFirstPath).at(kIndexStartTime);
    if (A_times0 != 0) {
      A_times0 -= 2;
    }
    uint32_t A_times1 = ic.times.at(kIndexFirstPath).at(kIndexEndTime) + 2;
    uint32_t B_times0 = ic.times.at(kIndexSecondPath).at(kIndexStartTime);
    if (B_times0 != 0) {
      B_times0 -= 2;
    }
    uint32_t B_times1 = ic.times.at(kIndexSecondPath).at(kIndexEndTime) + 2;
    // spdlog::debug("Range for A: [{},{}], B: [{},{}]", A_times0, A_times1,
    // B_times0, B_times1);
    bool success = [&]() -> bool {
      for (uint32_t t_A = A_times0; t_A <= A_times1; t_A++) {
        for (uint32_t t_B = B_times0; t_B <= B_times1; t_B++) {
          std::array<cplx, 2> A = {
              evolved_paths_.at(id_A).get_point(t_A).at(kIndexX),
              evolved_paths_.at(id_A).get_point(t_A + 1).at(kIndexX)};
          std::array<cplx, 2> B = {
              evolved_paths_.at(id_B).get_point(t_B).at(kIndexX),
              evolved_paths_.at(id_B).get_point(t_B + 1).at(kIndexX)};
          if (line_intersection(A, B, z)) {
            state_type new_state_A = evolved_paths_.at(id_A).get_point(t_A);
            state_type new_state_B = evolved_paths_.at(id_B).get_point(t_B);
            new_state_A.at(kIndexX) = z;
            new_state_B.at(kIndexX) = z;
            curve_->match_fiber(new_state_A);
            curve_->match_fiber(new_state_B);
            computed_intersection new_intersection;
            new_intersection.states = {new_state_A, new_state_B};
            new_intersection.ids = {id_A, id_B};
            new_intersection.times = {t_A, t_B};
            computed_intersections_.push_back(new_intersection);
            if (s.size() != 0) {
              s.append(",");
            }
            s.append(complex_to_string(z));
            found_intersections++;
            return true;  // Return from lambda; this just exits the double
                          // loop!
          }
        }
      }
      return false;
    }();
    if (!success) {
      spdlog::debug("No intersection found!");
      print_intersection(ic);
      spdlog::debug("Points for path A.");
      std::string s_A = "";
      std::string s_B = "";
      for (uint32_t t_A = A_times0; t_A <= A_times1 + 1; t_A++) {
        if (s_A.size() != 0) {
          s_A.append(",");
        }
        cplx A = evolved_paths_.at(id_A).get_point(t_A).at(kIndexX);
        spdlog::debug("{}", complex_to_string(A));
        s_A.append(complex_to_string(A));
      }
      spdlog::debug("Points for path B.");
      for (uint32_t t_B = B_times0; t_B <= B_times1 + 1; t_B++) {
        if (s_B.size() != 0) {
          s_B.append(",");
        }
        cplx B = evolved_paths_.at(id_B).get_point(t_B).at(kIndexX);
        spdlog::debug("{}", complex_to_string(B));
        s_B.append(complex_to_string(B));
      }
      std::string s_combined;
      /*s_combined = std::format("{}\n{}", s_A, s_B);
      save_string_to_file(
          std::format("data/intersection_data/unsuccessful_{}.csv",
                      unsuccess_id),
          s_combined);*/
      spdlog::debug("Missing intersection #{} found", unsuccess_id);
      unsuccess_id++;
    }
  }
  spdlog::debug("Intersection candidates :{}, Intersections found: {}",
                new_intersections_.size(), found_intersections);
  save_string_to_file("data/intersection_data/test.csv", s);
}

/* auto Network::determine_sign(const state_type& r, state_type& v) -> void {
  state_type dv;
  curve_->sw_differential(v, dv);
  cplx dx = v.at(kIndexX) - r.at(kIndexX);
  cplx next_dx = dv.at(kIndexX);
  if (next_dx.real() * dx.real() + next_dx.imag() * dx.imag() < 0) {
    cplx temp = v.at(kIndexY1);
    v.at(kIndexY1) = v.at(kIndexY2);
    v.at(kIndexY2) = temp;
  }
} */

auto Network::print_ramification_points() -> void {
  for (auto& r : ramification_points_) {
    std::cout << "The curve is ramified over x = " 
              << complex_to_string(r.at(kIndexX)) 
              << "with y = " 
              << complex_to_string(r.at(kIndexY))
              << std::endl;
  }
}
