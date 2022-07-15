#include "nets.hpp"

#include <sstream>
#include <fstream>

#include "Eigen/Dense"

// * Log has range (-i*pi,i*pi]

auto Network::start_path() -> void {
    for (auto &r : ramification_points_) {
        cplx b = r.at(kIndexX);
        cplx y = r.at(kIndexY);
        spdlog::debug("Numerical Check: H(b, y) = {}.", 
            complex_to_string(curve_->eval_H(b, y)));
        
        state_type start_state;
        start_state.at(kIndexX) = b;
        start_state.at(kIndexY1) = std::log(y);
        start_state.at(kIndexY2) = std::log(y);
        
        // Compute the expansion around a branch point.
        // TODO: I should move this to the SW_curve class and compute it generally (kappa = 1 / n! d^n H / dy^n * (1 / d H / dx))
        cplx kappa = -1.0 / 2 * curve_->eval_d2H_dy2(b, y) / curve_->eval_dH_dx(b,y);
        // spdlog::debug("kappa = {}", complex_to_string(kappa));
        
        cplx dx;
        cplx dy;
        for(uint32_t k = 0; k < 3; k++) {
            std::vector<state_type> v;
            std::vector<double> masses;
            state_type next_state;

            v.push_back(start_state);
            masses.push_back(0);

            dx = std::pow(std::pow(3.0 / 4 * std::pow(kappa, 1.0/2) * y * b * std::exp(J*theta_) * kInitialStepSize, 2.0), 1.0 / 3) * std::pow(kZeta3, k);  // ! This is only valid for exponential networks.
            dy = 3.0 / 4.0 * y * b * std::exp(J*theta_) * kInitialStepSize / dx;

            next_state.at(kIndexX) = b + dx;
            next_state.at(kIndexY1) = start_state.at(kIndexY1) + dy / y;
            next_state.at(kIndexY2) = start_state.at(kIndexY2) - dy / y;
            curve_->match_fiber(v.back());

            determine_sign(start_state, next_state);
            v.push_back(next_state);
            
            spdlog::debug("Numerical Check: \n |H(x, y_1)| = {}\n |H(x, y_2)| = {}.", 
                complex_to_string(std::abs(curve_->eval_H(v.back().at(kIndexX), std::exp(v.back().at(kIndexY1))))), 
                complex_to_string(std::abs(curve_->eval_H(v.back().at(kIndexX), std::exp(v.back().at(kIndexY2)))))); 

            cplx dlog_y1 = next_state.at(kIndexY1) - start_state.at(kIndexY1);
            cplx dlog_y2 = next_state.at(kIndexY2) - start_state.at(kIndexY2);
            state_type dv{dx, dlog_y1,dlog_y2};
            masses.push_back(compute_dm(next_state, dv));
            // First step is done.
  
            // Doing some runge kutta steps before going into the boost::odeint integration.
            for(uint32_t i = 0; i < kInitialSteps; i++) {
                ODE_runge_kutta_step(curve_, v, masses, kInitialStepSize, theta_);
            }

            Path path(v, masses, new_paths_.size());
            new_paths_.push_back(std::move(path));
            
            spdlog::debug("Path {} appended.", new_paths_.back().id_);
            spdlog::debug("Numerical Check: \n |H(x, y_1)| = {}\n |H(x, y_2)| = {}.", 
            complex_to_string(std::abs(curve_->eval_H(v.back().at(kIndexX), std::exp(v.back().at(kIndexY1))))), 
            complex_to_string(std::abs(curve_->eval_H(v.back().at(kIndexX), std::exp(v.back().at(kIndexY2))))));   
  
        }
        spdlog::debug("Paths are started.");
    }
}

auto Network::evolution_step() -> void {
    std::vector<future_type> futures;
    for(auto & new_path : new_paths_) {
        futures.emplace_back(new_path.integrate(curve_, theta_));
    }
    auto new_path_it = new_paths_.begin();
    for (auto &f : futures) {
        if (new_path_it == new_paths_.end()) {
            spdlog::debug("Iterator in future loop overshoots!");
        }
        new_path_it->update(std::ref(f));
        new_path_it->compute_map_points();
        for (auto pp_it = new_path_it->pp_vec.begin(); pp_it != std::prev(new_path_it->pp_vec.end()); ++pp_it) {
            if(!handle_new_intersections(map_.draw_line(*pp_it, *std::next(pp_it)), pp_it, new_path_it)) {
                spdlog::debug("Path {} crashed into other path", new_path_it->id_);
                break;
            }
        }

        new_path_it->save_data();
        evolved_paths_.push_back(std::move(*new_path_it));
        
        new_path_it++;
    }
    // map_.print_map_data();
    compute_intersection_points();
    new_paths_.clear();
}



//! This is somewhat misplaced.
auto neighbour_pixel(std::array<int32_t, 2> coord_arr1, std::array<int32_t, 2> coord_arr2) -> bool {
    if(std::max(std::abs(coord_arr1.at(kIndexCoordReal) - coord_arr2.at(kIndexCoordReal)), 
                std::abs(coord_arr1.at(kIndexCoordImag) - coord_arr2.at(kIndexCoordImag))) > 1) {
        return false;
    }
    else {
        return true;
    }
}

auto compare_states(state_type v1, state_type v2) -> bool {
  if (std::abs((v1.at(kIndexY1) - v1.at(kIndexY2)) - (v2.at(kIndexY2) - v2.at(kIndexY1))) < kFiberCompTolerance) {
    return true;
  }
  if (std::abs((v1.at(kIndexY1) - v1.at(kIndexY2)) - (v2.at(kIndexY1) - v2.at(kIndexY2))) < kFiberCompTolerance) {
    return true;
  }
  return false;
}

// This will either become a bool or and uint32_t to return some flags.
auto Network::handle_new_intersections(std::vector<std::vector<path_point>> intersection_candidates, std::vector<path_point>::iterator current_pp_it, std::vector<Path>::iterator current_path_it) -> bool {
    for (auto& ic : intersection_candidates) {
      spdlog::debug("Intersection candidate at ({},{}) where x = {}.", ic.front().coordinate_real, ic.front().coordinate_imag, complex_to_string(path_point_to_complex(ic.at(0))));
      for (auto& pp : ic) {
      // spdlog::debug("Found intersection of {} and {} found", pp.id, current_pp_it->id);
        bool candidate_inserted = false;
        std::array<int32_t, 2> coord_arr{pp.coordinate_real, pp.coordinate_imag};
        for (auto& prev_ic : new_intersections_) {
            if (pp.id != prev_ic.ids.at(kIndexFirstPath) || current_pp_it->id != prev_ic.ids.at(kIndexSecondPath)) {
                // spdlog::debug("Id check failed!");
                // print_intersection(prev_ic);
                // print_path_point(pp);
                // print_path_point(*current_pp_it);
                continue;
            }
            // spdlog::debug("Id check passed.");
            if (!neighbour_pixel(prev_ic.coordinates.back(), coord_arr)) {
                // spdlog::debug("Neighbour check failed.");
                // print_intersection(prev_ic);
                // print_path_point(pp);
                // print_path_point(*current_pp_it);
                continue;
            }
            // spdlog::debug("Neighbour check passed.");
            if ((prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) != pp.t.at(kIndexEndTime) + 1 && prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) != pp.t.at(kIndexEndTime)
                  && prev_ic.times.at(kIndexFirstPath).at(kIndexEndTime) != pp.t.at(kIndexStartTime) - 1 && prev_ic.times.at(kIndexFirstPath).at(kIndexEndTime) != pp.t.at(kIndexStartTime))
                  || prev_ic.times.at(kIndexSecondPath).at(kIndexEndTime) != current_pp_it->t.at(kIndexStartTime) -1) 
            {
                // spdlog::debug("Times check failed.");
                // print_intersection(prev_ic);
                // print_path_point(pp);
                // print_path_point(*current_pp_it);
                // continue;
            }
            // spdlog::debug("Time check passed.");
            prev_ic.coordinates.push_back(coord_arr);
            prev_ic.times.at(kIndexSecondPath).at(kIndexEndTime) = current_pp_it->t.at(kIndexEndTime);
            if (prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) == pp.t.at(kIndexEndTime) + 1 || prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) == pp.t.at(kIndexEndTime)) {
                prev_ic.times.at(kIndexFirstPath).at(kIndexStartTime) = pp.t.at(kIndexStartTime);
            } else {
                prev_ic.times.at(kIndexFirstPath).at(kIndexEndTime) = pp.t.at(kIndexEndTime);
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
            state_type v2 = current_path_it->get_point(current_pp_it->t.at(kIndexStartTime));
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
            spdlog::debug("Compare paths.");
            if (compare_states(v1, v2)) {
                spdlog::debug("\nOverlapping paths.\n");
                print_path_point(pp);
                print_path_point(*current_pp_it);
                return false;
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

auto Network::add_new_path() -> void {

}

auto save_string_to_file(std::string filename, std::string s) {
    std::fstream data_file;
    data_file.open(filename, std::ios::out);
    if(!data_file) {
      spdlog::debug("{} could not be created.", filename);
    } else {
        data_file << s;
    }
    spdlog::debug("Data saved to {}.", filename);
    data_file.close();
}



auto line_intersection(std::array<cplx,2> A, std::array<cplx,2> B, cplx& z) -> bool {
    spdlog::debug("Compute line intersection.");
    Eigen::Matrix2d M;
    M(0,0) = A.at(1).real() - A.at(0).real();
    M(1,0) = A.at(1).imag() - A.at(0).imag();
    M(0,1) = B.at(1).real() - B.at(0).imag();
    M(1,1) = B.at(1).imag() - B.at(0).imag();
    Eigen::Vector2d A0;
    Eigen::Vector2d B0;
    A0 << A.at(0).real(), A.at(0).imag();
    B0 << B.at(0).real(), B.at(0).imag();
    if (M.determinant() != 0) {
        Eigen::Matrix2d M_inv = M.inverse();
        Eigen::Vector2d t = M_inv * (B0 - A0);
        if ((t(0) >= 0 && t(0) <=1) && (t(1) >=0 && t(1) <= 1)) {
            z = (A0(0) + t(0) * M(0,0)) + J * (A0(1) + t(1) * M(1,0));
            spdlog::debug("Here is an intersection: {}!", complex_to_string(z));
            return true;
        }
    }
    return false;
}

auto Network::compute_intersection_points() -> void {
    spdlog::debug("Computing the precise intersections!");
    std::string s = "";
    for (auto& ic : new_intersections_) {
        // print_intersection(ic);
        cplx z;
        uint32_t id_A = ic.ids.at(kIndexFirstPath);
        uint32_t id_B = ic.ids.at(kIndexSecondPath);

        uint32_t A_times0 = ic.times.at(kIndexFirstPath).at(kIndexStartTime);
        if (A_times0 != 0) {A_times0--;}
        uint32_t A_times1 = ic.times.at(kIndexFirstPath).at(kIndexEndTime) + 1;
        uint32_t B_times0 = ic.times.at(kIndexSecondPath).at(kIndexStartTime);
        if (B_times0 != 0) {B_times0--;}
        uint32_t B_times1 = ic.times.at(kIndexSecondPath).at(kIndexEndTime) + 1;
        spdlog::debug("Range for A: [{},{}], B: [{},{}]", A_times0, A_times1, B_times0, B_times1); 
        [&] { 
            spdlog::debug("In lambda function.");
            for (uint32_t t_A = A_times0; t_A <= A_times1; t_A++) {
                for (uint32_t t_B = B_times0; t_B <= B_times1; t_B++) {
                  std::array<cplx,2> A = {evolved_paths_.at(id_A).get_point(t_A).at(kIndexX), evolved_paths_.at(id_A).get_point(t_A + 1).at(kIndexX)};
                  std::array<cplx,2> B = {evolved_paths_.at(id_B).get_point(t_B).at(kIndexX), evolved_paths_.at(id_B).get_point(t_B + 1).at(kIndexX)};
                  spdlog::debug("In double loop");
                  if (line_intersection(A, B, z)) {
                      if (s.size() != 0) {
                          s.append(",");
                      }
                      s.append(complex_to_string(z));
                      
                      spdlog::debug("Return from double loop");
                      return;  // Return from lambda; this just exits the double loop!
                  }
                }
            }
        }();
    }
    save_string_to_file("data/intersection_data/test.csv", s);
}

auto Network::determine_sign(const state_type &r, state_type &v) -> void {
    state_type dv;
    curve_->sw_differential(v, dv);
    cplx dx = v.at(kIndexX) - r.at(kIndexX);
    cplx next_dx = dv.at(kIndexX);
    if (next_dx.real() * dx.real() + next_dx.imag() * dx.imag() < 0) {
        cplx temp = v.at(kIndexY1);
        v.at(kIndexY1) = v.at(kIndexY2);
        v.at(kIndexY2) = temp;
    }
}

auto Network::print_ramification_points() -> void {
    for (auto& r : ramification_points_) {
        std::cout 
            << fmt::format("The curve is ramified over x = {} at y = {}.", 
            complex_to_string(r.at(kIndexX)), 
            complex_to_string(r.at(kIndexY))) 
            << std::endl;
    }
}

