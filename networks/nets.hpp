#ifndef NETS_HPP_
#define NETS_HPP_

#include <complex>
#include <cmath>
#include <vector>
#include <thread>
#include <memory>
#include <iostream>
#include <cstdint>
#include <ginac/ginac.h>
#include <spdlog/spdlog.h>
#include <fmt/core.h>

#include "path.hpp"
#include "maps.hpp"
#include "sw_curve.hpp"
#include "magic_numbers.h"
#include "ode_integrator.hpp"
#include "type_util.hpp"

class Network {
    public:
        Network(GiNaC::ex (*func)(const GiNaC::symbol&, const GiNaC::symbol&), double theta) 
            : theta_(theta), curve_(), ramification_points_(), new_paths_(), evolved_paths_() {
                curve_ = std::make_shared<SW_curve>(func);
                ramification_points_ = curve_->get_ramification_points();
                start_path();
            }

        // IO methods.
        auto print_ramification_points() -> void;
        auto evolution_step()-> void;
        auto draw_map()->void;
        auto save_intersections()->void;

    private:
        // Network parameter.
        double theta_;
        
        // SW_curve data.
        std::shared_ptr<SW_curve> curve_;
        std::vector<std::array<cplx,2>> ramification_points_;

        // Path stuff.
        auto start_path() -> void;
        auto determine_sign(const state_type &r, state_type &v) -> void;
        
        auto move_to_evolved(std::vector<Path>::iterator) -> void;
        auto add_new_path() -> void;
        std::vector<Path> new_paths_;
        std::vector<Path> evolved_paths_;

        // Map stuff.
        Map map_;

        // Intersection stuff.
        std::vector<intersection> new_intersections_;
        std::vector<state_type> computed_intersections_;
        auto handle_new_intersections(std::vector<std::vector<path_point>> intersection_candidates, std::vector<path_point>::iterator current_pp_it, std::vector<Path>::iterator current_path_it) -> bool;
        auto compute_intersection_points() -> void;
};

#endif  // NETS_HPP_
