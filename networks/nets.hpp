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
#include "sw_curve.hpp"
#include "magic_numbers.h"
#include "ode_integrator.hpp"
#include "type_util.hpp"
#include "format_util.hpp"

class Network {
    public:
        Network(GiNaC::ex (*func)(const GiNaC::symbol&, const GiNaC::symbol&), double theta) 
            : theta_(theta), curve_(), ramification_points_(), new_paths_(), evolved_paths_()  {
                curve_ = std::make_shared<SW_curve>(func);
                ramification_points_ = curve_->get_ramification_points();
                start_path();
            }

        // IO methods.
        auto print_ramification_points() -> void;
        auto evolution_step()-> void;

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
        std::vector<Path> new_paths_;
        std::vector<Path> evolved_paths_;    
};

#endif  // NETS_HPP_