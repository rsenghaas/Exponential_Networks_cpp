#include "nets.hpp"

// * Log has range (-i*pi,i*pi]

auto Network::start_path() -> void {
    for (auto &r : ramification_points_) {
        cplx b = r.at(kIndexX);
        cplx y = r.at(kIndexY);
        spdlog::debug("Numerical Check: H(b, y) = {}\n.", 
            complex_to_string(curve_->eval_H(b, y)));
        
        state_type start_state;
        start_state.at(kIndexX) = b;
        start_state.at(kIndexY1) = std::log(y);
        start_state.at(kIndexY2) = std::log(y);
        
        // Compute the expansion around a branch point.
        // TODO: I should move this to the SW_curve class and compute it generally (kappa = 1 / n! d^n H / dy^n * (1 / d H / dx))
        cplx kappa = -1.0 / 2 * curve_->eval_d2H_dy2(b, y) * curve_->eval_dH_dx(b,y);
        cplx dx = std::pow(3.0 / 4 * std::pow(kappa, 1.0/2) * y * b * std::exp(J*theta_) * kInitialStepSize, 2.0 / 3);  // ! This is only valid for exponential networks.

        cplx dy = 3.0 / 4 * (y * b * std::exp(J*theta_) * kInitialStepSize) / dx;
        for(uint32_t k = 0; k < 3; k++) {
            std::vector<state_type> v;
            std::vector<double> masses;
            v.push_back(start_state);
            masses.push_back(0);

            state_type next_state;
            
            spdlog::debug("(Zeta_3)^{} = {}.", k, complex_to_string(std::pow(kZeta3, k)));
            next_state.at(kIndexX) = b + std::pow(kZeta3, k) * dx;
            next_state.at(kIndexY1) = std::log(y + std::pow(kZeta3, k) * dy);
            next_state.at(kIndexY2) = std::log(y - std::pow(kZeta3, k) * dy);
            determine_sign(start_state, next_state);
            v.push_back(next_state);
            

            cplx dlog_y1 = next_state.at(kIndexY1) - start_state.at(kIndexY1);
            cplx dlog_y2 = next_state.at(kIndexY2) - start_state.at(kIndexY2);
            state_type dv{dx, dlog_y1,dlog_y2};
            masses.push_back(compute_dm(next_state, dv));
            
            
            for(uint32_t i = 0; i < kInitialEulerSteps; i++) {
                ODE_euler_step(curve_, v, masses);
            }
            curve_->match_fiber(std::prev(v.end()));
            Path path(v, masses, new_paths_.size());
            new_paths_.push_back(std::move(path));
            
            // new_paths_.back().print_data();
            std::vector<cplx> fiber = curve_->get_fiber(v.back().at(kIndexX));
            spdlog::debug("Fiber over x = {} is:", complex_to_string(v.back().at(kIndexX)));
            for (auto& f : fiber) {
                spdlog::debug("{}", complex_to_string(f));
            }
            spdlog::debug("Approximation: {} {}", complex_to_string(std::exp(v.back().at(kIndexY1))), complex_to_string(std::exp(v.back().at(kIndexY2))));
            
            spdlog::debug("Path {} appended.", new_paths_.back().id_);
            spdlog::debug("Numerical Check: \n H(x, y_1) = {}\n H(x, y_2) = {}.", 
            complex_to_string(std::abs(curve_->eval_H(v.back().at(kIndexX), std::exp(v.back().at(kIndexY1))))), 
            complex_to_string(std::abs(curve_->eval_H(v.back().at(kIndexX), std::exp(v.back().at(kIndexY2))))));   
        }
    }
}

auto Network::evolution_step() -> void {
    std::vector<future_type> futures;
    for(auto & new_path : new_paths_) {
        futures.emplace_back(new_path.integrate(curve_));
    }
    auto new_path_it = new_paths_.begin();
    for (auto &f : futures) {
        if (new_path_it == new_paths_.end()) {
            spdlog::debug("Iterator in future loop overshoots!");
        }
        new_path_it->update(std::ref(f));
        //! Only for Testing purposes!
        new_path_it->save_data();
        evolved_paths_.push_back(std::move(*new_path_it));
        new_path_it++;
        
    }
    new_paths_.clear();
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