#include "ode_integrator.hpp"

auto compute_dm(const state_type &v, const state_type &dv) -> double {
    double dm = std::abs((v.at(kIndexY2) - v.at(kIndexY1)) *dv.at(kIndexX) / v.at(kIndexX));
    return dm;
}

auto ODE_euler_step(const std::shared_ptr<SW_curve> curve, std::vector<state_type> &v, std::vector<double> &masses) -> void {
    state_type dv;
    curve->sw_differential(v.back(), dv);
    v.push_back(v.back() + kInitialStepSize * dv);
    masses.push_back(masses.back() + compute_dm(v.back(), dv));
}

auto ODE_integrator::integrate_ode() -> void {
    spdlog::debug("In ode integrator.");
    auto stepper = make_controlled(kOdeAbsError ,kOdeRelError ,runge_kutta_cash_karp54<state_type>());
    while (obs_.masses.back() < kCutoff) {
        integrate_adaptive(stepper, diff_, v0_, t0_, t0_ + kIntegratePeriod, kInitialStepSize, observer(obs_));
    }
}

auto ODE_integrator::print_v() -> void {
    for (uint32_t i = 0; i < obs_.states.size(); i++)
    {
        std::cout << "Step: " << i << "x: " << obs_.states.at(i).at(kIndexX) << " y1: " 
                << obs_.states.at(i).at(kIndexY1) << " y2: " << obs_.states.at(i).at(kIndexY2) << "with mass " << obs_.masses.at(i) << std::endl;
    }
}

auto ODE_integrator::get_times() -> std::vector<double> {
    return obs_.times;
}

auto ODE_integrator::get_v() -> std::vector<state_type> {
    return obs_.states;
}

auto ODE_integrator::get_masses() -> std::vector<double> {
    return obs_.masses;
}


void ODE_differential::operator() (const state_type &v, state_type &dv, const double) {  
        std::lock_guard<std::mutex> guard(curve_->sw_mutex);
        curve_->sw_differential(v, dv);
}