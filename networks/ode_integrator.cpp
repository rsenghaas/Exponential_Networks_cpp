#include "ode_integrator.hpp"

auto compute_dm(const state_type &v, const state_type &dv) -> double {
  double delta_m = std::abs((v.at(kIndexY2) - v.at(kIndexY1)) * dv.at(kIndexX) /
                            v.at(kIndexX));
  return delta_m;
}

auto rotate_dv(state_type &dv, const double theta) -> void {
  dv.at(kIndexX) *= std::exp(J * theta);
  dv.at(kIndexY1) *= std::exp(J * theta);
  dv.at(kIndexY2) *= std::exp(J * theta);
}

auto ODE_euler_step(const std::shared_ptr<SW_curve> &curve,
                    std::vector<state_type> &v, std::vector<double> &masses,
                    const double step_size, double theta) -> void {
  state_type dv;
  curve->sw_differential(v.back(), dv);
  rotate_dv(dv, theta);
  v.push_back(v.back() + dv);
  masses.push_back(masses.back() + compute_dm(v.back(), step_size * dv));
}

auto ODE_elliptic_euler_step(const std::shared_ptr<SW_curve> &curve,
                             std::vector<state_type> &v,
                             std::vector<double> &masses,
                             const double step_size, double theta) -> void {
  state_type dv;
  curve->elliptic_differential(v.back(), dv);
  rotate_dv(dv, theta);
  v.push_back(v.back() + step_size * dv);
  masses.push_back(masses.back() + compute_dm(v.back(), step_size * dv));
}

auto ODE_runge_kutta_step(const std::shared_ptr<SW_curve> &curve,
                          std::vector<state_type> &v,
                          std::vector<double> &masses, const double step_size,
                          double theta) -> void {
  state_type dv;
  state_type u = v.back();
  state_type k1;
  curve->sw_differential(u, k1);
  u = v.back() + step_size / 2.0 * k1;
  state_type k2;
  curve->sw_differential(u, k2);
  u = v.back() + step_size / 2.0 * k2;
  state_type k3;
  curve->sw_differential(u, k3);
  u = v.back() + step_size * k3;
  state_type k4;
  curve->sw_differential(u, k4);
  dv = step_size / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
  rotate_dv(dv, theta);
  v.push_back(v.back() + dv);
  masses.push_back(masses.back() + compute_dm(v.back(), dv));
}

auto ODE_integrator::integrate_ode(double cutoff) -> void {
  spdlog::debug("In ode integrator.");
  auto stepper = make_controlled(kOdeAbsError, kOdeRelError,
                                 runge_kutta_cash_karp54<state_type>());
  uint32_t loop_counter = 0;
  while (obs_.masses.back() < cutoff && obs_.states.size() < kMaxSteps) {
    if ((std::abs(v0_.at(kIndexX)) > kHighCutoffX ||
         std::abs(v0_.at(kIndexY1)) > kCutoffY ||
         std::abs(v0_.at(kIndexY2)) > kCutoffY ||
         std::abs(v0_.at(kIndexX)) < kLowCutoffX ||
         std::abs(v0_.at(kIndexY1) - v0_.at(kIndexY2)) < 1 / kCutoffY) &&
        loop_counter > 3) {
      break;
    }
    std::cout << obs_.states.size() << std::endl;
    integrate_adaptive(stepper, diff_, v0_, t0_, t0_ + kIntegratePeriod,
                       kIntegraterStepSize, observer(obs_));
    // integrate_const(runge_kutta4<state_type>(), diff_, v0_, t0_, t0_ +
    // kInitialStepSize, kInitialStepSize, observer(obs_));
    diff_.curve_->match_fiber(v0_);
    loop_counter++;
  }
}

auto ODE_integrator::print_v() -> void {
  for (uint32_t i = 0; i < obs_.states.size(); i++) {
    std::cout << "Step: " << i << "x: " << obs_.states.at(i).at(kIndexX)
              << " y1: " << obs_.states.at(i).at(kIndexY1)
              << " y2: " << obs_.states.at(i).at(kIndexY2) << "with mass "
              << obs_.masses.at(i) << std::endl;
  }
}

auto ODE_integrator::get_times() -> std::vector<double> { return obs_.times; }

auto ODE_integrator::get_v() -> std::vector<state_type> { return obs_.states; }

auto ODE_integrator::get_masses() -> std::vector<double> { return obs_.masses; }

void ODE_differential::operator()(const state_type &v, state_type &dv,
                                  const double time) {
  std::lock_guard<std::mutex> guard(curve_->sw_mutex);
  if (curve_->mode == "elliptic") {
    curve_->elliptic_differential(v, dv);
  } else {
    curve_->sw_differential(v, dv);
  }
  rotate_dv(dv, theta_);
}
