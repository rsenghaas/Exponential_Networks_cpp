#ifndef ODE_INTEGRATOR_HPP_
#define ODE_INTEGRATOR_HPP_

#include <spdlog/spdlog.h>

#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <vector>

#include "magic_numbers.h"
#include "sw_curve.hpp"

using namespace boost::numeric::odeint;

auto compute_dm(const state_type &v, const state_type &dv) -> double;

auto ODE_euler_step(const std::shared_ptr<SW_curve> curve,
                    std::vector<state_type> &v, std::vector<double> &masses,
                    const double step_size, double theta) -> void;
auto ODE_runge_kutta_step(const std::shared_ptr<SW_curve> curve,
                          std::vector<state_type> &v,
                          std::vector<double> &masses, const double step_size,
                          double theta) -> void;

struct observable {
  std::vector<double> times;
  std::vector<double> masses;
  std::vector<state_type> states;
};

class ODE_differential {
 public:
  ODE_differential(std::shared_ptr<SW_curve> curve, double theta)
      : theta_(theta), curve_(curve){};

  void operator()(const state_type &v, state_type &dv, const double time);

 private:
  double theta_;
  std::shared_ptr<SW_curve> curve_;
};

class ODE_integrator {
 public:
  ODE_integrator(ODE_differential diff, state_type v0)
      : diff_(diff), v0_(v0), t0_(0), obs_() {
    obs_.times.push_back(t0_);
    obs_.masses.push_back(0);
    obs_.states.push_back(v0);
  }

  auto integrate_ode(double cutoff) -> void;
  auto print_v() -> void;
  auto get_v() -> std::vector<state_type>;
  auto get_masses() -> std::vector<double>;
  auto get_times() -> std::vector<double>;

 private:
  ODE_differential diff_;
  state_type v0_;
  double t0_;
  observable obs_;
};

struct observer {
  observable &m_obs;

  observer(observable &obs) : m_obs(obs) {}

  void operator()(const state_type &v, const double t) {
    state_type dv;
    std::transform(std::begin(v), std::end(v), std::begin(m_obs.states.back()),
                   std::begin(dv),
                   [](std::complex<double> v1, std::complex<double> v2) {
                     return v1 - v2;
                   });
    m_obs.masses.push_back(m_obs.masses.back() + compute_dm(v, dv));
    m_obs.states.push_back(v);
    m_obs.times.push_back(t);
  }
};

#endif  // ODE_INTEGRATOR_HPP_
