#include "elliptic.hpp"

#include <arb.h>
#include <acb.h>
#include <acb_elliptic.h>

#include "ode_integrator.hpp"

auto H_c3(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex {
  return -x + GiNaC::pow(y, 2) + y;
}


/* auto arb_to_double(arb_t x) -> double {
  std::string arb_str = arb_get_str(x, 50, 0);
  size_t pos = arb_str.find(' ');
  if(pos != std::string::npos) {
    arb_str = arb_str.substr(1, pos);
  }
  double res = std::stold(arb_str);
  return res;
}

auto acb_to_cplx(acb_t z) -> cplx {
  arb_t re, im;
  arb_init(re);
  arb_init(im);
  acb_get_real(re, z);
  acb_get_imag(im, z);
  double re_d = arb_to_double(re);
  double im_d = arb_to_double(im);
  cplx res(re_d, im_d);
  return res;
} */

auto Elliptic::add_new_path(state_type start_point) -> void {
  std::vector<state_type> v;
  v.push_back(start_point);

  std::vector<double> masses;
  masses.push_back(0);

  Path path(v, masses, next_id_);
  paths_.push_back(std::move(path));
  next_id_++;
}

auto Elliptic::get_iterator_by_id(uint32_t id)
    -> std::vector<Path>::iterator {
  auto ret_it = paths_.begin();
  while (ret_it != paths_.end()) {
    if (ret_it->path_id_ == id) {
      break;
    }
    ret_it++;
  }
  return ret_it;
}
auto Elliptic::evolve_path(uint32_t path_id, double cutoff) -> void {
  auto path_it = get_iterator_by_id(path_id);
  future_type future = path_it->integrate(curve_, theta_, cutoff);
  path_it->update(std::ref(future));
}

auto Elliptic::custom_BPS() -> void {
  for(uint32_t i = 0; i < 10; i++) {
    cplx x = -0.27 + 0.01*J;
    spdlog::debug("Computing roots...");
    auto roots = curve_->get_fiber(x);
    cplx y1, y2;
    acb_t r, inv_p_y1, inv_p_y2, tau; 
    acb_init(r);
    acb_init(tau);
    acb_init(inv_p_y1);
    acb_init(inv_p_y2);
    acb_onei(tau);

    acb_set_d_d(r,  roots.at(0).real(), x.imag());
    acb_elliptic_inv_p(inv_p_y1, r,tau, 50);

    acb_set_d_d(r,  roots.at(1).real(), x.imag());
    acb_elliptic_inv_p(inv_p_y2, r, tau, 50);

    y1 = acb_to_cplx(inv_p_y1);
    y2 = acb_to_cplx(inv_p_y2);

    acb_clear(r);
    acb_clear(inv_p_y1);
    acb_clear(inv_p_y2);
    acb_clear(tau);
    std::cout << y1 << std::endl;
    std::cout << y2 << std::endl;
    std::vector<state_type> v;
    v.push_back(state_type{x, y1, y2});
    std::vector<double> masses;
    masses.push_back(0);
    state_type dv;
    print_state_type(v.back());
    curve_->elliptic_differential(v.back(), dv);
    print_state_type(dv);
    for (uint32_t i = 0; i < 200*kInitialSteps; i++) {
      ODE_elliptic_euler_step(curve_, v, masses, 40*kInitialStepSize, theta_);
    }
    Path path(v, masses, next_id_);
    paths_.push_back(std::move(path));
    next_id_++;
    save_data(i);
  }
}

auto Elliptic::save_data(uint32_t id) -> void {
  auto path_it = get_iterator_by_id(id);
  path_it->save_data();
}

