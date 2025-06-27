#ifndef CONI_HPP_
#define CONI_HPP_

#include <ginac/ginac.h>

#include "nets.hpp"
#include "sw_curve.hpp"

// TODO: Introduce dedicated type for this kind of function.
auto H_coni(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;
auto F_2_3(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;
auto H_trefoil(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;
auto F_fig_8(const GiNaC::symbol &y, const GiNaC::symbol &x) -> GiNaC::ex;
auto F_test(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;
auto F_2_5(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex;

class Coni : protected Network {
 public:
  explicit Coni(double theta) : Network(F_2_5, theta) {}

  auto custom_BPS(double cutoff) -> void;

 private:
  double theta_;
  auto custom_BPS_trifoil(double cutoff) -> void;
  auto custom_BPS_F() -> void;
  auto custom_BPS_fig_8(double cutoff) -> void;
  auto intersect_and_integrate(uint32_t k1, uint32_t k2, double cutoff) -> bool;
};

#endif  // CONI_HPP_
