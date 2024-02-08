#ifndef CONI_HPP_
#define CONI_HPP_

#include <ginac/ginac.h>

#include "nets.hpp"
#include "sw_curve.hpp"

// TODO: Introduce dedicated type for this kind of function.
auto H_local(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;

class Local: protected Network {
 public:
  explicit Local(double theta) : Network(H_local, theta) {}

  auto custom_BPS(double cutoff) -> void;

 private:
  double theta_;
  auto custom_BPS_local(double cutoff) -> void;
};

#endif  // CONI_HPP_
