#ifndef CONI_HPP_
#define CONI_HPP_

#include <ginac/ginac.h>

#include "spectral_nets.hpp"
#include "sw_curve.hpp"

// TODO: Introduce dedicated type for this kind of function.
auto H_su(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;


class PureYM : protected Spectral_Network {
 public:
  explicit PureYM(double theta) : Network(H_su, theta) {}

  auto custom_BPS(double cutoff) -> void;

 private:
  double theta_;
};

#endif  // CONI_HPP_
