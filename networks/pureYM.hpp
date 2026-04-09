#ifndef CONI_HPP_
#define CONI_HPP_

#include <ginac/ginac.h>

#include "spectral_nets.hpp"
#include "sw_curve.hpp"

// TODO: Introduce dedicated type for this kind of function.
auto H_su(const GiNaC::symbol &x, const GiNaC::symbol &y, cplx u) -> GiNaC::ex;


class PureYM : protected Spectral_Network {
 public:
  explicit PureYM(double theta, cplx u) : 
      Spectral_Network([u](const GiNaC::symbol &x, const GiNaC::symbol &y) 
              { return H_su(x,y, u); }, theta) {}

  auto custom_BPS(double cutoff) -> void;

 private:
  double theta_;
};

#endif  // CONI_HPP_
