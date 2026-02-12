#ifndef P2_HPP_
#define P2_HPP_

#include <ginac/ginac.h>

#include "nets.hpp"
#include "sw_curve.hpp"

auto H_p2(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;

class P2 : protected Network {
 public:
  explicit P2(double theta) : Network(H_p2, theta) {}

  auto custom_BPS(double cutoff) -> void;

 private:
  double theta_;
};

#endif // P2_HPP_
