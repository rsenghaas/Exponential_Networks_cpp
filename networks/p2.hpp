#ifndef P2_HPP_
#define P2_HPP_

#include <ginac/ginac.h>

#include "nets.hpp"
#include "sw_curve.hpp"

auto H_p2(const GiNaC::symbol &x, const GiNaC::symbol &y, cplx Q) -> GiNaC::ex;
auto H_1_1_2(const GiNaC::symbol& x, const GiNaC::symbol& y, cplx Q1, cplx Q2) -> GiNaC::ex;
auto custom_curve(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex;

class P2 : protected Network {
 public:
  explicit P2(double theta, cplx Q1, cplx Q2) : 
      Network([Q1,Q2](const GiNaC::symbol &x, 
              const GiNaC::symbol &y) {
              return custom_curve(x,y);
            // return H_1_1_2(x,y,Q1,Q2);
        }, theta) {}

  auto custom_BPS(double cutoff) -> void;

 private:
  double theta_;
};

#endif // P2_HPP_
