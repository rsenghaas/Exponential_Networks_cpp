#ifndef CONI_HPP_
#define CONI_HPP_ 

#include <ginac/ginac.h>

#include "nets.hpp"
#include "sw_curve.hpp"

auto H_coni(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;

class Coni : protected Network {
    public:
        explicit Coni(double theta) : Network(H_coni, theta) {}

        auto custom_BPS() -> void;

    private:
        double theta_;
};

#endif  // CONI_HPP_