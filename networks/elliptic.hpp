#ifndef ELLIPTIC_HPP_
#define ELLIPTIC_HPP_

#include <ginac/ginac.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <memory>
#include <vector>

#include "maps.hpp"
#include "path.hpp"
#include "sw_curve.hpp"
#include "type_util.hpp"
#include "arb_util.hpp"
#include "magic_numbers.h"

auto H_c3(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex;

class Elliptic {
  public:
    explicit Elliptic(double theta) : theta_(theta) {
      curve_ = std::make_shared<SW_curve>(H_c3, "elliptic");
    }
    auto custom_BPS() -> void;

  private:
    double theta_;
    uint32_t next_id_{0};
    
    // SW_curve data.
    std::shared_ptr<SW_curve> curve_;

    // Path stuff.
    std::vector<Path> paths_;
    auto add_new_path(state_type start_point) -> void;
    auto evolve_path(uint32_t path_id, double cutoff) -> void;
    auto get_iterator_by_id(uint32_t id) -> std::vector<Path>::iterator;

    // Saving.
    auto save_data(uint32_t id) -> void;

};

#endif  // ELLIPTIC_HPP_
