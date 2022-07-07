#ifndef PATH_HPP_
#define PATH_HPP_

#include <vector>
#include <complex>
#include <algorithm>
#include <memory>
#include <thread>
#include <future>
#include <spdlog/spdlog.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <fmt/core.h>

#include "magic_numbers.h"
#include "sw_curve.hpp"
#include "ode_integrator.hpp"
#include "format_util.hpp"

class Path {
  public:
    Path(std::vector<state_type> &v0, std::vector<double> &masses, uint32_t id) : id_(id), v_(), masses_()
    {
        v_.push_back(v0.front());
        masses_.push_back(masses.front());
        append_data(v0, masses);
    }  
    auto get_endpoint() -> state_type;
    auto update(future_type& future) -> void;
    auto integrate(std::shared_ptr<SW_curve> pCurve) -> future_type;
    

    auto print_data() -> void;
    auto save_data() -> void;
    uint32_t id_;
  private:
    auto append_data(std::vector<state_type>& v, std::vector<double>& masses) -> void;
    std::vector<state_type> v_; 
    std::vector<double> masses_;  
};

#endif  // PATH_HPP_