#ifndef PATH_HPP_
#define PATH_HPP_

#include <fmt/core.h>

#include <algorithm>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include "magic_numbers.h"
#include "ode_integrator.hpp"
#include "spdlog/spdlog.h"
#include "sw_curve.hpp"
#include "type_util.hpp"

class Path {
 public:
  Path(std::vector<state_type>& v0, std::vector<double>& masses,
       uint32_t path_id)
      : path_id_(path_id) {
    v_.push_back(v0.front());
    masses_.push_back(masses.front());
    append_data(v0, masses);
  }
  auto get_endpoint() -> state_type;
  auto update(future_type& future) -> void;
  auto add_single_point(state_type pt) -> void;
  auto integrate(std::shared_ptr<SW_curve> pCurve, double theta, double cutoff)
      -> future_type;

  auto compute_map_points(std::vector<path_point>& pp_vec,
                          std::vector<uint32_t>& index_vec) -> void;

  auto get_point(uint32_t t) -> state_type;
  auto get_length() -> uint32_t {return v_.size();}
  auto print_data() -> void;
  auto save_data() -> void;
  auto get_endmass() -> double {return masses_.back();}
  auto truncate(uint32_t t_start, uint32_t t_end) -> void;

  uint32_t path_id_;

 private:
  auto append_data(std::vector<state_type>& v, std::vector<double>& masses)
      -> void;
  std::vector<state_type> v_;
  std::vector<double> masses_;
};

#endif  // PATH_HPP_
