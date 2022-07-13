#ifndef MAPS_HPP_
#define MAPS_HPP_

#include <spdlog/spdlog.h>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdint>
#include <iostream>
#include <fmt/core.h>

#include "magic_numbers.h"
#include "type_util.hpp"

auto valid_coordinates(path_point pp) -> bool;
auto path_point_to_complex(path_point pp) -> cplx;

struct pixel {
  std::vector<path_point> pp_vec;
};

class Map {
  public:
    Map() : map_data_() {}
    
    auto draw_line(path_point pp1, path_point pp2) -> std::vector<std::vector<path_point>>;
    auto print_map_data() -> void;
    auto get_pixel_content(std::array<int32_t, 2> coordinates) -> std::vector<path_point>;
  private:
    std::map<std::array<int32_t, 2>, pixel> map_data_;
};

#endif  // MAPS_HPP_
