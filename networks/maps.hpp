#ifndef MAPS_HPP_
#define MAPS_HPP_

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <map>
#include <vector>

#include "magic_numbers.h"
#include "path.hpp"
#include "spdlog/spdlog.h"
#include "type_util.hpp"

auto valid_coordinates(path_point pp) -> bool;
auto path_point_to_complex(path_point pp) -> cplx;

struct pixel {
  std::vector<path_point> pp_vec;
};

class Map {
 public:
  auto get_intersections(const std::vector<path_point>& pp_vec)
      -> std::vector<intersection>;
  auto draw_line(path_point pp1, path_point pp2)
      -> std::vector<std::vector<path_point>>;
  auto print_map_data() -> void;
  auto get_pixel_content(std::array<int32_t, 2> coordinates)
      -> std::vector<path_point>;
  auto add_path(std::vector<path_point> pp_vec) -> void;

 private:
  // auto handle_intersection(const path_point& pp, const
  // std::vector<path_point>& map_pp_vec, std::vector<intersection>&
  // intersections) -> void;
  std::map<std::array<int32_t, 2>, std::vector<path_point>> map_data_;
};

class SinglePathMap {
 public:
  explicit SinglePathMap(std::vector<Path>::iterator path_it)
      : path_id(path_it->path_id_) {
    path_it->compute_map_points(pp_vec, index_vec);
  }

  auto get_self_intersections() -> void;
  uint32_t path_id;
  std::vector<path_point> pp_vec;
  std::vector<uint32_t> index_vec;
  std::vector<intersection> intersections;

 private:
  std::map<std::array<int32_t, 2>, std::vector<path_point>> map_data_;
};

#endif  // MAPS_HPP_
