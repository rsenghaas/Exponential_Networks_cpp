#include "maps.hpp"

// Note that we don't draw the first point, since this is the last point
// of the previous line and we don't want doubling.
//
/* auto Map::draw_line(path_point pp1, path_point pp2)
    -> std::vector<std::vector<path_point>> {
  std::vector<std::vector<path_point>> intersections = {};
  path_point next_pp;
  next_pp.id = pp1.id;
  next_pp.pp_vec_index = pp2.pp_vec_index;
  next_pp.coordinate_real = pp1.coordinate_real;
  next_pp.coordinate_imag = pp1.coordinate_imag;
  next_pp.t = {pp1.t.at(1), pp2.t.at(0)};
  int32_t start_coord_real = pp1.coordinate_real;
  int32_t start_coord_imag = pp1.coordinate_imag;
  int32_t end_coord_real = pp2.coordinate_real;
  int32_t end_coord_imag = pp2.coordinate_imag;

  // spdlog::debug("Enter connecting Pixel ({},{}) and ({}, {}).",
  // pp1.coordinate_real, pp1.coordinate_imag, pp2.coordinate_real,
  // pp2.coordinate_imag);
  double ghost_coord_real{0};
  double ghost_coord_imag{0};

  auto distance = static_cast<int32_t>(
      kLineStepsPerUnit *
      (std::pow(std::pow(end_coord_real - start_coord_real, 2) +
                    std::pow(end_coord_imag - start_coord_imag, 2),
                1.0 / 2)));

  for (int32_t counter = 0; counter <= distance; counter++) {
    ghost_coord_real = (1 - counter * 1.0 / distance) * (start_coord_real) +
                       counter * 1.0 / distance * (end_coord_real) +
                       1.0 /2;  // Draw line from middle of the pixel.
    ghost_coord_imag = (1 - counter * 1.0 / distance) * (start_coord_imag) +
                       counter * 1.0 / distance * (end_coord_imag) + 1.0 / 2;  //
    if (static_cast<int32_t>(ghost_coord_real + kNumOffset) !=
            next_pp.coordinate_real ||
        static_cast<int32_t>(ghost_coord_imag + kNumOffset) !=
            next_pp.coordinate_imag) {
      next_pp.coordinate_real =
          static_cast<int32_t>(ghost_coord_real + kNumOffset);
      next_pp.coordinate_imag =
          static_cast<int32_t>(ghost_coord_imag + kNumOffset);
      if (valid_coordinates(next_pp)) {
        auto map_pt = map_data_.find(std::array<int32_t, 2>{
            next_pp.coordinate_real, next_pp.coordinate_imag});
        if (map_pt != map_data_.end()) {
          intersections.push_back(map_pt->second);
        } else {
          map_data_.insert({std::array<int32_t, 2>{next_pp.coordinate_real,
                                                   next_pp.coordinate_imag},
                            {}});
          map_pt = map_data_.find(std::array<int32_t, 2>{
              next_pp.coordinate_real, next_pp.coordinate_imag});
        }
        map_pt->second.push_back(next_pp);
      }
      // spdlog::debug("({},{}) with ghost coordinates ({:.2f}, {:.2f})",
      // next_pp.coordinate_real, next_pp.coordinate_imag, ghost_coord_real,
      // ghost_coord_imag);
    }
  }
  // spdlog::debug("Exit Connecting pixel");
  return intersections;
} */

auto Map::get_pixel_content(std::array<int32_t, 2> coordinates)
    -> std::vector<path_point> {
  return map_data_[coordinates];
}

//WARN: This is defined twice.
// We will probably remove it from nets.cpp).
auto neighbour_pixel(std::array<int32_t, 2> coord_arr1,
                     std::array<int32_t, 2> coord_arr2) -> bool {
  return !(std::max(std::abs(coord_arr1.at(kIndexCoordReal) -
                             coord_arr2.at(kIndexCoordReal)),
                    std::abs(coord_arr1.at(kIndexCoordImag) -
                             coord_arr2.at(kIndexCoordImag))) > 1);
}

//WARN:: This defined twice.
auto valid_coordinates(path_point pp) -> bool {
  if (pp.coordinate_real < 0 || pp.coordinate_real >= kMapResolutionReal) {
    return false;
  }
  if (pp.coordinate_real < 0 || pp.coordinate_imag >= kMapResolutionImag) {
    return false;
  }
  return true;
}

auto handle_intersection(const path_point& pp, const std::vector<path_point>& map_pp_vec, std::vector<intersection>& intersections) -> void {
  std::array<int32_t, 2> coord_arr{pp.coordinate_real, pp.coordinate_imag};
  for(const auto& map_pp : map_pp_vec) {
    bool map_pp_inserted{false};
    for (auto& inter : intersections) {
        //NOTE: Check Ids.
        if(map_pp.id != inter.ids.at(kIndexFirstPath) || pp.id != inter.ids.at(kIndexSecondPath)) {
            continue;
        }
        //NOTE: Check if pixel are next to each other.
        if(!neighbour_pixel(coord_arr, inter.coordinates.back())) {
            // continue;
        }
        //WARN: This is maybe not very clean, since it could happen that we are close to an intersection in another path
        if(pp.t.at(kIndexEndTime) - inter.times.at(kIndexSecondPath).at(kIndexEndTime) > 1) {
            continue;
        }
        inter.coordinates.push_back(coord_arr);
        inter.times.at(kIndexSecondPath).at(kIndexEndTime) = pp.t.at(kIndexEndTime);
        if (map_pp.t.at(kIndexEndTime) > inter.times.at(kIndexFirstPath).at(kIndexEndTime))
        {
            inter.times.at(kIndexFirstPath).at(kIndexEndTime) = map_pp.t.at(kIndexEndTime);
        }
        if (map_pp.t.at(kIndexStartTime) < inter.times.at(kIndexFirstPath).at(kIndexStartTime)) {
            inter.times.at(kIndexFirstPath).at(kIndexStartTime) = map_pp.t.at(kIndexStartTime);
        }
        map_pp_inserted = true;
        break;
    }

    if(!map_pp_inserted) {
        intersection new_intersection;
        new_intersection.ids ={map_pp.id, pp.id};
        new_intersection.coordinates.push_back(coord_arr);
        new_intersection.times = {map_pp.t,pp.t};
        intersections.push_back(new_intersection);
    }
  }
}

auto Map::get_intersections(const std::vector<path_point>& pp_vec) -> std::vector<intersection> {
    std::vector<intersection> intersections;
    for (const auto& pp : pp_vec) {
        auto map_pt = map_data_.find(std::array<int32_t, 2>{pp.coordinate_real, pp.coordinate_imag});
        if (map_pt != map_data_.end()) {
            handle_intersection(pp, map_pt->second, intersections);
        }

    }
    return intersections;
}


/* auto Map::print_map_data() -> void {
  uint32_t counter = 0;
  for (auto px_it = map_data_.begin(); px_it != map_data_.end(); px_it++) {
    std::string line = fmt::format(
        "#{} At {} there are {} paths: ", counter,
        complex_to_string(path_point_to_complex(px_it->second.pp_vec.front())),
        px_it->second.pp_vec.size());
    for (auto pp_it = px_it->second.pp_vec.begin();
         pp_it != px_it->second.pp_vec.end(); pp_it++) {
      line.append(fmt::format(" {}", pp_it->id));
    }
    std::cout << line << std::endl;
  }
  counter++;
} */

auto path_point_to_complex(path_point pp) -> cplx {
  return (kMapRangeReal.at(kIndexLowerBound) +
          (pp.coordinate_real + 1.0 / 2) *
              (kMapRangeReal.at(kIndexUpperBound) -
               kMapRangeReal.at(kIndexLowerBound)) /
              kMapResolutionReal) +
         J * (kMapRangeImag.at(kIndexLowerBound) +
              (pp.coordinate_imag + 1.0 / 2) *
                  (kMapRangeImag.at(kIndexUpperBound) -
                   kMapRangeImag.at(kIndexLowerBound)) /
                  kMapResolutionImag);
}

auto SinglePathMap::get_self_intersections() -> void {
  spdlog::debug("In intersection check");
  for (auto& next_pp : pp_vec) {
    std::array<int32_t, 2> coord_array = {next_pp.coordinate_real, next_pp.coordinate_imag};
    // spdlog::debug("({},{})", coord_array.at(0), coord_array.at(1));
    auto map_pt = map_data_.find(coord_array);
    if (map_pt != map_data_.end()) {
        // spdlog::debug("Intersection added");
        handle_intersection(next_pp, map_pt->second, intersections);
    } else {
      map_data_.insert({coord_array,
                        {}});
      map_pt = map_data_.find(coord_array);
    }
        map_pt->second.push_back(next_pp);
  }
  spdlog::debug("Leaving intersection check.");
}

auto Map::add_path(std::vector<path_point> pp_vec) -> void {
    for (auto& next_pp : pp_vec) {
        std::array<int32_t, 2> coord_array = {next_pp.coordinate_real, next_pp.coordinate_imag};
    // spdlog::debug("({},{})", coord_array.at(0), coord_array.at(1));
    auto map_pt = map_data_.find(coord_array);
    if (map_pt == map_data_.end()) {
      map_data_.insert({coord_array, {}});
      map_pt = map_data_.find(coord_array);
    }
    map_pt->second.push_back(next_pp);
  }
}
