#include "type_util.hpp"

auto print_intersection(const intersection &intersect) -> void {
  spdlog::debug("Ids are {} and {}", intersect.ids.at(kIndexFirstPath),
                intersect.ids.at(kIndexSecondPath));
  spdlog::debug("Coordinate range is from ({},{}) to ({},{})",
                intersect.coordinates.front().at(kIndexCoordReal),
                intersect.coordinates.front().at(kIndexCoordImag),
                intersect.coordinates.back().at(kIndexCoordReal),
                intersect.coordinates.back().at(kIndexCoordImag));
  spdlog::debug("Time range for path {} is [{},{}]. ",
                intersect.ids.at(kIndexFirstPath),
                intersect.times.at(kIndexFirstPath).at(kIndexStartTime),
                intersect.times.at(kIndexFirstPath).at(kIndexEndTime));
  spdlog::debug("Time range for path {} is [{},{}].",
                intersect.ids.at(kIndexSecondPath),
                intersect.times.at(kIndexSecondPath).at(kIndexStartTime),
                intersect.times.at(kIndexSecondPath).at(kIndexEndTime));
}

auto print_path_point(const path_point &pp) -> void {
  spdlog::debug("Id is {}.", pp.id);
  spdlog::debug("Coordinate is ({},{}).", pp.coordinate_real,
                pp.coordinate_imag);
  spdlog::debug("Times range is [{},{}]", pp.t.at(kIndexStartTime),
                pp.t.at(kIndexEndTime));
}

auto print_state_type(const state_type &v) -> void {
  spdlog::debug("x = {}, y1 = {}, y2 = {}.", complex_to_string(v.at(kIndexX)),
                complex_to_string(v.at(kIndexY1)),
                complex_to_string(v.at(kIndexY2)));
}

auto complex_to_string(const std::complex<double> &z) -> std::string {
  std::string sign;
  sign = (z.imag() >= 0) ? "+" : "-";
  return fmt::format("{}{}{}j", z.real(), sign, std::abs(z.imag()));
}
