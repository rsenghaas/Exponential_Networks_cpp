#include "type_util.hpp"

#include <string>

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
  return std::to_string(z.real()) + sign + std::to_string(std::abs(z.imag())) +
         "j";
}

auto get_log_sheet(state_type &v) -> int32_t {
  cplx dv_y1 = v.at(kIndexY1) - std::log(std::exp(v.at(kIndexY1)));
  cplx dv_y2 = v.at(kIndexY2) - std::log(std::exp(v.at(kIndexY2)));
  int32_t k1 = static_cast<int32_t>(std::round((dv_y1 / (2 * pi * J)).real()));
  int32_t k2 = static_cast<int32_t>(std::round((dv_y2 / (2 * pi * J)).real()));
  return k1 - k2;
}

auto invert_state(state_type& v) -> void {
  cplx const y1_coord = v.at(kIndexY1);
  v.at(kIndexY1) = v.at(kIndexY2);
  v.at(kIndexY2) = y1_coord;
}

auto state_is_trivial(const state_type &v) -> bool {
    spdlog::debug("{}", std::abs(v.at(kIndexY1) - v.at(kIndexY2)));
    return (std::abs(v.at(kIndexY1) - v.at(kIndexY2)) < 1e-15);
}
