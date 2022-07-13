#include "maps.hpp"


// Note that we don draw the first point, since this is the last point
// of the previous line and we don't want doubling.
auto Map::draw_line(path_point pp1, path_point pp2) -> std::vector<std::vector<path_point>> {
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

    // spdlog::debug("Enter connecting Pixel ({},{}) and ({}, {}).", pp1.coordinate_real, pp1.coordinate_imag, pp2.coordinate_real, pp2.coordinate_imag);
    double ghost_coord_real;
    double ghost_coord_imag;

    int32_t distance = static_cast<int32_t>(kLineStepsPerUnit*(std::pow(std::pow(end_coord_real - start_coord_real,2) + std::pow(end_coord_imag - start_coord_imag,2) , 1.0 / 2)));

    for (int32_t counter = 0; counter <= distance; counter++) {
        ghost_coord_real = (1 - counter * 1.0 / distance) * (start_coord_real) + counter * 1.0 / distance * (end_coord_real) + 0.5;  // Draw line from middle of the pixel.
        ghost_coord_imag = (1 - counter * 1.0 / distance) * (start_coord_imag) + counter * 1.0 / distance * (end_coord_imag) + 0.5;  //
        if (static_cast<int32_t>(ghost_coord_real + kNumOffset) != next_pp.coordinate_real || static_cast<int32_t>(ghost_coord_imag + kNumOffset) != next_pp.coordinate_imag) {
            next_pp.coordinate_real = static_cast<int32_t>(ghost_coord_real + kNumOffset);
            next_pp.coordinate_imag = static_cast<int32_t>(ghost_coord_imag + kNumOffset);
            if (valid_coordinates(next_pp)) {
                auto map_pt = map_data_.find(std::array<int32_t,2>{next_pp.coordinate_real, next_pp.coordinate_imag});
                if (map_pt != map_data_.end()) {
                    intersections.push_back(map_pt->second.pp_vec);
                }
                else {
                    map_data_.insert({std::array<int32_t,2>{next_pp.coordinate_real, next_pp.coordinate_imag}, {}});
                    map_pt = map_data_.find(std::array<int32_t,2>{next_pp.coordinate_real, next_pp.coordinate_imag});
                }
                map_pt->second.pp_vec.push_back(next_pp);
            }
            // spdlog::debug("({},{}) with ghost coordinates ({:.2f}, {:.2f})", next_pp.coordinate_real, next_pp.coordinate_imag, ghost_coord_real, ghost_coord_imag);
        }
    }
    // spdlog::debug("Exit Connecting pixel");
    return intersections;
}

auto Map::get_pixel_content(std::array<int32_t, 2> coordinates) -> std::vector<path_point> {
    return map_data_[coordinates].pp_vec;
}


auto Map::print_map_data() -> void {
    uint32_t counter = 0;
    for (auto px_it = map_data_.begin(); px_it != map_data_.end(); px_it++) {
        std::string line = fmt::format("#{} At {} there are {} paths: ", counter, complex_to_string(path_point_to_complex(px_it->second.pp_vec.front())), px_it->second.pp_vec.size());
        for (auto pp_it = px_it->second.pp_vec.begin(); pp_it != px_it->second.pp_vec.end(); pp_it++) {
            line.append(fmt::format(" {}", pp_it->id));
        }
        std::cout << line << std::endl;
    }
    counter++;
}

auto path_point_to_complex(path_point pp) -> cplx {
    return (kMapRangeReal.at(kIndexLowerBound) + (pp.coordinate_real + 0.5) * (kMapRangeReal.at(kIndexUpperBound) - kMapRangeReal.at(kIndexLowerBound)) / kMapResolutionReal)
           + J * (kMapRangeImag.at(kIndexLowerBound) + (pp.coordinate_imag + 0.5) * (kMapRangeImag.at(kIndexUpperBound) - kMapRangeImag.at(kIndexLowerBound)) / kMapResolutionImag);
}

auto valid_coordinates(path_point pp) -> bool {
    if (pp.coordinate_real < 0 || pp.coordinate_real >= kMapResolutionReal) {
        return false;
    }
    if (pp.coordinate_real < 0 || pp.coordinate_imag >= kMapResolutionImag) {
        return false;
    }
    return true;
}


