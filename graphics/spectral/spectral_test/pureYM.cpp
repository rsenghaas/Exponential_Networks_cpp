#include "pureYM.hpp"

#include <spdlog/spdlog.h>

#include "ginac_util.hpp"
#include "type_util.hpp"


auto H_su(const GiNaC::symbol& x, const GiNaC::symbol& y, cplx u) -> GiNaC::ex { 
  return (y*y - x - 1 / x + 2*complex_to_ex(u)) * x;
}

auto PureYM::custom_BPS(double cutoff) -> void {
    spdlog::debug("In custom BPS");
    std::vector<size_t> inv = {}; //{2,8, 9,12,16,20,21};
    for (size_t &k : inv) {
        auto endpoint = new_paths_.at(k).get_endpoint();
        invert_state(endpoint);
        new_paths_.at(k).override_endpoint(endpoint);
    }
    initial_integration();
    auto path_it = get_iterator_by_id(new_paths_, 0);
    uint32_t k = 0;
    for ( ;k < new_paths_.size(); k += 1){
        path_it = get_iterator_by_id(new_paths_, k);
        auto endpoint = path_it->get_endpoint();
        if(std::abs(endpoint.at(kIndexX)) > 10) {
            // continue;
        }
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        path_it->truncate_mass(cutoff);
        path_it->truncate_x(50);
        save_data(path_it->path_id_);
    }
}
