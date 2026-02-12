#include "p2.hpp"

#include "ginac_util.hpp"
#include "type_util.hpp"

constexpr cplx kQ_coni =  11.0 / 8.0 - 3.0 / 4.0 * J;
// theta_1 = 0.31
// theta_2 = 0.047
// 

auto H_p2(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return (1 + y + x + complex_to_ex(kQ_coni) / ( x * y )) * x * y;
}

auto H_1_1_2(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  int f = 0;
  GiNaC::ex x_sub = x * GiNaC::pow(-y, f);
  return (1 + x_sub + complex_to_ex(kQ_coni) * (1 + y) * (1 + y) / ( GiNaC::pow(x_sub,2) * y )) * GiNaC::pow(x_sub,2) * y;
}

auto P2::custom_BPS(double cutoff) -> void {
    initial_integration();
    auto path_it = get_iterator_by_id(new_paths_, 0);
    for (uint32_t k = 0; k < new_paths_.size(); k += 1){
        auto path_it = get_iterator_by_id(new_paths_, k);
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        save_data(path_it->path_id_);
    }
}
