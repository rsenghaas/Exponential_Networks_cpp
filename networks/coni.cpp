#include "coni.hpp"

#include <spdlog/spdlog.h>

#include "type_util.hpp"
#include "ginac_util.hpp"

constexpr cplx kQ_coni = 5.0 + 0.0 * J;

auto H_coni(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return -1 + y + x*y - complex_to_ex(kQ_coni) * GiNaC::pow(y, 2) * x;
}

auto Coni::custom_BPS() -> void {
    spdlog::debug("Drawing state.");
    auto path_it = get_iterator_by_id(new_paths_, 3);
    evolve_path(path_it, kCutoff);
    save_data(3);
    path_it = get_iterator_by_id(new_paths_, 4);
    evolve_path(path_it, kCutoff);
    save_data(4);
    path_it = get_iterator_by_id(new_paths_, 5);
    evolve_path(path_it, kCutoff);
    save_data(5);
}
