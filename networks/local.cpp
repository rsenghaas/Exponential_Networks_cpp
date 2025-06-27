#include "local.hpp"

#include <spdlog/spdlog.h>

#include "ginac_util.hpp"
#include "type_util.hpp"

const cplx kQ_local = (-10.1 - 0.01*J);

auto H_local(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return -1 * complex_to_ex(kQ_local) * GiNaC::pow(x, 3) - y + x * y + GiNaC::pow(y, 2);
}

auto H_local_alt(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return -complex_to_ex(kQ_local) + y*  GiNaC::pow(x, 2) + x * GiNaC::pow(y, 2) + 1;
}

auto H_local_bm(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return y + GiNaC::pow(y, 2) + y * x - GiNaC::pow(x, 3) * complex_to_ex(kQ_local);
}


auto Local::custom_BPS_local(double cutoff) -> void {
  for (auto path_it = new_paths_.begin(); path_it != new_paths_.end();path_it++) {
       spdlog::debug("Path {}.", path_it->path_id_);
    try {
        evolve_path(path_it, cutoff);

    } catch (...) {
      continue;
    }
    save_data(path_it->path_id_);
  }
}

auto Local::custom_BPS(double cutoff) -> void { custom_BPS_local(cutoff); }
