#include "coni.hpp"

#include <spdlog/spdlog.h>

#include "ginac_util.hpp"
#include "type_util.hpp"

constexpr cplx kQ_coni =  11.0 / 8.0 - 3.0 / 4.0 * J;

auto H_coni(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return -1 + y + x * y - complex_to_ex(kQ_coni) * GiNaC::pow(y, 2) * x;
}

auto F_2_3(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
    int f = -1;
    GiNaC::ex x_sub = x * GiNaC::pow(-y, f);
    /* 
     GiNaC::ex aug = GiNaC::pow(x_sub, 2) * (y - 1)
     + x_sub * (
         GiNaC::pow(y, 4)
         - GiNaC::pow(y, 3) * complex_to_ex(kQ_coni)
         + 2 * GiNaC::pow(y, 2) * GiNaC::pow(complex_to_ex(kQ_coni), 2)
         - 2 * GiNaC::pow(y, 2) * complex_to_ex(kQ_coni)
         - y * GiNaC::pow(complex_to_ex(kQ_coni), 2)
         + GiNaC::pow(complex_to_ex(kQ_coni), 2)
     )
     + GiNaC::pow(y, 3) * GiNaC::pow(complex_to_ex(kQ_coni), 4)
     - GiNaC::pow(y, 4) * GiNaC::pow(complex_to_ex(kQ_coni), 3);*/

    GiNaC::ex aug = GiNaC::pow(complex_to_ex(kQ_coni), 6) * x_sub * GiNaC::pow(y, 4)
     + GiNaC::pow(complex_to_ex(kQ_coni), 2) * GiNaC::pow(y, 3) * 
     (1 - GiNaC::pow(complex_to_ex(kQ_coni), 2) * x_sub)
     + y * (1 - GiNaC::pow(complex_to_ex(kQ_coni), 2) * x_sub)
     + GiNaC::pow(complex_to_ex(kQ_coni), 2) * x_sub
     + GiNaC::pow(y, 2) * (
         GiNaC::pow(complex_to_ex(kQ_coni), 4) * (-GiNaC::pow(x_sub, 2))
         - 2 * GiNaC::pow(complex_to_ex(kQ_coni), 4) * x_sub
         + 2 * GiNaC::pow(complex_to_ex(kQ_coni), 2) * x_sub
         - GiNaC::pow(complex_to_ex(kQ_coni), 2)
         - 1); 

    /* return (-1 + complex_to_ex(kQ_coni) * y) +
         (GiNaC::pow(y, 3) - GiNaC::pow(y, 4) + 2 * GiNaC::pow(y, 5) -
          2 * complex_to_ex(kQ_coni) * GiNaC::pow(y, 5) -
          complex_to_ex(kQ_coni) * GiNaC::pow(y, 6) +
          complex_to_ex(kQ_coni * kQ_coni) * GiNaC::pow(y, 7)) *
             x +
         (-GiNaC::pow(y, 9) + GiNaC::pow(y, 10)) * GiNaC::pow(x, 2); */
    aug = aug * GiNaC::pow(y, -f);
    std::stringstream F;
    F << aug.expand().collect(y);
    spdlog::debug(F.str());
    return aug;
}

auto H_trefoil(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  return GiNaC::pow(y, 2) * GiNaC::pow(y - 1, 3) -
         x * (y - complex_to_ex(kQ_coni)) * (y - complex_to_ex(kQ_coni)) *
             (y - complex_to_ex(kQ_coni));
}

auto F_fig_8(const GiNaC::symbol& y, const GiNaC::symbol& x) -> GiNaC::ex {
  return (GiNaC::pow(y, 2) - complex_to_ex(kQ_coni) * GiNaC::pow(y, 3)) +
         (-1 + 2 * y - 2 * complex_to_ex(kQ_coni * kQ_coni) * GiNaC::pow(y, 4) +
          complex_to_ex(kQ_coni * kQ_coni) * GiNaC::pow(y, 5)) *
             x +
         (1 - 2 * complex_to_ex(kQ_coni) * y +
          2 * complex_to_ex(kQ_coni * kQ_coni) * GiNaC::pow(y, 4) -
          complex_to_ex(kQ_coni * kQ_coni * kQ_coni) * GiNaC::pow(y, 5)) *
             GiNaC::pow(x, 2) +
         (-complex_to_ex(kQ_coni * kQ_coni) * GiNaC::pow(y, 2) +
          complex_to_ex(kQ_coni * kQ_coni) * GiNaC::pow(y, 3)) *
             GiNaC::pow(x, 3);
}

const cplx kZeta5 = std::exp(2.0 * pi * J / 5.0);
auto F_test(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
  double N = 6;
  double epsilon = 1 / (64.0);
  const cplx zetaN = std::exp(2.0 * pi * J / (N + 13.5));
  GiNaC::ex x_pol = 1;
  for (uint32_t i = 0; i < N; i++) {
    x_pol *= (x + 1 -
              epsilon * (1 + std::pow(0.2, i)) *
                  GiNaC::pow(complex_to_ex(zetaN), i));
  }
  return (y - 4.0) * (y - 4.0) + x_pol;
}

std::vector<double> ratio_zero = {1, 4, 2, 1, 2,  8,  1, 32, 1, 1,
                                  2, 2, 1, 4, 64, 16, 1, 1,  1, 2,
                                  2, 1, 4, 4, 16, 4,  1, 1,  1, 32};

std::vector<double> ratio_pi_half = {2, 1, 1, 1, 2,  2, 16, 1,  1, 1,
                                     2, 2, 1, 4, 1,  1, 1,  16, 1, 2,
                                     2, 1, 4, 4, 16, 1, 1,  1,  1, 1};

auto Coni::intersect_and_integrate(uint32_t k1, uint32_t k2, double cutoff) -> bool {
    uint32_t k =  new_paths_.size();
    two_path_intersection_handler(k1, k2, false, false, 0, 0, false, false); 
    if (new_paths_.size() > k) {
        auto path_it = get_iterator_by_id(new_paths_, k);
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        save_data(path_it->path_id_);
        return true;
    }
    return false;
}

auto Coni::custom_BPS_trifoil(double cutoff) -> void {

    auto path_it = get_iterator_by_id(new_paths_, 0);
    for (uint32_t k = 0; k < new_paths_.size(); k += 1){
        auto path_it = get_iterator_by_id(new_paths_, k);
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        save_data(path_it->path_id_);
    }
    uint32_t k = new_paths_.size();
    two_path_intersection_handler(18, 21, true, false, 0, 0, false, false); 
    path_it = get_iterator_by_id(new_paths_, k);
    save_data(path_it->path_id_);
    evolve_path(path_it, cutoff);
    save_data(path_it->path_id_);
    save_data(21);
    save_data(18);
    k++; 
    two_path_intersection_handler(5, 23, false, true, 0, 0, false, false); 
    save_data(23);
    path_it = get_iterator_by_id(new_paths_, k);
    k++;

    two_path_intersection_handler(9, 21, true, false, 0, 0, false, false); 
    path_it = get_iterator_by_id(new_paths_, k);
    save_data(path_it->path_id_);
    evolve_path(path_it, cutoff);
    save_data(path_it->path_id_);
    save_data(25);
    save_data(9);
    k++;

    two_path_intersection_handler(10, 25, true, false, 0, 0, false, false); 
    save_data(10);
    two_path_intersection_handler(5, 25, false, true, 0, 0, false, false); 

    save_data(25); 
    path_it = get_iterator_by_id(new_paths_, 5);
    path_it->truncate(0, 500);
    two_path_intersection_handler(5, 19, true, true, 0, 0, false, false); 
    save_data(5);
    save_data(19);
    path_it = get_iterator_by_id(new_paths_, 3);
    path_it->truncate(0, 1350);
    save_data(3);
    path_it = get_iterator_by_id(new_paths_, 21);
    path_it->truncate(0, 1350);
    save_data(21);
    path_it = get_iterator_by_id(new_paths_, 5);






    uint32_t generation_start = 0;
    uint32_t rounds = 2;
    std::vector<uint32_t> v = {0,1,2,3,7,8,12,13,14,15, 16, 17, 18, 19};
    for (uint32_t l = 0; l < rounds; l++) {
        break;
        uint32_t generation_end = k;
        for(uint32_t k1 = 0; k1 < generation_end; k1++) {      
            for(uint32_t k2 = std::max(generation_start,k1 + 1); k2 < generation_end; k2++) 
            { 
                if(std::find(v.begin(), v.end(), k1) != v.end() ||
                    std::find(v.begin(), v.end(), k2) != v.end()) {
                    continue;
                }
                std::cout << k1 << k2 << std::endl;
                two_path_intersection_handler(k1, k2, false, false, 0, 0, false, false); 
                if (new_paths_.size() > k) {
                    path_it = get_iterator_by_id(new_paths_, k);
                    print_state_type(path_it->get_endpoint());
                    save_data(path_it->path_id_);
                    evolve_path(path_it, cutoff);
                    save_data(path_it->path_id_);
                    k++;
                }
            }
        }
        generation_start = generation_end;
    }

    /* two_path_intersection_handler(17, 23, false, true, 0, 0, false, false); 
    save_data(23);

    two_path_intersection_handler(17, 20, true, true, 0, 0, false, false); 
    save_data(17);
    save_data(20);*/

    // Central charge is close to 4 pi^2!!

    /* path_it = get_iterator_by_id(new_paths_, 17);
    auto state_A = path_it->get_endpoint();
    path_it = get_iterator_by_id(new_paths_, 24);
    state_A.at(kIndexX) = path_it->get_endpoint().at(kIndexX);
    curve_->match_fiber(state_A);
    auto state_B = path_it->get_endpoint();
    state_type next_state;
    intersect_states(state_A, state_B, next_state);
    print_state_type(next_state);
    path_it = get_iterator_by_id(new_paths_, 25);
    new_paths_.erase(path_it);
    path_it = get_iterator_by_id(new_paths_, 24);
    path_it->override_endpoint(next_state);
    evolve_path(path_it, cutoff);
    save_data(24); */
    
};

auto Coni::custom_BPS_fig_8(double cutoff) -> void {
    auto path_it = get_iterator_by_id(new_paths_, 1);
    evolve_path(path_it, cutoff);
    save_data(path_it->path_id_);
}

auto Coni::custom_BPS_F() -> void {
  spdlog::debug("Drawing state.");
  /*for(auto path_it = new_paths_.begin(); path_it != new_paths_.end();
  path_it++) { spdlog::debug("Path {}.", path_it-> path_id_); try {

    }
    catch(...) {
      continue;
    }
    save_data(path_it->path_id_);
  }*/
  cplx x_zero = 1e-30;
  auto vvec = curve_->get_fiber(x_zero);
  for(auto& v : vvec) {
    std::cout << v << std::endl;
  }
  cplx x_infty = 1e50 ;
  vvec = curve_->get_fiber(x_infty); 
  for(auto& v : vvec) {
    std::cout << v << std::endl;
  }

}

auto Coni::custom_BPS(double cutoff) -> void { custom_BPS_trifoil(kCutoff);}
