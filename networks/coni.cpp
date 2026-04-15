#include "coni.hpp"

#include <spdlog/spdlog.h>

#include "ginac_util.hpp"
#include "type_util.hpp"

const cplx kQ_coni = 11.0 / 8.0 * std::exp(2.0 * J * pi * 0.0);

auto H_coni(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex { 
  std::array<int32_t, 4> M = {1,1,0,1};
  auto x_sub = GiNaC::pow(x, M.at(0)) * GiNaC::pow(y, M.at(1));
  auto y_sub = GiNaC::pow(x, M.at(2)) * GiNaC::pow(y, M.at(3));
  int f = 0;
  x_sub = x_sub * GiNaC::pow(y_sub, f); 
  return 1 + x_sub + y_sub + complex_to_ex(kQ_coni) * x_sub * y_sub;
}

auto F_5_2_hyperbolic(const GiNaC::symbol& x, const GiNaC::symbol &y) -> GiNaC::ex {
    std::array<int32_t, 4> M = {1, 0,-1,1};
    auto x_sub = GiNaC::pow(x, M.at(0)) * GiNaC::pow(y, M.at(2));
    auto y_sub = GiNaC::pow(x, M.at(1)) * GiNaC::pow(y, M.at(3));
    auto Q = complex_to_ex(kQ_coni);

    auto poly =
    pow(Q,2)
    - 2*Q*x_sub
    + pow(x_sub,2)

    + pow(Q,3)*y_sub
    - 4*pow(Q,2)*x_sub*y_sub
    + 3*Q*pow(x_sub,2)*y_sub

    + 3*pow(Q,2)*x_sub*pow(y_sub,2)
    - 3*pow(Q,3)*x_sub*pow(y_sub,2)
    - 4*Q*pow(x_sub,2)*pow(y_sub,2)
    + 5*pow(Q,2)*pow(x_sub,2)*pow(y_sub,2)
    - Q*pow(x_sub,3)*pow(y_sub,2)

    + 2*pow(Q,3)*x_sub*pow(y_sub,3)
    - 3*pow(Q,2)*pow(x_sub,2)*pow(y_sub,3)
    + 3*pow(Q,3)*pow(x_sub,2)*pow(y_sub,3)
    - 2*pow(Q,2)*pow(x_sub,3)*pow(y_sub,3)

    - pow(Q,4)*x_sub*pow(y_sub,4)
    + 6*pow(Q,2)*pow(x_sub,2)*pow(y_sub,4)
    - 4*pow(Q,3)*pow(x_sub,2)*pow(y_sub,4)
    - pow(Q,3)*pow(x_sub,3)*pow(y_sub,4)

    - 2*pow(Q,4)*x_sub*pow(y_sub,5)
    - 3*pow(Q,3)*pow(x_sub,2)*pow(y_sub,5)
    + 3*pow(Q,4)*pow(x_sub,2)*pow(y_sub,5)
    + 2*pow(Q,3)*pow(x_sub,3)*pow(y_sub,5)

    - pow(Q,4)*x_sub*pow(y_sub,6)
    - 4*pow(Q,3)*pow(x_sub,2)*pow(y_sub,6)
    + 5*pow(Q,4)*pow(x_sub,2)*pow(y_sub,6)
    + 3*pow(Q,3)*pow(x_sub,3)*pow(y_sub,6)
    - 3*pow(Q,4)*pow(x_sub,3)*pow(y_sub,6)

    + 3*pow(Q,4)*pow(x_sub,2)*pow(y_sub,7)
    - 4*pow(Q,4)*pow(x_sub,3)*pow(y_sub,7)
    + pow(Q,4)*pow(x_sub,4)*pow(y_sub,7)

    + pow(Q,4)*pow(x_sub,2)*pow(y_sub,8)
    - 2*pow(Q,4)*pow(x_sub,3)*pow(y_sub,8)
    + pow(Q,4)*pow(x_sub,4)*pow(y_sub,8);
    return poly * y * y;
}

auto F_2_5(const GiNaC::symbol& x, const GiNaC::symbol &y) -> GiNaC::ex {
    std::array<int32_t, 4> M = {1, 0,-1,1};
    auto x_sub = GiNaC::pow(x, M.at(0)) * GiNaC::pow(y, M.at(2));
    auto y_sub = GiNaC::pow(x, M.at(1)) * GiNaC::pow(y, M.at(3));
    auto Q = complex_to_ex(kQ_coni);

    auto poly = 
    -0.99999* pow(Q,2)
    + x_sub
    - pow(Q,3)*y_sub
    + Q*x_sub*y_sub
    - 2*Q*x_sub*pow(y_sub,2)
    + 2*pow(Q,2)*x_sub*pow(y_sub,2)
    - 2*pow(Q,2)*x_sub*pow(y_sub,3)
    + 2*pow(Q,3)*x_sub*pow(y_sub,3)
    + pow(Q,2)*x_sub*pow(y_sub,4)
    - 4*pow(Q,3)*x_sub*pow(y_sub,4)
    + 3*pow(Q,4)*x_sub*pow(y_sub,4)
    + pow(Q,3)*x_sub*pow(y_sub,5)
    + pow(Q,4)*x_sub*pow(y_sub,5)
    - 2*pow(Q,3)*pow(x_sub,2)*pow(y_sub,5)
    + 2*pow(Q,4)*x_sub*pow(y_sub,6)
    - pow(Q,3)*pow(x_sub,2)*pow(y_sub,6)
    - pow(Q,4)*pow(x_sub,2)*pow(y_sub,6)
    - pow(Q,3)*pow(x_sub,2)*pow(y_sub,7)
    + 4*pow(Q,4)*pow(x_sub,2)*pow(y_sub,7)
    - 3*pow(Q,5)*pow(x_sub,2)*pow(y_sub,7)
    + 2*pow(Q,4)*pow(x_sub,2)*pow(y_sub,8)
    - 2*pow(Q,5)*pow(x_sub,2)*pow(y_sub,8)
    + 2*pow(Q,4)*pow(x_sub,2)*pow(y_sub,9)
    - 2*pow(Q,5)*pow(x_sub,2)*pow(y_sub,9)
    - pow(Q,5)*pow(x_sub,2)*pow(y_sub,10)
    + pow(Q,6)*pow(x_sub,3)*pow(y_sub,10)
    - pow(Q,5)*pow(x_sub,2)*pow(y_sub,11)
    + pow(Q,6)*pow(x_sub,3)*pow(y_sub,11);
    return poly * y;
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

auto A_super_3_1(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex 
{
    int f = -1;
    auto t = complex_to_ex(-1);
    auto a = complex_to_ex(kQ_coni);
    auto x_sub = y * GiNaC::pow(-x, f);
    auto y_sub = x;

    auto Poly = GiNaC::pow(a, 2) * GiNaC::pow(t, 4) * (x_sub - 1) 
        * GiNaC::pow(x_sub,3) 
    - a * (1 - GiNaC::pow(t,2) * x_sub
            + 2 * GiNaC::pow(t, 2) * (1 + a * t) * x_sub * x_sub 
            + a * GiNaC::pow(t,5) * GiNaC::pow(x_sub,3) 
            + GiNaC::pow(a,2) * GiNaC::pow(t,6) 
            * GiNaC::pow(x_sub,4) ) * y_sub
    + (1 + a * GiNaC::pow(t,3) * x_sub) * GiNaC::pow(y_sub,2);
    return Poly * GiNaC::pow(y_sub, 3);
}

auto F_trefoil_generic(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
    std::array<int32_t, 4> M = {1, 0,0,1};
    auto x_sub = GiNaC::pow(x, M.at(0)) * GiNaC::pow(y, M.at(2));
    auto y_sub = GiNaC::pow(x, M.at(1)) * GiNaC::pow(y, M.at(3));
    auto Q = complex_to_ex(kQ_coni);
    auto poly = 
       1 - x_sub - y_sub + x_sub * y_sub - 2 * x_sub * GiNaC::pow(y_sub,2) 
       + 1.99995 * Q * x_sub * GiNaC::pow(y_sub,2) 
       + Q * x_sub * GiNaC::pow(y_sub,3)
       - GiNaC::pow(x_sub,2) * GiNaC::pow(y_sub,3) 
       - Q * Q * x_sub * GiNaC::pow(y_sub,4) 
       + Q * GiNaC::pow(x_sub,2) * GiNaC::pow(y_sub,4);
    return (poly * y).expand();
} 

/*
auto F_2_5(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
    int f = -1;
    auto coni_ex = complex_to_ex(kQ_coni);
    GiNaC::ex x_sub = x * pow(-y, f); 
    GiNaC::ex poly = 
        -GiNaC::pow(coni_ex, 2) + x_sub 
        - GiNaC::pow(coni_ex, 3) * y 
        + coni_ex * x_sub * y
        - 2 * coni_ex * x_sub * GiNaC::pow(y, 2)
        + 2 * GiNaC::pow(coni_ex, 2) * x_sub * GiNaC::pow(y, 2)
        - 2 * GiNaC::pow(coni_ex, 2) * x_sub * GiNaC::pow(y, 3)
        + 2 * GiNaC::pow(coni_ex, 3) * x_sub * GiNaC::pow(y, 3)
        + GiNaC::pow(coni_ex, 2) * x_sub * GiNaC::pow(y, 4)
        - 4 * GiNaC::pow(coni_ex, 3) * x_sub * GiNaC::pow(y, 4)
        + 3 * GiNaC::pow(coni_ex, 4) * x_sub * GiNaC::pow(y, 4)
        + GiNaC::pow(coni_ex, 3) * x_sub * GiNaC::pow(y, 5)
        + GiNaC::pow(coni_ex, 4) * x_sub * GiNaC::pow(y, 5)
        - 2 * GiNaC::pow(coni_ex, 3) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 5)
        + 2 * GiNaC::pow(coni_ex, 4) * x_sub * GiNaC::pow(y, 6)
        - GiNaC::pow(coni_ex, 3) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 6)
        - GiNaC::pow(coni_ex, 4) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 6)
        - GiNaC::pow(coni_ex, 3) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 7)
        + 4 * GiNaC::pow(coni_ex, 4) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 7)
        - 3 * GiNaC::pow(coni_ex, 5) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 7)
        + 2 * GiNaC::pow(coni_ex, 4) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 8)
        - 2 * GiNaC::pow(coni_ex, 5) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 8)
        + 2 * GiNaC::pow(coni_ex, 4) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 9)
        - 2 * GiNaC::pow(coni_ex, 5) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 9)
        - GiNaC::pow(coni_ex, 5) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 10)
        + GiNaC::pow(coni_ex, 6) * GiNaC::pow(x_sub, 3) * GiNaC::pow(y, 10)
        - GiNaC::pow(coni_ex, 5) * GiNaC::pow(x_sub, 2) * GiNaC::pow(y, 11)
        + GiNaC::pow(coni_ex, 6) * GiNaC::pow(x_sub, 3) * GiNaC::pow(y, 11);    
    poly = poly * GiNaC::pow(y, -f);
    std::stringstream F;
    F << poly.expand().collect(y);
    spdlog::debug(F.str());
    return poly;
}
*/

auto H_trefoil(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
    std::array<int32_t, 4> M = {1, 0,0,1};
    auto x_sub = GiNaC::pow(x, M.at(0)) * GiNaC::pow(y, M.at(2));
    auto y_sub = GiNaC::pow(x, M.at(1)) * GiNaC::pow(y, M.at(3));
    auto Q = complex_to_ex(kQ_coni);
    auto poly = 
        1 - Q * y 
        + x * GiNaC::pow(y, 3) * (1 - y + 2 * GiNaC::pow(y,2) 
                - 2 * Q * GiNaC::pow(y,3) + GiNaC::pow(Q,2) * GiNaC::pow(y,4)) 
        - GiNaC::pow(x,2) * GiNaC::pow(y,9) * (1 - y);
    return poly;

}

auto F_fig_8(const GiNaC::symbol& y, const GiNaC::symbol& x) -> GiNaC::ex {
  std::array<int32_t, 4> M = {1, 0,-1,1};
  auto x_sub = GiNaC::pow(x, M.at(0)) * GiNaC::pow(y, M.at(2));
  auto y_sub = GiNaC::pow(x, M.at(1)) * GiNaC::pow(y, M.at(3));
  auto Q = complex_to_ex(kQ_coni);

  auto poly = x_sub - x_sub * x_sub 
      - 2 * Q * x_sub * y_sub + 2 * x_sub * x_sub * y_sub 
      - Q*Q *y_sub * y_sub 
      + GiNaC::pow(x_sub, 3) * y_sub * y_sub 
      + Q * Q * GiNaC::pow(y_sub,3) 
      - Q * GiNaC::pow(x_sub,3) * GiNaC::pow(y_sub, 3) 
      + 2 *  Q * Q * x_sub * GiNaC::pow(y_sub, 4)
      - 2 * Q * Q * GiNaC::pow(x_sub,2) * GiNaC::pow(y_sub,4) 
      - GiNaC::pow(Q,3) * x_sub * GiNaC::pow(y_sub, 5) 
      + 1 * Q * Q * x_sub * x_sub * GiNaC::pow(y_sub, 5);
  return poly * GiNaC::pow(y_sub, 2);
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
    initial_integration();
    auto path_it = get_iterator_by_id(new_paths_, 0);
    for (uint32_t k = 0; k < new_paths_.size(); k += 1){
        auto path_it = get_iterator_by_id(new_paths_, k);
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        save_data(path_it->path_id_);
    }
    uint32_t k = new_paths_.size();
    two_path_intersection_handler(18, 21, true, true, 0, 0, false, false); 
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

    two_path_intersection_handler(9, 21, true, true, 0, 0, false, false); 
    path_it = get_iterator_by_id(new_paths_, k);
    save_data(path_it->path_id_);
    evolve_path(path_it, cutoff);
    save_data(path_it->path_id_);
    save_data(25);
    save_data(9);
    save_data(21);
    k++;

    two_path_intersection_handler(10, 25, true, false, 0, 0, false, false); 
    save_data(10);
    path_it = get_iterator_by_id(new_paths_, 26);
    save_data(path_it->path_id_);
    evolve_path(path_it, cutoff);
    save_data(path_it->path_id_);
    save_data(26);


    two_path_intersection_handler(5, 25, false, false, 0, 0, false, false); 

    path_it = get_iterator_by_id(new_paths_, 5);
    path_it->truncate(0, 500);
    two_path_intersection_handler(5, 19, false, true, 0, 0, false, false); 
    
    path_it = get_iterator_by_id(new_paths_, 3);
    path_it->truncate(0, 750);

    two_path_intersection_handler(25, 3, true, true, 0, 0, false, false); 
    save_data(19);
    save_data(3);
    save_data(25); 

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
    cplx x = ramification_points_.at(2).at(kIndexX) + 0.005; 
    auto fiber = curve_->get_fiber(x);
    state_type v;
    for(uint32_t i = 0; i < 4; i++) {
        for(uint32_t j = 0; j < 4; j++) {
            if(i == j) {
                continue;
            }
            for(int32_t k = 0; k <= 0; k++) {
                v.at(kIndexX) = x;
                v.at(kIndexY1) = std::log(fiber.at(i)) + 2 * pi * k * J;
                v.at(kIndexY2) = std::log(fiber.at(j));
                uint32_t index = new_paths_.size();
                add_new_path(v);
                auto path_it = get_iterator_by_id(new_paths_, index);
                evolve_path(path_it, kCutoff);
                save_data(path_it->path_id_);
            }
        }
    }
}

auto Coni::evolve_all(double cutoff) -> void { 
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
    uint32_t k0 = 0;
    for (uint32_t l = 0; l < 2; l++) {
        for (uint32_t m = 0; m < k; m += 1){
            if(m >= 12 && m < 60) {
                continue;
            }
            for (uint32_t n = m + 1; n < k; n += 1){
                if(n >= 12 && n < 60) {
                    continue;
                }

                if (n >= k0) {
                    for (uint32_t q = 0; q < 1; q++) {
                        two_path_intersection_handler(m, n, false, false, 3, q, true, false);
                    }
                }
            }
        }
        k0 = k;
        for ( ;k < new_paths_.size(); k += 1) {
            path_it = get_iterator_by_id(new_paths_, k);
            if (path_it->get_endmass() > cutoff) {
                continue;
            }
            save_data(path_it->path_id_);
            evolve_path(path_it, cutoff);
            path_it->truncate_mass(cutoff);
            path_it->truncate_x(50);
            save_data(path_it->path_id_);
        }
    }
}

auto Coni::custom_BPS(double cutoff) -> void { 
    evolve_all(cutoff);
    // probe_curve();
}

/*
    // evolve_all(cutoff);
    cplx x = 15.5;
    auto fib = curve_->get_fiber(x);
    auto k = new_paths_.size();
    for(size_t j = 0;  j < fib.size(); ++j) {
        auto y = fib.at(j);
        state_type v = {x, std::log(y), std::log(y) + 2 * pi * J};
        add_new_path(v);
        auto path_it = get_iterator_by_id(new_paths_, k);
        evolve_path(path_it, cutoff);
        path_it->truncate_mass(4*pi*pi + 0.1);
        save_data(k);
        k++;
        auto target = std::exp(path_it->get_endpoint().at(kIndexY1));
        size_t nearest_idx = 0;
        double min_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < fib.size(); ++i) {
            double dist = std::abs(fib.at(i)- target); // Euclidean distance in complex plane
            if (dist < min_dist) {
                min_dist = dist;
                nearest_idx = i;
            }
        }
        spdlog::debug("Sheet {} -> Sheet {}.", j, nearest_idx);
        spdlog::debug("Endpoint of {}: {}, Sheet {}: {}.", 
                j, 
                complex_to_string(target), 
                nearest_idx, 
                complex_to_string(fib.at(nearest_idx))
            );
    }
    state_type v = {x, std::log(fib.at(1)), std::log(fib.at(2))};
    add_new_path(v);
    auto path_it = get_iterator_by_id(new_paths_, k);
    evolve_path(path_it, cutoff);
    save_data(k);
    k++;

    for(auto &y : fib) {
        std:: cout << complex_to_string(std::log(y)) << " ";
    }
    std::cout << std::endl;
}
*/
