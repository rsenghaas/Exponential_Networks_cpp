#include "p2.hpp"
#include "ginac_util.hpp"
#include "type_util.hpp"

constexpr cplx kQ_coni =  1.0 + 1.0 * J;
// Q = 11 / 8 - 3 / 4 J
// theta_1 = 0.31
// theta_2 = 0.047
// theta_3 = -0.346
//
// 
//
// CP1_1_4; Q = 1 + i
// t_1 = 0.11
// t_2 = 0.13
// t_3 = 0.28
// t_4 = -0.2
// t_5 = -0.28

auto H_p2(const GiNaC::symbol& x, const GiNaC::symbol& y, const cplx Q) -> GiNaC::ex {
  return (1 + y + x + complex_to_ex(Q) / ( x * y )) * x * y;
}

auto H_1_1_2(const GiNaC::symbol& x, const GiNaC::symbol& y, cplx Q1, cplx Q2) -> GiNaC::ex {
  int f = 0;
  int fy = 0;
  GiNaC::ex x_sub = x * GiNaC::pow(-y, f);
  GiNaC::ex y_sub = y * GiNaC::pow(-x, fy);
  return (1 + x_sub - 2 * complex_to_ex(Q1) * complex_to_ex(Q2) / (x_sub * x_sub)  + complex_to_ex(Q1) 
              * (1 + y_sub) * (1 + y_sub) / (x_sub * x_sub * y_sub)) * x_sub * x_sub * y_sub;
}

auto custom_curve(const GiNaC::symbol& x, const GiNaC::symbol& y) -> GiNaC::ex {
    // 0.4,0.67,-0.2, 0.2
    int f = 0;
    GiNaC::ex x_sub = x * GiNaC::pow(-y, f);
    /* return 1 * GiNaC::pow(x_sub, 2) * y + 1 * y + 
        GiNaC::pow(x_sub,3) * y + (-1.99 * 1 - 3) * x* y + y * y + 1;
        */ 

    // 0.51
    return (8 * GiNaC::pow(x_sub, 3) * y - 12 * GiNaC::pow(x_sub, 2) * y - 48 * x_sub * y + 27 * y * y + 25.99 * y + 27) / 27.0;
}

auto P2::custom_BPS(double cutoff) -> void {

    std::vector<size_t> inv = {};
    for (size_t &k : inv) {
        auto endpoint = new_paths_.at(k).get_endpoint();
        invert_state(endpoint);
        new_paths_.at(k).override_endpoint(endpoint);
    }
    initial_integration();
    auto path_it = get_iterator_by_id(new_paths_, 0); 
    for (uint32_t k = 0; k < new_paths_.size(); k += 1){
        path_it = get_iterator_by_id(new_paths_, k);
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        save_data(path_it->path_id_);
    }
    // two_path_intersection_handler(5, 8, true, true, 1, 0, false, false); 
    // path_it = get_iterator_by_id(new_paths_, 12);
    // save_data(path_it->path_id_);
    // evolve_path(path_it, cutoff);
    // save_data(12);
    // save_data(11);
}
