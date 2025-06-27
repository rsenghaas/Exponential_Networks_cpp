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

auto Coni::custom_BPS_trifoil(double cutoff) -> void {

    auto path_it = get_iterator_by_id(new_paths_, 0);
    for (uint k = 0; k < new_paths_.size(); k += 1){
        auto path_it = get_iterator_by_id(new_paths_, k);
        save_data(path_it->path_id_);
        evolve_path(path_it, cutoff);
        save_data(path_it->path_id_);
    }
    path_it = get_iterator_by_id(new_paths_, 9);
    path_it->truncate(0, 500);
    two_path_intersection_handler(5, 9, true, false, 0, 3, false, false);
    save_data(5);
    two_path_intersection_handler(7, 9, true, true, 0, 0, false, false);
    save_data(7);
    save_data(9);
    /*two_path_intersection_handler(15, 23, true, false, 0, 0, false, false); 
    save_data(15);
    save_data(24);
    two_path_intersection_handler(17, 23, false, true, 0, 0, false, false); 
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
  std::vector<uint32_t> relevant{0, 1,3,4,5,6,7};
  for (auto path_it = new_paths_.end(); path_it != new_paths_.begin();
       path_it--) {
    if (std::find(relevant.begin(), relevant.end(),
                  std::prev(path_it)->path_id_) == relevant.end()) {
      new_paths_.erase(std::prev(path_it));
    }
  }
  for (auto id_A = relevant.begin(); id_A != relevant.end(); id_A++) {
    auto path_it = get_iterator_by_id(new_paths_, *id_A);
    evolve_path(path_it, kCutoff * ratio_zero.at(path_it->path_id_));
    save_data(*id_A);
  }

  std::fstream data_file;
  std::string filename = "data/map_data/intersections.csv";
  std::string sep = ",";
  data_file.open(filename, std::ios::out);

  auto path_it = get_iterator_by_id(new_paths_, 13);
  path_it->truncate(0, 2000);
  save_data(13);
  path_it = get_iterator_by_id(new_paths_, 24);
  path_it->truncate(0, 1600);
  save_data(24);

  for (auto id_A = relevant.begin(); id_A != relevant.end(); id_A++) {
    auto path_it_A = get_iterator_by_id(new_paths_, *id_A);
    SinglePathMap path_A_map = SinglePathMap(path_it_A);
    Map two_path_map;
    two_path_map.add_path(path_A_map.pp_vec);
    for (auto id_B = id_A + 1; id_B != relevant.end(); id_B++) {
      auto path_it_B = get_iterator_by_id(new_paths_, *id_B);
      spdlog::debug("Intersections between {} and {}.", path_it_A->path_id_,
                    path_it_B->path_id_);
      SinglePathMap path_B_map = SinglePathMap(path_it_B);
      auto inters = two_path_map.get_intersections(path_B_map.pp_vec);

      for (auto inter_it = inters.begin(); inter_it != inters.end();
           inter_it++) {
        if (inter_it->times.at(0).at(0) == 0) {
          continue;
        }
        state_type next_state = {0, 0, 0};
        if (compute_intersection_points(*inter_it, path_it_A, path_it_B, 0,
                                        next_state)) {
          print_intersection(*inter_it);
          print_state_type(next_state);
          if ((*id_A == 7 && *id_B == 29)) {
            add_new_path(next_state);
            path_it_A->truncate(0, inter_it->times.at(0).at(1));
            save_data(*id_A);
            path_it_B->truncate(0, inter_it->times.at(1).at(1));
            save_data(*id_B);
            path_it = std::prev(new_paths_.end());
            evolve_path(path_it, 120);
            save_data(path_it->path_id_);
          }
          std::string output_line =
              complex_to_string(next_state.at(kIndexX)) + sep +
              complex_to_string(next_state.at(kIndexY1)) + sep +
              complex_to_string(next_state.at(kIndexY2)) + sep +
              std::to_string(path_it_A->path_id_) + sep +
              std::to_string(path_it_B->path_id_) + "\n";
          data_file << output_line;
        }
      }
    }
  }
  /* auto path_it_A = get_iterator_by_id(new_paths_, 13);
  SinglePathMap path_A_map = SinglePathMap(path_it_A);
  Map two_path_map;
  two_path_map.add_path(path_A_map.pp_vec);
  auto path_it_B = get_iterator_by_id(new_paths_, 30);
  spdlog::debug("Intersections between {} and {}.", path_it_A->path_id_,
  path_it_B->path_id_);

  SinglePathMap path_B_map = SinglePathMap(path_it_B);
  auto inters = two_path_map.get_intersections(path_B_map.pp_vec);
  for(auto inter_it = inters.begin(); inter_it != inters.end(); inter_it++) {
    state_type next_state = {0,0,0};
    if(compute_intersection_points(*inter_it, path_it_A, path_it_B, 0,
  next_state)) { print_intersection(*inter_it); print_state_type(next_state);
      std::string output_line =
            complex_to_string(next_state.at(kIndexX)) + sep +
            complex_to_string(next_state.at(kIndexY1)) + sep +
            complex_to_string(next_state.at(kIndexY2)) + sep +
            std::to_string(path_it_A->path_id_) + sep +
  std::to_string(path_it_B->path_id_) + "\n"; data_file << output_line;
    }
  } */
  data_file.close();
}

auto Coni::custom_BPS(double cutoff) -> void { custom_BPS_trifoil(cutoff); }
