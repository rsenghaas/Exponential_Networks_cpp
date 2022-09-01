#include <string>
#include <iostream>
#include <cstdint>

#include "spdlog/spdlog.h"
#include "adhm.hpp"
#include "files.hpp"
// #include "nets.hpp"

constexpr double kDefaultTheta = -0.01;

auto main(int argc, char *argv[]) -> int {
  create_directories();
  double theta{kDefaultTheta};
  std::cout << argc << std::endl;
  if(argc > 1) {
    std::string theta_string = std::string(argv[1]);
    theta = std::stod(theta_string);
    std::cout << "\u03D1 = " << theta << std::endl;
  }

  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Welcome to Networks!");
  // Network net(H_c3, theta);
  std::vector<uint32_t> pattern_vec = {1,2,3,2};
  uint32_t k = 0;
  for(auto& i : pattern_vec) {
      k += i;
  }
  cplx Z_total = kD4Cutoff*std::exp(J*kD4angle) + k*kD0Mass;
  double adhm_theta = std::arg(Z_total);
  std::cout << "\u03D1 = " << adhm_theta << ", M: " << std::abs(Z_total) << std::endl;
  ADHM adhm(adhm_theta);
  adhm.BPS_state(pattern_vec);

  spdlog::info("You made it through the network!");
  return 0;
}
