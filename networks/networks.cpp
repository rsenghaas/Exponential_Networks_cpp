#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <arb.h>
#include <acb.h>

#include "adhm.hpp"
#include "elliptic.hpp"
#include "files.hpp"
#include "spdlog/spdlog.h"
// #include "nets.hpp"

constexpr double kDefaultTheta = -0.01;

auto main(int argc, char* argv[]) -> int {
  create_directories();
  double theta{kDefaultTheta};
  std::cout << argc << std::endl;
  if (argc > 1) {
    std::string theta_string = std::string(argv[1]);
    theta = std::stod(theta_string);
    std::cout << "\u03D1 = " << theta << std::endl;
  }

  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Welcome to Networks!");
  // Network net(H_c3, theta);
  std::vector<uint32_t> pattern_vec{1,1,2};
  uint32_t k = 0;
  for (auto& i : pattern_vec) {
    k += i;
  }
  k = 2;
  const double epsilon = 0;
  double k_pert = k - epsilon;
  cplx Z_total = kD4Cutoff * std::exp(J * kD4angle) + k_pert * kD0Mass;

  double adhm_theta = std::arg(Z_total);
  std::cout << "\u03D1 = " << adhm_theta << ", M: " << std::abs(Z_total)
            << std::endl;
  // adhm_theta *= -1;
  // adhm_theta += 0.0008;
  Elliptic elliptic(0.0);
  // adhm.backwards(pattern_vec);
  // adhm.BPS_state(pattern_vec);
  elliptic.custom_BPS();
  spdlog::debug("{}", adhm_theta);
  spdlog::info("You made it through the network!");
  
    return 0;
}
