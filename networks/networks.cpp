#include <acb.h>
#include <arb.h>

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "adhm.hpp"
#include "elliptic.hpp"
#include "files.hpp"
#include "spdlog/spdlog.h"
// #include "nets.hpp"

constexpr double kDefaultTheta = -0.01003;

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
  std::vector<uint32_t> pattern_vec{1, 1, 2};
  uint32_t k = 0;
  for (auto& i : pattern_vec) {
    k += i;
  }
  k = 4;
  const double epsilon = 0;
  double k_pert = k + epsilon;
  cplx Z_total = kD4Cutoff * std::exp(J * kD4angle) + k_pert * kD0Mass;

  double adhm_theta = std::arg(Z_total);
  std::cout << "\u03D1 = " << adhm_theta << ", M: " << std::abs(Z_total)
            << std::endl;
  // adhm_theta *= -1;
  // adhm_theta += 0.0008;
  // double D0_D4_theta = -0.2;
  // ADHM cutoff(D0_D4_theta);
  // auto cutoff_point = cutoff.get_puncture_point(2500);
  // double D0_D4_mass = cutoff.get_puncture_mass();
  // spdlog::debug(D0_D4_mass);
  // Z_total = D0_D4_mass * std::exp(J * D0_D4_theta) + kD0Mass;
  // adhm_theta = std::arg(Z_total);
  ADHM adhm(adhm_theta);
  // adhm.backwards(pattern_vec);
  // adhm.BPS_state(pattern_vec);
  adhm.custom_BPS();
  spdlog::debug("{}", adhm_theta);
  spdlog::info("You made it through the network!");
}
