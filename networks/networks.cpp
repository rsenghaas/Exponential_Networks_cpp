#include <string>
#include <iostream>

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
  ADHM adhm(theta);
  adhm.BPS_state();

  spdlog::info("You made it through the network!");
  return 0;
}
