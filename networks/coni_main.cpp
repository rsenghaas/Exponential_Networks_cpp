#include <spdlog/spdlog.h>

#include "coni.hpp"
#include "files.hpp"
#include "magic_numbers.h"

constexpr double kDefaultTheta = std::numbers::pi / 2.0;

auto main(int argc, char* argv[]) -> int {
  create_directories();
  double theta{kDefaultTheta};
  if (argc > 1) {
    std::string theta_string = std::string(argv[1]);
    theta = std::stod(theta_string);
    std::cout << "\u03D1 = " << theta << std::endl;
  }
  double cutoff{kCutoff};
  if (argc > 2) {
    std::string cutoff_string = std::string(argv[2]);
    cutoff = std::stod(cutoff_string);
  }
  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Welcome to the Conifold Networks!");

  Coni coni(theta);
  coni.custom_BPS(cutoff);
  spdlog::info("You made it through the network!");
}
