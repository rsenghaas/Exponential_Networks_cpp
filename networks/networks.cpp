#include <spdlog/spdlog.h>

#include "adhm.hpp"
#include "files.hpp"
#include "nets.hpp"

auto main(int argc, char *argv[]) -> int {
  create_directories();

  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Welcome to Networks!");
  double theta = -0.04;
  // Network net(H_c3, theta);
  ADHM adhm(theta);
  adhm.BPS_state();

  spdlog::info("You made it through the network!");
  return 0;
}
