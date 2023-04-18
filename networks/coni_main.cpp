#include <spdlog/spdlog.h>

#include "coni.hpp"
#include "files.hpp"

constexpr double kDefaultTheta = -0.01;

auto main(int argc, char* argv[]) -> int {
    create_directories();
    double theta{kDefaultTheta};
    spdlog::set_level(spdlog::level::debug);
    spdlog::info("Welcome to the Conifold Networks!");

    Coni coni(theta);
    coni.custom_BPS();
    spdlog::info("You made it through the network!");
}

