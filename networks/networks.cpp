#include <spdlog/spdlog.h>
#include <iostream>
#include <complex>
#include <string>
#include <algorithm>
#include <ginac/ginac.h>
#include <memory>
#include <thread>
#include <future>
#include <chrono>
#include <mutex>

#include "files.hpp"
#include "sw_curve.hpp"
#include "magic_numbers.h"
#include "path.hpp"
#include "nets.hpp"

auto H_c3(const GiNaC::symbol &x, const GiNaC::symbol &y) -> GiNaC::ex {
    return -x + GiNaC::pow(y, 2) + y;
}

auto main(int argc, char *argv[]) -> int {
    create_directories();

    spdlog::set_level(spdlog::level::debug);
    spdlog::info("Welcome to Networks!");
    double theta = 0.01;
    Network net(H_c3, theta);
    net.print_ramification_points();
    net.evolution_step();

    spdlog::info("You made it through the network!");
    return 0;
}
