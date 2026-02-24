#include <spdlog/spdlog.h>

#include "p2.hpp"
#include "files.hpp"
#include "magic_numbers.h"
#include <cln/float.h>
#include <regex>

std::complex<double> parseComplex(const std::string& s)
{
    std::regex r(R"(([+-]?\d*\.?\d+)([+-]\d*\.?\d+)i)");
    std::smatch m;

    if(!std::regex_match(s, m, r))
        throw std::runtime_error("Bad format");

    return {std::stod(m[1]), std::stod(m[2])};
}


constexpr double kDefaultTheta = 0.1; // std::numbers::pi / 2.0 - 0.01;
// constexpr double kDefaultTheta = 0; // std::numbers::pi / 2.0 - 0.01;
// constexpr double kDefaultTheta = 0.06605624168680607;
// constexpr double kDefaultTheta =  0.12261412041038638;
const cplx kDefaultQ = 1.0 + 1.0 * J;
auto main(int argc, char* argv[]) -> int {
  cln::cl_inhibit_floating_point_underflow = true;
  create_directories();
  double theta{kDefaultTheta};
  if (argc > 1) {
    std::string theta_string = std::string(argv[1]);
    theta = std::stod(theta_string);
    std::cout << "\u03D1 = " << theta << std::endl;
  }
  double cutoff{100};
  cplx Q1{kDefaultQ};
  cplx Q2{0};
  if (argc > 2) {
    std::string cutoff_string = std::string(argv[2]);
    Q1 = parseComplex(cutoff_string);
  }
  if (argc > 3) {
    std::string cutoff_string = std::string(argv[3]);
    Q2 = parseComplex(cutoff_string);
  }

  spdlog::set_level(spdlog::level::debug);
  spdlog::info("Welcome to the Conifold Networks!");

  P2 p2(theta, Q1, Q2);
  p2.custom_BPS(cutoff);
  spdlog::info("You made it through the network!");
}
