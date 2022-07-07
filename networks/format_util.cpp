#include "format_util.hpp"

auto complex_to_string(const std::complex<double> &z) -> std::string {
    std::string sign;
    sign = (z.imag() > 0) ?  "+" : "-";
    return fmt::format("{}{}{}j", z.real(), sign, std::abs(z.imag()));
}