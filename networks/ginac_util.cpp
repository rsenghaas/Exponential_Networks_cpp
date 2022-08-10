#include "ginac_util.hpp"

auto numeric_to_complex(const GiNaC::numeric &z) -> cplx {
  return GiNaC::real(z).to_double() + J * GiNaC::imag(z).to_double();
}

auto complex_to_ex(const cplx &z) -> GiNaC::ex {
  GiNaC::ex z_ex = z.real() + GiNaC::I * z.imag();
  return z_ex;
}

auto sylvester_matrix(const GiNaC::ex &f, const GiNaC::ex &g,
                      const GiNaC::symbol &x) -> GiNaC::matrix {
  GiNaC::ex Poly_f = f.expand();
  GiNaC::ex Poly_g = g.expand();
  uint32_t deg_f = Poly_f.degree(x);
  uint32_t deg_g = Poly_g.degree(x);

  GiNaC::matrix s_matrix(deg_f + deg_g, deg_f + deg_g);

  for (uint32_t i = 0; i <= deg_f; i++) {
    for (uint32_t j = 0; j < deg_g; j++) {
      s_matrix(j, deg_f - i + j) = Poly_f.coeff(x, deg_f - i);
    }
  }
  for (uint32_t j = 0; j <= deg_g; j++) {
    for (uint32_t i = 0; i < deg_f; i++) {
      s_matrix(deg_g + i, deg_g - j + i) = Poly_g.coeff(x, deg_g - j);
    }
  }
  return s_matrix;
}

auto discriminant(const GiNaC::ex &f, const GiNaC::symbol &x) -> GiNaC::ex {
  GiNaC::matrix s_matrix = sylvester_matrix(f, f.diff(x), x);
  return GiNaC::determinant(s_matrix);
}

auto eval_ex_to_complex(const GiNaC::ex &f, const GiNaC::symbol &x,
                        const cplx &z) -> cplx {
  return numeric_to_complex(
      GiNaC::ex_to<GiNaC::numeric>(f.subs({x == complex_to_ex(z)})));
}

auto roots(const GiNaC::ex &f, const GiNaC::symbol &x) -> std::vector<cplx> {
  // spdlog::debug("Enter zero search.");
  GiNaC::ex Poly = f.expand();
  const uint32_t deg = Poly.degree(x);
  Poly = Poly / Poly.lcoeff(x);
  std::vector<cplx> approx_roots(deg);
  std::generate(std::begin(approx_roots), std::end(approx_roots),
                [i = 0]() mutable {
                  return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(
                      GiNaC::pow(complex_to_ex(kDummyRoot), i++)));
                });
  std::vector<double> norms(deg);
  double max_norm;
  uint32_t iteration_counter = 0;
  do {
    std::vector<cplx> new_approx_roots(deg);
    for (uint32_t i = 0; i < deg; i++) {
      cplx multi = 1;
      for (uint32_t j = 0; j < deg; j++) {
        if (j != i) {
          // std::cout  << "(" << i << "," << j << ")" << "->" <<
          // (approx_roots.at(j) - approx_roots.at(i)) << std::endl;
          multi *= (approx_roots.at(i) - approx_roots.at(j));
        }
      }
      cplx Poly_val = eval_ex_to_complex(Poly, x, approx_roots.at(i));
      new_approx_roots.at(i) = approx_roots.at(i) - Poly_val / multi;
      norms.at(i) =
          std::abs(eval_ex_to_complex(Poly, x, new_approx_roots.at(i)));
    }
    approx_roots = new_approx_roots;
    max_norm = *std::max_element(std::begin(norms), std::end(norms));
    iteration_counter++;
  } while (max_norm > kZerosPrecisions &&
           iteration_counter < kZerosMaxIterations);
  return approx_roots;
}
