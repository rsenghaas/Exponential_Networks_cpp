#include "ginac_util.hpp"

#include <iostream>
#include <sstream>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

auto numeric_to_complex(const GiNaC::numeric &z) -> cplx {
  double re = GiNaC::real(z).to_double();
  double im = GiNaC::imag(z).to_double();
  if (std::abs(re) < 1e-100){ re = 0.0; }
  if (std::abs(im) < 1e-100){im = 0.0; }
  return cplx(re,im);
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
  auto disc = GiNaC::determinant(s_matrix).expand();
  return disc;
}

auto eval_ex_to_complex(const GiNaC::ex &f, const GiNaC::symbol &x,
                        const cplx &z) -> cplx {
  return numeric_to_complex(
      GiNaC::ex_to<GiNaC::numeric>(f.subs({x == complex_to_ex(z)})));
}

auto roots(const GiNaC::ex &f, const GiNaC::symbol &x, bool initial_schur) -> std::vector<cplx> {
  // spdlog::debug("Enter zero search.");
  GiNaC::ex Poly = f.expand(); 
  int min_deg = std::numeric_limits<int>::max();
  for (const auto& term : Poly) {
        int power = term.degree(x);
        if (term.coeff(x, power).is_zero())
            continue; // skip zero coefficients
        min_deg = std::min(min_deg, power);
  }
  Poly = (Poly * GiNaC::pow(x, -min_deg)).expand();
  uint32_t deg = Poly.degree(x);
  GiNaC::ex result = 0;

  for (uint32_t i = 0; i <= deg; ++i) {
    GiNaC::numeric coeff = GiNaC::ex_to<GiNaC::numeric>(Poly.coeff(x, i).evalf());
    if (!coeff.is_zero()) 
    {
        result += coeff * GiNaC::pow(x, i);
    }
    else {
        spdlog::debug("Dropped summand {} x^{}", complex_to_string(numeric_to_complex(coeff)), i);
    }
  }
  result = result.expand();
  deg = result.degree(x);
  Poly = result / result.lcoeff(x);
  std::vector<cplx> approx_roots(deg);
  if (initial_schur) {
    approx_roots = schur_roots(Poly,x);
  }
  else {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(-1.0,1.0);
      std::generate(std::begin(approx_roots), std::end(approx_roots),
                [i = 0, &gen, &dis]() mutable {
                  return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(
                      GiNaC::pow(complex_to_ex(
                            dis(gen) + dis(gen) * 1 *J
                              ), i++)));
                }); 
  }
  std::vector<double> norms(deg);
  double max_norm;
  int max_iterations = initial_schur ? kZerosMaxIterations : 10 * kZerosMaxIterations;
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
      if (std::abs(multi) < 1e-25) {
        multi = 1e-25;  
      }
      cplx Poly_val = eval_ex_to_complex(Poly, x, approx_roots.at(i));
      new_approx_roots.at(i) = approx_roots.at(i) - Poly_val / multi;

      norms.at(i) =
          std::abs(eval_ex_to_complex(Poly, x, new_approx_roots.at(i)));
    }
    approx_roots = new_approx_roots;
    max_norm = *std::max_element(std::begin(norms), std::end(norms));
    iteration_counter++;
  } while (iteration_counter < max_iterations / 10 || max_norm > kZerosPrecisions &&
           iteration_counter < max_iterations);
  /* for (auto& r : approx_roots) {
    std::cout << complex_to_string(std::log(r)) << "," << std::abs(eval_ex_to_complex(Poly, x, r)) << "   ";
  }

  if (iteration_counter >= kZerosMaxIterations && initial_schur) {
    // approx_roots = schur_roots(f,x);
    for (auto& r : approx_roots) {
        std::cout << complex_to_string(std::log(r)) << "," << std::abs(eval_ex_to_complex(Poly, x, r)) << "   ";
     } 
  } */
  /*for (uint32_t i = 0; i < approx_roots.size(); i++) {
    approx_roots.at(i) = refine_root_newton(f, x, approx_roots.at(i));
  }*/
  std::sort(approx_roots.begin(), approx_roots.end(), [](const std::complex<double>& a, const std::complex<double>& b) {
        return a.real() < b.real();
  });
  return approx_roots;
}

auto companion_matrix(const std::vector<cplx>& coefficients) -> Eigen::MatrixXcd {
    uint32_t n = coefficients.size();
    Eigen::MatrixXcd M = Eigen::MatrixXcd::Zero(n - 1, n - 1);
    // std::cout << std::endl << coefficients[n - 1] << std::endl;
    for (uint32_t i = 0; i < n - 1; ++i) {
        M(i, n - 2) = -coefficients[i] / coefficients[n - 1];
    }

    // Fill the subdiagonal with ones
    for (uint32_t i = 0; i < n - 2; ++i) {
        M(i + 1, i) = 1;
    }   
    // std::cout << M << std::endl;
    return M;
}

auto schur_roots(const GiNaC::ex &f, const GiNaC::symbol &x) -> std::vector<cplx> {
  // spdlog::debug("Enter zero search.");
  GiNaC::ex Poly = f.expand();
  int min_deg = std::numeric_limits<int>::max();
  for (const auto& term : Poly) {
        int power = term.degree(x);
        if (term.coeff(x, power).is_zero())
            continue; // skip zero coefficients
        min_deg = std::min(min_deg, power);
  }
  Poly = (Poly * GiNaC::pow(x, -min_deg)).expand();
  std::stringstream F;
  F << Poly.expand().collect(x);
  const uint32_t deg = Poly.degree(x);
  // Poly = Poly / Poly.lcoeff(x);
  std::vector<cplx> coefficients;
  for (uint32_t i = 0; i <= deg; i++) {
        GiNaC::ex coeff = Poly.coeff(x, i);
        coefficients.push_back(eval_ex_to_complex(coeff, x, 0)); 
        // std::cout<< coefficients.at(i);
  }
  std::vector<cplx> approx_roots;

  auto M = companion_matrix(coefficients);
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(M);

  // Get the eigenvalues (roots)
  Eigen::VectorXcd eigenvalues = solver.eigenvalues();

  // Convert eigenvalues to a vector of complex numbers
  // std::cout << std::endl;
  for (auto& e : eigenvalues) {
    approx_roots.push_back(e);
  }
  return approx_roots;
}
auto refine_root_newton(const GiNaC::ex& f, const GiNaC::symbol& x, cplx z0) -> cplx {
    GiNaC::ex df = f.diff(x);
    for (uint32_t i = 0; i <  kZerosMaxIterations; ++i) {
        cplx fz = eval_ex_to_complex(f, x, z0);
        cplx dfz = eval_ex_to_complex(df, x, z0);
        if (std::abs(dfz) < 1e-10) break;  // avoid division by near-zero
        cplx z1 = z0 - fz / dfz;
        if (std::abs(z1 - z0) < 1e-12) return z1;  // convergence
        z0 = z1;
    }
    return z0;
}
