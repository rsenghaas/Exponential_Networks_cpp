#include "arb_util.hpp"

auto arb_to_double(const arb_t x) -> double {
  std::string arb_str = arb_get_str(x, 50, 0);
  size_t pos = arb_str.find(' ');
  if (pos != std::string::npos) {
    arb_str = arb_str.substr(1, pos);
  }
  double res = std::stold(arb_str);
  return res;
}

auto acb_to_cplx(const acb_t z) -> cplx {
  arb_t re, im;
  arb_init(re);
  arb_init(im);
  acb_get_real(re, z);
  acb_get_imag(im, z);
  double re_d = arb_to_double(re);
  double im_d = arb_to_double(im);
  cplx res(re_d, im_d);
  arb_clear(re);
  arb_clear(im);
  return res;
}

auto acb_abs_diff_small(const acb_t a,
                        const acb_t b,
                        double eps) -> bool
{
    slong prec = 64;
    acb_t diff;
    acb_init(diff);
    acb_sub(diff, a, b, prec);

    arb_t absval;
    arb_init(absval);
    acb_abs(absval, diff, prec);

    arb_t eps_arb;
    arb_init(eps_arb);
    arb_set_d(eps_arb, eps);

    bool result = arb_lt(absval, eps_arb);

    arb_clear(eps_arb);

    arb_clear(absval);
    acb_clear(diff);
    return result;
}

bool acb_abs_too_small_or_large(const acb_t a,
                                double bound)
{
    slong prec = 64;
    arb_t absval;
    arb_init(absval);
    acb_abs(absval, a, prec);

    arb_t arb_bound;
    arb_init(arb_bound);
    arb_set_d(arb_bound,bound);

    bool too_small = arb_lt(absval, arb_bound);
    arb_inv(arb_bound, arb_bound, prec);
    bool too_large = arb_gt(absval, arb_bound);

    arb_clear(arb_bound);
    arb_clear(absval);
    return too_small || too_large;
}
