#include "arb_util.hpp"

auto arb_to_double(arb_t x) -> double {
  std::string arb_str = arb_get_str(x, 50, 0);
  size_t pos = arb_str.find(' ');
  if(pos != std::string::npos) {
    arb_str = arb_str.substr(1, pos);
  }
  double res = std::stold(arb_str);
  return res;
}

auto acb_to_cplx(acb_t z) -> cplx {
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


