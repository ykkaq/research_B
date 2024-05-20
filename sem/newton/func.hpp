#include "init.hpp"

void pointer(void) {
  double (*g)(double) = pt::f;
  double (*dg)(double) = pt::df;

  double xo = 2.0;
  double xn = 3.0;

  while (dabs(xn - xo) > 1.0e-10) {
    xo = xn;
    xn = pt::newton(g, dg, xo);
  }

  std::cout << std::scientific << xn << std::endl;
}
