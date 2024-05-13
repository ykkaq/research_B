#include <iostream>

#include "interval.hpp"

template <typename T>
T add(T a, T b) {
  return a + b;
}

int main(void) {
  element<double> A, B;
  A.inf = 10.0;
  A.sup = 10.0;

  interval<double> a, b, c;
  a.setting(1.0f, 2.0f);
  b.setting(2.0f, 3.0f);

  c = -a - b;
  // c = a.add(b);
  c.disp();
  c = sqrt(c);
  c.disp();
}