#include <iostream>

#include "interval.hpp"

int main(void) {
  interval a, b, c;
  a.setting(1.0, 2.0);
  b.setting(2.0, 3.0);

  c = a + b;
  // c = a.add(b);
  c.disp();
  c = sqrt(c);
  c.disp();

  interval d;
  d = a - b;
  d.disp();
}
