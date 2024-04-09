#include <iostream>

#include "interval.hpp"

int main(void) {
  interval a, b, c;
  a.setting(1.0, 2.0);
  b.setting(2.0, 3.0);

  //c = a.add(b);
  c = a+b;
  c.disp();
}