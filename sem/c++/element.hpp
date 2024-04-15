#include <cfenv>
#include <cmath>
#include <iostream>

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

struct element {
  double inf;
  double sup;
};