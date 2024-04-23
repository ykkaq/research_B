#include <cfenv>
#include <cmath>
#include <iostream>

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template <typename T>
struct element {
  T inf;
  T sup;
};