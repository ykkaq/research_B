#include <iostream>
#include <cfenv>
#include <cmath>

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template <typename _T>
struct element{
    _T inf;
    _T sup;
};