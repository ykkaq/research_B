#include <cfenv>
#include <cmath>
#include <iostream>

#include "element.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template <typename _T>
struct iadd : public element<_T> {
  void add(const iadd<_T>& b) {
    volatile _T ainf = this->inf, asup = this->sup;
    volatile _T binf = b.inf, bsup = b.sup;
    volatile _T rinf, rsup;

    fesetround(FE_DOWNWARD);
    rinf = ainf + binf;

    fesetround(FE_UPWARD);
    rsup = asup + bsup;

    fesetround(FE_TONEAREST);
    this->sup = rsup;
    this->inf = rinf;
  }
};