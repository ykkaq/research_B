#include <cfenv>
#include <cmath>
#include <iostream>

#include "add.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template <typename T>

struct isubs : public iadd<T> {
  void subs(const iadd<T>& b) {
    volatile T ainf = this->inf, asup = this->sup;
    volatile T binf = b.inf, bsup = b.sup;
    volatile T rinf, rsup;

    fesetround(FE_DOWNWARD);
    rinf = ainf - bsup;

    fesetround(FE_UPWARD);
    rsup = asup - binf;

    fesetround(FE_TONEAREST);
    this->sup = rsup;
    this->inf = rinf;
  }
};