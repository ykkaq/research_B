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
  void subs_aminusb(const isubs<T>& b) {
    volatile double ainf = this->inf, asup = this->sup;
    volatile double binf = b.inf, bsup = b.sup;
    volatile double rinf, rsup;

    fesetround(FE_DOWNWARD);
    rinf = ainf - bsup;

    fesetround(FE_UPWARD);
    rsup = asup - binf;

    fesetround(FE_TONEAREST);
    this->sup = rsup;
    this->inf = rinf;
  }

  void subs_bminusa(const isubs < T & b) {
    volatile double ainf = this->inf, asup = this->sup;
    volatile double binf = b.inf, bsup = b.sup;
    volatile double rinf, rsup;

    fesetround(FE_DOWNWARD);
    rinf = binf - asup;

    fesetround(FE_UPWARD);
    rsup = bsup - ainf;

    fesetround(FE_TONEAREST);
    this->sup = rsup;
    this->inf = rinf;
  }

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