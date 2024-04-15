#include <cfenv>
#include <cmath>
#include <iostream>

#include "element.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

struct isubs : public virtual element {
  void subs(const iadd& b) {
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
};