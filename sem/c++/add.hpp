#include <cfenv>
#include <cmath>
#include <iostream>

#include "element.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

struct iadd : public virtual element {
  void add(const iadd& b) {
    volatile double ainf = this->inf, asup = this->sup;
    volatile double binf = b.inf, bsup = b.sup;
    volatile double rinf, rsup;

    fesetround(FE_DOWNWARD);
    rinf = ainf + binf;

    fesetround(FE_UPWARD);
    rsup = asup + bsup;

    fesetround(FE_TONEAREST);
    this->sup = rsup;
    this->inf = rinf;
  }
};