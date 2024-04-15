#include <cfenv>
#include <cmath>
#include <iostream>

#include "add.hpp"
#include "subs.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

class interval : protected iadd, protected isubs {
 public:
  void setting(const double& i, const double& s) {
    this->inf = i;
    this->sup = s;
  }

  void disp(void) const {
    std::cout << "[" << this->inf << "," << this->sup << "]" << std::endl;
  }

  friend interval operator+(const interval& a, const interval& b) {
    interval res = a;
    res.add(b);
    return res;
  }

  friend interval operator-(const interval& a, const interval& b) {
    interval res = a;
    res.subs(b);
    return res;
  }

  // 野良関数をclass内に書くにはfriend
  //     → メソッド関数ではない!! * a.xxxx(  )
  friend interval sqrt(const interval& a) {
    interval res;
    volatile double asup = a.sup, ainf = a.inf;
    volatile double rsup, rinf;
    fesetround(FE_UPWARD);
    rsup = std::sqrt(asup);
    fesetround(FE_DOWNWARD);
    rinf = std::sqrt(ainf);
    fesetround(FE_TONEAREST);
    res.sup = rsup;
    res.inf = rinf;
    return res;
  }
};