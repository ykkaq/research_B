#include <cfenv>
#include <cmath>
#include <iostream>

#include "subs.hpp"

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

template <typename T>

class interval : protected isubs<T> {
 public:
  void setting(const T& i, const T& s) {
    this->inf = i;
    this->sup = s;
  }

  void disp(void) const {
    std::cout << "[" << this->inf << "," << this->sup << "]" << std::endl;
  }

  friend interval<T> operator+(const interval<T>& a, const interval<T>& b) {
    interval<T> res = a;
    res.add(b);
    return res;
  }

  friend interval<T> operator-(const interval)

      friend interval<T> operator+(const interval<T>& a, const interval<T>& b) {
    interval res = a;
    res.add(b);
    return res;
  }

  friend interval<T> operator-(const interval<T>& a, const interval<T>& b) {
    interval res = a;
    res.subs(b);
    return res;
  }

  // 野良関数をclass内に書くにはfriend
  //     → メソッド関数ではない!! * a.xxxx(  )
  friend interval<T> sqrt(const interval<T>& a) {
    interval<T> res;
    volatile T asup = a.sup, ainf = a.inf;
    volatile T rsup, rinf;
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

struct F {
  friend double operator()(double a) { return a * a; }
};