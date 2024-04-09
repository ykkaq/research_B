#include <cfenv>
#include <iostream>
#include <cmath>

// fesetround(FE_UPWARD);
// fesetround(FE_DOWNWARD);
// fesetround(FE_TOWARDZERO);
// fesetround(FE_TONEAREST);

class interval {
 private:
  double inf;
  double sup;

 public:
  void setting(const double& i, const double& s) {
    this->inf = i;
    this->sup = s;
  }

  void disp() const {
    std::cout << "[" << this->inf << "," << this->sup << "]" << std::endl;
  }

  interval add(const interval& b) const {
    volatile double ainf = this->inf, asup = this->sup, binf = b.inf,
                    bsup = b.sup, rinf, rsup;
    interval res;

    fesetround(FE_DOWNWARD);
    rinf = ainf + binf;

    fesetround(FE_UPWARD);
    rsup = asup + bsup;

    fesetround(FE_TONEAREST);
    res.sup = rsup;
    res.inf = rinf;

    return res;
  }

  friend interval operator+(const interval& a, const interval& b) {
    interval res;
    res = a.add(b);
    return res;
  }

  friend interval sqrt(const interval& a) {
    interval res;
    fesetround(FE_UPWARD);
    res.sup = std::sqrt(a.sup);
    fesetround(FE_DOWNWARD);
    res.inf = std::sqrt(a.inf);
    fesetround(FE_TONEAREST);
    return res;
  }
};
