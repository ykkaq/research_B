
// pointer
double dabs(double x) { return x < 0 ? -x : x; }

namespace pt {

double f(double x) { return x * x * x - 2; }

double df(double x) { return 3 * x * x; }

double newton(double (*f)(double), double (*df)(double), double x0) {
  return x0 - f(x0) / df(x0);
}

}  // namespace pt

// object
namespace obj {
// 方針：
// - newton関数に，fとdfをオブジェクトとして渡したい．

struct Object {
  double operator()(double xn) { return xn - pt::f(xn) / pt::df(xn); }
};

}  // namespace obj