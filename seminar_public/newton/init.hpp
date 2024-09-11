
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
struct f {
  double operator()(double x) { return x * x + 2 * x - 4; }
};

struct df {
  double operator()(double x) { return 2 * x + 2; }
};

template <typename F, typename DF>
double newton(F g, DF dg, double x) {
  return x - g(x) / dg(x);
}
}  // namespace obj

namespace vt {
class base {
 public:
  double virtual f(double x) { return 1; }

  double virtual df(double x) { return 0; }

  double newton(double x) { return x - f(x) / df(x); }
};

class advanced : public base {
 public:
  double f(double x) override { return x * x - 3 * x - 5; }

  double df(double x) override { return 2 * x - 3; }
};
}  // namespace vt