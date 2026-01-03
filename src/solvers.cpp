#include <vector>
#include "systems.h"
#include "solvers.h"

std::vector<double> operator+(const std::vector<double>& a,
                              const std::vector<double>& b) {
    if (a.size() != b.size())
        throw std::invalid_argument("Size mismatch");

    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] + b[i];
    return r;
}

std::vector<double> operator*(double s, const std::vector<double>& v) {
    std::vector<double> r(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        r[i] = v[i] * s;
    return r;
}

// for opposite ordering
std::vector<double> operator*(const std::vector<double>& v, double s) {
    return s * v;
}

// TODO Overload /

void Solver::append_state(std::vector<std::vector<double>>& x, const std::vector<double>& x_new,
    const size_t& dim) {
    for (size_t i = 0; i < dim; i++) {
        x[i].push_back(x_new[i]);
    }
}

std::vector<std::vector<double>> ExplicitOneStep::solve(const std::vector<double>& x0, double t0, double te, double h,
    const System& sys) {
    size_t dim = x0.size();
    std::vector<std::vector<double>> x_out(dim);
    std::vector<double> t;
    append_state(x_out, x0, dim);
    t.push_back(t0);

    std::vector<double> x = x0;

    for (double i = t0+h; i < te; i+=h) {
        x = step(x, h, sys);
        append_state(x_out, x, dim);
        t.push_back(i);
    }
    std::vector<std::vector<double>> output;
    output.push_back(t);
    for (size_t i = 0; i < dim; i++) {
        output.push_back(x_out[i]);
    }
    return output;
}

// Euler Cauchy:
std::vector<double>& EulerCauchy::step(std::vector<double>& x, double h, const System& sys) {
    std::vector<double> f = sys(x);
    x = x + h * f;
    return x;
}

// ModifiedEulerCauchy
std::vector<double>& ModifiedEulerCauchy::step(std::vector<double>& x, double h, const System& sys) {
    std::vector<double> f = sys(x);
    std::vector<double> x_pred = x + h*0.5*f;
    x = x + h*sys(x_pred);
    return x;
}

// Heun
std::vector<double>& Heun::step(std::vector<double>& x, double h, const System& sys) {
    std::vector<double> f = sys(x);
    std::vector<double> x_pred = x + h*f;
    x = x + h * 0.5 * (f + sys(x_pred));
    return x;
}

// RungeKutta:
std::vector<double>& RungeKutta::step(std::vector<double>& x, double h, const System& sys) {
    std::vector<double> k1 = sys(x);
    std::vector<double> k2 = sys(x + h*0.5*k1);
    std::vector<double> k3 = sys(x + h*0.5*k2);
    std::vector<double> k4 = sys(x + h*k3);
    x = x + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * (h / 6.0) ;
    return x;
}


