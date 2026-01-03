#include <cmath>
#include <iostream>
#include <stdexcept>
#include <stacktrace>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "systems.h"


UnderdampedHarmonicOscillator::UnderdampedHarmonicOscillator(double gamma, double omega0, std::vector<double> x0): gamma(gamma),
omega0(omega0), x0(x0) {
    //gamma = d / m;
    //omega = std::sqrt(c/m);
    //zeta = d / (std::sqrt(4.0*m*c));
    //pybind11::print("zeta: ", zeta);
    if(gamma >= omega0) {
        throw std::runtime_error("Error: gamma larger or equal to omega0");
    }

};

std::vector<double> UnderdampedHarmonicOscillator::operator()(std::vector<double> x) const {
    std::vector<double> output(x.size());
    output[0] = x[1];
    output[1] = -2 * gamma * x[1] - omega0 * omega0 * x[0];
    //output[1] = -omega * omega * x[0] -gamma * x[1];

    return output;
};

double UnderdampedHarmonicOscillator::exact_sol(double t) const {
    double Om = std::sqrt(omega0 * omega0 - gamma * gamma);
    double oscillating = x0[0]*std::cos(Om*t) + (gamma*x0[0]+x0[1])*std::sin(Om*t)/Om;
    //pybind11::print("oscillating: ", oscillating);
    return std::exp(-gamma*t) * oscillating;
    //return std::exp(-omega*t*zeta) * std::cos(omega*t*std::sqrt(1 - (zeta * zeta)));
};


