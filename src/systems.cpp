#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stacktrace>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "systems.h"


// UnderdampadHarmonicOscillator
UnderdampedHarmonicOscillator::UnderdampedHarmonicOscillator(double gamma, double omega0, std::vector<double> x0): gamma(gamma),
omega0(omega0), x0(x0) {
    if(gamma >= omega0) {
        std::stringstream ss;
        ss << "Error: gamma: " << gamma << " larger or equal to omega0: " << omega0;
        throw std::invalid_argument(ss.str());
    }
    if(gamma < 0) throw std::invalid_argument("Error: gamma must not be < 0");
    if(omega0 <= 0) throw std::invalid_argument("Error: omega0 must be > 0");
    if(x0.size() != 2) throw std::invalid_argument("Error: x0 must have two elements");
};

std::vector<double> UnderdampedHarmonicOscillator::operator()(std::vector<double> x) const {
    if(x0.size() != 2) throw std::invalid_argument("Error: x0 must have two elements");

    std::vector<double> output(x.size());
    output[0] = x[1];
    output[1] = -2 * gamma * x[1] - omega0 * omega0 * x[0];

    return output;
};

double UnderdampedHarmonicOscillator::exact_sol(double t) const {
    if(t < 0) throw std::invalid_argument("Error: t must not be < 0");
    double Om = std::sqrt(omega0 * omega0 - gamma * gamma);
    double oscillating = x0[0]*std::cos(Om*t) + (gamma*x0[0]+x0[1])*std::sin(Om*t)/Om;
    return std::exp(-gamma*t) * oscillating;
};


