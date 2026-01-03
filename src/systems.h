#pragma once
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <stacktrace>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

class System{
    public:
        virtual std::vector<double> operator()(std::vector<double> x) const = 0;
        virtual double exact_sol(double t) const = 0;
};

class UnderdampedHarmonicOscillator: public System {
    /**
    Assumes ODE: m*x_dd * + d * x_d + c * x = 0

    Assumes ODE: x_dd * + 2 * gamma * x_d + w_0² * x = 0
    **/
    // x0 = {1.0, 0.0}
    public:
        //UnderdampedHarmonicOscillator(double m, double d, double c) {
        UnderdampedHarmonicOscillator(double gamma, double omega0, std::vector<double> x0);
        std::vector<double> operator()(std::vector<double> x) const override;
        double exact_sol(double t) const override;
    private:
        //double gamma, omega, zeta;
        double gamma, omega0;
        std::vector<double> x0;
};