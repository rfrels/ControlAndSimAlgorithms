#pragma once
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <stacktrace>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

/**
 * @brief Abstract base class for systems of differential equations.
 */
class System{
    public:
        /**
         * @brief Computes the derivative(s) of the system at a given state.
         *
         * @param x The current state vector
         *
         * @return A vector of derivatives dx/dt, with the same dimension as x
         *
         * @note Must be implemented by derived classes to provide the specific
         *       dynamics of the system
         */
        virtual std::vector<double> operator()(std::vector<double> x) const = 0;

        /**
         * @brief Computes the analytical solution of the system at a given time.
         *
         * @param t The time at which to evaluate the exact solution
         *
         * @return The analytical solution value at time t
         *
         * @note Must be implemented by derived classes. For systems without a known
         *       analytical solution, this method must throw an exception
         */
        virtual double exact_sol(double t) const = 0;
};


/**
 * @brief System modeling an underdamped harmonic oscillator.
 *
 * This class represents the second-order differential equation:
 * x'' + 2*gamma*x' + omega0^2*x = 0
 *
 * The behavior depends on the damping ratio zeta = gamma / omega0:
 * - zeta < 1: Underdamped (oscillatory with exponential decay)
 * - zeta = 1: Critically damped (fastest non-oscillatory decay)
 * - zeta > 1: Overdamped (slow non-oscillatory decay)
 *
 * This implementation focuses on the underdamped case where the solution exhibits
 * characteristic oscillations with exponential envelope decay.
 */
class UnderdampedHarmonicOscillator: public System {
    public:
        /**
         * @brief Constructs an underdamped harmonic oscillator system.
         *
         * @param gamma The damping coefficient (must be non-negative). Controls the rate
         *              of exponential decay of oscillations. Units: [1/time]
         * @param omega0 The natural (undamped) angular frequency (must be positive).
         *               Determines the oscillation frequency in the absence of damping.
         *               Units: [rad/time] or [1/time]
         * @param x0 The initial state vector [x(0), x'(0)], where x(0) is the initial
         *           displacement and x'(0) is the initial velocity
         *
         * @throws std::invalid_argument if gamma >= omega0
         * @throws std::invalid_argument if omega0 <= 0 or gamma < 0
         */
        UnderdampedHarmonicOscillator(double gamma, double omega0, std::vector<double> x0);

        /**
         * @brief Computes the time derivatives of the system state.
         *
         * @param x The state vector [x, x'], where x is displacement and x' is velocity
         *
         * @return A vector [x', x''] containing the displacement derivative (velocity)
         *         and velocity derivative (acceleration)
         *
         * @throws std::invalid_argument if x does not have exactly 2 elements
         */
        std::vector<double> operator()(std::vector<double> x) const override;

        /**
         * @brief Computes the analytical solution for displacement at a given time.
         *
         * @param t The time at which to evaluate the solution (must be non-negative)
         *
         * @return The analytical displacement x(t) at time t
         *
         * @throws std::invalid_argument if t is < 0
         */
        double exact_sol(double t) const override;
    private:
        double gamma, omega0;
        std::vector<double> x0;
};