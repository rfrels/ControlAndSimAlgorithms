#pragma once
#include <vector>
#include "systems.h"


/**
 * @brief Abstract base class for numerical ODE solvers.
 *
 * This class defines the interface for solving systems of first-order ordinary
 * differential equations using various numerical integration methods.
 */
class Solver {
    public:
        /**
         * @brief Solves a system of ODEs using the implemented numerical method.
         *
         * @param x0 Initial state vector at time t0
         * @param t0 Initial time
         * @param te Final time (might not be reached exactly, if it doesn't line up with the step size)
         * @param h Step size for numerical integration
         * @param sys The system of ODEs to solve
         *
         * @return A 2D vector where the first row contains time values and
         *         subsequent rows contain the state vector at each time step.
         *         Format: [[t0, t1, t2, ...], [x0, x1, x2, ...], ...]
         *
         * @note Pure virtual function - must be implemented by derived classes
         */
        virtual std::vector<std::vector<double>> solve(const std::vector<double>& x0,
        double t0, double te, double h, const System& sys) = 0;
        static void append_state(std::vector<std::vector<double>>& x, const std::vector<double>& x_new,
            const size_t& dim);
};

class ExplicitOneStep: public Solver {
    public:
        virtual std::vector<std::vector<double>> solve(const std::vector<double>& x0,
            double t0, double te, double h, const System& sys) override;
    private:
        virtual std::vector<double>& step(std::vector<double>& x, double h, const System& sys) = 0;
};

class EulerCauchy: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};

class RungeKutta: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};