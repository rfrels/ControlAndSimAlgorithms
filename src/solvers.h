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

        /**
         * @brief Appends a new state vector to the solution matrix.
         *
         * @param x The solution matrix to which the new state is appended
         * @param x_new The new state vector to append
         * @param dim The dimensionality (size) of the state vector
         */
        static void append_state(std::vector<std::vector<double>>& x, const std::vector<double>& x_new,
            const size_t& dim);
};

/**
 * @brief Abstract base class for explicit one-step ODE solvers.
 */
class ExplicitOneStep: public Solver {
    public:
        /**
         * @brief Implements the standard integration loop for explicit one-step methods.
         *
         * This method iteratively calls the step() function to advance the solution
         * from t0 to te using fixed step size h. It manages the overall integration
         * process while delegating individual step calculations to derived classes.
         *
         * @param x0 Initial state vector at time t0
         * @param t0 Initial time
         * @param te Final time
         * @param h Step size for numerical integration
         * @param sys The system of ODEs to solve
         *
         * @return Solution matrix with time values and state vectors at each step
         *
         * @note The simulation will stop right before te
         */
        virtual std::vector<std::vector<double>> solve(const std::vector<double>& x0,
            double t0, double te, double h, const System& sys) override;
    private:
        /**
         * @brief Pure virtual function that performs a single integration step.
         *
         * This method advances the state vector x by one step of size h according
         * to the specific numerical method implemented by the derived class.
         *
         * @param x The current state vector, which is updated in-place to the next state
         * @param h The step size
         * @param sys The system of ODEs providing derivative information
         *
         * @return Reference to the updated state vector x
         *
         * @note Implemented by specific method classes (Euler, Runge-Kutta, etc.)
         */
        virtual std::vector<double>& step(std::vector<double>& x, double h, const System& sys) = 0;
};

class EulerCauchy: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};

class ModifiedEulerCauchy: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};

class Heun: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};

class Simpson: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};

class RungeKutta: public ExplicitOneStep {
    private:
        std::vector<double>& step(std::vector<double>& x, double h, const System& sys) override;
};