#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>
#include <memory>
#include "systems.h"
#include "solvers.h"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

struct SystemType {
  enum Value {
    UHO
  };
  static std::unique_ptr<System> get_system(int systemtype, std::vector<double> x0) {
    if (systemtype == SystemType::UHO) return std::make_unique<UnderdampedHarmonicOscillator>(0.1, 2.0, x0);
    else return nullptr;
  }
};

struct SimulationType {
  enum Value {
    EC
  };
  static std::unique_ptr<Solver> get_simulation(int simulationtype) {
    if (simulationtype == SimulationType::EC) return std::make_unique<EulerCauchy>();
    else return nullptr;
  }
};

struct ControlType {
  enum Value {
    PID
  };
};



std::vector<std::vector<double>> solve_exact_and_sim(double te, double h, std::vector<double> x0,
    int systemtype, int simulationtype) {
    std::unique_ptr<System> system_ptr = SystemType::get_system(systemtype, x0);
    std::unique_ptr<Solver> simulation_ptr = SimulationType::get_simulation(simulationtype);
    std::vector<std::vector<double>> output;

    // Simulated Solution
    output = simulation_ptr->solve(x0, 0., te, h, *system_ptr);
    // Exact Solution
    std::vector<double> exact_solution(output[0].size()); // number of included time steps
    for (size_t i=0; i < exact_solution.size(); i++) {
        exact_solution[i] = system_ptr->exact_sol(output[0][i]);
    }
    output.push_back(exact_solution);
    pybind11::print("C++ done");
    // output[0]: t, output[1]: sim x1, ... output[n] = sim xn, output[n+1]: exact x1, for n system dimension
    return output;
}

//std::vector<std::vector<double>> solve_sim(double te, double h, std::vector<double> x0,
//    int systemtype, int simulationtype) {

//std::vector<std::vector<double>> solve_sim_and_control(double te, double h, std::vector<double> x0,
//    int systemtype, int simulationtype, int controltype) {}

namespace py = pybind11;

PYBIND11_MODULE(control_and_sim_algorithms, module) {
  py::module submodule_system = module.def_submodule("SystemType", "Type enumerator");
  py::enum_<SystemType::Value>(submodule_system, "SystemType")
    .value("UHO", SystemType::UHO)
    .export_values();

  py::module submodule_simulation = module.def_submodule("SimulationType", "Type enumerator");
  py::enum_<SimulationType::Value>(submodule_simulation, "SimulationType")
    .value("EC", SimulationType::EC)
    .export_values();

    module.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: python_example

        .. autosummary::
           :toctree: _generate


    )pbdoc";

    module.def("solve_exact_and_sim", &solve_exact_and_sim, R"pbdoc(
        Add something to a list

        Some other explanation about the add function.
    )pbdoc");

#ifdef VERSION_INFO
    module.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    module.attr("__version__") = "dev";
#endif
}
