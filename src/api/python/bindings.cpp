#include "src/api/suplex_solver.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <stdexcept>

namespace py = pybind11;

PYBIND11_MODULE(pysuplex, module) {
    module.doc() = "Python bindings for the Suplex optimization solver";
    py::enum_<suplex::SolverStatus>(module, "Status")
        .value("NOT_STARTED", suplex::SolverStatus::NOT_STARTED)
        .value("OPTIMAL", suplex::SolverStatus::OPTIMAL)
        .value("INFEASIBLE", suplex::SolverStatus::INFEASIBLE)
        .value("UNBOUNDED", suplex::SolverStatus::UNBOUNDED)
        .value("INF_OR_UNBD", suplex::SolverStatus::INF_OR_UNBD)
        .value("ITERATION_LIMIT", suplex::SolverStatus::ITERATION_LIMIT)
        .value("TIME_LIMIT", suplex::SolverStatus::TIME_LIMIT)
        .value("NUMERICAL_ERROR", suplex::SolverStatus::NUMERICAL_ERROR)
        .value("USER_INTERRUPT", suplex::SolverStatus::USER_INTERRUPT);
    py::enum_<suplex::Algorithm>(module, "Algorithm")
        .value("AUTO", suplex::Algorithm::AUTO)
        .value("PRIMAL_SIMPLEX", suplex::Algorithm::PRIMAL_SIMPLEX)
        .value("DUAL_SIMPLEX", suplex::Algorithm::DUAL_SIMPLEX)
        .value("IPM", suplex::Algorithm::IPM);
    py::class_<suplex::Suplex>(module, "Solver")
        .def(py::init<>())
        .def("read_mps", [](suplex::Suplex& self, const std::string& filename) {
            if (!self.read_mps(filename)) throw std::runtime_error(self.last_error());
        })
        .def("read_lp", [](suplex::Suplex& self, const std::string& filename) {
            if (!self.read_lp(filename)) throw std::runtime_error(self.last_error());
        })
        .def("set_algorithm", &suplex::Suplex::set_algorithm)
        .def("set_presolve", &suplex::Suplex::set_presolve)
        .def("set_gpu", &suplex::Suplex::set_gpu)
        .def("set_time_limit", &suplex::Suplex::set_time_limit)
        .def("set_threads", &suplex::Suplex::set_threads)
        .def("set_mip_gap", &suplex::Suplex::set_mip_gap)
        .def("solve", &suplex::Suplex::solve, py::call_guard<py::gil_scoped_release>())
        .def_property_readonly("status", [](const suplex::Suplex& self) { return self.solution().status; })
        .def_property_readonly("objective", [](const suplex::Suplex& self) { return self.solution().objective_value; })
        .def_property_readonly("solution", [](const suplex::Suplex& self) {
            const auto& values = self.solution().primal_values;
            py::array_t<double> array(values.size());
            std::copy(values.begin(), values.end(), array.mutable_data());
            return array;
        })
        .def_property_readonly("dual_values", [](const suplex::Suplex& self) {
            const auto& values = self.solution().dual_values;
            py::array_t<double> array(values.size());
            std::copy(values.begin(), values.end(), array.mutable_data());
            return array;
        })
        .def_property_readonly("reduced_costs", [](const suplex::Suplex& self) {
            const auto& values = self.solution().reduced_costs;
            py::array_t<double> array(values.size());
            std::copy(values.begin(), values.end(), array.mutable_data());
            return array;
        })
        .def_property_readonly("solve_time", [](const suplex::Suplex& self) { return self.solution().solve_time_seconds; })
        .def_property_readonly("iterations", [](const suplex::Suplex& self) { return self.solution().num_iterations; })
        .def_property_readonly("nodes", [](const suplex::Suplex& self) { return self.solution().num_nodes; })
        .def_property_readonly("mip_gap", [](const suplex::Suplex& self) { return self.solution().mip_gap; });
    module.attr("__version__") = "0.1.0";
}
