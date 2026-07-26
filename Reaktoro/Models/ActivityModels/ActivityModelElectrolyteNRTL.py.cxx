// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// pybind11 includes
#include <Reaktoro/pybind11.hxx>

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/ActivityModelElectrolyteNRTL.hpp>
using namespace Reaktoro;

void exportActivityModelElectrolyteNRTL(py::module& m)
{
    py::class_<ActivityModelParamsElectrolyteNRTL>(m, "ActivityModelParamsElectrolyteNRTL")
        .def(py::init<>())
        .def_readwrite("nrtl", &ActivityModelParamsElectrolyteNRTL::nrtl, "Molecule-molecule NRTL parameters per species pair.")
        .def_readwrite("gmelc", &ActivityModelParamsElectrolyteNRTL::gmelc, "Electrolyte-molecule pair parameters per (cation, anion, molecule).")
        .def_readwrite("gmelcee", &ActivityModelParamsElectrolyteNRTL::gmelcee, "Electrolyte-electrolyte pair parameters (shared common ion).")
        .def_readwrite("dielectric", &ActivityModelParamsElectrolyteNRTL::dielectric, "Pure-solvent dielectric-constant coefficients (CPDIEC).")
        .def_readwrite("bornradius", &ActivityModelParamsElectrolyteNRTL::bornradius, "Born radius (m) per ionic species.")
        ;

    m.def("ActivityModelElectrolyteNRTL", py::overload_cast<>(ActivityModelElectrolyteNRTL));
    m.def("ActivityModelElectrolyteNRTL", py::overload_cast<ActivityModelParamsElectrolyteNRTL const&>(ActivityModelElectrolyteNRTL));
    m.def("ActivityModelElectrolyteNRTL", py::overload_cast<Params const&>(ActivityModelElectrolyteNRTL));
}
