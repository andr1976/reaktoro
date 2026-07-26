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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Models/ActivityModels/ActivityModelElectrolyteNRTL.hpp>
using namespace Reaktoro;

namespace {

// Raw composition (mole numbers, un-normalized) at the test state.
Map<String, double> comp = {
    {"H2O",0.78}, {"CO2",0.002}, {"PZ",0.01}, {"AMP",0.05}, {"PZCOOH",0.005},
    {"PZH+",0.02}, {"PZH2+2",0.001}, {"AMPH+",0.03}, {"PZCOO-",0.02}, {"PZ(COO)2-2",0.005},
    {"H3O+",1e-6}, {"OH-",1e-6}, {"HCO3-",0.03}, {"CO3-2",0.005},
};

// Reference mole-fraction ln gamma* (unsymmetric) from the validated Python engine at T=313.15 K.
// NOTE: the C++ model reports ln_g on the molal scale for solutes (adds ln(x_w)); water is symmetric.
Map<String, double> expected_mf = {
    {"H2O",-0.178605}, {"CO2",-2.235330}, {"PZ",-1.112748}, {"AMP",1.078538}, {"PZCOOH",0.386692},
    {"PZH+",0.140750}, {"PZH2+2",-5.100344}, {"AMPH+",1.736854}, {"PZCOO-",7.653608}, {"PZ(COO)2-2",0.878029},
    {"H3O+",-1.132182}, {"OH-",-0.598759}, {"HCO3-",0.003643}, {"CO3-2",-1.299555},
};

} // namespace

TEST_CASE("Testing ActivityModelElectrolyteNRTL", "[ActivityModelElectrolyteNRTL]")
{
    auto const db = Database::embedded("ElectrolyteNRTL.v2024.yaml");
    auto const params = Params::embedded("ElectrolyteNRTL.v2024.yaml");
    auto const species = db.species();

    ActivityModel fn = ActivityModelElectrolyteNRTL(params)(species);
    ActivityProps props = ActivityProps::create(species.size());

    auto const N = species.size();

    // Build mole-fraction vector in the database species order
    ArrayXr x(N);
    for(auto i = 0; i < N; ++i)
        x[i] = comp.count(species[i].name()) ? comp.at(species[i].name()) : 0.0;
    x = x / x.sum();

    real const T = 313.15; // K
    real const P = 1.0e5;  // Pa

    fn(props, {T, P, x});

    // index of water and ln(x_w) for the molal-scale conversion applied to solutes
    Index iw = 0; for(auto i = 0; i < N; ++i) if(species[i].name() == "H2O") iw = i;
    auto const ln_xw = std::log(x[iw].val());

    for(auto i = 0; i < N; ++i)
    {
        auto const name = species[i].name();
        auto const expected = expected_mf.at(name) + (name == "H2O" ? 0.0 : ln_xw);
        INFO("species = " << name << "  ln_g(C++) = " << props.ln_g[i].val() << "  expected = " << expected);
        CHECK( props.ln_g[i].val() == Approx(expected).margin(5e-3) );
    }
}
