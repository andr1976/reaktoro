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

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ActivityModel.hpp>
#include <Reaktoro/Core/Params.hpp>

namespace Reaktoro {

/// The parameters in the electrolyte NRTL activity model for aqueous electrolyte solutions.
/// Implements the Aspen Plus ELECNRTL model (Chen-Britt-Boston-Evans 1982; Chen & Evans 1986):
/// Gᴱ = Gᴱ(PDH, long range) + Gᴱ(Born) + Gᴱ(local composition, NRTL).
/// Temperature laws (Tref = 298.15 K):
///   molecule-molecule:   τ_ij = a + b/T + e·lnT + f·T ;  α_ij = c + d·(T-273.15)
///   electrolyte-molecule: τ = C + D/T + E·[(Tref-T)/T + ln(T/Tref)] ;  α = N (GMELCN)
struct ActivityModelParamsElectrolyteNRTL
{
    /// Molecule-molecule NRTL parameters, per species pair (i, j).
    /// coeffs = [ a_ij, a_ji, b_ij, b_ji, alpha, dalpha, e_ij, e_ji, f_ij, f_ji ]
    /// (alpha = NRTL/3 nonrandomness, symmetric; dalpha = NRTL/4 its T-slope).
    Tuples<String, String, Vec<real>> nrtl;

    /// Electrolyte-molecule pair parameters, per (cation, anion, molecule) triple.
    /// coeffs = [ C_caB, D_caB, E_caB, C_Bca, D_Bca, E_Bca, alpha ]  (GMELCC/GMELCD/GMELCE/GMELCN)
    Tuples<String, String, String, Vec<real>> gmelc;

    /// Electrolyte-electrolyte pair parameters (shared common ion), per (elec1, elec2) written
    /// as cation-anion strings; coeffs = [ C_12, D_12, E_12, C_21, D_21, E_21, alpha ].
    /// Optional — absent (τ=0, α=0.2) for the AMP-PZ-CO2-H2O system.
    Tuples<String, String, Vec<real>> gmelcee;

    /// Pure-solvent dielectric-constant coefficients (CPDIEC), per solvent species.
    /// coeffs = [ A, B, C ] with eps_B(T) = A + B·(1/T - 1/C). Water uses a built-in correlation.
    Tuples<String, Vec<real>> dielectric;

    /// Born radius (m) per ionic species. Default 3.0e-10 m if omitted.
    Tuples<String, real> bornradius;

    /// Liquid molar volume (m3/mol) per non-water solvent (species with dielectric data), used for the
    /// mixed-solvent density in the Pitzer-Debye-Hückel term; coeffs = [ V0, dV/dT ] with
    /// V(T) = V0 + dV/dT·(T - 298.15). Water uses Mw/997 kg/m3. Solvents without an entry use water's
    /// mass density.
    Tuples<String, Vec<real>> molarvolume;
};

/// Return the electrolyte NRTL activity model for aqueous electrolyte phases.
auto ActivityModelElectrolyteNRTL() -> ActivityModelGenerator;

/// Return the electrolyte NRTL activity model with given parameters.
auto ActivityModelElectrolyteNRTL(ActivityModelParamsElectrolyteNRTL const& params) -> ActivityModelGenerator;

/// Return the electrolyte NRTL activity model with parameters in a Params object.
auto ActivityModelElectrolyteNRTL(Params const& params) -> ActivityModelGenerator;

} // namespace Reaktoro
