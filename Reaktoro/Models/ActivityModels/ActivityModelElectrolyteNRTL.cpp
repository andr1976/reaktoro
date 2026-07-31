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

#include "ActivityModelElectrolyteNRTL.hpp"

// C++ includes
#include <map>
#include <array>
#include <cmath>

// Reaktoro includes
#include <Reaktoro/Common/Constants.hpp>
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Common/Matrix.hpp>
#include <Reaktoro/Core/Embedded.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Serialization/Models/ActivityModels.hpp>

namespace Reaktoro {
namespace {

const double QE   = 1.602176634e-19;   // elementary charge (C)
const double NA   = 6.02214076e23;     // Avogadro
const double KB   = 1.380649e-23;      // Boltzmann (J/K)
const double EPS0 = 8.8541878128e-12;  // vacuum permittivity (F/m)
const double RHO  = 14.9;              // PDH closest-approach parameter
const double TREF = 298.15;            // reference temperature (K)

/// Water dielectric constant correlation (adequate for the PDH/Born terms).
inline auto epsWater(real const& T) -> real
{
    auto const t = T - 298.15;
    return 78.54*(1.0 - 4.579e-3*t + 1.19e-5*t*t);
}

} // anonymous namespace

auto createActivityModelElectrolyteNRTL(SpeciesList const& species, ActivityModelParamsElectrolyteNRTL const& params) -> ActivityModel
{
    // Create the aqueous solution and gather fixed structural data
    AqueousMixture solution(species);
    AqueousMixtureState aqstate;

    const auto N = species.size();
    const ArrayXd z = solution.charges();
    const Index iw = solution.indexWater();
    const auto Mw = species[iw].molarMass(); // kg/mol

    errorif(iw >= N, "ActivityModelElectrolyteNRTL expects an H2O species in the aqueous phase.");

    // Classify species by charge sign
    Indices icat, ian, imol;
    for(auto i = 0; i < N; ++i)
    {
        if(z[i] > 0.0) icat.push_back(i);
        else if(z[i] < 0.0) ian.push_back(i);
        else imol.push_back(i);
    }

    // --- Molecule-molecule NRTL parameter matrices (constants; index by species) ---
    MatrixXd Amm = MatrixXd::Zero(N, N), Bmm = MatrixXd::Zero(N, N);
    MatrixXd Emm = MatrixXd::Zero(N, N), Fmm = MatrixXd::Zero(N, N);
    MatrixXd almm = MatrixXd::Constant(N, N, 0.2), dalmm = MatrixXd::Zero(N, N);
    for(auto const& [si, sj, c] : params.nrtl)
    {
        auto i = species.findWithName(si);
        auto j = species.findWithName(sj);
        if(i >= N || j >= N) continue;
        errorif(c.size() < 6, "ActivityModelElectrolyteNRTL: nrtl entry needs >=6 coeffs.");
        auto v = [&](int k){ return c.size() > k ? c[k].val() : 0.0; };
        Amm(i,j)=v(0); Amm(j,i)=v(1); Bmm(i,j)=v(2); Bmm(j,i)=v(3);
        almm(i,j)=almm(j,i)=v(4); dalmm(i,j)=dalmm(j,i)=v(5);
        Emm(i,j)=v(6); Emm(j,i)=v(7); Fmm(i,j)=v(8); Fmm(j,i)=v(9);
    }

    // --- Electrolyte-molecule pair params: (cation,anion,molecule) -> [C_caB,D_caB,E_caB,C_Bca,D_Bca,E_Bca,alpha] ---
    std::map<std::array<int,3>, std::array<double,7>> gm;
    for(auto const& [sc, sa, sB, c] : params.gmelc)
    {
        auto ic = species.findWithName(sc), ia = species.findWithName(sa), iB = species.findWithName(sB);
        if(ic >= N || ia >= N || iB >= N) continue;
        std::array<double,7> a{}; for(int k = 0; k < 7 && k < (int)c.size(); ++k) a[k] = c[k].val();
        gm[{(int)ic,(int)ia,(int)iB}] = a;
    }

    // --- Dielectric coefficients per solvent, and Born radii ---
    std::map<int, std::array<double,3>> diel;
    for(auto const& [s, c] : params.dielectric)
    {
        auto i = species.findWithName(s);
        if(i >= N) continue;
        std::array<double,3> a{0,0,298.15}; for(int k = 0; k < 3 && k < (int)c.size(); ++k) a[k] = c[k].val();
        diel[(int)i] = a;
    }
    ArrayXd rborn = ArrayXd::Constant(N, 3.0e-10);
    for(auto const& [s, r] : params.bornradius)
    {
        auto i = species.findWithName(s);
        if(i < N) rborn[i] = r.val();
    }

    // Helpers (capture the constant matrices/maps)
    auto tau_mm = [=](int around, int center, real const& T) -> real {
        return Amm(around,center) + Bmm(around,center)/T + Emm(around,center)*log(T) + Fmm(around,center)*T;
    };
    auto alpha_mm = [=](int i, int j, real const& T) -> real {
        return almm(i,j) + dalmm(i,j)*(T - 273.15);
    };
    auto has_ca = [=](int c, int a, int B){ return gm.find({c,a,B}) != gm.end(); };
    auto tau_caB = [=](int c, int a, int B, real const& T) -> real {
        auto it = gm.find({c,a,B}); if(it==gm.end()) return real(0.0);
        auto const& p = it->second; return p[0] + p[1]/T + p[2]*((TREF-T)/T + log(T/TREF));
    };
    auto tau_Bca = [=](int c, int a, int B, real const& T) -> real {
        auto it = gm.find({c,a,B}); if(it==gm.end()) return real(0.0);
        auto const& p = it->second; return p[3] + p[4]/T + p[5]*((TREF-T)/T + log(T/TREF));
    };
    auto alpha_ca = [=](int c, int a, int B) -> double {
        auto it = gm.find({c,a,B}); return it==gm.end() ? 0.2 : it->second[6];
    };
    auto epsSolvent = [=](int i, real const& T) -> real {
        if(i == (int)iw) return epsWater(T);
        auto it = diel.find(i); if(it==diel.end()) return epsWater(T);
        auto const& d = it->second; return d[0] + d[1]*(1.0/T - 1.0/d[2]);
    };

    auto aqstateptr = std::make_shared<AqueousMixtureState>();
    auto aqsolutionptr = std::make_shared<AqueousMixture>(solution);

    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        auto const& T = args.T;
        auto const& P = args.P;
        auto const& x = args.x;
        auto const RT = universalGasConstant*T;

        auto const& aqstate = *aqstateptr = solution.state(T, P, x);
        auto const& m = aqstate.m;
        auto const& ln_m = log(m);
        auto const& Ie = aqstate.Ie; // effective ionic strength (molal)

        props.som = StateOfMatter::Liquid;
        props.extra["AqueousMixtureState"] = aqstateptr;
        props.extra["AqueousMixture"] = aqsolutionptr;

        auto const xw = x[iw];
        auto const ln_xw = log(xw);

        // Effective mole fractions X_j = x_j * C_j (C=|z| ions, 1 molecules)
        ArrayXr X(N);
        for(auto i = 0; i < N; ++i) X[i] = (z[i] != 0.0 ? std::abs(z[i]) : 1.0) * x[i];

        // Local-composition ln gamma (symmetric), computed at a given composition array XX.
        // Returns ArrayXr of ln gamma^lc for all species.
        auto lc = [&](ArrayXr const& XX) -> ArrayXr
        {
            real sXc = 0.0, sXa = 0.0;
            for(auto c : icat) sXc += XX[c];
            for(auto a : ian) sXa += XX[a];
            ArrayXr out = ArrayXr::Zero(N);
            if(sXc <= 0.0 || sXa <= 0.0) return out; // no electrolyte -> lc contribution handled by mol-mol only below

            // aggregate G_cB, tau_cB (per cation c, molecule B) and G_aB, tau_aB
            auto G_cB = [&](int c, int B) -> real {
                real s = 0.0; for(auto a : ian) s += XX[a]*exp(-alpha_ca(c,a,B)*tau_caB(c,a,B,T)); return s/sXa; };
            auto a_cB = [&](int c, int B) -> real {
                real s = 0.0; for(auto a : ian) s += XX[a]*alpha_ca(c,a,B); return s/sXa; };
            auto tau_cB = [&](int c, int B) -> real { return -log(G_cB(c,B))/a_cB(c,B); };
            auto G_aB = [&](int a, int B) -> real {
                real s = 0.0; for(auto c : icat) s += XX[c]*exp(-alpha_ca(c,a,B)*tau_caB(c,a,B,T)); return s/sXc; };
            auto a_aB = [&](int a, int B) -> real {
                real s = 0.0; for(auto c : icat) s += XX[c]*alpha_ca(c,a,B); return s/sXc; };
            auto tau_aB = [&](int a, int B) -> real { return -log(G_aB(a,B))/a_aB(a,B); };

            // per-species neighbour G,tau in each cell type
            // molecule-B cell: neighbour k
            auto molcell = [&](int B, int k) -> std::pair<real,real> {
                if(z[k]==0.0){ real t=tau_mm(k,B,T); return {exp(-alpha_mm(k,B,T)*t), t}; }
                if(z[k]>0.0){ real g=G_cB(k,B); return {g, tau_cB(k,B)}; }
                real g=G_aB(k,B); return {g, tau_aB(k,B)};
            };

            // Assemble term1 (own-cell avg) + term2 (contribution to all cells) per species.
            // cells: molecule B (weight X_B); cation c ref anion a' (weight X_c*X_a'/sXa);
            //        anion a ref cation c' (weight X_a*X_c'/sXc).
            ArrayXr t1 = ArrayXr::Zero(N), t2 = ArrayXr::Zero(N);

            // --- molecule cells ---
            for(auto B : imol) { if(XX[B]<=0.0) continue;
                real D=0.0, Rn=0.0;
                std::vector<std::pair<real,real>> gt(N);
                for(auto k = 0; k < N; ++k){ if(XX[k]<=0.0){ gt[k]={0.0,0.0}; continue; } gt[k]=molcell(B,k); D+=XX[k]*gt[k].first; Rn+=XX[k]*gt[k].first*gt[k].second; }
                real R = D>0.0 ? Rn/D : real(0.0);
                t1[B] += R; // B's own cell is a molecule cell (handled below via central match)
                for(auto i = 0; i < N; ++i) if(XX[i]>0.0 && D>0.0) t2[i] += XX[B]*(gt[i].first/D)*(gt[i].second - R);
            }
            // fix term1 for molecules: it's the R of their own molecule cell (added above only if B==i). Undo the generic add:
            // (the loop above added R to t1[B] for each molecule cell B, i.e. exactly own-cell for molecules) -> OK.

            // --- cation cells (c, a') ---
            for(auto c : icat){ if(XX[c]<=0.0) continue;
                for(auto ap : ian){ if(XX[ap]<=0.0) continue;
                    real w = XX[c]*(XX[ap]/sXa);
                    real D=0.0, Rn=0.0; std::vector<std::pair<real,real>> gt(N);
                    for(auto k = 0; k < N; ++k){
                        if(XX[k]<=0.0 || z[k]>0.0){ gt[k]={0.0,0.0}; continue; } // like-ion (cations) excluded
                        real g, t;
                        if(z[k]==0.0){ t = tau_cB(c,k) - tau_caB(c,ap,k,T) + tau_Bca(c,ap,k,T); g = exp(-alpha_ca(c,ap,k)*t); }
                        else { t = 0.0; g = 1.0; } // anion neighbour: e-e default
                        gt[k]={g,t}; D+=XX[k]*g; Rn+=XX[k]*g*t;
                    }
                    real R = D>0.0 ? Rn/D : real(0.0);
                    t1[c] += (XX[ap]/sXa)*R;
                    for(auto i = 0; i < N; ++i) if(XX[i]>0.0 && z[i]<=0.0 && D>0.0) t2[i] += w*(gt[i].first/D)*(gt[i].second - R);
                }
            }
            // --- anion cells (a, c') ---
            for(auto a : ian){ if(XX[a]<=0.0) continue;
                for(auto cp : icat){ if(XX[cp]<=0.0) continue;
                    real w = XX[a]*(XX[cp]/sXc);
                    real D=0.0, Rn=0.0; std::vector<std::pair<real,real>> gt(N);
                    for(auto k = 0; k < N; ++k){
                        if(XX[k]<=0.0 || z[k]<0.0){ gt[k]={0.0,0.0}; continue; } // like-ion (anions) excluded
                        real g, t;
                        if(z[k]==0.0){ t = tau_aB(a,k) - tau_caB(cp,a,k,T) + tau_Bca(cp,a,k,T); g = exp(-alpha_ca(cp,a,k)*t); }
                        else { t = 0.0; g = 1.0; }
                        gt[k]={g,t}; D+=XX[k]*g; Rn+=XX[k]*g*t;
                    }
                    real R = D>0.0 ? Rn/D : real(0.0);
                    t1[a] += (XX[cp]/sXc)*R;
                    for(auto i = 0; i < N; ++i) if(XX[i]>0.0 && z[i]>=0.0 && D>0.0) t2[i] += w*(gt[i].first/D)*(gt[i].second - R);
                }
            }
            for(auto i = 0; i < N; ++i) out[i] = (z[i]!=0.0 ? std::abs(z[i]) : 1.0) * (t1[i] + t2[i]);
            return out;
        };

        // ln gamma^lc at current composition (symmetric)
        ArrayXr lcg = lc(X);

        // SOLVENTS = water + molecular solvents carrying CPDIEC dielectric data (i.e. AMP/PZ);
        // SOLUTES = Henry molecules (CO2), ions, carbamate. (Same test the Born term uses.)
        auto const isSolvent = [&](int i){ return (i == (int)iw) || diel.count(i) > 0; };

        // Unsymmetric reference for SOLUTES: ln gamma_i^inf at infinite dilution in the MIXED SOLVENT
        // (water + amines). Send only the SOLUTES to trace while KEEPING the solvent composition, so the
        // reference is the standard Aspen mixed-solvent infinite dilution — NOT the earlier "ratio-
        // preserving trace of all non-water", which also diluted the amine solvents and made lcinf
        // loading-dependent (so no constant standard-state shift could undo it).
        real solutesum = 0.0; for(auto i = 0; i < N; ++i) if(!isSolvent(i)) solutesum += x[i];
        ArrayXr lcinf = ArrayXr::Zero(N);
        {
            real s = solutesum > 0.0 ? real(1e-9)/solutesum : real(0.0);
            ArrayXr xref(N); real sumref = 0.0;
            for(auto i = 0; i < N; ++i){ xref[i] = isSolvent(i) ? x[i] : x[i]*s; sumref += xref[i]; }
            for(auto i = 0; i < N; ++i) xref[i] /= sumref;               // renormalize; solvent ratios preserved
            ArrayXr Xref(N);
            for(auto i = 0; i < N; ++i) Xref[i] = (z[i]!=0.0 ? std::abs(z[i]) : 1.0) * xref[i];
            lcinf = lc(Xref);
        }

        // --- Pitzer-Debye-Huckel (mole-fraction basis, eq 6) ---
        real Ix = 0.0; for(auto i = 0; i < N; ++i) Ix += x[i]*z[i]*z[i]; Ix *= 0.5;
        ArrayXr ln_gPDH = ArrayXr::Zero(N);
        if(Ix.val() > 0.0)
        {
            auto epsw = epsWater(T);
            real dsolv = 997.0; // kg/m3 (water-dominated)
            real Aphi = (1.0/3.0)*sqrt(2*3.14159265358979323846*NA*dsolv)*pow(QE*QE/(4*3.14159265358979323846*EPS0*epsw*KB*T), 1.5);
            real sI = sqrt(Ix);
            real pref = std::sqrt(1000.0/(Mw*1000.0)) * Aphi; // Mw*1000 = g/mol solvent MW
            for(auto i = 0; i < N; ++i){ real z2 = z[i]*z[i];
                ln_gPDH[i] = -pref*((2*z2/RHO)*log(1+RHO*sI) + (z2*sI - 2*pow(Ix,1.5))/(1+RHO*sI)); }
        }

        // --- Born term (eq 8), mixed-solvent dielectric ---
        // eps_mix = solvent-mole-fraction-weighted over {H2O + species with CPDIEC data}
        real epsmix, epsw2 = epsWater(T);
        {
            real num = 0.0, den = 0.0;
            for(auto i = 0; i < N; ++i){ bool solv = (i==(int)iw) || diel.count(i); if(solv && x[i]>0.0){ num += x[i]*epsSolvent(i,T); den += x[i]; } }
            epsmix = den>0.0 ? num/den : epsw2;
        }
        ArrayXr ln_gBorn = ArrayXr::Zero(N);
        {
            real bpref = (QE*QE/(2*4*3.14159265358979323846*EPS0*KB*T))*(1.0/epsmix - 1.0/epsw2)*1e-2;
            for(auto i = 0; i < N; ++i) if(z[i]!=0.0) ln_gBorn[i] = bpref*z[i]*z[i]/rborn[i];
        }

        // --- combine ---
        // SOLVENTS (water + AMP/PZ): symmetric lc + PDH + Born (Born=0 for neutral), MOLE-FRACTION scale.
        // SOLUTES (CO2 + ions + carbamate): unsymmetric lc (lcg - lcinf) + PDH + Born, MOLAL scale.
        ArrayXr ln_g(N);
        for(auto i = 0; i < N; ++i)
        {
            if(isSolvent(i)) ln_g[i] = lcg[i] + ln_gPDH[i] + ln_gBorn[i];
            else             ln_g[i] = (lcg[i] - lcinf[i]) + ln_gPDH[i] + ln_gBorn[i] + ln_xw;
        }
        props.ln_g = ln_g;

        props.ln_a = props.ln_g + ln_m;                                     // molal scale (solutes)
        for(auto i = 0; i < N; ++i)
            if(isSolvent(i)) props.ln_a[i] = props.ln_g[i] + log(x[i]);      // mole-fraction scale (solvents)
    };

    return fn;
}

auto ActivityModelElectrolyteNRTL() -> ActivityModelGenerator
{
    return ActivityModelElectrolyteNRTL(Params::embedded("ElectrolyteNRTL.v2024.yaml"));
}

auto ActivityModelElectrolyteNRTL(ActivityModelParamsElectrolyteNRTL const& params) -> ActivityModelGenerator
{
    return [=](SpeciesList const& species) { return createActivityModelElectrolyteNRTL(species, params); };
}

auto ActivityModelElectrolyteNRTL(Params const& params) -> ActivityModelGenerator
{
    auto const& data = params.data();
    errorif(!data.exists("ActivityModelParams"), "Expecting ElectrolyteNRTL activity model parameters with an `ActivityModelParams` section.");
    errorif(!data.at("ActivityModelParams").exists("ElectrolyteNRTL"), "Expecting an `ElectrolyteNRTL` subsection under `ActivityModelParams`.");
    ActivityModelParamsElectrolyteNRTL enrtlparams =
        data["ActivityModelParams"]["ElectrolyteNRTL"].as<ActivityModelParamsElectrolyteNRTL>();
    return ActivityModelElectrolyteNRTL(enrtlparams);
}

} // namespace Reaktoro
