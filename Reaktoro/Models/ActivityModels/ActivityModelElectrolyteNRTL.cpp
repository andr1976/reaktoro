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
#include <algorithm>
#include <vector>

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

    // Dense copy of the molecule-ion-pair parameters for the per-call precomputation in lc(): index
    // t = (ci*Na + ai)*Nm + mi over the positions of the cation, anion and molecule in icat/ian/imol.
    const int Nc = static_cast<int>(icat.size()), Na = static_cast<int>(ian.size()), Nm = static_cast<int>(imol.size());
    std::vector<int> pos(N, -1);
    for(int k = 0; k < Nc; ++k) pos[icat[k]] = k;
    for(int k = 0; k < Na; ++k) pos[ian[k]] = k;
    for(int k = 0; k < Nm; ++k) pos[imol[k]] = k;
    std::vector<std::array<double,7>> tpar(Nc*Na*Nm);
    std::vector<char> thas(Nc*Na*Nm, 0);
    for(auto const& [key, a] : gm)
    {
        const int t = (pos[key[0]]*Na + pos[key[1]])*Nm + pos[key[2]];
        if(pos[key[0]] < 0 || pos[key[1]] < 0 || pos[key[2]] < 0 || z[key[0]] <= 0.0 || z[key[1]] >= 0.0 || z[key[2]] != 0.0) continue;
        tpar[t] = a; thas[t] = 1;
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
    std::map<int, std::array<double,2>> vmol;   // solvent liquid molar volume [V0, dV/dT] (m3/mol)
    for(auto const& [s, c] : params.molarvolume)
    {
        auto i = species.findWithName(s);
        if(i >= N) continue;
        vmol[(int)i] = { c.size() > 0 ? c[0].val() : 0.0, c.size() > 1 ? c[1].val() : 0.0 };
    }
    ArrayXd mw(N); for(auto i = 0; i < N; ++i) mw[i] = species[i].molarMass();   // kg/mol
    ArrayXd rborn = ArrayXd::Constant(N, 3.0e-10);
    for(auto const& [s, r] : params.bornradius)
    {
        auto i = species.findWithName(s);
        if(i < N) rborn[i] = r.val();
    }

    // --- Optional Aspen mixed-solvent Henry rule (params.henry / params.henryderived; absent -> no groups) ---
    struct HenrySolvent { int A; std::array<double,4> lnH; double vc23; };
    struct HenryGroup { int i; int iW; std::vector<HenrySolvent> solv; std::vector<std::pair<int,double>> derived; };
    std::vector<HenryGroup> henry;
    for(auto const& [si, sa, c] : params.henry)
    {
        auto i = species.findWithName(si), a = species.findWithName(sa);
        if(i >= N || a >= N) continue;
        errorif(c.size() < 5, "ActivityModelElectrolyteNRTL: henry entry needs [A, B, C, D, Vc].");
        auto it = std::find_if(henry.begin(), henry.end(), [&](HenryGroup const& g){ return g.i == (int)i; });
        if(it == henry.end()) { henry.push_back(HenryGroup{(int)i, -1, {}, {}}); it = henry.end() - 1; }
        it->solv.push_back(HenrySolvent{(int)a, {c[0].val(), c[1].val(), c[2].val(), c[3].val()}, std::pow(c[4].val(), 2.0/3.0)});
        if(a == iw) it->iW = (int)a;
    }
    for(auto& g : henry)
    {
        errorif(g.iW < 0, "ActivityModelElectrolyteNRTL: the henry rule of ", species[g.i].name(), " needs water among its solvents.");
        for(auto const& [si, sk, nu] : params.henryderived)
        {
            auto i = species.findWithName(si), k = species.findWithName(sk);
            if((int)i == g.i && k < N) g.derived.emplace_back((int)k, nu.val());
        }
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
        // Temperature-dependent interaction parameters, computed once per call (they were recomputed, with map
        // lookups, inside the cell loops of lc for every neighbour). Same expressions as tau_caB/tau_Bca/alpha_ca
        // and tau_mm/alpha_mm; absent triples have tau = 0, alpha = 0.2.
        const int NT = Nc*Na*Nm;
        std::vector<real> tCA(NT), tBC(NT), Gca(NT);
        std::vector<double> alca(NT);
        {
            const real ft = (TREF - T)/T + log(T/TREF);
            for(int t = 0; t < NT; ++t)
            {
                if(thas[t]) {
                    auto const& q = tpar[t];
                    tCA[t] = q[0] + q[1]/T + q[2]*ft;
                    tBC[t] = q[3] + q[4]/T + q[5]*ft;
                    alca[t] = q[6];
                } else { tCA[t] = 0.0; tBC[t] = 0.0; alca[t] = 0.2; }
                Gca[t] = exp(-alca[t]*tCA[t]);
            }
        }
        std::vector<real> tMM(Nm*Nm), GMM(Nm*Nm);   // [k*Nm + B]: neighbour molecule k around central molecule B
        for(int kk = 0; kk < Nm; ++kk) for(int bb = 0; bb < Nm; ++bb) {
            const int k = imol[kk], B = imol[bb];
            tMM[kk*Nm + bb] = tau_mm(k, B, T);
            GMM[kk*Nm + bb] = exp(-alpha_mm(k, B, T)*tMM[kk*Nm + bb]);
        }
        auto tix = [&](int c, int a, int B) { return (pos[c]*Na + pos[a])*Nm + pos[B]; };

        // Local-composition ln gamma (symmetric), computed at a given composition array XX.
        std::vector<real> GcB(Nc*Nm), tcB(Nc*Nm), GaB(Na*Nm), taB(Na*Nm), gg(N), gtau(N);
        auto lc = [&](ArrayXr const& XX) -> ArrayXr
        {
            real sXc = 0.0, sXa = 0.0;
            for(auto c : icat) sXc += XX[c];
            for(auto a : ian) sXa += XX[a];
            ArrayXr out = ArrayXr::Zero(N);
            if(sXc <= 0.0 || sXa <= 0.0) return out; // no electrolyte -> lc contribution handled by mol-mol only below

            // aggregate G_cB, tau_cB (per cation c, molecule B) and G_aB, tau_aB
            for(int ci = 0; ci < Nc; ++ci) for(int mi = 0; mi < Nm; ++mi) {
                real sg = 0.0, sa = 0.0;
                for(int ai = 0; ai < Na; ++ai) { const int t = (ci*Na + ai)*Nm + mi; sg += XX[ian[ai]]*Gca[t]; sa += XX[ian[ai]]*alca[t]; }
                const real g = sg/sXa, al = sa/sXa;
                GcB[ci*Nm + mi] = g; tcB[ci*Nm + mi] = -log(g)/al;
            }
            for(int ai = 0; ai < Na; ++ai) for(int mi = 0; mi < Nm; ++mi) {
                real sg = 0.0, sa = 0.0;
                for(int ci = 0; ci < Nc; ++ci) { const int t = (ci*Na + ai)*Nm + mi; sg += XX[icat[ci]]*Gca[t]; sa += XX[icat[ci]]*alca[t]; }
                const real g = sg/sXc, al = sa/sXc;
                GaB[ai*Nm + mi] = g; taB[ai*Nm + mi] = -log(g)/al;
            }

            // Assemble term1 (own-cell avg) + term2 (contribution to all cells) per species.
            // cells: molecule B (weight X_B); cation c ref anion a' (weight X_c*X_a'/sXa);
            //        anion a ref cation c' (weight X_a*X_c'/sXc).
            ArrayXr t1 = ArrayXr::Zero(N), t2 = ArrayXr::Zero(N);

            // --- molecule cells: every species is a neighbour ---
            for(int bb = 0; bb < Nm; ++bb) { const int B = imol[bb]; if(XX[B]<=0.0) continue;
                real D=0.0, Rn=0.0;
                for(auto k = 0; k < N; ++k){
                    if(XX[k]<=0.0){ gg[k]=0.0; gtau[k]=0.0; continue; }
                    if(z[k]==0.0){ gg[k]=GMM[pos[k]*Nm + bb]; gtau[k]=tMM[pos[k]*Nm + bb]; }
                    else if(z[k]>0.0){ gg[k]=GcB[pos[k]*Nm + bb]; gtau[k]=tcB[pos[k]*Nm + bb]; }
                    else { gg[k]=GaB[pos[k]*Nm + bb]; gtau[k]=taB[pos[k]*Nm + bb]; }
                    D+=XX[k]*gg[k]; Rn+=XX[k]*gg[k]*gtau[k];
                }
                real R = D>0.0 ? Rn/D : real(0.0);
                t1[B] += R;
                for(auto i = 0; i < N; ++i) if(XX[i]>0.0 && D>0.0) t2[i] += XX[B]*(gg[i]/D)*(gtau[i] - R);
            }

            // --- cation cells (c, a'): molecules and anions are neighbours (like-ion repulsion) ---
            for(auto c : icat){ if(XX[c]<=0.0) continue;
                for(auto ap : ian){ if(XX[ap]<=0.0) continue;
                    real w = XX[c]*(XX[ap]/sXa);
                    real D=0.0, Rn=0.0;
                    for(auto k = 0; k < N; ++k){
                        if(XX[k]<=0.0 || z[k]>0.0){ gg[k]=0.0; gtau[k]=0.0; continue; }
                        if(z[k]==0.0){ const int t = tix(c,ap,k); gtau[k] = tcB[pos[c]*Nm + pos[k]] - tCA[t] + tBC[t]; gg[k] = exp(-alca[t]*gtau[k]); }
                        else { gtau[k] = 0.0; gg[k] = 1.0; } // anion neighbour: e-e default
                        D+=XX[k]*gg[k]; Rn+=XX[k]*gg[k]*gtau[k];
                    }
                    real R = D>0.0 ? Rn/D : real(0.0);
                    t1[c] += (XX[ap]/sXa)*R;
                    for(auto i = 0; i < N; ++i) if(XX[i]>0.0 && z[i]<=0.0 && D>0.0) t2[i] += w*(gg[i]/D)*(gtau[i] - R);
                }
            }
            // --- anion cells (a, c'): molecules and cations are neighbours ---
            for(auto a : ian){ if(XX[a]<=0.0) continue;
                for(auto cp : icat){ if(XX[cp]<=0.0) continue;
                    real w = XX[a]*(XX[cp]/sXc);
                    real D=0.0, Rn=0.0;
                    for(auto k = 0; k < N; ++k){
                        if(XX[k]<=0.0 || z[k]<0.0){ gg[k]=0.0; gtau[k]=0.0; continue; }
                        if(z[k]==0.0){ const int t = tix(cp,a,k); gtau[k] = taB[pos[a]*Nm + pos[k]] - tCA[t] + tBC[t]; gg[k] = exp(-alca[t]*gtau[k]); }
                        else { gtau[k] = 0.0; gg[k] = 1.0; }
                        D+=XX[k]*gg[k]; Rn+=XX[k]*gg[k]*gtau[k];
                    }
                    real R = D>0.0 ? Rn/D : real(0.0);
                    t1[a] += (XX[cp]/sXc)*R;
                    for(auto i = 0; i < N; ++i) if(XX[i]>0.0 && z[i]>=0.0 && D>0.0) t2[i] += w*(gg[i]/D)*(gtau[i] - R);
                }
            }
            for(auto i = 0; i < N; ++i) out[i] = (z[i]!=0.0 ? std::abs(z[i]) : 1.0) * (t1[i] + t2[i]);
            return out;
        };

        // ln gamma^lc at current composition (symmetric)
        ArrayXr lcg = lc(X);

        // Unsymmetric reference for ions & aqueous solutes: ln gamma_i^inf at ratio-preserving trace in water
        real nonw = 0.0; for(auto i = 0; i < N; ++i) if(i != (int)iw) nonw += x[i];
        ArrayXr lcinf = ArrayXr::Zero(N);
        if(nonw > 0.0)
        {
            real s = 1e-9/nonw;
            ArrayXr xref(N); real sumref = 0.0;
            for(auto i = 0; i < N; ++i){ xref[i] = (i==(int)iw) ? real(0.0) : x[i]*s; sumref += xref[i]; }
            xref[iw] = 1.0 - sumref;
            ArrayXr Xref(N);
            for(auto i = 0; i < N; ++i) Xref[i] = (z[i]!=0.0 ? std::abs(z[i]) : 1.0) * xref[i];
            lcinf = lc(Xref);
        }

        // --- Mixed-solvent properties over the solvents {H2O + species with CPDIEC data} ---
        // eps_s: mass-fraction-weighted dielectric constant (Aspen); M_s: mole-averaged molar mass;
        // d_s: mass density from ideal mixing of the solvent liquid molar volumes. Pure water gives
        // eps_w, Mw and 997 kg/m3 (the previous water-only values), so aqueous results are unchanged.
        const real epsw = epsWater(T);
        real epsmix = epsw, Msolv = Mw, dsolv = 997.0;
        {
            real nx = 0.0, mass = 0.0, epsm = 0.0, vol = 0.0;
            for(auto i = 0; i < N; ++i)
            {
                const bool solv = (i==(int)iw) || diel.count(i);
                if(!solv || !(x[i] > 0.0)) continue;
                const real mi = x[i]*mw[i];
                nx += x[i]; mass += mi; epsm += mi*epsSolvent(i,T);
                auto it = vmol.find(i);
                vol += (i==(int)iw || it==vmol.end()) ? mi/997.0 : x[i]*(it->second[0] + it->second[1]*(T - 298.15));
            }
            if(mass > 0.0) { epsmix = epsm/mass; Msolv = mass/nx; dsolv = mass/vol; }
        }

        // --- Pitzer-Debye-Huckel (mole-fraction basis, eq 6), mixed-solvent eps, density and molar mass ---
        real Ix = 0.0; for(auto i = 0; i < N; ++i) Ix += x[i]*z[i]*z[i]; Ix *= 0.5;
        ArrayXr ln_gPDH = ArrayXr::Zero(N);
        if(Ix.val() > 0.0)
        {
            real Aphi = (1.0/3.0)*sqrt(2*3.14159265358979323846*NA*dsolv)*pow(QE*QE/(4*3.14159265358979323846*EPS0*epsmix*KB*T), 1.5);
            real sI = sqrt(Ix);
            real pref = sqrt(1.0/Msolv) * Aphi;   // (1000/M_s[g/mol])^1/2 = (1/M_s[kg/mol])^1/2
            for(auto i = 0; i < N; ++i){ real z2 = z[i]*z[i];
                ln_gPDH[i] = -pref*((2*z2/RHO)*log(1+RHO*sI) + (z2*sI - 2*pow(Ix,1.5))/(1+RHO*sI)); }
        }

        // --- Born term (eq 8): Aspen's 1e-2 converts its cgs form with r in m; in SI (4·pi·eps0) it is absent ---
        ArrayXr ln_gBorn = ArrayXr::Zero(N);
        {
            real bpref = (QE*QE/(2*4*3.14159265358979323846*EPS0*KB*T))*(1.0/epsmix - 1.0/epsw);
            for(auto i = 0; i < N; ++i) if(z[i]!=0.0) ln_gBorn[i] = bpref*z[i]*z[i]/rborn[i];
        }

        // --- combine ---
        // water: symmetric lc + PDH; solutes: unsymmetric lc (lcg - lcinf) + PDH + Born, molal-scale (+ln_xw)
        ArrayXr ln_g(N);
        for(auto i = 0; i < N; ++i)
        {
            if(i == (int)iw) ln_g[i] = lcg[i] + ln_gPDH[i];
            else ln_g[i] = (lcg[i] - lcinf[i]) + ln_gPDH[i] + ln_gBorn[i] + ln_xw;
        }

        // --- optional mixed-solvent Henry rule (see ActivityModelParamsElectrolyteNRTL::henry) ---
        for(auto const& g : henry)
        {
            real sw = 0.0;
            for(auto const& s : g.solv) sw += x[s.A]*s.vc23;
            if(!(sw > 0.0)) continue;
            auto tau = [&](int a, int b) { return tau_mm(a, b, T); };
            auto lnH = [&](HenrySolvent const& s) { return s.lnH[0] + s.lnH[1]/T + s.lnH[2]*log(T) + s.lnH[3]*T; };
            real ssum = 0.0, refW = 0.0;
            for(auto const& s : g.solv)
            {
                const real tia = tau(g.i, s.A);
                const real lng_inf = tau(s.A, g.i) + tia*exp(-alpha_mm(g.i, s.A, T)*tia);   // binary NRTL, i in A
                ssum += x[s.A]*s.vc23/sw*(lnH(s) - lng_inf);
                if(s.A == g.iW) refW = lnH(s) - lng_inf;
            }
            // ln gamma_inf of i in the salt-free solvent mixture (multicomponent NRTL over {solvents, i}, x_i -> 0)
            const int ns = static_cast<int>(g.solv.size());
            std::vector<int> idx; for(auto const& s : g.solv) idx.push_back(s.A); idx.push_back(g.i);
            std::vector<real> xs(ns + 1, real(0.0)); real xt = 0.0;
            for(int k = 0; k < ns; ++k) { xs[k] = x[idx[k]]; xt += xs[k]; }
            for(int k = 0; k < ns; ++k) xs[k] /= xt;
            auto G = [&](int a, int b) { return exp(-alpha_mm(idx[a], idx[b], T)*tau(idx[a], idx[b])); };
            auto t = [&](int a, int b) { return tau(idx[a], idx[b]); };
            const int ii = ns;
            real num = 0.0, den = 0.0;
            for(int j = 0; j <= ns; ++j) { num += xs[j]*t(j, ii)*G(j, ii); den += xs[j]*G(j, ii); }
            real lnginf_mix = num/den;
            for(int j = 0; j <= ns; ++j)
            {
                if(!(xs[j] > 0.0)) continue;
                real dj = 0.0, nj = 0.0;
                for(int k = 0; k <= ns; ++k) { dj += xs[k]*G(k, j); nj += xs[k]*t(k, j)*G(k, j); }
                lnginf_mix += xs[j]*G(ii, j)/dj*(t(ii, j) - nj/dj);
            }
            HenrySolvent const* sW = nullptr; for(auto const& s : g.solv) if(s.A == g.iW) sW = &s;
            const real delta = ssum - refW;
            const real delta2 = ssum + lnginf_mix - lnH(*sW);
            ln_g[g.i] += delta;
            for(auto const& [k, nu] : g.derived) ln_g[k] += nu*delta2;
        }
        props.ln_g = ln_g;

        props.ln_a = props.ln_g + ln_m;
        props.ln_a[iw] = props.ln_g[iw] + ln_xw;
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
