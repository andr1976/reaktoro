# Reaktoro fork for the CAPE-OPEN amine packages (andr1976/reaktoro)

This fork carries the changes to Reaktoro v2.13 that the CAPE-OPEN property packages in
`andr1976/CapeOpenReaktoro` (ReaktoroShim.dll) depend on. It is **not** proposed upstream.

| Branch | Role |
|---|---|
| `capeopen-production` | **The only branch that may be built into `rkt-install`** (the Reaktoro.dll shipped with the shims). |
| `ExtendedNRTL` | History of the first eNRTL commits (3c5cefc7, c033e1cc); superseded by `capeopen-production`. |
| `fix-gas-ideal-entropy` | Sits on top of the symmetric-amine EXPERIMENT (850e60e9). **Never install from it**: on 2026-08-09 … 10-11 a build from this branch shipped and broke CESAR-1 CO2 VLE by ×44–671. |
| `symmetric-amine-lcinf`, `enrtl-mixed-solvent-fix` | Research branches (symmetric amine reference, re-fit harness). Build in a separate directory only. |

Base: `f587235e` (upstream main, "Minor changes in CondaAware.cmake"). All production changes are the seven
commits below.

---

## 1. Changes on `capeopen-production`

| # | Commit | Date | Area | Files |
|---|---|---|---|---|
| 1 | `3c5cefc7` | 2026-07-26 | New electrolyte-NRTL activity model + AMP/PZ parameter set | ActivityModelElectrolyteNRTL.{hpp,cpp,py.cxx,test.cxx}; ActivityModels.hpp/.py.cxx; Serialization/Models/ActivityModels.{hpp,cpp}; embedded/params + embedded/databases/reaktoro `ElectrolyteNRTL.v2024.yaml` |
| 2 | `c033e1cc` | 2026-07-27 | Database: neutral-amine Hr, CO2(aq) from the paper's Henry law | ElectrolyteNRTL.v2024.yaml (database) |
| 3 | `48196742` | 2026-08-09 | Cubic EOS: Gibbs/entropy of the gas phase | ActivityModelCubicEOS.cpp |
| 4 | `9ff15a45` | 2026-10-09 | Cubic EOS: opt-in root selection | Support/CubicEOS.{hpp,cpp}; ActivityModelCubicEOS.hpp |
| 5 | `5c996d5d` | 2026-10-12 | eNRTL: Born units, mixed-solvent PDH/Born, Aspen non-water defaults | ActivityModelElectrolyteNRTL.{hpp,cpp}; Serialization; params yaml |
| 6 | `8fa7631d` | 2026-10-12 | Database: PZ zwitterion on Aspen's pure-liquid basis | ElectrolyteNRTL.v2024.yaml (database) |
| 7 | `d7daad0e` | 2026-10-12 | eNRTL: per-call precomputation (5.3× faster, bit-identical) | ActivityModelElectrolyteNRTL.cpp |

### 1.1 `ActivityModelElectrolyteNRTL` (commits 1, 5, 7)

Aspen Plus ELECNRTL (Chen et al. 1982; Chen & Evans 1986; Mock et al. 1986), mixed-solvent form:

ln γ = local composition (NRTL) + Pitzer–Debye–Hückel (long range) + Born.

**Public API (new):**

- `struct ActivityModelParamsElectrolyteNRTL { nrtl, gmelc, gmelcee, dielectric, bornradius, molarvolume }`
- `ActivityModelElectrolyteNRTL()` — uses the embedded `ElectrolyteNRTL.v2024.yaml` params.
- `ActivityModelElectrolyteNRTL(ActivityModelParamsElectrolyteNRTL const&)`
- `ActivityModelElectrolyteNRTL(Params const&)` — reads `ActivityModelParams: ElectrolyteNRTL:`.
- Python bindings for all three; YAML encode/decode for the params struct.

**Parameter schema (YAML, `ActivityModelParams: ElectrolyteNRTL:`):**

| Key | Shape | Meaning |
|---|---|---|
| `nrtl: {i: {j: [a_ij, a_ji, b_ij, b_ji, α, dα, e_ij, e_ji, f_ij, f_ji]}}` | ≥ 6 values | τ_ij = a + b/T + e ln T + f T; α = α + dα (T − 273.15). Aspen NRTL-1 order differs (see §3). |
| `gmelc: {cation: {anion: {molecule: [C_caB, D_caB, E_caB, C_Bca, D_Bca, E_Bca, α]}}}` | 7 values | τ = C + D/T + E[(Tref − T)/T + ln(T/Tref)], Tref = 298.15 K (GMELCC/D/E/N). |
| `gmelcee` | parsed | **Parsed and serialised but not used**: ion-pair–ion-pair τ = 0, G = 1 is hard-coded. |
| `dielectric: {solvent: [A, B, C]}` | 3 values | ε(T) = A + B (1/T − 1/C) (CPDIEC). **A species listed here becomes a solvent** (mixed-solvent PDH/Born). **An `H2O` entry is ignored**: water always uses the built-in quadratic `epsWater`. |
| `bornradius: {ion: r}` | m | Default 3.0e-10 m. |
| `molarvolume: {solvent: [V0, dV/dT]}` | m³/mol | Mixed-solvent density for PDH; solvents without an entry count with water's mass density. |

**Behaviour that is fixed in code (not parameters):**

1. **Reference states.** Water: symmetric (pure liquid). **Every other species** (ions *and* neutral
   molecules, including amines and CO2): unsymmetric, ln γ* = ln γ^lc − ln γ^lc,∞, where the infinite-dilution
   value is evaluated at a ratio-preserving trace (1e-9) of all non-water species **in pure water**; ln γ is
   reported on the molal scale (+ ln x_w).
2. **Absent interaction parameters.** Missing `nrtl` pairs: τ = 0, α = 0.2. Missing `gmelc` triples:
   τ = 0, α = 0.2. Aspen instead uses water 8/−4 and non-water 10/−2; **write Aspen defaults explicitly** in
   the yaml (commit 5 did this for the 11 non-water triples of the AMP/PZ set).
3. **Mixed solvent (commit 5).** Solvents = H2O + species with a `dielectric` entry. ε_mix is
   mass-fraction weighted over solvents; PDH uses ε_mix, the mole-averaged solvent molar mass and the density
   from `molarvolume` (ideal mixing). Pure water reproduces the earlier water-only values (997 kg/m³, Mw).
4. **Born term (commit 5).** (e²/(2·4πε0·kT))·(1/ε_mix − 1/ε_w)·z²/r. The 1e-2 factor in Aspen's documented
   form is a cgs conversion and was removed; zero in pure water.
5. **Performance (commit 7).** τ/G tables for (cation, anion, molecule) and molecule–molecule pairs are
   evaluated once per call into dense arrays. ln γ is bit-identical to the previous implementation on 203
   states; property update 5.3× faster.

### 1.2 Database and parameters `ElectrolyteNRTL.v2024.yaml` (commits 1, 2, 5, 6)

AMP–PZ–CO2–H2O transcribed from Yi et al. (2024, Sep. Purif. Technol. 350, 127924), open-access Aspen input.

- Species: H2O, CO2, PZ, AMP, PZCOOH (zwitterion HPZCOO), PZH+, PZH2+2, AMPH+, PZCOO-, PZ(COO)2-2, H3O+, OH-,
  HCO3-, CO3-2 + gases CO2(g), H2O(g). Standard states from DGAQFM/DHAQFM (ions), DGFORM + Antoine (amines).
- Commit 2: neutral-amine Hr re-derived (PZ −74.6, AMP −343.0 kJ/mol); CO2(aq) Gr/Hr/Cp refit so the
  CO2(aq) ↔ CO2(g) transfer reproduces the paper's Henry correlation (25–120 °C).
- Commit 5: 11 non-water gmelc triples absent from Yi's file set to Aspen's non-water default (C = −2/10);
  `dielectric` for AMP and PZ; `molarvolume` for AMP and PZ.
- Commit 6: PZCOOH Gr −278.596 → −279.0447 kJ/mol (= + RT0 ln γ∞, ln γ∞ = −0.181) because Aspen treats the
  zwitterion on a pure-liquid (symmetric) basis and this model on the aqueous (unsymmetric) one.

Result (2026-10-12): ln γ vs Aspen ELECNRTL within 0.015 (lean) / 0.04 (rich) per z²; CO2 exact; Aspen's lean
and rich streams re-speciate within ±8 % per species.

### 1.3 Cubic EOS (commits 3, 4) — affects every package (AMP/PZ, HPC, biogas)

- **Commit 3 — gas entropy/Gibbs.** `ActivityModelCubicEOS` set `Gx = Gres` (residual only). It now sets
  `Gx = Gres + RT(Σ x ln x + ln P/bar)` so that G = Σ n μ and dS/d ln P = −R. Chemical potentials and
  equilibria are unchanged; P-H/P-S flashes, compressor and expander work were wrong before (isentropic
  compression ratio 3: +20 K → +119 K).
- **Commit 4 — root selection.** `CubicEOS::EquationModel::rootstate` (default `Unspecified` = upstream
  behaviour, minimum residual Gibbs). `Gas` forces Zmax, `Liquid` forces Zmin when there are three real roots.
  `ActivityModelCubicEOS(bipgen, eqmodel)` is exported so a custom equation model can be used. The shims set
  `Gas` for the vapour phase (supersaturated water vapour otherwise took the liquid root, −46 kJ/mol in H).

---

## 2. What must not break — regression tests

Reference values are the production results of 2026-10-14 (this branch at d7daad0e, the deployed shim). They
reproduce the techref validation; the July 2026 numbers (p_CO2 0.0921/1.0058/12.13 kPa at loading
0.30/0.519/0.70, 40 °C) predate commits 5–6 and are history (current: 0.1250/1.0552/11.25 kPa).

| Suite | Where | What it checks |
|---|---|---|
| Catch2 `[capeopen]` | `ActivityModelElectrolyteNRTL.test.cxx` (generated), `ActivityModelCubicEOS.test.cxx` (new test case) | ln γ golden values at 10 states (Catch2 original state, Aspen S1/S2-LEAN/S3-RICH, 6 random CESAR-1), tolerance 1e-9; pure-water and infinite-dilution limits; absent gmelc triple ≡ code default; H2O dielectric entry ignored; params YAML round trip; standard Gibbs energies (PZCOOH −279.0447 kJ/mol); cubic EOS Gx = RT Σx ln a, dS/d ln P = −R, rootstate Gas/Liquid/default |
| Catch2 upstream `[ActivityModelCubicEOS]` | existing upstream test case | the EOS changes keep upstream results |
| pytest (CapeOpenReaktoro `tests/`) | `test_fork_enrtl.py`, `test_fork_cubic.py` | the same ln γ golden values through reaktoro4py, plus accuracy against Aspen ELECNRTL (ions ≤ 0.05 z², CO2 and water ≤ 0.01) and mixed-solvent sensitivity |
| pytest (CapeOpenReaktoro `tests/`) | `test_shim_cesar1.py` | the deployed ReaktoroShim.dll: p_CO2 (27 points, 40–120 °C, loading 0.05–0.7), amine fugacities, water fugacity (IAPWS check), liquid and gas enthalpy, water-vapour root (latent heat), density/viscosity, N2/O2 Henry, TVF/PVF/TP flashes |

Golden values are generated by `CapeOpenReaktoro/tests/make_golden.py` (it also writes the Catch2 file with
`--fork`). Regenerate them only after a deliberate, validated change, and commit golden files and code together.
A 0.01 change in one gmelc τ fails 35 of the 46 shim tests.

**How to run**

- Locally (test executables are blocked by WDAC): `tests\catch2uild_catch2_dll.bat` then
  `python tests\catch2un_catch2.py "[capeopen]"` (Catch2 inside a DLL, linked to `rkt-install`), and
  `python -m pytest tests` in the conda `reaktoro` env. Both use the deployed Reaktoro.dll.
- CI: `.github/workflows/capeopen-production.yml` builds this branch with MSVC and the production pins and runs
  `reaktoro-cpptests "[capeopen],[ActivityModelCubicEOS]"` on every push to `capeopen-production`.

## 3. Guidance for adding the MEA package (or any new eNRTL system)

1. **New files, never edit v2024.** Add `ElectrolyteNRTL-MEA.<year>.yaml` (params + database) and load it
   explicitly; `ActivityModelElectrolyteNRTL()` with no arguments must keep returning the AMP/PZ set.
2. **Anything that changes ln γ for an existing parameter set must be opt-in** through a new parameter key
   with a default equal to today's behaviour. The symmetric-amine EXPERIMENT (850e60e9) was unconditional
   and broke CESAR-1. Candidates that MEA needs:
   - MEA on the **pure-liquid (symmetric) reference** (Aspen: MEA is a solvent; pKa 9.61 vs 10.20 otherwise):
     e.g. a `solvents:` / `symmetric:` list, default = H2O only.
   - **Aspen mixed-solvent Henry rule for CO2** (ln(H_i/γ_i∞) = Σ_A w_A ln(H_iA/γ_iA∞), w_A ∝ x_A V_cA^(2/3)):
     new optional key (e.g. `henry:` with V_c per solvent); absent → current CO2 treatment.
   - **Water dielectric from CPDIEC** (Aspen 78.54/31989.38/298.15) instead of `epsWater`: only if an `H2O`
     `dielectric` entry is present *and* an explicit switch is set — today an `H2O` entry is silently ignored,
     and the AMP/PZ file must keep the quadratic.
3. **Write Aspen defaults explicitly**: the code default for a missing triple is τ = 0, α = 0.2. In the yaml
   order `[C_caB, D_caB, E_caB, C_Bca, D_Bca, E_Bca, α]` (ion pair → molecule first) Aspen's water default is
   `[-4, 0, 0, 8, 0, 0, 0.2]` and its non-water default `[-2, 0, 0, 10, 0, 0, 0.2]`.
4. **Parameter order differs from Aspen NRTL-1**: yaml `[a_ij, a_ji, b_ij, b_ji, α, dα, e_ij, e_ji, f_ij,
   f_ji]`; Aspen BPVAL lists `a_ij b_ij c_ij d_ij e_ij f_ij` per ordered pair.
5. A species with a `dielectric` entry becomes a solvent in PDH/Born: give MEA a `molarvolume` too.
6. Before installing any rebuild: branch check (`git branch --show-current` = `capeopen-production`), then the
   §2 checks. Build MEA research work in its own build directory.

---

## 4. Build, install and verify

- Build dir `C:\Users\AndersAndreasen\rb2` (VS 2022 x64, conda env `reaktoro` for dependencies, Eigen 3.4,
  `/LIBPATH` for yaml-cpp), install prefix `C:\Users\AndersAndreasen\rkt-install`; recipe in the
  CapeOpenReaktoro memory/HANDOVER notes.
- `rb2` also builds `reaktoro4py` (`rb2\Reaktoro\Release\reaktoro4py.cp312-win_amd64.pyd`); with
  `os.add_dll_directory` on `rb2\Reaktoro\Release` and the conda `Library\bin` it runs in the conda Python —
  the only way to execute the model directly on this machine.
- The upstream workflows (`linux.yml`, `osx.yml`, `windows.yml`) trigger on `main` and pull requests;
  `capeopen-production.yml` runs the CAPE-OPEN tests for this branch (§2).
- Note: `rb2\python\packageuild\libeaktoro` holds a stale copy of Reaktoro.dll (2026-10-11); the tests load
  the pyd from `rb2\Reaktoro\Release` against the deployed DLL instead.
