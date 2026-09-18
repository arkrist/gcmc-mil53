# CO2 in MIL-53(Al) lp by rigid-framework GCMC (RASPA2): working notes

Each physical or technical choice is recorded here with its reason. Items marked
**OPEN** still need a decision from the project owner.

---

## Phase 1: installation and sanity check

### Installation
- `raspa2 2.0.50` from conda-forge, native **osx-arm64** build (Apple M2), in env
  `gcmc-mil53` (`environment.yml`). Package manager: `micromamba` (Homebrew). No
  source build was needed.
  ```
  brew install micromamba
  micromamba env create -f environment.yml
  export RASPA_DIR=~/micromamba/envs/gcmc-mil53   # simulate reads $RASPA_DIR/share/raspa/...
  ```
- The conda package ships `bin/simulate` and `share/raspa/{forcefield,molecules,structures}`
  but **no `examples/` directory**. The examples, together with their reference
  outputs, were taken from a shallow clone of the tag matching the binary:
  `git clone --depth 1 --branch v2.0.50 https://github.com/iRASPA/RASPA2 ~/src/RASPA2-2.0.50`
  (outside this repo; read only).

### What RASPA's "+/-" means
RASPA splits the production run into 5 blocks and prints `avg +/- 2.776 * s/sqrt(5)`
(`src/statistics.h`, `ERROR_CONFIDENCE_INTERVAL_95`): a **95 % confidence
half-width**, not a standard error. `parse_raspa_output.py` keeps it as `*_err95`
and also gives `*_sem = err95/2.776`. Every error bar in this project is the
95 % CI unless a figure says otherwise.

### Choice of sanity-check examples
v2.0.50 has no single-component CO2 or CH4 isotherm in IRMOF-1 with a shipped
reference output. Chosen instead (decision of 2026-09-18):

| example | conditions | what it tests |
|---|---|---|
| `Basic/7_Adsorption_of_Methane_in_MFI` | CH4, 300 K, 10 and 100 kPa, MFI 2x2x2 | swap/reinsertion GCMC machinery, conversion to mol/kg, an isotherm with 2 points |
| `Basic/8_Adsorption_of_CO2_in_CU-BTC` | CO2, 323 K, 1 MPa, Cu-BTC 1x1x1 | the chain MIL-53 needs (see below) |

**Pass criterion:** |L_ours - L_ref| < sqrt(e_ours^2 + e_ref^2), with L the absolute
loading in molecules/uc and e the 95 % errors. The two builds use different random
streams (2.0.45 clang-12 vs 2.0.50 conda-forge arm64), so only statistical
agreement is expected. If a point fails, nothing is rerun silently: the values,
the insertion acceptance and the equilibration length are reported first.

**What Cu-BTC validates for MIL-53:** a rigid 3-site charged CO2 (group move,
`RIGID_BOND`), framework partial charges read from the CIF
(`UseChargesFromCIFFile yes`), Ewald summation for framework-CO2 and CO2-CO2
electrostatics, and rotation moves. MIL-53(Al) relies on each of these. It does
not validate the PR fugacity (the example sets `FugacityCoefficient 1.0`, an
ideal gas), the TraPPE parameters, UFF, or tail corrections (its force field is
`ExampleMOFsForceField`: shifted LJ, no tail corrections, and a different CO2
model, q_C = +0.6512, C=O 1.149 A).

Caveat on the examples' **excess** loading: both inputs hard-code
`HeliumVoidFraction 0.29`. That is plausible for MFI, but for Cu-BTC it looks like
a copied placeholder. Our run uses the same value, so the comparison is
consistent, but the Cu-BTC excess number carries no physical meaning. The pass
test is therefore made on the absolute loading.

### Results
- MFI CH4: **PASS 2/2** (`results/sanity_mfi_ch4.{csv,png}`)

  | p [kPa] | ours [molec/uc] | reference [molec/uc] | delta | tolerance | acc. insertion ours / ref |
  |---|---|---|---|---|---|
  | 10  | 0.3485 +/- 0.0086 | 0.3482 +/- 0.0103 | +0.0003 | 0.0134 | 33.0 % / 33.2 % |
  | 100 | 2.863 +/- 0.042   | 2.898 +/- 0.033   | -0.034  | 0.053  | 29.9 % / 29.9 % |

  Wall time on the M2: 413 s for both points (reference: 691 s).
- Cu-BTC CO2: running (reference took 3.3 h).

---

## Phase 2: MIL-53(Al) lp setup

### Structure (**OPEN**: waiting for `structures/MIL-53_Al_lp.cif`, CoRE MOF DDEC set)
**The phase is identified by its cell, not by a CSD refcode** (decision of 2026-09-19).
The accepted phase is MIL-53(Al) lp (high-temperature, empty) from Loiseau et al.,
Chem. Eur. J. 2004, 10, 1373: orthorhombic Imma (No. 74), a = 6.608, b = 16.675,
c = 12.813 A, V ~ 1412 A^3, Al4C32H20O20 = 76 atoms, 832.4 g/mol per cell, pore
diameter 8.5 A, Langmuir area 1590 m^2/g.
**Rejected**: the as-synthesised Pnma cell (a = 17.129, b = 6.628, c = 12.182 A,
terephthalic acid in the pores) and the monoclinic Cc / P21/c narrow-pore hydrated forms.

`scripts/check_cif.py` applies this as its **first gate**: sorted axis lengths
within 1 % and angles within 0.5 deg of 90. On a mismatch it stops with exit status
1 before looking at anything else (tested on the lp cell with permuted axes: PASS;
on the as-synthesised, monoclinic and MIL-53(Cr) ht cells: FAIL). Then it
reports the cell, symmetry, atom count and composition,
charges (present or absent, net charge, per-element range), the unit-cell mass,
and whether every atom label has LJ parameters. Expected for a correct P1 lp
cell: 4 formula units Al(OH)(O2C-C6H4-CO2), i.e. Al4 O20 C32 H20 = **76 atoms**,
M_uc = 4 x 208.10 = **832.4 g/mol**, neutral.

### Force field (`forcefield/`)
| element | choice | reason |
|---|---|---|
| CO2 | TraPPE (Potoff & Siepmann, AIChE J. 2001, 47, 1676): rigid linear, C=O 1.16 A; C: eps/kB 27.0 K, sigma 2.80 A, q +0.70 e; O: 79.0 K, 3.05 A, -0.35 e | fitted to the CO2 vapour-liquid equilibrium, so the bulk fluid (the reservoir) is right; its point quadrupole makes the framework charges matter. Rigid: bending and stretching barely change adsorption at 303 K and would require CBMC growth. |
| framework LJ | UFF (Rappe et al., JACS 1992, 114, 10024), converted from (D, x) with eps = D and sigma = x / 2^(1/6): Al 254.126 K / 4.0082 A, O 30.193 / 3.1181, C 52.838 / 3.4309, H 22.142 / 2.5711 | generic and transferable, covers every element, and is the standard choice in MOF screening. Carboxylate O and hydroxyl O share one UFF LJ type (O_2/O_3/O_R all have D = 0.060, x = 3.500). |
| cross terms | Lorentz-Berthelot (arithmetic sigma, geometric eps) | standard with UFF + TraPPE. RASPA echoes the full pair table in the output header (check it there). |
| LJ cutoff | **12.0 A** (changed from 12.8 A on 2026-09-19, see Replication), **truncated, not shifted**, with **analytic tail corrections** | tail corrections assume g(r) = 1 beyond the cutoff and are only consistent with an unshifted potential. RASPA confirms: "TailCorrections are used / All potentials are unshifted". |
| electrostatics | Ewald, real-space cutoff 12.0 A, `EwaldPrecision 1e-6` | long-range CO2 quadrupole-framework interactions cannot be truncated. RASPA derives alpha (0.248 A^-1 in the smoke test, done at r_c = 12.8 A) and k-vectors from the precision. |
| framework charges | taken from the CIF (`UseChargesFromCIFFile yes`) | pending the CIF check. If they are absent: **stop and decide** (EQeq via RASPA is the fallback, see below). |

`.def` files are **line-positional** (RASPA skips a fixed number of header lines),
so explanations belong here, not in extra comment lines. An extra comment line in
`CO2.def` shifted Tc/Pc/omega and crashed the smoke test.

**Safety net.** For a pair with no LJ parameters, RASPA prints `WARNING: THERE ARE
ATOM-PAIRS WITH NO VDW INTERACTION` and **sets it to zero without stopping**.
`check_cif.py` fails if a label (digits stripped, as with
`RemoveAtomNumberCodeFromLabel yes`) is missing from the force field,
`run_isotherm.sh` aborts if the warning appears, and the parser flags it
(`missing_vdw_pairs`).

### Fugacity
RASPA2 has no `ComputeFugacityCoefficient` keyword. Peng-Robinson is the
hard-coded default EOS (`input.c`: `EquationOfState=PENG_ROBINSON`), and the
fugacity coefficient is computed whenever `FugacityCoefficient` is **absent**
(default -1). So the template omits it. PR uses Tc = 304.1282 K,
Pc = 7.3773 MPa, omega = 0.22394 from `CO2.def`. Check: phi(303 K, 1 MPa) = 0.949.

### Cutoff and replication
The minimum-image convention needs every **perpendicular** box width > 2 r_c. For a
triclinic cell, width_i = V / |a_j x a_k|, and n_i is the smallest integer with
n_i * width_i > 2 r_c (`make_inputs.py`, printed by `check_cif.py`).

With r_c = 12.8 A, two cells along c give 2 x 12.813 = 25.626 A against 25.6 A
required: a 0.03 A margin, rejected as too thin. **Decision: r_c = 12.0 A** for
both the LJ and the real-space Ewald parts, with LJ tail corrections. This is
standard for MOF GCMC; the tail correction accounts for the truncated LJ beyond
r_c, and Ewald accounts for the full electrostatics at any r_c. For the Imma
cell this gives:

| axis | cell [A] | n | box width [A] | margin over 24.0 A |
|---|---|---|---|---|
| a | 6.608  | 4 | 26.43 | 2.43 |
| b | 16.675 | 2 | 33.35 | 9.35 |
| c | 12.813 | 2 | 25.63 | 1.63 |

So the box is **4 x 2 x 2 = 16 cells**, 1216 framework atoms. It is recomputed
from the actual CIF, whose axes may be permuted.

### Moves and run length
Translation : rotation : reinsertion : swap = 0.5 : 0.5 : 0.5 : 1.0, so 20/20/20/40 %
of the moves. Swap (insertion/deletion, CBMC with 10 trial positions) is what
equilibrates N, so it gets the largest share. Reinsertion moves a molecule to a
random position, which helps it leave a pore. 10,000 initialization cycles,
50,000 production cycles. One RASPA cycle = max(20, N) moves.

### Excess loading and helium void fraction (method decided 2026-09-19; value pending the CIF)
Excess = absolute - rho_bulk(PR) * V_pore, with V_pore = theta_He x cell volume.
theta_He is **method-dependent**, so it is recorded with its parameters:

- method: RASPA Widom insertion of a helium probe, theta_He = <exp(-beta U_He-fw)>,
  T = **298 K**, rigid framework, 500,000 cycles (`runs/helium_void_fraction/`)
- probe: single LJ site, eps/kB = **10.9 K**, sigma = **2.64 A**, no charge
- framework: the same UFF LJ parameters as the CO2 runs, Lorentz-Berthelot, r_c = 12.0 A,
  truncated with tail corrections, no electrostatics
- **theta_He = (to fill in after the run) +/- (95 % CI)**

### Comparison with the literature (rule)
- **Absolute is compared with absolute and excess with excess.** Wherever a
  reference reports absolute loading (simulation papers usually do), the comparison
  uses our absolute loading. Experimental papers report excess, which is compared
  with our excess, with the caveat that it depends on our theta_He.
- Every reference point in `reference/*.csv` carries its **convention**
  (absolute/excess) and **temperature** (schema in `reference/README.md`).
- **We never convert a reference point between absolute and excess** unless the
  reference itself states the pore volume or void fraction it used. Without it,
  the point is plotted only against the matching column.

### Resume logic
`run_isotherm.sh` skips a point whose output contains "Simulation finished". An
unfinished point is rerun from scratch; its partial output is kept as
`Output.incomplete.<timestamp>`. Every result therefore comes from one
uninterrupted run.

---

## What a rigid-framework GCMC can and cannot reproduce for MIL-53

**What it gives.** The adsorption isotherm of CO2 in the **lp structure, held
fixed**. This is the lp branch: the loading the large-pore form would have at each
pressure if it stayed open. It tests the force field (UFF + DDEC charges + TraPPE)
for CO2-lp interactions, and it is directly comparable with other rigid-lp
simulations and with experimental data in the pressure range where the real
material is in the lp form.

**What it cannot give, by construction.** MIL-53(Al) breathes. Experimentally, CO2
adsorption drives the empty lp form to the narrow-pore (np) form at low pressure,
and the np form reopens to lp at higher pressure, which gives a stepped isotherm
(Bourrelly et al., JACS 2005, 127, 13519; check the step pressures there for the
temperature used). With the cell and atoms frozen in the lp geometry, the
simulation has no np state to go to. **The np-lp step cannot appear.** Any
disagreement in the pressure range where the real solid is np reflects the
missing structural transition, not the force field.

**What would be needed to capture the transition.**
1. **Osmotic ensemble** (N_host, mu_CO2, sigma, T). The stable phase at each
   pressure minimises the osmotic potential, which combines, for each rigid phase,
   the free-energy difference of the empty hosts, dF_host(np-lp) (from DFT or
   experiment), and the grand potential computed from that phase's rigid-GCMC
   isotherm. So one rigid np isotherm, the lp one from this work, and dF_host give
   a first estimate of the step pressures (Coudert et al., JACS 2008, 130, 14294).
2. **Hybrid GCMC/MD with a flexible framework.** GCMC insertion/deletion moves
   alternate with MD or cell-volume moves (NPT, or the full osmotic ensemble with
   volume moves). This needs a flexible MIL-53 force field that reproduces the
   np and lp minima and the barrier between them (e.g. the flexible MIL-53 models of
   the Maurin group; Ghoufi & Maurin, J. Phys. Chem. C 2010, 114, 6496). RASPA
   supports flexible frameworks and volume moves (see its
   `examples/Advanced/5_Adsorption_of_CO2_in_Flexible_IRMOF-1_Osmotic`), and so
   does LAMMPS (`fix gcmc` plus `fix npt`).

*Citations in this section are from memory; check them before the presentation.*
