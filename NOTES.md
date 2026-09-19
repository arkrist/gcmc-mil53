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

### Structure: DECIDED 2026-09-19. `structures/MIL-53_Al_lp.cif` (SABVUN, CoRE MOF 2014 DDEC)
**Accepted:** `SABVUN_clean.cif` from the CoRE MOF 2014 DDEC set (Nazarian, Camp &
Sholl, Chem. Mater. 2016, 28, 785; Zenodo 3986573). The lattice decides: the implied
conventional cell 6.6085 / 12.8130 / 16.6750 A, V = 1412.0 A^3, is identical to
Loiseau's lp values. DDEC charges are in the same file as the coordinates, the
channels are empty and the mu-OH hydrogens are present.

**Refcode correction.** The refcode first given for the lp form, SABWAU, was
**wrong**. The CoRE 2024 index shows SABWAU01 is a monoclinic P21/c **narrow-pore**
form (LCD 2.8 A), while **SABVUN** is the lp form: its lattice matches, and the
CoRE 2024 index lists it under the Loiseau 2004 DOI, 10.1002/chem.200305413. The
`_srcid 'CAKYAQ_clean'` in the file's metadata line is a **database artefact**. Cell,
composition and DOI govern.

**Files and checksums**
- `structures/original/SABVUN_clean.cif`: the downloaded original, **untouched**
  (read-only). md5 `7b79c613f3c5da8c31a051df553575ff`, sha256
  `9cc2086f8654414508382fedd2e77874a636d63e082c4d5376186a5cd0d6b8b8`. It comes from
  `core-mof-1.0-ddec.tar` (md5 `a2e3578003968739640b6faeed361299`, identical to
  the md5 published by Zenodo), member `./re-labeled/SABVUN_clean.cif`.
- `structures/SABVUN_clean_primitive.cif`: a copy whose **only** change is that the
  Python-dict metadata line (not CIF) is commented out (diff-verified).
- `structures/MIL-53_Al_lp.cif`: the **file used**. An exact change of basis to the
  conventional orthorhombic cell (`scripts/primitive_to_conventional.py`), chosen
  because it matches the convention of every paper we compare against, gives the
  4 x 2 x 2 replication the cutoff was sized for, and is cheaper (16 conventional
  cells instead of 22.5).

**Change of basis and its verification** (run on the re-read output file; the
agreed rule was to fall back to the primitive cell on any failure, never to patch):
C = M P with M = [[-1,0,0],[0,-1,1],[1,1,1]] (rows = conventional a, b, c).

| check | result |
|---|---|
| det(M) = 2 | 2.000000 PASS |
| 76 atoms, no losses, no duplicates | 76 = 2 x 38, each primitive atom exactly twice; shortest distance 0.8615 A (the O-H) PASS |
| composition, molar mass | Al4C32H20O20, 832.415 g/mol PASS |
| net charge (tolerance 1e-4 e) | +2e-6 e; every atom keeps its source charge PASS |
| crystal density | 0.978966 g/cm^3 for both files PASS |
| nearest-neighbour distances | all 76 per-atom neighbour lists within 5 A identical to the source atom's, max delta 2e-8 A; NN histogram identical PASS |
| resulting cell | 6.6085 / 16.6750 / 12.8130 A; angles within 1.7e-4 deg of 90, inherited from the rounding in the primitive file (b = 11.02159403 vs c = 11.02160000), kept as is |

`check_cif.py structures/MIL-53_Al_lp.cif`: CELL PASS (0.01/0.00/0.00 %),
COMPOSITION PASS, charges present (net +2e-6 e; Al +1.847, O -0.65 (carboxylate)
/ -1.12 (hydroxyl), H +0.11 (aromatic) / +0.477 (hydroxyl)), FF CHECK PASS,
replication 4 x 2 x 2. RASPA reads 1216 framework atoms and 13318.65 g/mol
(= 16 x 832.415), with no missing-VDW warning.

**Known systematic: mu-OH geometry and parameters.** O-H = 0.86 A as published: the
short X-ray-like distance, where a real O-H bond is about 0.97 A. It is **kept
deliberately**, because the DDEC charges were computed on this geometry and
consistency between charges and geometry matters more than a realistic bond length.
Why it matters *here*: Bourrelly 2005 concludes that CO2 first sorbs onto the
**hydroxyl groups**, at **0.75 CO2 per structural OH** before the step, so the mu-OH
is the key adsorption site in this material. A short O-H (H closer to O, so a
smaller effective OH dipole and a different H position) and generic UFF parameters
for H (eps 22.1 K, sigma 2.57 A) act on exactly the interaction that matters most.
First quantification (Widom, infinite dilution, 500 cycles): switching the
framework charges off lowers K_H only from 1.91e-4 to 1.66e-4 mol/kg/Pa (-15 %) and
<U_gh> from -22.9 to -22.0 kJ/mol. At infinite dilution the binding is dominated by
UFF dispersion in the channel, not by the OH electrostatics (see the Henry section).


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

### Structure provenance (search of 2026-09-19)
Search order as instructed: CoRE MOF (GitHub mirrors, then the Zenodo DDEC record),
QMOF, simulation SI, Ghoufi's repository. No CSD access. Every CIF was screened by
cell against the Loiseau lp cell (sorted axes, 1 %). Files in a primitive/reduced
setting were compared through their implied conventional cell (lattice-vector
search, analysis only, file untouched). 23,845 CoRE 2014/2019 CIFs + 2,932 DDEC
CIFs + 2,665 CoRE 2024 SI CIFs were scanned.

| # | file | source (exact) | SG in file | cell in file [A, deg] | atoms / composition / M | charges | channels | verdict |
|---|---|---|---|---|---|---|---|---|
| **1** | `re-labeled/SABVUN_clean.cif` (sha256 `9cc2086f...d0d6b8b8`) | **CoRE MOF 2014 DDEC Database**, Nazarian, Camp & Sholl, Chem. Mater. 2016, 28, 785; Zenodo record 3986573 (DOI 10.5281/zenodo.3986573, published 2016-01-07), `https://zenodo.org/records/3986573/files/core-mof-1.0-ddec.tar`, md5 `a2e3578003968739640b6faeed361299` (matches Zenodo) | P 1, **reduced primitive setting** of the I-centred lattice | 6.6085 11.0216 11.0216 / 98.308 107.445 107.445, V 705.98. **Implied conventional: 6.6085 12.8130 16.6750, 90/90/90, V 1412.0 = Loiseau lp exactly** | 38 = Al2C16H10O10 per primitive cell (x2 = Al4C32H20O20, 76 atoms, 832.4 g/mol) | **present, DDEC** (file header: "VASP DFT with PBE;DDEC"), net +0.000001 e; Al +1.847, O(H) -1.12, H(O) +0.477 | empty: all 38 atoms form one bonded network | **lp. Best candidate.** The CoRE 2024 index gives SABVUN the DOI 10.1002/chem.200305413 (Loiseau 2004) |
| 2 | `re-labeled/WAYMEQ_clean.cif` | same DDEC tar | P 1, reduced primitive | 6.612 11.034 11.034 / 70.4 72.6 72.6; implied conv. 6.6125 12.7143 16.7807 (<= 0.8 %) | **36** = Al2C16H8O10: **mu-OH H missing** | DDEC | empty | lp lattice, **rejected** (incomplete) |
| 3 | `WAYMEQ_SL.cif` | CoRE MOF 2019 ASR (public, v1.1.4 per package README), GitHub `coudertlab/CoRE-MOF` commit `a931fb5809696bbd43fbdec54341ee95b1d323db`, `src/CoRE_MOF/data/2019-ASR.tar.xz` | P1 | 16.7807 6.6125 12.7143 / 90 90 90, V 1410.8 | 76, Al4C32H20O20, 832.415 | **absent** | empty | lp, **not usable** (no charges); H added by CoRE curators (header `WAYMEQ_clean_h`, Materials Studio 2015) |
| 4 | `cm503311x_aloh300K_clean.cif` (+ 150-500 K series) | same CoRE 2019 ASR/FSR; SI of Chem. Mater., DOI prefix 10.1021/cm503311x | P 1 | 6.6295 16.7579 12.7936 / 90 90 90 | **56 = Al4C32O20: no H at all** | absent | - | lp cell, **rejected** (no H, no charges) |
| 5 | `CR/2019[Al][bpq]3[ASR]3.cif` (refcode `S0885715619000460sup002`) | **CoRE MOF 2024**, GitHub `Chung-Research-Group/CoRE-MOF-Tools` commit `3f9fff76136b83135195013fa88d313478b23118`, `CoREMOF/data/SI/CR.zip`; structure from SI of Powder Diffraction, DOI 10.1017/S0885715619000460 | P 1, primitive | 11.0559 x3 / 145.09 109.00 81.62; implied conv. 6.632 12.840 16.736, V 1425 (<= 0.4 %) | 38 = Al2C16H10O10 | **PACMAN v1.1** (ML surrogate of DDEC6, not DDEC) | - | lp, not Loiseau; **backup only** |
| - | SABVOH, SABVOH01 | CoRE 2014 / 2019 / DDEC | Pnma or P1 | 17.129 6.628 12.182 | 76 | - | - | **as-synthesised, rejected** |
| - | SABWAU01 | CoRE 2024 index only (CSD-derived CIFs are not on GitHub) | P21/c | LCD 2.8 A, VF 0.40 | 152 | - | - | **narrow-pore, rejected** |
| - | EGELUY, EGELUY01 ("MIL-53ht"), EHALOP, QONQEQ, WAYMIU, WAYMOA, QONQAM | CoRE 2014/2019/DDEC | P1 / reduced | conventional deviations 1.6-8.6 % | 76 or 38 | - | - | not the Loiseau lp cell, rejected |

Other sources checked:
- **QMOF** (`Andrew-S-Rosen/QMOF`): the GitHub repo holds tools only; the structures are
  on Figshare. Not searched, since candidate 1 already exists.
- **Simulation-paper SI**: not needed after candidate 1.
- **A. Ghoufi, `aghoufi/DUT-49-Cu--OMD`** (single commit `cd684dfebb`, 2022-11-23):
  two files, `FIELD-DUT49Cu` and `FIELD-MIL53Cr`, both **DL_POLY FIELD files, with no
  coordinates** (no CONFIG, no CIF). `FIELD-MIL53Cr` is titled "Mil53_Cr_Lp": 2432
  framework atoms (= 32 cells x 76), fully **flexible** (2304 bonds, 3584 angles,
  1536 dihedrals), charges Cr +1.418, O(H) -0.730, H(O) +0.299, O(carboxylate) -0.566.
  The Cr LJ parameters (eps 1.297 kJ/mol = 156.0 K, sigma 3.911 A) equal DREIDING's
  **Al** values. Guests: CO2 (q_C +0.6512, C=O 1.162 A, harmonic bend) and N2.
  **No Al structure and no lp geometry usable for us** (Cr, and no coordinates).
  **Worth knowing before meeting A. Ghoufi:** his MIL-53(Cr) force field uses the
  **DREIDING Al** Lennard-Jones parameters for Cr, and a CO2 model with
  **q_C = +0.6512**, the same charge as RASPA's built-in `ExampleDefinitions/CO2.def`
  (Garcia-Sanchez type), **not TraPPE** (q_C = +0.70), which is what we use.

**Caveats on candidate 1**, all RESOLVED 2026-09-19 (see the Structure section above): (1) option (b) chosen, (2) header commented out in a copy, (3) O-H kept and recorded as a systematic:
1. **Setting.** Primitive 38-atom cell. Options: (a) use it **as-is** (RASPA
   handles triclinic cells; replication 5 x 3 x 3 = 45 primitive = 22.5
   conventional cells, about 1.4x the cost), or (b) an **exact change of basis** to the
   conventional orthorhombic cell (a = p2+p3, etc.), mapping the same atoms, the same
   charges and all interatomic distances unchanged, 76 atoms, replication 4 x 2 x 2.
   Neither changes the physics.
2. The **first line is not CIF** but a Python-dict metadata line, and its `_srcid`
   reads **`CAKYAQ_clean`**, not SABVUN (apparently a metadata artefact of the
   database; the cell, composition and DOI all identify SABVUN). It must be
   commented out for RASPA.
3. **mu-OH geometry**: O-H = 0.86 A (the X-ray-like short distance, vs about 0.97 A
   for a real O-H). This is the experimental geometry, on which the DDEC charges were
   computed. Kept as is and noted.

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

### Excess loading and helium void fraction (DONE 2026-09-19)
Excess = absolute - rho_bulk(PR) * V_pore, with V_pore = theta_He x cell volume.
theta_He is **method-dependent**, so it is recorded with its parameters:

- method: RASPA Widom insertion of a helium probe, theta_He = <exp(-beta U_He-fw)>,
  T = **298 K**, rigid framework, 500,000 cycles (`runs/helium_void_fraction/`)
- probe: single LJ site, eps/kB = **10.9 K**, sigma = **2.64 A**, no charge
- framework: the same UFF LJ parameters as the CO2 runs, Lorentz-Berthelot, r_c = 12.0 A,
  truncated with tail corrections, no electrostatics
- structure: `structures/MIL-53_Al_lp.cif`, 4 x 2 x 2 cells, `runs/helium_void_fraction/` (2791 s)
- **theta_He = 0.7115 +/- 0.0004** (95 % CI; blocks 0.7111-0.7118)
- The corresponding pore volume is theta_He x V_uc / m_uc = 0.7115 x 1411.96 A^3 / 832.415 g/mol
  = **0.727 cm^3/g**.
- This is **not** the geometric void fraction. CoRE 2024 lists AV_VF = 0.587 for
  SABVUN (probe-accessible volume). The Widom value <exp(-beta U)> gives weights > 1
  to weakly attractive regions near the walls, so it is larger. This definitional
  gap is why the method and probe are stated whenever the value is quoted.

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
- **Pressure window for Bourrelly 2005: P >= 9 bar only** (the fully open lp region,
  see the last section). Points at 5-9 bar are in the breathing step and are shown
  greyed out, labelled "not compared". The figure caption states the restriction.
- Temperature: reference at 304 K, simulation at 303 K. The 1 K offset is noted and
  not corrected.
- **Convention of the Bourrelly points (corrected 2026-09-19): EXCESS, assumed.**
  Bourrelly 2005 reports n^a from **manometry**, which measures excess, and states
  no conversion to absolute amounts. The CSV says: "excess (assumed; source reports
  n^a from manometry, convention not stated)". Therefore:
  - **The primary comparison is our EXCESS isotherm** vs the points (P >= 9 bar).
  - Our **absolute** isotherm is drawn as a **shaded band** (from excess to absolute)
    so the systematic between the two conventions stays visible. With our theta_He
    (V_pore = 0.727 cm^3/g) and rho_bulk ~ 67 kg/m^3 (CO2, 303 K, 30 bar), it is about
    **1.1 mmol/g (~11 %) at 30 bar**. The earlier 0.8 mmol/g / 8 % figure assumed
    V_pore ~ 0.55 cm^3/g and is superseded. It uses **our** theta_He, never a number
    from the source; RASPA's exact PR-based value is in the results CSV.
  - The reference points are **never converted**.
  - Caption and NOTES state this.
- **Units:** every loading is reported in **mol/kg and molecules per conventional
  unit cell** (Al4(OH)4(bdc)4, 832.4 g/mol; 1 molecule/uc = 1.2013 mol/kg), both in
  the CSV and on the figure's second y-axis. The reference CSV carries both columns,
  consistent with 832.4 g/mol (5.30 mmol/g <-> 4.41/uc). Order-of-magnitude sanity
  check (different metal, not a validation): Ghoufi & Maurin 2010, MIL-53(Cr), report
  3.0 CO2/uc at 4.7 bar in the np form and a rigid-lp isotherm reaching roughly 7-8
  CO2/uc by 15 bar. If our Al lp result came out at ~2 or ~20 CO2/uc, something
  would be wrong.

### Low-pressure test: Henry constant vs Coudert's K_lp (Phase 3; approved 2026-09-19)
Target: K_lp ~ 2.6e-5 mol kg^-1 Pa^-1 (Coudert 2008, Langmuir fit, no error bar).
Three independent numbers are reported side by side:
1. **Weighted linear fit through zero**, n_abs = K_H p, over the lowest pressures that
   are still linear (weights 1/sem^2, uncertainty from the covariance). If the 0.1 bar
   point already shows curvature, **ask before adding lower pressures**.
2. **Langmuir fit** of the whole computed isotherm, K = q_sat b (the same functional
   form Coudert used), uncertainty from the covariance.
3. **Direct Widom test-particle insertion** of CO2 in the empty rigid framework
   (`runs/henry_widom_CO2`, 20,000 cycles, same force field, charges, Ewald and T).

**Caution on comparing them.** (1) and (3) estimate the true zero-loading Henry
constant. (2) and Coudert's K_lp are Langmuir parameters, i.e. an effective slope of
a fit over the whole isotherm. On an energetically heterogeneous surface the true
Henry constant exceeds the Langmuir K. So the like-for-like comparison with Coudert
is (2), and (1)/(3) vs (2) measures how non-Langmuir our isotherm is.

**Result (3), final, `runs/henry_widom_CO2`, 20,000 cycles, 303 K: K_H(Widom) =
1.933e-4 +/- 0.007e-4 mol/kg/Pa** (95 % CI; blocks 1.924-1.937e-4),
<U_gh> - <U_h> = -22.85 +/- 0.03 kJ/mol, **7.4x Coudert's K_lp**. (500-cycle test:
1.91e-4 +/- 0.11e-4.) With the framework charges off, 1.66e-4, so
the difference is mostly UFF dispersion, not the mu-OH electrostatics.
<U_gh> - <U_h> = -22.9 kJ/mol. The full run and the isotherm-based numbers (1, 2) will
show whether this is a strong-site Henry regime (large K_H, smaller Langmuir K) or a
real over-binding of the force field.

### Resume logic
`run_isotherm.sh` skips a point whose output contains "Simulation finished". An
unfinished point is rerun from scratch; its partial output is kept as
`Output.incomplete.<timestamp>`. Every result therefore comes from one
uninterrupted run.

---

## What a rigid-framework GCMC can and cannot reproduce for MIL-53

(Sources read by the project owner on 2026-09-19: Bourrelly 2005, Coudert 2008,
Ghoufi & Maurin 2010.)

**The experiment.** Bourrelly et al., JACS 2005, 127, 13519, Fig. 2: CO2 and CH4
isotherms at **304 K up to 30 bar** for **both** MIL-53(Al) (top panel) and
MIL-53(Cr) (bottom panel). CO2 shows a step at **~6 bar** in both metals. The Al
panel is our experimental reference (`reference/bourrelly2005_MIL53Al_CO2_304K.csv`,
digitised by the project owner, not by us).

**What happens physically** (Coudert et al., JACS 2008, 130, 14294, Fig. 5b,
arXiv:1904.09588, which analyses exactly this system). At 304 K the empty framework
is **lp**. CO2 first **closes** it, lp -> np near **0.3 bar**, because np has the
higher affinity (Langmuir Henry constants K_np ~ 9.0e-5 vs **K_lp ~ 2.6e-5
mol kg^-1 Pa^-1**). CO2 then **reopens** it, np -> lp near **6 bar**, because lp has
the larger pore volume. This is a double transition, case "c" of Coudert's
taxonomy. The free-energy difference between the two empty structures is
dF(lp -> np) ~ **2.5 kJ/mol per unit cell**. The predicted low-pressure transition
at 0.3 bar was confirmed by microcalorimetry at 0.25 bar. The reopening step spans
**5-9 bar**: np below 5 bar, fully open lp above 9 bar.

**What our simulation is.** With the lp structure held rigid, we compute the
**"virtual" rigid-host lp isotherm**, the blue dashed curve of Coudert Fig. 5b.
The np-lp breathing cannot appear, **by construction**: the simulation has no np
state to go to. Consequences for the comparison:
- It should agree with experiment **only for P >= 9 bar**, where the real solid is
  fully lp. The comparison with the Bourrelly points is therefore restricted to
  P >= 9 bar, and every figure caption says so.
- Below 5 bar the real solid is np. Our lp loading there is not expected to match
  and is not a force-field test.
- **Quantitative low-pressure test instead:** the initial slope (Henry constant) of
  our isotherm is compared with K_lp ~ 2.6e-5 mol kg^-1 Pa^-1 from Coudert's
  Langmuir fit, with the uncertainty of our slope (Phase 3). This tests the
  CO2-lp interaction where the lp assumption holds by definition. Note: Coudert
  reports K_lp without an error bar, and it comes from a Langmuir fit, so "agreement"
  will be judged against our own uncertainty plus a stated tolerance, not as exact.

**What would be needed to capture the transition.**
1. **Osmotic ensemble** (N_host, mu_CO2, sigma, T), as in Coudert 2008. The stable
   phase at each pressure minimises the osmotic potential, which combines, for each
   rigid phase, dF_host(lp -> np) and the grand potential from that phase's
   rigid-GCMC isotherm. Needed on top of this work: a rigid np isotherm plus
   dF_host (~2.5 kJ/mol/uc). This yields the step pressures.
2. **Hybrid GCMC/MD (HOMC) with a flexible framework**, with a force field that has
   both np and lp minima. Ghoufi & Maurin, J. Phys. Chem. C 2010, 114, 6496
   (DOI 10.1021/jp911484g): their HOMC scheme **captured the lp -> np transition at
   0.35 bar but failed to capture the np -> lp reopening**, which they attributed to
   the highly ordered orientation of CO2 in the np channels. They needed a
   **phase-mixture model** to reproduce the full isotherm. So even a flexible
   simulation does not trivially give the second step.
   Their settings, useful as a sanity reference: **32 unit cells**, Ewald
   electrostatics, van der Waals truncated at **12 A**, **Harris-Yung rigid CO2**,
   **300 K**.
   (For comparison, the DL_POLY `FIELD-MIL53Cr` published by A. Ghoufi, see
   "Structure provenance", describes 2432 atoms = 32 cells x 76, with a flexible
   framework and a CO2 model with q_C = +0.6512, i.e. not the rigid HY model of 2010.)
