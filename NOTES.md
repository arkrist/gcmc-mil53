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
- Cu-BTC CO2, 323 K, 1 MPa: **PASS 1/1** (`results/sanity_cubtc_co2.{csv,png}`)

  | quantity | ours | reference |
  |---|---|---|
  | absolute loading [molec/uc] | 113.40 +/- 0.80 | 113.68 +/- 0.34 |
  | absolute loading [mol/kg] | 11.717 +/- 0.082 | 11.747 +/- 0.035 |
  | delta / tolerance | -0.29 / 0.87 molec/uc | |
  | acceptance: insertion / deletion / reinsertion | 21.9 / 21.9 / 8.6 % | 21.7 / 21.8 / 8.5 % |
  | acceptance: translation / rotation | 49.8 / 49.8 % | 51.5 / 50.8 % |

  Equilibration: N ~ 116 molecules already at init cycle 5000 of 25,000; production
  fluctuates between 97 and 116. Wall time on the M2: 9450 s (reference: 11,838 s),
  run partly alongside other jobs.
  This passes the chain MIL-53 needs: rigid 3-site charged CO2, framework charges
  from the CIF, Ewald, rotation moves.

  **Do not reuse this example's excess loading.** Its input hard-codes
  `HeliumVoidFraction 0.29`, which appears to be a copied placeholder: a Widom-He
  value for Cu-BTC is far larger. RASPA uses it silently to compute excess = absolute
  - rho_bulk * theta_He * V, so the example's "excess" number (11.59 mol/kg here)
  is not physical. For our system theta_He is computed, not assumed (0.7115, Phase 2).

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
  **Note on this force field:** his MIL-53(Cr) model uses the
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

### Temperature: 304 K (changed from 303 K on 2026-09-20)
Bourrelly 2005 measured at 304 K and Coudert 2008 analyses at 304 K. Everything --
the isotherm and the direct-insertion Henry run -- is at **304 K**, so no avoidable
mismatch is carried. (The helium void fraction stays at its conventional 298 K; it
is a property of the empty framework and the probe, not of the adsorption
temperature.)

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

### What every result carries (CSV and figures)
- **Both conventions, always**: absolute AND excess loading, each with its 95 % error
  bar, in **mol/kg AND molecules per conventional unit cell** (832.4 g/mol;
  1 molecule/uc = 1.2013 mol/kg). The figure carries a second y-axis for molecules/uc.
- **theta_He = 0.7115 (Widom He, 298 K, eps/k 10.9 K, sigma 2.64 A) is stated on the
  figure and in the CSV header**, because excess depends on it.
- The absolute-excess difference (~11 % at 30 bar) is drawn as a **shaded band
  between the two curves**, not buried in a column: the reference points are excess
  (assumed), so the band shows what the convention choice is worth.
- **Per point, the number of ACCEPTED insertions and deletions** (not only the
  acceptance rate), since that count is what the statistics depend on. Rates below
  1 % are flagged, with their counts beside them.
- **Relative error on loading per point**: err95/loading. At the three low pressures
  (0.01, 0.02, 0.05 bar) this is reported explicitly; **if it exceeds 10 % at
  0.01 bar, that is stated and the point is drawn as an open symbol with its error
  bar**, not plotted as if it were solid.

### Low-pressure comparison, REFRAMED 2026-09-20
**Coudert's K_lp = 2.6e-5 mol kg^-1 Pa^-1 is not a Henry constant.** It is the
initial slope of a **Langmuir function fitted to the 9-30 bar branch** of the
experimental isotherm and extrapolated back to zero pressure. Comparing our
directly computed K_H against it is not like-for-like, and the factor 7.4 seen in
Phase 2 is partly an artefact of that extrapolation. Phase 3 therefore does this:

1. **Fit a Langmuir isotherm to the digitised experimental points with P >= 9 bar**
   (`reference/bourrelly2005_MIL53Al_CO2_304K.csv`, 12 of the 13 points; the 7.37 bar
   point is inside the step). Report K and N_max. **Validation of the procedure:**
   if this K reproduces Coudert's 2.6e-5, our fitting procedure matches his and the
   comparison is meaningful; if it does not, the discrepancy is in the procedure and
   is reported as such before anything is said about the force field.
2. **Fit the same Langmuir form to our simulated isotherm over the same window**
   (P >= 9 bar) and compare K and N_max **pairwise** (experiment vs simulation),
   with uncertainties from the fit covariance.
3. **Plot three curves together**: our simulated isotherm, the experimental points
   for P >= 9 bar, and Coudert's virtual rigid-lp Langmuir curve. **The low-pressure
   divergence between them is the result, not a failure**: it is where the real
   material is np and the rigid-lp model cannot follow.
4. **The direct-insertion K_H stays a separate, clearly labelled number.** It is the
   true zero-coverage limit **of our force field** (Widom test insertion, no fit, no
   extrapolation) and is the right quantity to quote when discussing the force field
   itself, never as "the Henry constant of MIL-53(Al) lp" from experiment.

**Values in hand.**
- **K_H(Widom, 304 K) = 1.869e-4 +/- 0.016e-4 mol/kg/Pa**, <U_gh>-<U_h> = -22.84 +/-
  0.04 kJ/mol (`runs/henry_widom_CO2/`, 20,000 cycles, 1166 s). This is the true
  zero-coverage limit **of our force field**, no fit, no extrapolation.
- **Langmuir fit to the digitised experimental points, P >= 9 bar (procedure
  validation, step 1): K = 2.606e-5 +/- 0.060e-5 mol/kg/Pa vs Coudert's 2.6e-5,
  ratio 1.00** (N_max = 12.18 +/- 0.07 mol/kg = 10.14 molec/uc, b = 0.214 +/- 0.006
  bar^-1, chi2_red = 0.49, 12 points). The digitisation and the fitting procedure
  therefore reproduce Coudert's published number, and the pairwise comparison is
  meaningful. It also means this fitted curve *is* his virtual rigid-lp Langmuir
  curve, so his N_max is not needed to draw it.
- Ratio K_H(Widom, force field) / K(Langmuir, experimental lp branch) = **7.2**. The
  two are different quantities (zero-coverage limit vs fitted initial slope over
  9-30 bar); the pairwise Langmuir-vs-Langmuir comparison is the one that tests the
  force field, and it comes with the isotherm.

The 303 K run
(1.933e-4 +/- 0.007e-4 mol/kg/Pa, <U_gh>-<U_h> = -22.85 kJ/mol) is kept in
`runs/henry_widom_CO2_303K/` for the record. The linear fit through the origin on
our own low-pressure points is reported too, as the sampled counterpart of K_H.

### Force-field limitations
The combination used here -- **UFF** Lennard-Jones for the framework, **DDEC**
charges, **TraPPE** CO2, Lorentz-Berthelot -- is a **generic** combination for MOF
screening. It was **not parameterised for MIL-53(Al)**, and nothing in it was fitted
to CO2 adsorption in this material.

What we measured (Phase 2, Widom at infinite dilution): switching the framework
charges off changes K_H by only about 15 % (1.93e-4 -> 1.66e-4 mol/kg/Pa) and
<U_gh>-<U_h> from -22.9 to -22.0 kJ/mol. **The high zero-coverage affinity comes from
the Lennard-Jones term, not from the mu-OH electrostatics.**

Context: the Maurin group derived **system-specific** force fields for MIL-53 rather
than using generic ones, and Ghoufi & Maurin 2010 combined **Harris-Yung** CO2 with
their own framework force field (32 unit cells, Ewald, van der Waals truncated at
12 A, 300 K). Our generic parameters are a deliberate, cheaper choice, and its
consequences must be stated with the result.

**Expected signature.** Bourrelly 2005 finds that CO2 sorbs first onto the hydroxyl
groups (0.75 CO2 per structural OH before the step), so the **mu-OH is the key
adsorption site** in this material. If our simulated isotherm overestimates
low-pressure uptake relative to the virtual rigid-lp curve, that is exactly the
expected signature of a generic force field in a material whose key site is the
mu-OH group -- reinforced by the short published O-H distance (0.86 A) kept for
consistency with the DDEC charges. This is a statement about the model, not a bug to
be tuned away: no parameter will be adjusted to improve the agreement.

### Timing test (2026-09-19, idle M2, sequential, `runs/timing/`)
1000 cycles per point (500 init + 500 prod). Estimated with
`scripts/estimate_walltime.py` (cost per cycle ~ max(20, N); production half of the
test; see the docstring for the assumptions).

| p [bar] | test [s] | N/uc at end | acc. insertion | s/cycle (prod) | full 60k cycles [h] |
|---|---|---|---|---|---|
| 0.1 | 38  | 2.1 | 45.8 % | 0.039 | 0.6 (lower bound, N still rising) |
| 0.5 | 124 | 5.4 | 18.2 % | 0.131 | 2.2 |
| 1   | 160 | 7.5 | 8.1 %  | 0.171 | 2.8 |
| 2   | 174 | 7.6 | 4.2 %  | 0.185 | 3.1 |
| 5   | 189 | 8.6 | 1.9 %  | 0.205 | 3.4 |
| 10  | 194 | 9.3 | 1.3 %  | 0.208 | 3.5 |
| 20  | 200 | 9.6 | **0.76 %** | 0.217 | 3.6 |
| 30  | 200 | 9.4 | **0.59 %** | 0.216 | 3.6 |
| 50  | 204 | 9.8 | **0.56 %** | 0.221 | 3.7 |

**Full isotherm, sequential: ~27 h.** Four points in parallel on the four performance
cores would take roughly 8 h wall time (to be measured, not assumed).

**Acceptance flag (< 1 %): 20, 30 and 50 bar.** The lp channel is near saturation
there (~9.5 CO2/uc), so insertions rarely find room. 50,000 cycles x ~150
molecules x 40 % swaps x 0.6 % still gives ~2e4 accepted exchanges per point,
probably enough, but it will be checked with the loading-vs-cycle traces in Phase 3.
Remedies if needed: more CBMC trial positions, or CFCMC.

**Loadings after 1000 cycles are NOT results** (not equilibrated or converged at low
p; 500 production cycles). Order of magnitude only: ~7 CO2/uc at 1 bar, ~9.3 at
10 bar, ~9.6 at 50 bar. This is the same order as Ghoufi & Maurin's rigid-lp Cr
value (7-8/uc by 15 bar), not 2 or 20. The early plateau is consistent with the
large Widom K_H.

### Resume logic
`run_isotherm.sh` skips a point whose output contains "Simulation finished". An
unfinished point is rerun from scratch; its partial output is kept as
`Output.incomplete.<timestamp>`. Every result therefore comes from one
uninterrupted run.

---

## Phase 3 results (2026-09-20)

### The run
12 points, 304 K, 0.01-50 bar, `runs/production/`, driver with 4 concurrent jobs:
**batch wall time 43,677 s = 12.1 h** (44.6 h of CPU). All 12 finished, none failed,
no missing-VDW warnings.

**The 1000-cycle timing estimator was ~50 % low: 8 h predicted, 12.1 h actual.**
Two reasons, both in the estimator's own docstring as assumptions: (i) it takes the
cost per MC move as constant, whereas it grows with loading (more guest-guest pairs
and a larger Ewald guest contribution), so the high-pressure points, which dominate
the total, are underestimated; (ii) four concurrent jobs share memory bandwidth and
last-level cache on the M2, so four points do not run 4x faster than one. Per-point
comparison: 50 bar was estimated at 3.7 h and took 4.6 h. **Treat the estimate as a
lower bound; multiply by ~1.5 for planning.**

### Convergence
- Relative error on absolute loading: **0.3-1.6 % at every point**, including 1.05 %
  at 0.01 bar (200,000 production cycles there, 10x the rest). No point needed the
  "> 10 %, draw as open symbol" treatment.
- Loading vs cycle (`results/convergence_loading_vs_cycle.png`): at 50 bar the box
  fills within ~2000 cycles and then fluctuates by +/- 0.3 molecules/uc around
  9.99; at 0.1 bar around 1.69; at 0.01 bar the instantaneous count swings between 0
  and 0.5 molecules/uc (2.5 molecules in the box), which the long run averages to
  1 %. Initialisation (10,000 / 20,000 cycles) is amply long everywhere.
- **Insertion acceptance below 1 % at 20, 30 and 50 bar** (0.67 / 0.51 / 0.37 %), as
  the timing test predicted, but the **accepted counts** are 10,305 / 7,998 / 5,993
  insertions (and the same number of deletions). Combined with the flat traces and
  the 0.4-0.6 % relative errors, the sampling is adequate; no rerun is needed.

### The isotherm
`results/isotherm_MIL53_lp_CO2_304K.{csv,png}`: absolute and excess, mol/kg and
molecules/uc, theta_He = 0.7115 stated, accepted counts and relative errors per
point. Absolute rises from 0.19 mol/kg (0.16 molec/uc) at 0.01 bar to
12.01 +/- 0.07 mol/kg (9.99 molec/uc) at 50 bar.

**Our excess loading passes through a maximum near 20 bar** (10.66, 10.89, 10.73,
9.93 mol/kg at 10, 20, 30, 50 bar) because the absolute isotherm saturates while
rho_bulk * V_pore keeps growing. The experimental points still rise at 29 bar.

### Langmuir comparison on the lp branch (P >= 9 bar)
| fit | N_max [mol/kg] | b [bar^-1] | K [mol/kg/Pa] | chi2_red |
|---|---|---|---|---|
| **experiment** (Bourrelly, 12 pts, excess) | 12.183 +/- 0.069 | 0.2139 +/- 0.0061 | **2.606e-5 +/- 0.060e-5** | 0.49 |
| **simulation, absolute** (4 pts) | 12.248 +/- 0.039 | 0.840 +/- 0.045 | **1.029e-4 +/- 0.052e-4** | 2.73 |
| simulation, absolute, 9-30 bar (3 pts) | 12.21 +/- 0.03 | - | 1.069e-4 +/- 0.039e-4 | 1.0 |
| simulation, **excess** (4 pts) | - | - | **FIT REFUSED** | - |

- **Procedure validated (step 1):** our Langmuir fit to the digitised experimental
  points reproduces Coudert's K_lp = 2.6e-5 exactly (ratio 1.00). The digitisation
  and the fitting are sound, so the pairwise comparison is meaningful. This fitted
  curve *is* Coudert's virtual rigid-lp curve.
- **N_max: sim/exp = 1.005** (12.25 vs 12.18 mol/kg; 10.20 vs 10.14 molecules/uc).
  The saturation capacity of the rigid lp framework is reproduced within 0.5 %.
- **K: sim/exp = 3.95** (1.03e-4 vs 2.61e-5). The simulated isotherm rises much
  faster, i.e. the generic force field over-binds at low coverage - the signature
  anticipated in "Force-field limitations".
- **The excess fit is refused, not fudged.** `langmuir.py` checks monotonicity, the
  parameter bound and dK < K, and prints the reason instead of a number. **This is a
  deliberate design choice** (decision of 2026-09-20): a script that says why it will
  not fit is worth more than one that returns K = 0.32 mol/kg/Pa with an error bar
  2400x the value, which is what the unguarded fit did here.
- **CONVENTION MISMATCH, resolved as option (a) (decision of 2026-09-20).** The main
  table and figures compare our **absolute** Langmuir fit with the experiment's
  **excess** fit, and say so in the caption and in this footnote. Our own
  absolute-excess gap is **2.8 % at 10 bar, 5.6 % at 20 bar, 8.8 % at 30 bar and
  17.3 % at 50 bar**.
- **(b) excess-aware fit, for robustness:** fitting
  n_exc(p) = Langmuir_abs(p) - rho_bulk(p) V_pore / M to our excess data gives
  **K = 1.029e-4 mol/kg/Pa, a change of -0.00 %** versus the absolute fit, and the
  same N_max. It is identical *by construction*, not by coincidence: RASPA's
  (absolute - excess) equals rho_bulk V_pore / M to within 1e-6 mol/kg at every
  point (verified numerically), so adding that term back to the excess data
  reconstructs the absolute isotherm exactly. **The conclusion K_sim/K_exp ~ 4 is
  therefore independent of the excess/absolute convention**, which is worth saying
  out loud -- but it is an algebraic identity, not independent evidence.
- Direct-insertion **K_H = 1.869e-4 +/- 0.016e-4 mol/kg/Pa** stays separate: the
  zero-coverage limit of the force field, 7.2x the experimental Langmuir K and 1.8x
  our own fitted Langmuir K - the gap between the last two measures how
  non-Langmuir (energetically heterogeneous) the simulated isotherm is.

### Comparison with experiment, point by point (excess, P >= 9 bar)
| p [bar] | ours, excess [mol/kg] | Bourrelly (nearest p) |
|---|---|---|
| 10 | 10.66 +/- 0.06 | 8.09 at 9.46 bar |
| 20 | 10.89 +/- 0.05 | 9.86 at 20.04 bar |
| 30 | 10.73 +/- 0.05 | 10.51 at 29.03 bar |
| 50 | 9.93 +/- 0.07 | (no data above 29 bar) |
Agreement improves with pressure: ~30 % high at 10 bar, ~10 % at 20 bar, ~2 % at
30 bar. This is consistent with over-binding at low coverage plus a correct
saturation capacity.

## Phase 4: cross-code check against LAMMPS `fix gcmc` (2026-09-21)

Scope, as agreed: **three pressures (10, 20, 30 bar)**, same structure, same force
field, same temperature, absolute loading, compared point by point against RASPA.
A cross-code check, not a second isotherm. **Nothing is tuned to make the two
codes agree**; if they disagreed beyond tolerance, the single-point energies would
say where.

LAMMPS 2025.07.22 from conda-forge (MC, RIGID, KSPACE, MOLECULE packages present).

### Protocol decisions (and why)
| item | choice | reason |
|---|---|---|
| reservoir | `pressure` + `fugacity_coeff`, with **RASPA's own Peng-Robinson phi** (0.949190 / 0.899667 / 0.851271) | mu = kT ln(phi P Lambda^3/kT): passing phi explicitly makes the two reservoirs identical instead of merely similar |
| electrostatics | `kspace_style ewald 1.0e-6`, `full_energy` (LAMMPS applies it automatically with kspace and tail corrections) | every trial move costs a full energy evaluation; this is what makes the runs expensive |
| rigid CO2 | molecule template, **MC only, no time integration**; NOT `fix rigid`/`fix shake` | with shake/rigid the docs allow exchange moves only (M = 0), so translations would have to come from MD. With no integrator the geometry cannot change and `fix gcmc` does exchanges + rigid-body translations/rotations, exactly as RASPA does |
| intramolecular exclusions | bonds in the template + `bond_style zero` + `special_bonds lj/coul 0 0 0` | this also corrects the Ewald reciprocal term. `neigh_modify exclude` does NOT correct kspace and would leave the 1.16 A C-O pair in the reciprocal sum. Verified: one isolated CO2 gives -0.0005 kcal/mol |
| tail corrections | `pair_modify tail yes` | matches RASPA's truncated + analytic tail; verified numerically below |
| framework | rigid, `neigh_modify exclude group framework framework` | the host-host term is constant and cancels in every MC energy difference; excluding it only makes the runs cheaper |

### Step 1: single-point energies on ONE identical configuration
Framework (4x2x2, 1216 atoms) + **8 CO2** whose coordinates were taken from a RASPA
restart file (12 decimals). RASPA prints its decomposition directly; in LAMMPS each
cross term is **A - B - C** (framework+CO2, framework only, CO2 only), which is exact
for the pair terms and for the Ewald reciprocal term (a quadratic form in the
charges) and removes the host-host energy RASPA never computes.
Reproduce with `python scripts/crosscheck_energy.py`.

| term | RASPA [K] | LAMMPS [K] | difference | rel. |
|---|---|---|---|---|
| host-guest LJ (no tail) | -17916.900 | -17916.899 | +0.002 K | 1e-7 |
| host-guest Coulomb, real space | 386.035 | 19673.446 | (split differs) | - |
| host-guest Coulomb, reciprocal | 3.709 | -19283.579 | (split differs) | - |
| **host-guest Coulomb, total** | **389.744** | **389.868** | +0.124 K | 3e-4 |
| guest-guest LJ (no tail) | -500.852 | -500.852 | +0.000 K | 1e-8 |
| guest-guest Coulomb, total | 522.263 | 521.842 | -0.421 K | 8e-4 |
| tail correction (host-guest + guest-guest) | -646.319 | -646.319 | 0.000 K | 1e-8 |
| **TOTAL interaction energy** | **-18152.064** | **-18152.360** | **-0.295 K** | **2e-5** |

**The two codes agree to 0.002 % on the total interaction energy.** Two conventions
had to be handled, both verified rather than assumed:
1. **LAMMPS `E_vdwl` includes the tail correction** (`E_tail` is printed separately
   for information); RASPA reports VDW without it. Verified by re-running with
   `pair_modify tail no`: the difference is exactly `E_tail`. Subtracting it is what
   turns a 3.6 % apparent disagreement into 1e-7.
2. **The real-space/reciprocal split is convention-dependent** and differs wildly
   (19673 vs 386 K in real space), because the codes choose the Ewald convergence
   parameter by different criteria: RASPA alpha = 0.265058 A^-1 with kvec 7x9x7;
   LAMMPS G = 0.275223 A^-1 with kmax1d = 10 (1438 vectors), from a relative force
   accuracy of 1.05e-6. **Only the sums are physical**, and the sums agree.

### Step 2, first attempt (from an empty pore): NOT EQUILIBRATED, discarded
10/20/30 bar, 1500 equilibration + 7000 production steps, each step 20 exchange +
20 translation/rotation attempts (340,000 trial moves per point), 8 h per point on
2 threads with three points in parallel (`runs/crosscheck/gcmc/`).

| p [bar] | LAMMPS mean [molec/uc] | RASPA | delta | tolerance | verdict |
|---|---|---|---|---|---|
| 10 | 7.74 +/- 0.66 | 9.127 +/- 0.053 | -1.39 | 0.67 | not equilibrated |
| 20 | 8.46 +/- 0.71 | 9.604 +/- 0.038 | -1.15 | 0.71 | not equilibrated |
| 30 | 8.39 +/- 0.56 | 9.800 +/- 0.042 | -1.41 | 0.56 | not equilibrated |

**The loading was still rising at the last step** (10 bar: 6.4 at step 1000, 7.8 at
5000, 8.25 at the end; same at 20 and 30 bar). The runs were still filling the pore,
so the means are biased low and the block error bars describe a **trend, not noise**:
last block minus first block = +1.28 / +1.45 / +1.10 molec/uc, larger than the error
bars themselves. **This is a sampling failure, not a code disagreement** -- the
single-point energies (step 1) already showed the two force fields agree to 0.002 %.

Two consequences, both kept:
1. `crosscheck_gcmc.py` now **flags drift** (last block vs first block against the
   error bar) and refuses to call such a point a pass or a failure, in the same
   spirit as the Langmuir refusal logic.
2. Filling an almost-full pore from empty with single-position insertions is the
   wrong way to spend the cycles. LAMMPS must recompute the full Ewald sum for every
   trial move (`full_energy` is mandatory with kspace) whereas RASPA updates it
   incrementally and inserts with CBMC using 10 trial positions -- an inherent
   efficiency difference of the two implementations, not of the physics.

### Step 2, second attempt: four run families, and the verdict
Because a run started from one side keeps memory of where it started, the loading
was approached from **both** sides at every pressure (`scripts/crosscheck_summary.py`,
`results/crosscheck_gcmc_summary.csv`):

| family | 10 bar | 20 bar | 30 bar |
|---|---|---|---|
| from an empty pore (first attempt) | 7.74 (rising) | 8.46 (rising) | 8.39 (rising) |
| continued from those | 8.67 (rising) | 9.17 | 9.14 (rising) |
| seeded from RASPA's configuration | 9.076 +/- 0.110 | 9.521 +/- 0.042 | 9.423 |
| seeded from above | 9.285 (falling) | 10.080 (falling) | 9.907 (falling) |
| **RASPA** | **9.127 +/- 0.053** | **9.604 +/- 0.038** | **9.800 +/- 0.042** |

**Verdict: consistent, not confirmed.** RASPA's loading lies inside the two-sided
LAMMPS bracket at all three pressures (midpoints 9.18 +/- 0.10, 9.80 +/- 0.28,
9.67 +/- 0.24 vs 9.127, 9.604, 9.800). The Phase-1 criterion cannot be applied to a
single LAMMPS run, because its block error bar is not the true uncertainty:

- **sd(N) is 0.23-0.49x RASPA's** (RASPA 3.0-3.3 molecules, LAMMPS 0.7-1.5).
  <dN^2> is a physical property of the grand-canonical ensemble, so the two codes
  must reproduce it; a too-narrow distribution means N is nearly frozen within a run
  and the samples are correlated.
- **Runs retain their starting configuration**: from below they end low, from above
  they end high, and the sides were still converging when the runs ended. The spread
  between them is the honest uncertainty, not the block error bar of either.

**Cause, established not guessed.** With kspace and tail corrections LAMMPS requires
`full_energy`, so every trial move costs a **full Ewald evaluation of the whole
system**, and insertions are tried at a single random position; RASPA updates Ewald
incrementally and inserts with CBMC (10 trial positions). Near saturation
(~0.4 % acceptance) that decides everything: 12 RASPA points cost 12 h, while the
LAMMPS runs cost ~70 h of wall time for 3 pressures. This is an implementation
difference, not physics -- and step 1 had already shown the force fields agree to
0.002 %.

Both failure modes are now **permanent guards** in the analysis scripts: a drift test
(last block vs first block against the error bar) and a fluctuation test (sd(N)
against RASPA's). Neither lets an unconverged run be reported as a code
disagreement. Full write-up: `docs/technical_note.md` (section 6).

## What a rigid-framework GCMC can and cannot reproduce for MIL-53
(rewritten 2026-09-20 against our own numbers)

**What it reproduces.** The **saturation capacity of the open form, to 0.5 %**:
N_max = 12.248 +/- 0.039 mol/kg (10.195 molecules/uc) from our isotherm against
12.183 +/- 0.069 mol/kg (10.141 molecules/uc) from the experimental lp branch. The
pore geometry of the lp structure, the charge model and the GCMC machinery are
therefore sound: the model holds the right amount of CO2 when the channel is full.

**What it cannot do, by construction.** The np-lp step cannot appear: the cell and
the atoms are frozen in the lp geometry, so there is no np state to go to. The
model also says nothing about **P < 9 bar**, where the real material is np (it
closes near 0.3 bar and reopens near 6 bar, Coudert 2008). Our isotherm is the
**virtual rigid-lp branch** in Coudert's sense -- the loading the lp form would
have if it stayed open -- and only its P >= 9 bar part is comparable with
experiment.

**What remains at P >= 9 bar is a force-field effect, not a flexibility effect.**
Over the window where the real solid is fully lp, the deviation falls
monotonically towards 1: simulation/experiment = **1.32 at 10 bar, 1.17 at 20 bar,
1.12 at 30 bar** on absolute loading (1.28 / 1.10 / 1.02 on excess, the reference's
convention); see `results/deviation_vs_pressure.png`. Flexibility is not the
explanation, because these pressures are above the reopening. The explanation is
low-coverage over-binding by the generic force field:
- **K_sim / K_exp = 3.95** (1.029e-4 vs 2.606e-5 mol/kg/Pa, same Langmuir form,
  same window, same procedure), while N_max agrees to 0.5 %. The simulated isotherm
  rises too steeply and saturates too early, then matches once the channel is full.
- **K_H(Widom) = 1.869e-4 mol/kg/Pa is 1.8x our own fitted Langmuir K**, so the
  simulated isotherm is more energetically heterogeneous than a single-site
  Langmuir: a few strong sites dominate at low coverage.
- Switching the framework charges off moves K_H by only ~15 %, so **the excess sits
  in the Lennard-Jones term, not in the mu-OH electrostatics**.

**What would be needed next.**
1. A **system-specific framework force field** for MIL-53, as the Maurin group
   derived, instead of generic UFF + DDEC + TraPPE. That addresses the
   low-coverage over-binding, which is the deviation we actually see.
2. A **flexible-framework treatment** to reach the transition itself: the osmotic
   ensemble (rigid np isotherm + dF_host, Coudert 2008) or hybrid GCMC/MD
   (Ghoufi & Maurin 2010, whose HOMC captured lp -> np at 0.35 bar but not the
   reopening, and needed a phase-mixture model for the full isotherm).
These are separate problems: (1) is about the interaction potential on the lp
branch, (2) is about which structure the host adopts.

