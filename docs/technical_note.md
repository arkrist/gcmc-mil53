# Technical note — GCMC of CO₂ in rigid large-pore MIL-53(Al)

Version 2, 22 September 2026. Technical record of the workflow, its validation and its limits.

## 1. System

MIL-53(Al), Al(OH)(bdc), large-pore (lp) form. Chains of corner-sharing AlO₄(OH)₂ octahedra linked by terephthalate, forming one-dimensional lozenge channels of about 8.5 Å free diameter.

| Quantity | Value |
|---|---|
| Space group | Imma (No. 74) |
| Cell | a = 6.608, b = 16.675, c = 12.813 Å; V ≈ 1412 Å³ |
| Content | Al₄C₃₂H₂₀O₂₀, 76 atoms, 832.4 g mol⁻¹, 2 channels per cell |
| Supercell | 4 × 2 × 2 = 16 cells, 1216 framework atoms, 26.4 × 33.4 × 25.6 Å |
| Boundary conditions | periodic in all three directions; framework rigid |

Structure file: `SABVUN_clean.cif`, CoRE MOF 2014 DDEC set (Nazarian, Camp & Sholl 2016), Zenodo record 3986573, md5 verified. Stored there as a reduced primitive cell (6.6085 / 11.0216 / 11.0216 Å, 98.3 / 107.4 / 107.4°) describing the same lattice; an exact change of basis (determinant 2) gives the 76-atom conventional cell. Atom count, composition, molar mass, net charge, density and nearest-neighbour distances were checked to be unchanged by the transformation.

`scripts/check_cif.py` verifies the cell first (1 % per axis, 0.5° on angles, axes compared after sorting) and only then composition, µ-OH hydrogens, partial charges and force-field coverage. It exits on the first failure.

## 2. Force field

| Item | Choice | Reason |
|---|---|---|
| CO₂ | TraPPE, rigid, 3 sites | Fitted to CO₂ vapour–liquid equilibria, so the bulk fluid is correct at 304 K, 0.13 K below the critical temperature; the site charges reproduce the quadrupole, which drives the interaction with the µ-OH groups |
| Framework LJ | UFF | Generic, covers Al, standard baseline in MOF screening; not parameterised for MIL-53 |
| Framework charges | DDEC, from the same CIF as the coordinates | Needed because of the CO₂ quadrupole and the polar µ-OH groups |
| Mixing | Lorentz–Berthelot | Standard |
| Cut-off | 12.0 Å, tail corrections; Ewald 1e-6 | 12.8 Å would leave a 0.03 Å margin along c after replication |

Known systematic: the µ-OH bond length in the CIF is 0.86 Å, the X-ray value rather than the ≈ 0.97 Å of a real O–H bond. It was kept because the DDEC charges were computed on that geometry. This matters here specifically, since CO₂ binds first at the hydroxyl groups.

## 3. Method

Grand canonical Monte Carlo (μVT): temperature, cell volume and adsorbate chemical potential imposed, number of molecules fluctuating. The reservoir is set by pressure through the fugacity f = φp, with φ from Peng–Robinson (0.949, 0.900, 0.851 at 10, 20, 30 bar). Moves: translation, rotation, reinsertion, and swap with configurational-bias insertion at 10 trial positions. 10 000 initialisation cycles discarded, 50 000 production cycles averaged (220 000 at 0.01 bar). Error bars are 95 % confidence intervals over five blocks.

Absolute loading is the total amount inside the pore volume; excess subtracts ρ_bulk(p)·V_pore. Both are reported throughout, with the convention of each reference dataset stated.

Auxiliary quantities:

- Helium void fraction by test insertions at 298 K (ε/k = 10.9 K, σ = 2.64 Å, same framework and cut-off, 500 000 cycles): **θ_He = 0.7115 ± 0.0004**, pore volume 0.727 cm³ g⁻¹. The CoRE geometric value for the same structure is 0.587; the definitions differ.
- Henry coefficient by Widom test insertions: **K_H = 1.869 × 10⁻⁴ mol kg⁻¹ Pa⁻¹** at 304 K, ΔU = −22.85 kJ mol⁻¹. Switching off framework charges lowers it only to 1.66 × 10⁻⁴.

## 4. Validation against the code's own references

| Check | Conditions | This work | Reference | Δ / tolerance |
|---|---|---|---|---|
| CH₄ in MFI | 300 K, 10 kPa | 0.3485 ± 0.0086 | 0.3482 ± 0.0103 | +0.0003 / 0.013 |
| CH₄ in MFI | 300 K, 100 kPa | 2.863 ± 0.042 | 2.898 ± 0.033 | −0.034 / 0.053 |
| CO₂ in Cu-BTC | 323 K, 1 MPa, charged CO₂, Ewald, rotations | 113.40 ± 0.80 | 113.68 ± 0.34 | −0.29 / 0.87 |

Loadings in molecules per unit cell. The random-number stream differs between builds, so a check passes when |Δ| < √(e₁² + e₂²) with the 95 % intervals. The Cu-BTC case exercises what MIL-53 needs: rigid charged CO₂, framework charges, Ewald summation, rotation moves. Its excess loading was not used — the bundled example hard-codes a helium void fraction of 0.29, apparently copied from another example.

Three properties of RASPA established during this stage: there is no `ComputeFugacityCoefficient` keyword (Peng–Robinson is the default when `FugacityCoefficient` is absent); the bundled CO₂ model is not TraPPE (charges ±0.6512, C–O 1.149 Å), so TraPPE was written explicitly; missing Lennard-Jones parameters are set to zero with only a warning, so the structure checker and driver abort in that case.

## 5. Isotherm

Twelve pressures from 0.01 to 50 bar at 304 K, 12.1 h wall time on four cores (44.6 h CPU), no failures. Relative errors 0.3–1.6 % at every point, including 1.05 % at 0.01 bar. Insertion acceptance falls to 0.67, 0.51 and 0.37 % at 20, 30 and 50 bar, with 10 305, 7 998 and 5 993 accepted insertions respectively and flat loading traces.

Absolute loading runs from 0.19 mol kg⁻¹ at 0.01 bar to 12.01 ± 0.07 mol kg⁻¹ (9.99 molecules per cell) at 50 bar. Excess passes through a maximum near 20 bar, as it must once absolute loading saturates. The absolute–excess gap is 2.8 % at 10 bar, 5.6 % at 20, 8.8 % at 30, 17.3 % at 50.

Reference data: CO₂ on MIL-53(Al) at 304 K, digitised from Figure 2 of Bourrelly et al. 2005, 13 points from 7.4 to 29 bar, ±0.1 bar and ±0.05 mmol g⁻¹. The paper reports nᵃ from manometry without stating a conversion, so the data are treated as excess (assumed) and the convention is carried in the CSV and every caption. Only p ≥ 9 bar is used.

### Langmuir fits over the lp branch

| Fit | N_max (mol kg⁻¹) | K (mol kg⁻¹ Pa⁻¹) | χ²_red |
|---|---|---|---|
| Experiment (Bourrelly 2005, 12 points, excess) | 12.183 ± 0.069 | 2.606 ± 0.060 × 10⁻⁵ | 0.49 |
| Simulation (absolute, 4 points) | 12.248 ± 0.039 | 1.029 ± 0.052 × 10⁻⁴ | 2.73 |
| Ratio | 1.005 | 3.95 | — |

The fit to the digitised experimental points returns the K_lp of Coudert et al. (2008) to a ratio of 1.00, which validates the procedure and identifies the fitted experimental curve with their virtual rigid-lp branch. An excess-aware fit of the simulation returns the same K, because RASPA's absolute − excess equals ρ_bulk·V_pore/M exactly; that is an identity, not independent evidence. A direct Langmuir fit of the excess curve is refused by the fitting script, since the curve is non-monotonic.

Point by point, simulation / experiment is 1.32, 1.17, 1.12 at 10, 20, 30 bar on the absolute convention and 1.28, 1.10, 1.02 on excess.

## 6. Cross-code check with LAMMPS

Verdict: **consistent, not confirmed.**

**Energies.** On one identical configuration the two codes agree to 0.002 % on the total interaction energy (−18152.06 vs −18152.36 K); LJ and tail match to 1e-7, Coulomb totals to 3e-4. Two conventions had to be established rather than assumed: LAMMPS's `E_vdwl` already contains the tail correction (subtracting it produced an apparent 3.6 % disagreement), and the real/reciprocal Ewald split is not comparable between codes — 19673 vs 386 K — because the Ewald parameter is chosen differently (RASPA α = 0.2651 Å⁻¹, 7×9×7 vectors; LAMMPS G = 0.2752 Å⁻¹, 1438 vectors). Only the sums are physical, and they agree. The force field and its implementation are therefore code-independent, and any difference in loading is sampling.

**Loadings**, molecules per unit cell:

| p | LAMMPS from below | from above | midpoint | RASPA |
|---|---|---|---|---|
| 10 bar | 9.076 | 9.285 | 9.18 ± 0.10 | 9.127 ± 0.053 |
| 20 bar | 9.521 | 10.080 | 9.80 ± 0.28 | 9.604 ± 0.038 |
| 30 bar | 9.423 | 9.907 | 9.67 ± 0.24 | 9.800 ± 0.042 |

RASPA lies inside the two-sided LAMMPS bracket at every pressure. The single-run block error bar is not the true uncertainty here: sd(N) from LAMMPS is 0.23–0.49× the RASPA value, and ⟨δN²⟩ is a physical property of the grand canonical ensemble that both codes must reproduce; runs also retain memory of their starting configuration, ending low from below and high from above. The spread between the two sides is the honest uncertainty.

Cause of the sampling gap, established rather than guessed: with kspace and tail corrections, `full_energy` is mandatory in LAMMPS, so every trial move costs a full Ewald evaluation and insertions are attempted at a single position, while RASPA updates Ewald incrementally and uses configurational bias with 10 trial positions. At ~0.4 % acceptance this decides everything: 12 RASPA points cost 12 h, three LAMMPS pressures cost ~70 h.

Nothing was tuned to improve agreement. A drift test and a fluctuation test are now permanent guards in the scripts, either of which blocks an unconverged run from being reported as a code disagreement.

## 7. What this establishes, and what it does not

- The workflow is correct: two reference calculations reproduced, and two codes agreeing to 0.002 % on energies.
- Structure, charge model and pore geometry are sound: the saturation capacity of the open framework is reproduced to 0.5 % (10.20 vs 10.14 molecules per cell).
- The generic UFF/DDEC/TraPPE force field over-binds CO₂ at low coverage: Langmuir K 3.95× the experimental lp value, Widom K_H 7.2×, and simulated low-pressure points above the single-site Langmuir fit, so the computed isotherm is more heterogeneous than Langmuir. The excess sits in the Lennard-Jones term, not in the µ-OH electrostatics.
- **Not addressed, by construction:** the narrow-pore/large-pore transition. A rigid lp framework yields only the virtual rigid-lp branch; below 9 bar the real material is in the np form. Reaching the transition requires the osmotic-ensemble construction (a rigid np isotherm plus ΔF_host) or a flexible framework with volume moves — hybrid GCMC/MD.
- Other known systematics: the 0.86 Å O–H distance; the assumed excess convention of the reference data; the digitisation uncertainty; and a 1 000-cycle timing estimator that underestimated wall time by ≈ 50 %, since the per-move cost grows with loading and parallel jobs share memory bandwidth.

## 8. Reproduce

```
micromamba create -n gcmc-mil53 -f environment.yml
micromamba activate gcmc-mil53
bash scripts/run_example.sh mfi_ch4        # sanity checks
bash scripts/run_example.sh cubtc_co2
python scripts/fetch_structure.py          # Zenodo 3986573, md5 checked
python scripts/check_cif.py structures/MIL-53_Al_lp.cif
bash scripts/run_helium.sh                 # void fraction
bash scripts/run_widom.sh                  # Henry coefficient
python scripts/make_inputs.py
bash scripts/run_isotherm.sh --jobs 4      # ~12 h on 4 cores, resumable
python scripts/isotherm.py                 # CSV and figures
python scripts/langmuir.py                 # fits
```

Note: `simulate -v` misreports "RASPA 2.0.41" while the package and the run headers report 2.0.50.

Every decision, including those corrected along the way, is recorded in `NOTES.md`.

## References

- Bourrelly, S. et al. *J. Am. Chem. Soc.* **2005**, 127, 13519. [10.1021/ja054668v](https://doi.org/10.1021/ja054668v)
- Coudert, F.-X. et al. *J. Am. Chem. Soc.* **2008**, 130, 14294. Postprint: [arXiv:1904.09588](https://arxiv.org/abs/1904.09588)
- Ghoufi, A.; Maurin, G. *J. Phys. Chem. C* **2010**, 114, 6496. [10.1021/jp911484g](https://doi.org/10.1021/jp911484g)
- Loiseau, T. et al. *Chem. Eur. J.* **2004**, 10, 1373. [10.1002/chem.200305413](https://doi.org/10.1002/chem.200305413)
- Nazarian, D.; Camp, J. S.; Sholl, D. S. *Chem. Mater.* **2016**, 28, 785. Data: Zenodo 3986573
- Ramsahye, N. A. et al. *Adsorption* **2007**. [10.1007/s10450-007-9025-5](https://doi.org/10.1007/s10450-007-9025-5)
- Dubbeldam, D. et al. RASPA. *Mol. Simul.* **2016**, 42, 81
- Frenkel, D.; Smit, B. *Understanding Molecular Simulation*, 3rd ed., 2023
