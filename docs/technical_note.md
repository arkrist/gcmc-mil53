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

The narrow-pore half needs two further environments, kept separate so the validated RASPA
environment is untouched:

```
micromamba create -n zeopp -c conda-forge zeopp-lsmo
micromamba create -n mace -c conda-forge python=3.11 pip
micromamba run -n mace pip install mace-torch ase torch-dftd
```

```
python scripts/build_np_co2_cell.py          # Al framework into the Serre 2007 CO2 cell
python scripts/mace_relax.py --validate      # MACE-MP-0 on both empty cells, ~50 min
python scripts/mace_relax.py --cif structures/derived/MIL-53_Al_np_Serre2007_CrCO2np.cif \
       --fixed-cell --dispersion --label npCO2cell_Serre --fmax 0.03
python scripts/traj_to_cif.py logs/mace_npCO2cell_Serre.traj \
       structures/derived/MIL-53_Al_np_Serre2007_CrCO2np.cif \
       structures/MIL-53_Al_np_CO2cell.cif   # re-attaches labels and charges
network -ha -res out.res structures/MIL-53_Al_np_CO2cell.cif          # Zeo++, env "zeopp"
network -ha -volpo 1.65 1.65 50000 out.volpo structures/MIL-53_Al_np_CO2cell.cif
python scripts/make_inputs.py --cif structures/MIL-53_Al_np_CO2cell.cif --helium-run
python scripts/make_inputs.py --cif structures/MIL-53_Al_np_CO2cell.cif --widom-run --cycles 100000
python scripts/make_inputs.py --cif structures/MIL-53_Al_np_CO2cell.cif --tag np_co2cell \
       --helium-vf 0.280486 --cycles 50000 --init 20000
bash scripts/run_isotherm.sh np_co2cell 4    # resumable; see the cycle-count note below
python scripts/isotherm.py --tag np_co2cell --theta-he 0.280486 --z 2 --no-reference \
       --phase np --fit-window 0 --out-prefix isotherm_MIL53_npCO2cell_CO2_304K
python scripts/insertion_energy.py           # insertion-energy distributions
python scripts/osmotic_numeric.py            # equation 8 numerically; validates itself first
python scripts/dual_site.py                  # dual-site cross-check
python scripts/figure_osmotic.py             # the three-panel figure
scripts/status.sh                            # state of every run
```

Cycle counts differ by point and are recorded per point in the CSV (`cycles_init`,
`cycles_prod`). The four lowest pressures ran 20 000 + 200 000 cycles; the rest 20 000 +
50 000, which was enough for the same relative precision because loading is higher. The
0.5 bar point needed 100 000 initialisation cycles to clear the equilibration guard. The
20, 30 and 50 bar points affect only the capacity figure, not the construction.

Note: `simulate -v` misreports "RASPA 2.0.41" while the package and the run headers report 2.0.50.

Every decision, including those corrected along the way, is recorded in `NOTES.md`.

## 9. The narrow-pore phase: a structure problem, not a force-field one

The narrow-pore structure was taken from the SI of *Powder Diffraction* 2019 (DOI 10.1017/S0885715619000460) through CoRE MOF 2024, the sibling of the lp file, with PACMAN-DDEC6 charges. Identity was fixed by the lattice: in the setting a′ = a + c the cell is 19.499 / 7.617 / 6.569 Å, β = 104.22°, against Loiseau's lt cell 19.51 / 7.61 / 6.58, β = 104.2°. Volume 945.8 Å³ per formula cell, mass 832.415 g/mol, identical to lp.

**In that cell the rigid framework adsorbs essentially nothing**: 0.214 CO₂ per formula cell at 10 bar, still linear, against a literature np capacity near 3.0. The Widom Henry constant is 2.863 × 10⁻⁷ mol kg⁻¹ Pa⁻¹, **326 times below** the experimental 9.0 × 10⁻⁵. That is far too large to be a force-field error, since the same force field gets the lp saturation capacity right to 0.5 %.

Three measurements identify the cause as the structure.

**The framework is thermodynamically invisible to CO₂.** The Widom Rosenbluth factor is 1.058 and the excess chemical potential −0.14 kJ mol⁻¹ — zero within a tenth of a kJ — while the adsorption energy is −25.96 kJ mol⁻¹, *deeper* than lp's −22.83. CO₂ binds more strongly where it fits; there is almost nowhere it fits.

**The distribution of insertion energies says the same.** Sampling 200 000 random rigid-CO₂ insertions per cell with the production force field and a full Ewald treatment of the guest–host term, favourable sites are 224 times rarer in np than in lp (0.068 % against 15.2 % below zero) and 1320 times rarer below −20 kJ mol⁻¹, while the deepest wells are comparable (−26.2 against −31.6 kJ mol⁻¹).

**The pore is closed.** Zeo++ gives a pore-limiting diameter of 2.52 Å and a largest included sphere of 2.83 Å against CO₂'s 3.30 Å kinetic diameter, and **zero percolating accessible volume at every probe size down to helium's** — only 65.9 Å³ of isolated pockets.

### The CO₂-loaded cell

Serre et al. (2007) measured the np cell of MIL-53(**Cr**) under 1 bar CO₂ by in-situ synchrotron diffraction: C2/c, a = 19.713(1), b = 8.310(1), c = 6.806(1) Å, β = 105.85(1)°, V = 1072.5(1) Å³, against 19.685 / 7.849 / 6.782, β = 104.90°, V = 1012.8 Å³ for the hydrated form. **The opening is almost entirely along b, 7.849 → 8.310 Å**, with a and c fixed to 0.1–0.4 % — a shear of the lozenge, not a dilation, which is why isotropic scaling was rejected. Dundar et al. (2017) used that same Cr cell for MIL-53(Al).

Our Al framework was placed in it by the volume-preserving basis change [a + c, b, −a] (determinant exactly +1), keeping fractional coordinates, atom order, labels and charges, with b doubled to preserve the genuine superstructure. The transplant stretches every bond by 2–9 %, so internal coordinates were relaxed at fixed cell with MACE-MP-0 + D3 (54 steps, f_max 0.028 eV Å⁻¹). Bond lengths returned to within 1.5 % of the original and onto the independently measured lp values; only the two X–H bonds changed materially, both *towards* physical values that X-ray systematically underestimates.

MACE-MP-0 was validated first on the two empty cells: with D3 it gives 1442.3 and 947.0 Å³ against 1412.0 and 945.8 measured (+2.14 %, +0.13 %); without dispersion it tracks the dispersion-free PBE reference of Stavitski et al. instead (+1.71 %, −1.19 %). **That validation licenses only what it tests** — basin geometry at fixed or near-fixed cell. A local relaxation started from the experimental geometry stays in that basin whatever the potential, so it says nothing about the relative stability of np and lp, and no phase-stability or transition prediction rests on it.

### What changes

| | original np cell | CO₂-loaded cell |
|---|---|---|
| V per formula cell | 945.8 Å³ | 1072.5 Å³ |
| pore-limiting diameter | 2.52 Å | 3.03 Å |
| percolating volume, He probe | 0 | 485 Å³ (68 % of lp's) |
| helium void fraction | 0.1210 | 0.2805 |
| Widom K_H | 2.863 × 10⁻⁷ | **5.627 × 10⁻⁴** |

**Changing the structure alone raises the Henry constant by a factor of 1965**, with the same force field, charges, cut-off and Ewald treatment. And the sign of the error inverts: against Coudert's experimental K_np the branch goes from **326× too small to 6.3× too large** — the same size and sign as the error the identical force field makes on the lp branch (7.2×). **The discrepancy was the structure; what remains is the low-coverage over-binding already characterised in section 5.**

## 10. The osmotic construction, on the isotherms themselves

Coudert's equation 11 is the analytical specialisation of equation 8 for a single-site Langmuir isotherm. Neither of our simulated branches is one: weighted χ²_red is 318 for the np branch and 1368 for lp. For lp the failure is not site heterogeneity — a dense search over dual-site parameters improves the raw χ² by a factor 0.99, i.e. not at all, and the residual reaches +81σ near 0.5 bar, the signature of cooperative filling. **There is no analytic Langmuir form, of any number of sites, that describes the lp branch.**

Equation 8 was therefore integrated numerically on the measured points: monotone PCHIP interpolation in log p, an analytic Henry segment below the lowest point anchored on the measured Widom K_H, and the Peng–Robinson molar volume. This also disposes of N_max, which had been the most sensitive input: the integral runs from zero to the transition pressure, both transitions lie below 10 bar, and points above never enter.

**Validation.** Synthetic points generated from Coudert's experimental parameters on our own pressure grid, pushed through this machinery, reproduce his analytical transitions to +0.04 % (0.273 bar) and +0.27 % (3.983 against 3.972 bar).

**Result.** The double transition appears from simulation alone, where the Langmuir route gave none at any ΔF_host: 0.060 and 0.167 bar at ΔF = 2.5 kJ mol⁻¹, against 0.25–0.3 and 5–6 bar measured. Both are roughly twentyfold low. A dual-site construction on the np branch agrees to 20–30 %, so this is not an artefact of the interpolation.

### Decomposing the discrepancy

**Corresponding states.** Stretching each branch's pressure axis by its own measured K_sim/K_exp (7.19 for lp, 6.25 for np) puts the closing transition at **0.290 bar against 0.27 measured, 7 % high**, from twentyfold low. **The closing discrepancy is entirely Henry-regime over-binding.** The reopening moves to 1.66 bar and stops there, **3.3 times below** the measured 5–6.

**Inverse ΔF_host** — one free parameter fitted to one transition, the other predicted. Coudert's 2.5 kJ mol⁻¹ was extracted from experimental isotherms, so it is not force-field-consistent. With the unscaled simulated branches **no value of ΔF_host reproduces the experimental closing pressure at all**: the window opens at small ΔF and annihilates near 2.8 kJ mol⁻¹ with the closing never exceeding 0.097 bar. Under corresponding states, ΔF_host fitted to the closing gives **2.40–2.45 kJ mol⁻¹**, within 4 % of Coudert's independently extracted 2.5, and predicts a reopening of 1.6–1.7 bar.

![Osmotic construction](../results/osmotic_summary.png)

### The residual factor of three

Three candidate explanations, none tested here, in the order we would test them.

1. **Hysteresis.** The measured 5–6 bar step is on the *adsorption* branch, which brackets the thermodynamic transition from above; Coudert notes that a single branch only brackets it. Our construction predicts the equilibrium pressure, which should lie below the adsorption step. This alone could account for a substantial part of the residual with no force-field error at all.
2. **A rigid host cannot expand with loading.** Salles et al. (2008) show the np phase expanding progressively above about 4 CO₂ per unit cell. Ours is frozen at the 1 bar geometry, so beyond that loading it must under-hold, which moves the reopening down. A flexible-framework treatment would test this directly.
3. **High-loading error is not the Henry-regime error.** The rescaling applies one factor per branch, measured at zero coverage. Neither branch is Langmuir, so that factor need not describe the error at the loadings that set the reopening.

### Capacity of the virtual rigid host

Reported, and deliberately not used in the construction: the rigid host held at the CO₂-loaded geometry takes **3.9 CO₂ per formula cell** (dual-site fit 3.907, observed 3.889 at 50 bar), against ~3.0 for the real flexible np phase and a pore-volume ceiling of 4.56 at liquid CO₂ density. A rigid host cannot relax away from a geometry that is already open, so it keeps filling where the real solid would have transformed. The 30 and 50 bar points fail the frozen-particle-number guard — sd(N)/√⟨N⟩ of 0.189 and 0.145 against 0.240 and 0.247 for the lp control, with 808 and 494 accepted insertions — so their error bars understate the uncertainty; the plateau is flat enough that 3.9 is safe to two significant figures.

## References

- Bourrelly, S. et al. *J. Am. Chem. Soc.* **2005**, 127, 13519. [10.1021/ja054668v](https://doi.org/10.1021/ja054668v)
- Coudert, F.-X. et al. *J. Am. Chem. Soc.* **2008**, 130, 14294. Postprint: [arXiv:1904.09588](https://arxiv.org/abs/1904.09588)
- Ghoufi, A.; Maurin, G. *J. Phys. Chem. C* **2010**, 114, 6496. [10.1021/jp911484g](https://doi.org/10.1021/jp911484g)
- Loiseau, T. et al. *Chem. Eur. J.* **2004**, 10, 1373. [10.1002/chem.200305413](https://doi.org/10.1002/chem.200305413)
- Nazarian, D.; Camp, J. S.; Sholl, D. S. *Chem. Mater.* **2016**, 28, 785. Data: Zenodo 3986573
- Ramsahye, N. A. et al. *Adsorption* **2007**. [10.1007/s10450-007-9025-5](https://doi.org/10.1007/s10450-007-9025-5)
- Serre, C. et al. *Adv. Mater.* **2007**, 19, 2246 — in-situ synchrotron XRD of CO₂-induced breathing; the np(Cr) cell under 1 bar CO₂ used here
- Salles, F. et al. *Angew. Chem. Int. Ed.* **2008**, 47, 8487. [10.1002/anie.200803067](https://doi.org/10.1002/anie.200803067)
- Dundar, E. et al. *J. Mol. Model.* **2017** (volume/pages not verified here). [10.1007/s00894-017-3281-4](https://doi.org/10.1007/s00894-017-3281-4) — precedent for the Cr CO₂-np cell applied to Al
- Stavitski, E. et al. *Langmuir* **2011**, 27, 3970 — plain-PBE MIL-53(Al) cells, used as a dispersion-free reference
- Liu, Y. et al. *J. Am. Chem. Soc.* **2008**, 130, 11813 — empty low-temperature np MIL-53(Al)
- Batatia, I. et al. MACE-MP-0, **2023**. [arXiv:2401.00096](https://arxiv.org/abs/2401.00096)
- Willems, T. F. et al. Zeo++. *Microporous Mesoporous Mater.* **2012**, 149, 134
- Dubbeldam, D. et al. RASPA. *Mol. Simul.* **2016**, 42, 81
- Frenkel, D.; Smit, B. *Understanding Molecular Simulation*, 3rd ed., 2023
