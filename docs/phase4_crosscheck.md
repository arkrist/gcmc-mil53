# Phase 4: cross-code check, RASPA2 vs LAMMPS `fix gcmc`

CO₂ in rigid MIL-53(Al) lp, 304 K, same structure (4 × 2 × 2, 1216 framework atoms),
same force field (UFF + DDEC charges + TraPPE CO₂, Lorentz–Berthelot, 12.0 Å cut-off,
tail corrections, Ewald). Scope: **three pressures, a cross-code check, not a second
isotherm**. LAMMPS 2025.07.22 (conda-forge, MC/RIGID/KSPACE/MOLECULE).

**Nothing was tuned to make the two codes agree.** Every disagreement below is
reported with the diagnostic that explains it.

## Step 1 — single-point energies on one identical configuration

Framework + **8 CO₂** whose coordinates come from a RASPA restart file (12 decimals).
In LAMMPS every cross term is **A − B − C** (framework+CO₂, framework, CO₂), which is
exact for the pair terms and for the Ewald reciprocal term (quadratic in the charges)
and removes the host–host energy RASPA never computes for a rigid host.
Reproduce: `python scripts/crosscheck_energy.py`.

| term | RASPA [K] | LAMMPS [K] | rel. diff |
|---|---|---|---|
| host–guest LJ (no tail) | −17916.900 | −17916.899 | 1e−7 |
| host–guest Coulomb, real space | 386.035 | 19673.446 | *split differs, see below* |
| host–guest Coulomb, reciprocal | 3.709 | −19283.579 | *split differs* |
| **host–guest Coulomb, total** | **389.744** | **389.868** | 3e−4 |
| guest–guest LJ (no tail) | −500.852 | −500.852 | 1e−8 |
| guest–guest Coulomb, total | 522.263 | 521.842 | 8e−4 |
| tail correction (host–guest + guest–guest) | −646.319 | −646.319 | 1e−8 |
| **TOTAL interaction energy** | **−18152.064** | **−18152.360** | **2e−5** |

**The force fields are identical to 0.002 %.** Two conventions had to be found rather
than assumed, and both would have been silent traps:

1. **LAMMPS `E_vdwl` already includes the tail correction** (`E_tail` is printed
   separately for information); RASPA reports VDW without it. Verified with
   `pair_modify tail no`. Subtracting it turns an apparent 3.6 % disagreement into 1e−7.
2. **The real/reciprocal split is not comparable** — 19673 vs 386 K in real space —
   because the codes choose the Ewald convergence parameter by different criteria.
   Only the sums are physical, and the sums agree.

| | RASPA | LAMMPS |
|---|---|---|
| Ewald parameter | α = 0.265058 Å⁻¹ | G = 0.275229 Å⁻¹ |
| k-vectors | 7 × 9 × 7 | kmax1d = 10, 1438 vectors |
| criterion | `EwaldPrecision 1e-6` | relative force accuracy 1.05e−6 |

## Step 2 — loadings at 10, 20, 30 bar

Reproduce: `python scripts/crosscheck_summary.py` (families below),
`python scripts/crosscheck_gcmc.py --dir runs/crosscheck/seeded --equil 500`.

| family | start | 10 bar | 20 bar | 30 bar |
|---|---|---|---|---|
| from an empty pore | 0 | 7.74 ↑ | 8.46 ↑ | 8.39 ↑ |
| continued from those | below | 8.67 ↑ | 9.17 | 9.14 ↑ |
| seeded from RASPA's configuration | ≈ RASPA | 9.08 ± 0.11 | 9.52 ± 0.04 | 9.42 |
| seeded from above | above | 9.29 ↓ | 10.08 ↓ | 9.91 ↓ |
| **RASPA** | | **9.127 ± 0.053** | **9.604 ± 0.038** | **9.800 ± 0.042** |

(molecules per unit cell; ↑ / ↓ = still drifting at the end of the run.)

**Verdict: consistent, not confirmed.** RASPA's value lies inside the two-sided
LAMMPS bracket at every pressure:

| p [bar] | LAMMPS from below | from above | bracket midpoint | RASPA | verdict |
|---|---|---|---|---|---|
| 10 | 9.076 | 9.285 | 9.18 ± 0.10 | 9.127 | inside |
| 20 | 9.521 | 10.080 | 9.80 ± 0.28 | 9.604 | inside |
| 30 | 9.423 | 9.907 | 9.67 ± 0.24 | 9.800 | inside |

The Phase-1 criterion |Δ| < √(e₁²+e₂²) **cannot be applied to a single LAMMPS run
here**, because that run's block error bar is not the true uncertainty:

- **The particle number is nearly frozen.** ⟨δN²⟩ is a physical property of the
  grand-canonical ensemble, so both codes must reproduce it. LAMMPS gives
  sd(N) = 0.23–0.49 × RASPA's (RASPA: 3.0–3.3 molecules; LAMMPS: 0.7–1.5), i.e. its
  samples are strongly correlated within a run.
- **Runs keep memory of their starting configuration**: from below they end low,
  from above they end high, and the two sides were still converging when the runs
  ended. The spread between them, not the block error bar, is the honest uncertainty.

**Why, concretely.** With kspace and tail corrections LAMMPS requires `full_energy`,
so **every trial move costs a full Ewald evaluation over the whole system**, and
insertions are attempted at a single random position. RASPA updates the Ewald sum
incrementally and inserts with CBMC using 10 trial positions. Near saturation
(~0.4 % acceptance) this is decisive: 12 RASPA points took 12 h, while 10 LAMMPS runs
of 340 k trial moves each took about 70 h of wall time for 3 pressures. This is an
implementation difference, not a difference in the physics.

## Protocol (and why each choice)

| item | choice | reason |
|---|---|---|
| reservoir | `pressure` + `fugacity_coeff` with **RASPA's own Peng–Robinson φ** (0.949190 / 0.899667 / 0.851271) | μ = kT ln(φPΛ³/kT); passing φ makes the reservoirs identical rather than similar |
| rigid CO₂ | molecule template, **MC only, no time integration** | `fix rigid`/`fix shake` restrict `fix gcmc` to exchange moves only (M = 0), which would force translations to come from MD. With no integrator the geometry cannot change and `fix gcmc` does exchanges + rigid-body translations/rotations, as RASPA does |
| intramolecular exclusions | template bonds + `bond_style zero` + `special_bonds lj/coul 0 0 0` | also corrects the Ewald reciprocal term; `neigh_modify exclude` does not. Verified: one isolated CO₂ gives −0.0005 kcal/mol |
| framework | rigid, `neigh_modify exclude group framework framework` | host–host is constant and cancels in every MC energy difference; excluding it only makes the runs cheaper |
| tail corrections | `pair_modify tail yes` | matches RASPA; verified numerically in step 1 |

## What this establishes

- **The force field and its implementation are code-independent** (0.002 % on energies).
  Any remaining difference in loading is sampling, not physics.
- **The loadings are consistent between the codes** within the bracket, at all three
  pressures — but LAMMPS cannot reach RASPA's ±0.05 molecules/uc precision at a
  reasonable cost with this protocol, so this is a consistency check, not a
  validation to that precision.
- **Two guards now live in the analysis scripts**: a drift test (last block vs first
  block against the error bar) and a fluctuation test (sd(N) against RASPA's). Both
  refuse to report an unconverged run as a code disagreement.
