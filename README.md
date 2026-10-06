# CO₂ breathing in MIL-53(Al) from rigid GCMC

## What I did

- CO₂ isotherms at 304 K with RASPA2, in the rigid large-pore (lp) form and in two narrow-pore (np) structures.
- Checked the setup against RASPA reference cases, against experiment (Bourrelly et al. 2005) and against LAMMPS.
- Used the osmotic construction of Coudert et al. (2008) to predict the lp→np and np→lp transitions from the two rigid isotherms.

## Main results

1. **The np structure matters more than the force field.** In the hydrated np structure CO₂ hardly fits, and the Henry constant is about 300 times below experiment. In the np cell measured with CO₂ inside (Serre et al. 2007), it is about 6 times above — the same over-binding the force field gives for lp.
2. **With the right np structure, the two rigid isotherms give the double transition**: the framework closes, then reopens.
3. **The pressures are too low**, because the generic force field over-binds CO₂ at low loading. After correcting for it, the closing comes out at 0.29 bar (experiment: 0.27 bar), and the free-energy difference between the empty phases at about 2.4 kJ/mol (Coudert: 2.5). The reopening stays about three times too low; possible reasons are in NOTES.md.

![Osmotic construction](results/osmotic_summary.png)

## Question

Can two rigid GCMC isotherms, one per phase, predict the pressures at which MIL-53(Al) closes and reopens under CO₂? They can, and what limits the answer is the structure used for the narrow-pore phase and the low-coverage over-binding of a generic force field.

## What did not work

- The first narrow-pore structure was the hydrated one. Its pore-limiting diameter is 2.52 Å against CO₂ at 3.30 Å, so almost no CO₂ fits and the Henry constant came out about 300 times below experiment.
- A single-site Langmuir could not describe either isotherm, so the construction was done by numerical integration of the isotherms instead.
- A truncated Coulomb sum gave wrong insertion energies. It failed the check against RASPA and was discarded, not patched.
- The point at 0.5 bar was under-equilibrated. It was rerun with 100,000 initialisation cycles instead of 20,000.
- The drift guard taken from the earlier cross-check fired on about half of all converged points. It was recalibrated before being used.
- A bond-length table silently dropped every O–H bond, because its lookup key was not sorted. That is the bond the narrow-pore binding site depends on.

## What I learned

- For this system the structure of the CO₂-loaded phase mattered more than the force field. Changing it moved the Henry constant by three orders of magnitude, with everything else fixed.
- A fit form chosen for convenience can hide the physics. When the form failed I integrated the data directly rather than keep the form.
- Every loading needs its convention written next to it. Absolute and excess differ by about 11 % at 30 bar here, and the reference data are excess.
- A check is only useful once I have tested that it can fail. One of mine looked strict until I measured how often it fires on data with nothing wrong.

## Reproduce

From `docs/technical_note.md`:

```
micromamba create -n gcmc-mil53 -f environment.yml      # environment
bash scripts/run_isotherm.sh production 4               # lp isotherm, ~12 h on 4 cores
bash scripts/run_isotherm.sh np_co2cell 4               # np isotherm, ~23 h on 4 cores
python scripts/osmotic_numeric.py                       # transitions, seconds
python scripts/figure_osmotic.py                        # the figure above, seconds
```

Both isotherm runs are resumable and skip points that are already finished. The input files come from `scripts/make_inputs.py`; the technical note gives the full sequence.

## Repository

- `scripts/` — setup, runs and analysis
- `results/` — tables and figures
- `NOTES.md` — every decision, including the mistakes I found and corrected
- `docs/technical_note.md` — full technical description and how to reproduce everything

## References

- Bourrelly et al., J. Am. Chem. Soc. 2005, 127, 13519
- Serre et al., Adv. Mater. 2007, 19, 2246
- Coudert et al., J. Am. Chem. Soc. 2008, 130, 14294 (arXiv:1904.09588)
- Salles et al., Angew. Chem. Int. Ed. 2008, 47, 8487
- Loiseau et al., Chem. Eur. J. 2004, 10, 1373
