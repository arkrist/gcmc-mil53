# Literature reference points

One CSV per source, digitised by hand. **Never add invented points.**

The digitised file (`bourrelly2005_MIL53Al_CO2_304K.csv`) uses
these columns, which are the schema from now on:
`source, system, adsorbate, T_K, convention, pressure_bar, loading_mmol_per_g,
loading_molec_per_uc, notes`. `loading_molec_per_uc` is per **conventional** cell
(832.4 g/mol). The earlier generic schema below is kept for other sources.

Generic schema (earlier):

| column | meaning |
|---|---|
| `source` | short citation, e.g. `Bourrelly2005_JACS127_13519_Fig3` |
| `T_K` | temperature of that isotherm [K] |
| `pressure` | numeric value |
| `pressure_unit` | `Pa`, `kPa`, `MPa`, `bar` or `atm` |
| `loading` | numeric value |
| `loading_unit` | `mol/kg`, `mmol/g`, `molecules/uc` or `cm3STP/g` |
| `convention` | `absolute` or `excess`, **as stated by the source** |
| `kind` | `experiment` or `simulation` |

Optional: `void_fraction_stated` (only if the source gives it), `notes`.

Rules (NOTES.md, "Comparison with the literature"): absolute is compared with
absolute and excess with excess. A point is never converted between the two
unless the source itself states the pore volume or void fraction.
