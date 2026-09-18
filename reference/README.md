# Literature reference points

One CSV per source, digitised by hand. **Never add invented points.**

Required columns:

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
