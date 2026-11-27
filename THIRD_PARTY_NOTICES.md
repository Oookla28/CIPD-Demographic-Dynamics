# Third-Party Notices

This repository relies only on pre-installed open-source Python packages in the local runtime and does not redistribute third-party GPL source code.

- **NumPy 2.1.3** — BSD-3-Clause License
- **SciPy 1.15.2** (embedding **HiGHS**) — BSD-3-Clause License; HiGHS under MIT License
- **Pandas 2.2.3** — BSD-3-Clause License
- **Matplotlib 3.10.1** — Matplotlib License (PSF-compatible)
- **Seaborn 0.13.2** — BSD-3-Clause License
- **SymPy 1.13.3** — BSD-3-Clause License

External reference repositories consulted in literature review (not vendored):
- `bcallaway11/did` (GPL-2), `asheshrambachan/HonestDiD` (MIT), `mrc-ide/ccmpp.tmb` (MIT), `PPgp/bayesTFR` (GPL-3), Ellison et al. (2026) Zenodo archive (`15179793`).


## Level-3 sprint additions (2026-11-26)

### External statistical software (invoked as separate processes; not vendored)
| Package | Version | License | Use |
|---|---|---|---|
| did (Callaway & Sant'Anna) | 2.1.2 | GPL-2 | E05 staggered DiD baseline |
| HonestDiD | 0.2.8 | GPL-3 | E05 Delta^SD(M) robust CIs |
| Rglpk | 0.6-5.1 | GPL-2 | LP interface (compiled against GLPK 5.0) |
| GLPK | 5.0 | GPL-3 | LP solver for Rglpk |
| ccmpp.tmb | git 2f14f12 | MIT | E02/E04 cohort-component baseline |
| bayesTFR | git 2d949b8 | GPL-2 | E04 Bayesian TFR baseline |
| rstan / StanHeaders | CRAN | GPL-3 / BSD-3 | HMC engine for Ellison et al. (2026) replication |
| TMB | CRAN | GPL-2 | ccmpp.tmb engine |
| nimble | CRAN | BSD-3 | optional MCMC baseline |
| wpp2024, wpp2024extra | git HEAD (PPgp) | GPL-3 | WPP2024 data mirror |
| PyFixest | >=0.24 | MIT | E05 FE event study |
| RapidOCR-onnxruntime | PyPI | Apache-2.0 | Census table OCR |
| Ellison et al. (2026) replication code | Zenodo 15179793 | per Zenodo record | E02/E04 GAM+HMC baseline |

All GPL-licensed components are used as separate executables/scripts and are **not** linked
or vendored into this MIT-licensed repository. Their sources and exact versions/commits are
recorded in `provenance/code_sources.yaml`.
