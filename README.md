# AutoLatch — Coupled Thermal Activation and Autocatalytic Cure Kinetics

**Author:** Leon Sandler, Independent Researcher
**Preprint server:** ChemRxiv
**Target journal:** Macromolecular Theory and Simulations (Wiley)
**Article type:** Theoretical / computational study — primary evidence is mathematical and numerical, not experimental.

> **Version 2** of this model. The original draft (retained in `manuscript/AutoLatch_npj_ThermSciEng_Manuscript.docx`) used a prescribed release-availability function and an isothermal assumption, and was submitted to npj Thermal Science and Engineering. Peer-review feedback on that draft — thermal transients/exotherm neglected, the release function lacking physical justification, an arbitrarily scaled control comparison, a structural (not predictive) induction result, and a single dimensionless ratio being insufficient to organize the system — motivated this substantially revised v2, now targeting ChemRxiv (preprint) and Wiley's *Macromolecular Theory and Simulations* (peer review).

> **Version 3** (this release) leaves every v2 number unchanged and adds four analyses that test the v2 conclusions: a capsule-population model in which the gate width emerges from a spread of melting temperatures (w_eff ≈ 0.55σ); a storage-stability objective with five tuned control families (Arrhenius, same energy balance) to test the Pareto claim; a design rule for the admissible melting-point spread; and Sobol variance attribution. The v2 constant-rate first-order control is superseded for the non-dominance claim. A status-of-claims table and two further falsifiable predictions are added. See `results_v3.json`.

## Abstract (v2, extended in v3 manuscript)

Delayed transformation of a reactive material into a solid phase is limited, in most thermal-management approaches, by a basic trade-off between transport stability and transformation speed. We formulate and computationally analyze a coupled zero-dimensional model in which (i) a lumped thermal energy balance governs the fluid element's temperature, including exothermic self-heating from the reaction itself; (ii) release of a reaction-triggering species follows a physically motivated first-order relaxation toward a temperature-dependent, Arrhenius-rate-controlled saturation level, rather than a prescribed availability function; and (iii) the triggered reaction follows Kamal–Sourour autocatalytic kinetics. Four dimensionless groups — Λ=0.911 (release-to-reaction timescale ratio), Θ=0.454 (thermal-to-reaction timescale ratio), ψ=10.08 (reduced exothermicity), and Ar=18.9 (activation-sensitivity number) — organize the parameter space. An objectively matched (root-solved, not arbitrarily scaled) comparison against a first-order control shows the coupled autocatalytic system is Pareto-optimal: none of 60 sampled first-order rate constants simultaneously matches or beats it on both premature conversion and completion time. A two-parameter regime map over (Θ, Λ) identifies thermal-limited, release-limited, and reaction-limited regions, with peak exothermic overshoots up to 204°C above setpoint in the thermal-limited region. Because complete suppression of conversion below the nominal activation threshold is a structural consequence of the idealized (zero-width) gate rather than an independent prediction, the falsifiable claim is reframed around a finite-width thermal gate: increasing gate width from 0.5°C to 15°C is predicted to reduce premature conversion (α at t=25 min falling from 0.073 to 0.039) at the cost of slower completion (t90 rising from 34.1 to 40.4 min). Vectorized Monte Carlo propagation of joint uncertainty across nine kinetic, release, and thermal parameters (N=20,000) indicates 69.9% of draws satisfy both a minimum induction time and a 75-minute completion target, with reaction-kinetics parameters (Ea2, Spearman r=0.69) dominating predicted variance. Three falsifiable predictions, each paired with a standard dynamic-DSC or rheo-DSC validation protocol, are stated for future experimental testing.

## Repository contents

```
manuscript/   AutoLatch_ChemRxiv_WileyMTS_Manuscript_v3.docx (CC BY 4.0, current v3; build_v3.py reads results_v2/v3.json)
              AutoLatch_ChemRxiv_WileyMTS_Manuscript.docx   (CC BY 4.0, v2)
              AutoLatch_npj_ThermSciEng_Manuscript.docx     (CC BY 4.0, superseded v1, kept for record)
code/         theoretical_study.py                          (MIT; v2 results, unchanged)
              extended_core.py, extended_analyses.py        (v3: population, control families, design rule, Sobol)
              test_extended.py                              (7 checks: v2 reference, energy bound, Sobol on Ishigami)
figures/      generated on run — f0, t1-t9 (v2) and t10-t13 (v3)
```

## Reproducing the figures

```bash
pip install numpy scipy matplotlib
python code/theoretical_study.py          # v2 figures and results_v2.json
python code/extended_analyses.py          # v3 figures and results_v3.json (about 6 minutes)
python code/test_extended.py
```

All parameters are literature-informed illustrative values, not fitted to any measured system — the script says so inline everywhere it matters. Every number and figure referenced in the manuscript is genuinely computed by this script, not asserted after the fact.

## Citation

If you use this model or code, please cite both the paper and the software record (concept DOIs shown below always resolve to the latest version):

- Paper (Zenodo preprint, concept DOI): [10.5281/zenodo.22073390](https://doi.org/10.5281/zenodo.22073390) — v3: [10.5281/zenodo.23201205](https://doi.org/10.5281/zenodo.23201205)
- Code (Zenodo software, concept DOI): [10.5281/zenodo.22073392](https://doi.org/10.5281/zenodo.22073392) — v3: [10.5281/zenodo.23201204](https://doi.org/10.5281/zenodo.23201204)

## License

- Manuscript text and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — see `manuscript/LICENSE`
- Code: [MIT](LICENSE)
