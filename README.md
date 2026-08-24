# AutoLatch — Thermally Gated Autocatalytic Kinetics in Flowing Media

**Author:** Leon Sandler, Independent Researcher
**Target journal:** npj Thermal Science and Engineering (Nature Portfolio, open access)
**Article type:** Theoretical / computational study — primary evidence is mathematical and numerical, not experimental.

## Abstract

Delayed transformation of a flowing reactive material into a solid phase is limited, in most existing thermal-management approaches, by a basic trade-off: reaction mechanisms fast enough to complete transformation in a practical window also tend to permit some conversion during transport, because conversion begins gradually from the moment of mixing under monotonic heat accumulation rather than being withheld until a defined thermal threshold is crossed. We formulate and computationally analyze a coupled model in which a discrete thermal gate, H(T − T_melt), controls release of a reaction-triggering species into an autocatalytic (Kamal–Sourour) cure reaction, and ask whether the coupling is predicted to produce a sharper, more controllable separation between a stable transport regime and a rapidly solidifying regime than either mechanism achieves alone. A full 2×2 computational comparison across gated/ungated and autocatalytic/first-order kinetics predicts that the coupled system uniquely minimizes premature conversion during a nominal 25-minute transport period (α ≈ 0.09) while achieving the sharpest transition of the four cases (S = 1.84); the conventional (ungated, first-order) baseline is predicted to reach α ≈ 0.88 with S ≈ 20.6 under the same conditions. A dimensionless coupling ratio, Λ = τ_release/τ_cure — analogous to a Damköhler number for this thermally gated reactive-transport problem — organizes the parameter space and reveals that predicted transition sharpness and absolute transformation time can diverge at large Λ, meaning a favorable relative sharpness metric alone does not guarantee the transformation completes within a practical absolute time. Monte Carlo propagation of joint parameter uncertainty indicates the reaction-kinetics parameters, not the thermal trigger, dominate predicted performance variability under the illustrative ranges used. A non-ideal (finite-width) thermal gate is also analyzed, predicting that a broader activation range reduces premature conversion at the cost of slower completion — a quantifiable trade-off relevant to real phase-change materials with a melting range rather than a sharp transition point. Three falsifiable predictions, each paired with a standard dynamic-DSC or rheo-DSC validation protocol, are stated for future experimental testing.

## Repository contents

```
manuscript/   AutoLatch_npj_ThermSciEng_Manuscript.docx  (CC BY 4.0)
code/         theoretical_study.py                       (MIT)
figures/      generated on run — every Result 1-9 figure in the manuscript
```

## Reproducing the figures

```bash
pip install numpy scipy matplotlib
python code/theoretical_study.py
```

All parameters are literature-informed illustrative values, not fitted to any measured system — the script says so inline everywhere it matters. Every number and figure referenced in the manuscript is genuinely computed by this script, not asserted after the fact.

## Citation

If you use this model or code, please cite both the paper and the software record:

- Paper (Zenodo preprint DOI): [10.5281/zenodo.22073391](https://doi.org/10.5281/zenodo.22073391)
- Code (Zenodo software DOI): [10.5281/zenodo.22073393](https://doi.org/10.5281/zenodo.22073393)

## License

- Manuscript text and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) — see `manuscript/LICENSE`
- Code: [MIT](LICENSE)
