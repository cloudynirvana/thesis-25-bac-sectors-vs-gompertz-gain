# Aging-like BAC sectors and load×gain Gompertz regimes on one shared subsystem toy

**Thesis #25.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #13 (Gompertz-like hazard from load×gain) and Thesis #18 (Bounded Adaptive Coherence, grounded-Laplacian λ_min sectors).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

On one shared subsystem toy, do aging-like Bounded Adaptive Coherence sectors coincide with load×gain regimes that produce Gompertz-like hazard, or do the two objects carve different regions?

They carve different regions. The shared object is the five-index weight matrix of Thesis #18. One reading is that deposit's sector classifier. The other reading, on the same paths, is a damage equation whose coupling is those weights, scored by the load×gain residual test of Thesis #13. At damage-rate gain zero the exponential map is an exact Gompertz hazard on every class (slope 0.0800 per toy year), and only global decay is aging-like. The additive map of the same damage fails the χ² budget on every class, including the aging-like path (residual 9354 against 100.75). At gain 0.01 the aging-like path fails (residual 260.9) and the uneven path, classed neither, still passes (residual 79.53). On a 13×13 rate grid at gain 0.02 the two sets are disjoint.

No number is taken from either deposit's results file. Thesis #13's six-node slope 0.0840 and residual 5776 are not results of this toy. This deposit does not claim rejuvenation or cure.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/joint_toy.py` | Shared toy: BAC labels and Gompertz residuals (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Eigenvalue paths, hazards, residuals, rate plane |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/joint_toy.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Aging-like BAC sectors and load×gain Gompertz regimes on one shared subsystem toy [Internet]. Thesis #25 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-25-bac-sectors-vs-gompertz-gain

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
