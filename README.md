# A Six-Week Heat-Accumulation Window for Weekly EHEC Notifications in Korea, 2005–2024

Reproducible analysis code for the manuscript submitted to *Environmental Health and Preventive Medicine* (EHPM).

**Authors:** Seongdae Kim, Byung Chul Chun

## Overview
Twenty years (2005–2024) of weekly national enterohaemorrhagic *Escherichia coli* (EHEC) notifications in the Republic of Korea were paired with weekly climate covariates and analysed with a negative-binomial generalized additive mixed model (NB-GAMM) and distributed-lag formulations. Residual autocorrelation is modelled with a first-order autoregressive [AR(1)] structure plus a quarterly harmonic.

> **Methodological note.** `mgcv::bam(family = nb(), rho = ...)` silently ignores the AR(1) term for generalized (non-Gaussian) models. Proper NB + AR(1) here uses `gamm() + corAR1`. All reported estimates are autocorrelation-adjusted.

## Headline results (autocorrelation-adjusted)
- Single-week temperature effect significant at lag 3 (IRR ≈ 1.05 per 1 °C, p ≈ 0.010).
- Cumulative temperature effect: +9.2 % per 1 °C at the five-week window (significant) and a peak of +8.6 % at six weeks (p ≈ 0.055, borderline).
- Population offset, family choice (NB best), period stratification, and outbreak treatment do not change the conclusion; worst-case season–temperature concurvity = 0.976.

## Files
- `EHEC_GAM_통합full_v2_260615_reproducible.R` — main reproducible pipeline (Tables 1–3, S1–S6; diagnostics).
- `figures/` — figure scripts (Figure 1–2, S1–S4) and outputs.

## Data
The analysis uses aggregated weekly counts with **no personally identifiable information**:
- Weekly EHEC notifications — Korea Disease Control and Prevention Agency (KDCA) Infectious Disease Portal.
- Climate (temperature, humidity, precipitation, wind) — Korea Meteorological Administration (KMA) ASOS, Open MET Data Portal.
- PM10 — AirKorea / National Institute of Environmental Research (NIER).

Set `BASE_IV` at the top of the script to your local data directory containing the preprocessed CSVs.

## Reproduce
```r
# R 4.5.x
install.packages(c("dplyr","mgcv","nlme","MASS","dlnm","splines","ggplot2","patchwork"))
# edit BASE_IV, then:
Rscript EHEC_GAM_통합full_v2_260615_reproducible.R
```
The script prints Table 1–3, Supplementary Tables S1–S6 values, and the autocorrelation diagnostic; the companion figure script produces Figures 1–2 and S1–S4.

## Citation
Kim S, Chun BC. A six-week heat-accumulation window for weekly EHEC notifications in Korea, 2005–2024. *Environ Health Prev Med.* (under review).
Code archived at Zenodo: https://doi.org/10.5281/zenodo.20695800 (concept DOI)

## License
MIT License.

## Contact
Byung Chul Chun — chun@korea.ac.kr
