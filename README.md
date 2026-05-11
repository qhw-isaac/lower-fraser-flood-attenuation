---
output:
  html_document: default
  pdf_document: default
---
# Lower Fraser Flood Attenuation & Benefiting Areas

Regional adaptation of the national flood-prevention assessment by Duarte et al. (2024, *Ecosystem Services*), scoped to the **Metro Vancouver Regional District (MVRD)** and **Fraser Valley Regional District (FVRD)** for the **UBC Sustainability Scholars 2026 project #2026-030**.

This R pipeline (1) maps potential runoff retention by upstream natural ecosystems, (2) identifies the downstream built-up, agricultural, and population areas that benefit from that retention, and (3) ranks sub-basins by realized benefit so that conservation and capital-project investments can be prioritized in support of *Metro 2050* and the *BC Flood Strategy*.

---

## 1. Approach

Duarte's national framework is preserved: SCS Curve Number runoff modelling (USDA TR-55) with slope adjustment (Huang et al. 2006), a natural-vegetation counterfactual to isolate retention attributable to ecosystems, sub-basin aggregation with distance-decay-weighted downstream demand, and a realized-benefit ranking.