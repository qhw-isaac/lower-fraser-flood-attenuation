# Lower Fraser Flood Attenuation & Benefiting Areas

Regional adaptation of Duarte et al. (2024, *Ecosystem Services*), scoped to **Metro Vancouver (MVRD)** and **Fraser Valley (FVRD)** for UBC Sustainability Scholars 2026 project #2026-030.

---

## Code Architecture

```
R/
├── 00_setup.R              libs, paths, CRS/grid helpers, data_path()
├── 01_aoi.R                outcome + downstream area of interests (AOIs) (MVRD ∪ FVRD)
├── 02_subbasins_fwa.R      BC FWA assessment watersheds, code-derived topology
├── 03_lulc.R               Land use land cover raster (NALCMS, AAFC)
├── 04_soils_hysog.R        HSG raster (HYSOGs250m)
├── 05_dem_slope.R          Slope raster (GLO-30 DEM)
├── 06_precipitation.R      wettest-month precipitation (PCIC PRISM)
├── 07_floodplains.R        100-yr floodplain mask (Mohanty)
│                           
├── 08_curve_numbers.R      CN baseline + barren counterfactual
├── 09_runoff_retention.R   SCS-CN -> Potential Runoff Retention (PRR) per pixel
│                           
```

Every script sources `00_setup.R`, reads raw data via `data_path("layer_id")` from `data_sources.csv`, and writes to `data/processed/` (or `output/` for 99).

---

## Data

- [`data_sources.csv`](data_sources.csv) — data directory (`id`, `local_path`, …)
- [`lookup/`](lookup/) — CN tables, class codes, crop vulnerability scores