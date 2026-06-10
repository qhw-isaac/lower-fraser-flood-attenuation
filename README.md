# Lower Fraser Flood Attenuation & Benefiting Areas

Regional adaptation of Duarte et al. (2024, *Ecosystem Services*), scoped to **Metro Vancouver (MVRD)** and **Fraser Valley (FVRD)** for UBC Sustainability Scholars 2026 project #2026-030.

The R pipeline maps potential runoff retention by upstream natural ecosystems, identifies downstream built-up and agricultural areas that benefit, and ranks sub-basins by realised benefit (provision × demand).

---

## Approach

Duarte's framework is preserved: SCS Curve Number runoff (USDA TR-55) with Huang slope adjustment, a natural-vegetation counterfactual, sub-basin aggregation, distance-decay-weighted downstream demand, and percentile ranking.

Regional inputs: 30 m NALCMS + AAFC ACI, BC Soil Survey (HYSOGs gap-fill), Copernicus GLO-30, PCIC PRISM wettest-month precipitation, Mohanty 100-yr floodplain, HydroBASINS L12, MV RNIN patches/corridors inside MVRD. CRS: EPSG:3005; decay half-life: 20 km; upstream cap: 100 km flow distance.

---

## Pipeline order

Scripts run in numeric order. **Sub-basins come before rasters** so the working grid is clipped to the hydrologically refined upstream AOI — not a regional-district placeholder that would miss Fraser/Thompson basins draining into the Lower Mainland from outside MVRD ∪ FVRD ∪ SLRD.

```
R/
├── 00_setup.R              libs, paths, CRS/grid helpers, data_path()
├── 01_aoi.R                outcome + downstream AOIs (MVRD ∪ FVRD); district refs
├── 02_subbasins.R          HydroBASINS L12, topology, 01_aoi_upstream.gpkg
│                           ── harmonize inputs (Duarte 01_harmonize) ──
├── 03_lulc.R               NALCMS + AAFC + RNIN → LULC + natural mask
├── 04_soils.R              HSG raster (BC Soil Survey + HYSOGs)
├── 05_dem_slope.R          GLO-30 DEM + slope (degrees)
├── 06_precipitation.R      wettest-month P (mm)
├── 07_floodplains.R        100-yr floodplain mask (Mohanty)
│                           ── service provision (Duarte 02_spa) ──
├── 08_curve_numbers.R      CN baseline + barren counterfactual, slope-adj
├── 09_runoff_retention.R   SCS-CN → PRR (mm) per pixel
│                           ── service-benefiting area (Duarte 03_sba) ──
├── 10_demand.R             built-up + crop area in floodplain, per sub-basin
├── 11_routing_decay.R      flow-distance decay → TDA per sub-basin
│                           ── realised benefit (Duarte 04) ──
├── 12_realized_benefit.R   PRR × TDA → RI + percentile intervals
├── 13_population.R         benefiting population by RI interval
├── 14_interactive_map.R    web data for the shared-watersheds interactive map
└── 99_figures_tables.R     publication figures + summary tables
```

The **interactive map** (`output/interactive_map/index.html`) lets you click a
sub-basin to see its upstream ecosystems, click municipalities to find their
upstream protectors, and surface the watersheds that protect *multiple*
municipalities — the cross-jurisdiction coordination view. Built by
`14_interactive_map.R`; double-click `index.html` to open. It stays in sync with
the model: re-run `R/14_interactive_map.R` after any pipeline change to refresh
`data.js` (see [`output/interactive_map/README.md`](output/interactive_map/README.md)
for the input→dashboard contract).

Every script sources `00_setup.R`, reads raw data via `data_path("layer_id")` from `data_sources.csv`, and writes to `data/processed/` (or `output/` for 99).

---

## AOIs

| AOI | Built by | Extent |
|---|---|---|
| `outcome_aoi` | `01_aoi.R` | MVRD ∪ FVRD |
| `downstream_aoi` | `01_aoi.R` | MVRD ∪ FVRD (demand / floodplain tally) |
| `upstream_aoi` | `02_subbasins.R` | HydroBASINS L12 ≤ 100 km upstream of downstream AOI (+ 5 km buffer) |

Raster scripts (03–07) mask to `upstream_aoi`. Demand and population scripts use `downstream_aoi`.

---

## Run

```r
source("packages.R")   # one-time

# Drop datasets into data/raw/<local_path> per data_sources.csv, then:
for (f in sort(list.files("R", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))) {
  message("→ ", f)
  source(f)
}
```

Set `FLOOD_SCENARIO=wettest_month` (default) before `99_figures_tables.R` if you want a different precip scenario folder under `data/processed/runoff/`.

QA previews land in `output/figures/qa/`. Final figures and tables in `output/figures/` and `output/tables/`.

---

## Data

- [`data_sources.csv`](data_sources.csv) — machine-readable registry (`id`, `local_path`, licence, …)
- [`lookup/`](lookup/) — CN tables, class codes, crop vulnerability scores

Download each layer by hand into `data/raw/<local_path>`.
