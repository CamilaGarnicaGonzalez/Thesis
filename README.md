# Drought Characterization in South American River Basins

This repository contains part of a research pipeline for characterizing meteorological and hydrological droughts across six major South American river basins, using standardized climate indices and CMIP6 climate projections.

## Overview

The pipeline covers three stages:

1. **Index calculation** — computation of standardized drought indices (SPEI at 6- and 12-month accumulation periods, and SSI at 12 months) from climate model outputs.
2. **Drought characterization** — event-based analysis using Run Theory (Yevjevich, 1967) to extract drought frequency, duration, and severity for each basin, climate model, and scenario.
3. **Exposure analysis** — assessment of population/land exposure to drought conditions based on the characterized events.

## Study area

Six major river basins in South America:

- Amazon
- Orinoco
- Rio Negro
- Rio Parana
- Sao Francisco
- Tocantins

## Data and scenarios

- **Indices:** SPEI-6, SPEI-12, SSI-12
- **Climate models (CMIP6):** MPI-ESM1-2-HR, UKESM1-0-LL, GFDL-ESM4, IPSL-CM6A-LR, MRI-ESM2-0
- **Scenarios:** Historical (1985–2014), SSP1-2.6 (2041–2070), SSP5-8.5 (2041–2070)

## Methodology — drought event characterization

Drought events are identified using **Run Theory** (Yevjevich, 1967):

- **Event start:** index value below −1 for at least 2 consecutive months
- **Event end:** index value returns to ≥ 0

For each event, three metrics are computed per pixel:

| Metric | Definition |
|---|---|
| Frequency | Number of drought events in the time series |
| Duration | Median event length (months) |
| Severity | Median cumulative absolute index value during an event |

A 5-model ensemble median is calculated for each basin, index, and scenario to summarize results across climate models.

## Repository structure
1_calculo_indices/          - Index calculation scripts (R)
  spei-12/
  spei-6/
  ssi-12/
2_caracterizacion_sequia/   - Drought event characterization (Python)
  spei-12/
  spei-6/
  ssi-12/
3_analisis_exposicion/      - Exposure analysis notebooks
  spei-12/
  spei-6/

## Notes
- This repository shares the analysis code only. Raw and processed climate data (.nc files) are not included due to size and are managed separately on a HPC cluster.

## References

- Yevjevich, V. (1967). An objective approach to definitions and investigations of continental hydrologic droughts. Hydrology Papers, Colorado State University.
- Spinoni, J., et al. (2014). World drought frequency, duration, and severity for 1951-2010. International Journal of Climatology.
- Tabari, H., et al. (2021). Global change impacts on drought.
