# ============================================================================
# SPEI calculation — Amazon basin
# Scenarios: Historical (1985-2014), SSP126 (2041-2070), SSP585 (2041-2070)
# Baseline: Common 1965-2070 (106 years)
# Distribution: log-Logistic (Vicente-Serrano et al. 2010)
# 5 climatic models
# ============================================================================

library(SPEI)
library(ncdf4)

models        <- c("MPI-ESM1-2-HR", "UKESM1-0-LL", "GFDL-ESM4", "IPSL-CM6A-LR", "MRI-ESM2-0")
basin         <- "amazon"
path_D_base   <- "/data/brussel/vo/000/bvo00033/vsc11346/SPEI_R/1.archivos_D_python/"
path_out_base <- "/data/brussel/vo/000/bvo00033/vsc11346/SPEI_R/3.outputs_SPEI_R/4.baseline_106_12ac_OF/Amazon/"

scenarios <- list(
    historical = list(
        path_D_hist    = paste0(path_D_base, "1.baseline_50/Amazon/historical/"),
        path_out       = paste0(path_out_base, "historical/"),
        start_yr       = 1965,
        time_start     = as.Date("1985-01-31"),
        time_ref_start = c(1965, 1),
        time_ref_end   = c(2014, 12),
        is_future      = FALSE
    ),
    ssp126 = list(
        path_D_hist    = paste0(path_D_base, "1.baseline_50/Amazon/historical/"),
        path_D_gap     = paste0(path_D_base, "2.baseline_120/Amazon/ssp126/"),
        path_D_fut     = paste0(path_D_base, "1.baseline_50/Amazon/ssp126/"),
        path_out       = paste0(path_out_base, "ssp126/"),
        start_yr       = 1965,
        time_start     = as.Date("2041-01-31"),
        time_ref_start = c(1965, 1),
        time_ref_end   = c(2070, 12),
        is_future      = TRUE
    ),
    ssp585 = list(
        path_D_hist    = paste0(path_D_base, "1.baseline_50/Amazon/historical/"),
        path_D_gap     = paste0(path_D_base, "2.baseline_120/Amazon/ssp585/"),
        path_D_fut     = paste0(path_D_base, "1.baseline_50/Amazon/ssp585/"),
        path_out       = paste0(path_out_base, "ssp585/"),
        start_yr       = 1965,
        time_start     = as.Date("2041-01-31"),
        time_ref_start = c(1965, 1),
        time_ref_end   = c(2070, 12),
        is_future      = TRUE
    )
)

# ============================================================================
# LOOP — SCENARIOS x MODELS
# ============================================================================
for (scenario_name in names(scenarios)) {
    cfg <- scenarios[[scenario_name]]
    dir.create(cfg$path_out, showWarnings = FALSE, recursive = TRUE)

    for (model in models) {
        print(paste("Processing:", scenario_name, "-", model))

        # --- READ HISTORICAL D (1965-2014) ---
        nc_hist <- nc_open(paste0(cfg$path_D_hist, "D_", basin, "_historical_", model, ".nc"))
        D_hist  <- ncvar_get(nc_hist, "D")
        lat     <- ncvar_get(nc_hist, "lat")
        lon     <- ncvar_get(nc_hist, "lon")
        nc_close(nc_hist)

        n_lon  <- length(lon)
        n_lat  <- length(lat)
        n_hist <- dim(D_hist)[3]  # 600 meses

        if (cfg$is_future) {

            # --- READ GAP D (2015-2020) ---
            nc_gap <- nc_open(paste0(cfg$path_D_gap, "D_", basin, "_", scenario_name,
                                     "_", model, "_gap_2015_2020.nc"))
            D_gap  <- ncvar_get(nc_gap, "D")
            nc_close(nc_gap)
            n_gap  <- dim(D_gap)[3]   # 72 meses

            # --- READ FUTURE D (2021-2070) ---
            nc_fut <- nc_open(paste0(cfg$path_D_fut, "D_", basin, "_", scenario_name,
                                     "_", model, ".nc"))
            D_fut  <- ncvar_get(nc_fut, "D")
            nc_close(nc_fut)
            n_fut  <- dim(D_fut)[3]   # 600 meses

            n_total <- n_hist + n_gap + n_fut  # 1272 meses
            n_out   <- 360                     # 2041-2070

            print(paste("D_hist dims:", paste(dim(D_hist), collapse = " x ")))
            print(paste("D_gap  dims:", paste(dim(D_gap),  collapse = " x ")))
            print(paste("D_fut  dims:", paste(dim(D_fut),  collapse = " x ")))
            print(paste("n_total:", n_total))

        } else {
            n_total <- n_hist   # 600 meses
            n_out   <- 360      # 1985-2014
            print(paste("D dims:", paste(dim(D_hist), collapse = " x ")))
        }

        # --- INSTANTIATE OUTPUT ARRAY ---
        spei_grid <- array(NA_real_, dim = c(n_lon, n_lat, n_out))

        # --- LOOP PIXEL BY PIXEL ---
        for (i in 1:n_lon) {
            for (j in 1:n_lat) {

                if (cfg$is_future) {
                    # Concatenar hist + gap + fut
                    d_combined <- c(D_hist[i, j, ], D_gap[i, j, ], D_fut[i, j, ])
                } else {
                    d_combined <- D_hist[i, j, ]
                }

                D_ts <- ts(d_combined, start = c(cfg$start_yr, 1), frequency = 12)

                result <- tryCatch({
                    spei(D_ts,
                         scale        = 12,
                         ref.start    = cfg$time_ref_start,
                         ref.end      = cfg$time_ref_end,
                         distribution = "log-Logistic"
                    )$fitted
                }, error = function(e) {
                    rep(NA_real_, n_total)
                })

                if (length(result) != n_total) result <- rep(NA_real_, n_total)

                if (cfg$is_future) {
                    # Extraer 2041-2070 = posiciones 913-1272
                    spei_grid[i, j, ] <- as.numeric(result)[913:1272]
                } else {
                    # Extraer 1985-2014 = posiciones 241-600
                    spei_grid[i, j, ] <- as.numeric(result)[241:600]
                }
            }
            if (i %% 5 == 0) print(paste("  lon", i, "of", n_lon, "done"))
        }

        print(paste("SPEI completed:", scenario_name, model))

        # --- REPLACE NA/NaN/Inf ---
        spei_grid[is.nan(spei_grid)]      <- -9999
        spei_grid[is.infinite(spei_grid)] <- -9999
        spei_grid[is.na(spei_grid)]       <- -9999

        # --- SAVE NetCDF ---
        output_file <- paste0(cfg$path_out, "spei12_", basin, "_",
                              scenario_name, "_", model, ".nc")

        time_origin <- if (cfg$is_future) "2041-01-01" else "1985-01-01"
        time_values <- as.numeric(
            as.Date(seq(cfg$time_start, by = "month", length.out = n_out)) -
            as.Date(time_origin)
        )

        lon_dim  <- ncdim_def("lon",  "degrees_east",  lon)
        lat_dim  <- ncdim_def("lat",  "degrees_north", lat)
        time_dim <- ncdim_def("time", paste("days since", time_origin), time_values)

        spei_var <- ncvar_def(
            name     = "SPEI_12",
            units    = "-",
            dim      = list(lon_dim, lat_dim, time_dim),
            missval  = -9999,
            longname = "SPEI-12 (R SPEI package, log-Logistic, ub-pwm, baseline 1965-2070)"
        )

        nc_out <- nc_create(output_file, spei_var)
        ncvar_put(nc_out, spei_var, spei_grid)
        nc_close(nc_out)

        print(paste("NetCDF saved:", output_file))
    }
}

print("All scenarios completed — Amazon, 5 models.")