# ============================================================================
# SSI-12 calculation — Orinoco basin
# Scenarios: Historical (1985-2014), SSP126 (2041-2070), SSP585 (2041-2070)
# Baseline: Historical 1965-2014 (50 yrs) / SSP 1965-2070 (106 yrs)
# Distribution: Gamma (Vicente-Serrano et al. 2010)
# Accumulation: 12 months
# 5 climatic models
# ============================================================================

library(SPEI)
library(ncdf4)

models        <- c("MPI-ESM1-2-HR", "UKESM1-0-LL", "GFDL-ESM4", "IPSL-CM6A-LR", "MRI-ESM2-0")
basin         <- "orinoco"

path_flo_base <- "/data/brussel/vo/000/bvo00033/vsc11346/SSI_R/2.Codigos_SSI/3.Outputs_flo_out/Orinoco/"
path_gap_base <- "/data/brussel/vo/000/bvo00033/vsc11346/SSI_R/0.Gap_filling_2015-2020/2.Outputs_gap_filling/Orinoco/"
path_out_base <- "/data/brussel/vo/000/bvo00033/vsc11346/SSI_R/3.Outputs_SSI/2.baseline_105/Orinoco/"

# ============================================================================
# SCENARIO CONFIGURATION
# ============================================================================

scenarios <- list(
    historical = list(
        path_flo_hist  = paste0(path_flo_base, "historical/"),
        path_out       = paste0(path_out_base, "historical/"),
        start_yr       = 1965,
        time_start     = as.Date("1985-01-31"),
        time_origin    = "1985-01-01",
        time_ref_start = c(1965, 1),
        time_ref_end   = c(2014, 12),
        n_baseline     = 600,
        n_analysis     = 360,
        skip_months    = 240,
        is_future      = FALSE
    ),
    ssp126 = list(
        path_flo_hist  = paste0(path_flo_base, "historical/"),
        path_flo_gap   = paste0(path_gap_base, "ssp126/"),
        path_flo_ssp   = paste0(path_flo_base, "ssp126/"),
        path_out       = paste0(path_out_base, "ssp126/"),
        start_yr       = 1965,
        time_start     = as.Date("2041-01-31"),
        time_origin    = "2041-01-01",
        time_ref_start = c(1965, 1),
        time_ref_end   = c(2070, 12),
        n_baseline     = 1272,
        n_analysis     = 360,
        skip_months    = 912,
        is_future      = TRUE
    ),
    ssp585 = list(
        path_flo_hist  = paste0(path_flo_base, "historical/"),
        path_flo_gap   = paste0(path_gap_base, "ssp585/"),
        path_flo_ssp   = paste0(path_flo_base, "ssp585/"),
        path_out       = paste0(path_out_base, "ssp585/"),
        start_yr       = 1965,
        time_start     = as.Date("2041-01-31"),
        time_origin    = "2041-01-01",
        time_ref_start = c(1965, 1),
        time_ref_end   = c(2070, 12),
        n_baseline     = 1272,
        n_analysis     = 360,
        skip_months    = 912,
        is_future      = TRUE
    )
)

# ============================================================================
# LOOP — scenarios x models
# ============================================================================

for (scenario_name in names(scenarios)) {
    cfg <- scenarios[[scenario_name]]
    dir.create(cfg$path_out, showWarnings = FALSE, recursive = TRUE)

    print(paste(rep("=", 50), collapse=""))
    print(paste("SCENARIO:", scenario_name))
    print(paste(rep("=", 50), collapse=""))

    for (model in models) {
        print(paste("Processing:", scenario_name, "-", model))

        # ----- Read historical flo_out (1965-2014) -----------------------
        nc_hist <- nc_open(paste0(cfg$path_flo_hist,
                           "flo_masked_", basin, "_historical_", model, ".nc"))
        flo_hist <- ncvar_get(nc_hist, "flo_out")
        lat      <- ncvar_get(nc_hist, "lat")
        lon      <- ncvar_get(nc_hist, "lon")
        nc_close(nc_hist)

        flo_hist <- aperm(flo_hist, c(3, 2, 1))  # → (lon, lat, time)
        flo_hist[flo_hist < -9000] <- NA
        n_hist <- dim(flo_hist)[3]  # 600

        n_lon <- length(lon)
        n_lat <- length(lat)

        if (cfg$is_future) {

            # ----- Read gap flo_out (2015-2019) --------------------------
            nc_gap <- nc_open(paste0(cfg$path_flo_gap,
                              "flo_gap_", basin, "_", scenario_name,
                              "_", model, "_2015_2019.nc"))
            flo_gap <- ncvar_get(nc_gap, "flo_out")
            nc_close(nc_gap)

            flo_gap <- aperm(flo_gap, c(3, 2, 1))  # → (lon, lat, time)
            flo_gap[flo_gap < -9000] <- NA
            n_gap <- dim(flo_gap)[3]  # 60

            # ----- Read SSP flo_out (2020-2070) --------------------------
            nc_ssp <- nc_open(paste0(cfg$path_flo_ssp,
                              "flo_masked_", basin, "_", scenario_name,
                              "_", model, ".nc"))
            flo_ssp <- ncvar_get(nc_ssp, "flo_out")
            nc_close(nc_ssp)

            flo_ssp <- aperm(flo_ssp, c(3, 2, 1))  # → (lon, lat, time)
            flo_ssp[flo_ssp < -9000] <- NA
            n_ssp <- dim(flo_ssp)[3]  # 612

            n_total <- n_hist + n_gap + n_ssp  # 1272
            print(paste("flo_hist:", paste(dim(flo_hist), collapse=" x ")))
            print(paste("flo_gap :", paste(dim(flo_gap),  collapse=" x ")))
            print(paste("flo_ssp :", paste(dim(flo_ssp),  collapse=" x ")))
            print(paste("n_total :", n_total))

        } else {
            n_total <- n_hist  # 600
            print(paste("flo_hist:", paste(dim(flo_hist), collapse=" x ")))
        }

        # ----- SSI grid --------------------------------------------------
        ssi_grid <- array(NA_real_, dim = c(n_lon, n_lat, cfg$n_analysis))

        for (i in 1:n_lon) {
            for (j in 1:n_lat) {

                if (cfg$is_future) {
                    flo_combined <- c(flo_hist[i, j, ],
                                      flo_gap[i, j, ],
                                      flo_ssp[i, j, ])
                } else {
                    flo_combined <- flo_hist[i, j, ]
                }

                # Skip non-river cells
                if (all(is.na(flo_combined))) next

                flo_ts <- ts(flo_combined,
                             start     = c(cfg$start_yr, 1),
                             frequency = 12)

                result <- tryCatch({
                    spi(flo_ts,
                        scale        = 12,
                        ref.start    = cfg$time_ref_start,
                        ref.end      = cfg$time_ref_end,
                        distribution = "Gamma",
                        na.rm        = TRUE
                    )$fitted
                }, error = function(e) {
                    rep(NA_real_, n_total)
                })

                if (length(result) != n_total) result <- rep(NA_real_, n_total)

                ssi_grid[i, j, ] <- as.numeric(result)[
                    (cfg$skip_months + 1):(cfg$skip_months + cfg$n_analysis)
                ]
            }
            if (i %% 5 == 0) print(paste("  lon", i, "of", n_lon, "done"))
        }

        print(paste("SSI completed:", scenario_name, "-", model))

        # ----- Replace non-finite values ---------------------------------
        ssi_grid[is.nan(ssi_grid)]      <- -9999
        ssi_grid[is.infinite(ssi_grid)] <- -9999
        ssi_grid[is.na(ssi_grid)]       <- -9999

        # ----- Write NetCDF ----------------------------------------------
        output_file <- paste0(cfg$path_out,
                       "ssi12_", basin, "_", scenario_name, "_", model, ".nc")

        time_values <- as.numeric(
            as.Date(seq(cfg$time_start, by = "month",
                        length.out = cfg$n_analysis)) -
            as.Date(cfg$time_origin)
        )

        lon_dim  <- ncdim_def("lon",  "degrees_east",  lon)
        lat_dim  <- ncdim_def("lat",  "degrees_north", lat)
        time_dim <- ncdim_def("time", paste("days since", cfg$time_origin),
                              time_values)

        ssi_var  <- ncvar_def(
            name     = "SSI_12",
            units    = "-",
            dim      = list(lon_dim, lat_dim, time_dim),
            missval  = -9999,
            longname = paste0(
                "SSI-12 (R SPEI package, spi(), Gamma, na.rm=TRUE, ",
                "baseline ", cfg$time_ref_start[1], "-", cfg$time_ref_end[1],
                ", mask Strahler s6)"
            )
        )

        nc_out <- nc_create(output_file, ssi_var)
        ncvar_put(nc_out, ssi_var, ssi_grid)
        nc_close(nc_out)

        print(paste("NetCDF saved:", output_file))
    }
}

print("Done — Orinoco historical + SSP126 + SSP585 SSI-12.")