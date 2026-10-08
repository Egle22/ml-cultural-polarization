# =============================================================================
# scripts/05_recovery/03_geography.R
# Purpose : Build the labelled geographic data used by the "geo" predictor
#           blocks of 04_recovery_all_k.R (town size, urban/rural, state from
#           the interview coordinates), with missing-value diagnostics, and
#           draw choropleths of the dominant cluster by state.
# Thesis  : not reported (geography robustness of Section 4.3)
# Input   : data/processed/geo_raw.rds (raw geographic variables, same rows as
#           the analysis sample), data/processed/imputed_labeled.rds,
#           results/objects/clusters/all_clusters_aligned.rds
# Output  : results/objects/Geo/geo_usa_labeled.rds
#           results/figures/05_recovery/geo/*.png
# Note    : requires internet access on the first run (tigris downloads the
#           state boundaries; cached afterwards). The spatial join uses the
#           UNSHIFTED state boundaries; tigris::shift_geometry() (which moves
#           Alaska and Hawaii) is applied only to the plotting object.
#           This script never writes all_clusters_aligned.rds.
# =============================================================================

source(here::here("R", "00_setup.R"))
library(sf)
library(tigris); options(tigris_use_cache = TRUE)

objects_dir <- ensure_dir(DIR_OBJ, "Geo")
plots_dir   <- ensure_dir(DIR_FIG, "05_recovery", "geo")

# -- 1. Data ---------------------------------------------------------------------
geo_raw         <- readRDS(file.path(DIR_DATA, "geo_raw.rds"))
clusters        <- readRDS(PATH_CLUSTERS)
imputed_labeled <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds"))

stopifnot(nrow(geo_raw) == 2579)
stopifnot(nrow(clusters)        == nrow(geo_raw))
stopifnot(nrow(imputed_labeled) == nrow(geo_raw))
stopifnot("Q223" %in% names(imputed_labeled))
stopifnot(identical(as.integer(clusters$D_INTERVIEW), as.integer(geo_raw$D_INTERVIEW)))

# Party voted for (Q223) attached by position: imputed_labeled has the same
# rows, in the same order, as geo_raw and the cluster file.
q223_df <- data.frame(D_INTERVIEW = as.integer(geo_raw$D_INTERVIEW),
                      Q223        = imputed_labeled$Q223)
cat("Q223 (imputed, labelled):\n"); print(table(q223_df$Q223, useNA = "ifany"))

# -- 2. Missing-code handling (-1, -2, -3, -5 -> NA) -------------------------------
to_na <- function(x) { x <- suppressWarnings(as.numeric(x)); x[x %in% c(-1,-2,-3,-5)] <- NA; x }

# -- 3. Coordinates and sanity checks ----------------------------------------------
lon <- to_na(geo_raw$O1_LONGITUDE)
lat <- to_na(geo_raw$O2_LATITUDE)

cat("Non-missing coordinates:", sum(!is.na(lon) & !is.na(lat)), "/", nrow(geo_raw), "\n")
if (sum(!is.na(lon) & !is.na(lat)) == 0L)
  stop("All coordinates are NA: reverse geocoding impossible, use N_REGION_ISO.")

cat("Range lon:", paste(round(range(lon, na.rm = TRUE), 2), collapse = " / "), "\n")
cat("Range lat:", paste(round(range(lat, na.rm = TRUE), 2), collapse = " / "), "\n")

# Sign fix: U.S. longitudes are negative. If they are positive, negate them.
if (median(lon, na.rm = TRUE) > 0) { lon <- -abs(lon); cat("[fix] longitudes negated\n") }

# Flag points outside the U.S. bounding box (Alaska/Hawaii included -> wide limits)
in_box <- lon >= -180 & lon <= -66 & lat >= 18 & lat <= 72
cat("Points outside the U.S. bounding box:", sum(!in_box, na.rm = TRUE), "\n")

cat("lon == -999:", sum(lon == -999, na.rm = TRUE), "\n")
cat("lat == -999:", sum(lat == -999, na.rm = TRUE), "\n")
cat("points outside the box that are NOT -999:",
    sum(!in_box & lon != -999 & lat != -999, na.rm = TRUE), "\n")

# Distribution of the negative (missing sentinel) values for each geo variable
geo_vars <- c("G_TOWNSIZE","G_TOWNSIZE2","H_URBRURAL","H_SETTLEMENT",
              "N_REGION_ISO","N_TOWN","O1_LONGITUDE","O2_LATITUDE")

for (v in geo_vars) {
  x <- suppressWarnings(as.numeric(geo_raw[[v]]))
  neg <- x[x < 0 & !is.na(x)]
  cat("\n---", v, "---\n")
  cat("  n negative:", length(neg), "/", length(x), "\n")
  if (length(neg) > 0) print(table(neg))
}

# -- 4. State boundaries -----------------------------------------------------------
# Join object: UNSHIFTED boundaries (true coordinates).
states_sf <- tigris::states(cb = TRUE, resolution = "20m", year = 2021) |>
  dplyr::select(region_name = NAME, region_iso = STUSPS) |>
  sf::st_transform(4326)

# Plotting object only: Alaska and Hawaii moved below the contiguous states.
states_plot <- tigris::shift_geometry(states_sf) |>
  sf::st_transform(4326)

# -- 5. Point-in-polygon join ------------------------------------------------------
pts <- data.frame(D_INTERVIEW = as.integer(geo_raw$D_INTERVIEW),
                  lon = lon, lat = lat, row_id = seq_len(nrow(geo_raw)))

pts_valid <- pts[!is.na(pts$lon) & !is.na(pts$lat), ]
pts_sf <- sf::st_as_sf(pts_valid, coords = c("lon","lat"), crs = 4326, remove = FALSE)

joined <- sf::st_join(pts_sf, states_sf, join = sf::st_intersects) |>
  sf::st_drop_geometry()

cat("Points not assigned to any state:",
    sum(is.na(joined$region_name)), "/", nrow(joined), "\n")

# -- 6. Labelled geographic data frame ---------------------------------------------
geo <- data.frame(
  D_INTERVIEW = as.integer(geo_raw$D_INTERVIEW),
  town_size   = factor(to_na(geo_raw$G_TOWNSIZE), levels = 1:8,
                       labels = c("Under 2,000","2,000-5,000","5,000-10,000",
                                  "10,000-20,000","20,000-50,000","50,000-100,000",
                                  "100,000-500,000","500,000 and more")),
  urban_rural = factor(to_na(geo_raw$H_URBRURAL), levels = 1:2,
                       labels = c("Urban","Rural")),
  longitude   = lon,
  latitude    = lat,
  stringsAsFactors = FALSE
)
geo <- dplyr::left_join(geo,
                        joined[, c("D_INTERVIEW","region_name","region_iso")],
                        by = "D_INTERVIEW")

cat("\nRegion coverage:\n"); print(sort(table(geo$region_name), decreasing = TRUE))

saveRDS(geo, file.path(objects_dir, "geo_usa_labeled.rds"))

# -- 7. Missingness diagnostics: do the NAs fall on the same rows? ------------------
lon_raw <- suppressWarnings(as.numeric(geo_raw$O1_LONGITUDE))
lat_raw <- suppressWarnings(as.numeric(geo_raw$O2_LATITUDE))

miss_town   <- is.na(geo$town_size)                       # G_TOWNSIZE (-5 etc.)
miss_region <- is.na(geo$region_name)                     # failed join (-999)
miss_coord  <- lon_raw <= -900 | lat_raw <= -900          # coordinate sentinel

cat("Total interviews:", nrow(geo_raw), "\n\n")
cat("town_size missing :", sum(miss_town),   "\n")
cat("region missing    :", sum(miss_region), "\n")
cat("coord = -999       :", sum(miss_coord), "\n\n")

cat("Coord -999 that also have region NA:", sum(miss_coord & miss_region), "\n")
cat("Coord -999 but with a region        :", sum(miss_coord & !miss_region), "\n")
cat("Region NA but with valid coordinates:", sum(!miss_coord & miss_region), "\n\n")

patterns <- data.frame(
  town   = ifelse(miss_town,   "town_NA",   "town_ok"),
  region = ifelse(miss_region, "region_NA", "region_ok")
)
cat("Joint patterns town_size x region:\n")
print(table(patterns$town, patterns$region))

cat("\ntown_size NA with VALID coordinates:", sum(miss_town & !miss_coord), "\n")
cat("town_size NA with coord -999       :", sum(miss_town &  miss_coord), "\n")

# -- 8. Merge with the clusters (by D_INTERVIEW, not by position) --------------------
df <- geo |>
  dplyr::inner_join(clusters, by = "D_INTERVIEW") |>
  dplyr::inner_join(q223_df,  by = "D_INTERVIEW")
stopifnot(nrow(df) == nrow(clusters))          # fail fast: no row lost

# -- 9. Choropleth of the dominant cluster by state ------------------------------------
MIN_N        <- 20   # states with fewer respondents are flagged "n < 20"
small_states <- c("Rhode Island","Delaware","Connecticut","New Jersey",
                  "New Hampshire","Vermont","Massachusetts","Maryland",
                  "District of Columbia")

# One map for column `col`: dominant category, its share, and state size.
make_choropleth <- function(col, palette = NULL, legend_name = "Cluster",
                            title = col, min_n = MIN_N) {

  thr_lab <- sprintf("n < %d", min_n)

  # dominant category + share + n_total by state
  sm <- df |>
    dplyr::filter(!is.na(region_name), !is.na(.data[[col]])) |>
    dplyr::transmute(region_name, cat = as.character(.data[[col]])) |>
    dplyr::add_count(region_name, name = "n_total") |>
    dplyr::count(region_name, cat, n_total, name = "n") |>
    dplyr::group_by(region_name) |>
    dplyr::mutate(share = n / sum(n)) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      cat   = ifelse(n_total < min_n, thr_lab, cat),
      share = ifelse(n_total < min_n, NA_real_, share)
    )

  # map + centroids + labels (shifted boundaries: plotting only)
  mm <- dplyr::left_join(states_plot, sm, by = "region_name")
  xy <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(mm)))
  mm$lab_x <- xy[, 1]; mm$lab_y <- xy[, 2]
  mm$lab   <- ifelse(is.na(mm$share), "", paste0(round(mm$share * 100), "%"))
  mm$is_small <- mm$region_name %in% small_states & !is.na(mm$share)

  # palette: the one passed, otherwise automatic; always add "n < min"
  levs <- setdiff(sort(unique(stats::na.omit(mm$cat))), thr_lab)
  pal  <- if (is.null(palette)) setNames(scales::hue_pal()(length(levs)), levs)
          else palette
  pal[thr_lab] <- "grey85"

  # legend table with the size of each category
  cnt <- df |>
    dplyr::filter(!is.na(.data[[col]])) |>
    dplyr::count(cat = as.character(.data[[col]]), name = "n")
  legend_text <- paste(sprintf("%s: n = %d", cnt$cat, cnt$n), collapse = "\n")

  ggplot(mm) +
    geom_sf(aes(fill = cat), color = "white", linewidth = 0.2) +
    geom_text(data = subset(mm, !is_small & lab != ""),
              aes(lab_x, lab_y, label = lab), size = 2.3, color = "black") +
    ggrepel::geom_text_repel(
      data = subset(mm, is_small),
      aes(lab_x, lab_y, label = lab), size = 2.3, color = "black",
      min.segment.length = 0, segment.size = 0.3, segment.color = "grey50",
      box.padding = 0.5, force = 3, nudge_x = 6, direction = "y", hjust = 0) +
    annotate("label", x = -Inf, y = -Inf, hjust = 0, vjust = 0,
             label = legend_text, size = 2.6, label.size = 0.3,
             fill = "white", alpha = 0.8) +
    scale_fill_manual(values = pal, na.value = "grey70", name = legend_name) +
    labs(title = title) +
    theme_void()
}

# Target columns: KM-PC2, KM-PC5, KM-UMAP (k = 2..6)
cluster_cols <- grep("pc2|pc5|umap", names(df), value = TRUE, ignore.case = TRUE)
cat("Columns found:", length(cluster_cols), "\n"); print(cluster_cols)

for (col in cluster_cols) {
  p <- make_choropleth(col, palette = CLUSTER_COLORS,
                       legend_name = "Cluster",
                       title = paste0("Dominant cluster by state — ", col))
  fname <- gsub("[^A-Za-z0-9]+", "_", col)
  ggsave(file.path(plots_dir, paste0("geo_choropleth_", fname, ".png")),
         p, width = 9, height = 6, dpi = 300)
}

# Q223 separately (automatic palette: parties are not clusters)
p_q223 <- make_choropleth("Q223", palette = NULL,
                          legend_name = "Party (Q223)",
                          title = "Dominant party by state — Q223")
ggsave(file.path(plots_dir, "geo_choropleth_Q223.png"),
       p_q223, width = 9, height = 6, dpi = 300)

cat("Saved", length(cluster_cols) + 1, "maps in", plots_dir, "\n")
