# =============================================================================
# scripts/01_data/08_descriptive_stats.R
# Purpose : Descriptive statistics of the analysis sample (before clustering):
#           variable inventory, age x sex, ethnic composition vs the 2017 U.S.
#           population, religiosity, per-question distributions.
# Thesis  : Section 2.3, Figure 2.1 (a) ethnic composition, (b) age x sex
# Input   : data/processed/imputed_labeled.rds
# Output  : results/figures/01_data/{age_sex_distribution,ethnic_composition,
#           religious_identification}.png
#           results/tables/01_data/{variable_inventory,category_counts,
#           distribution_all_questions}.csv
# =============================================================================

source(here::here("R", "00_setup.R"))

PLOTS_DIR <- ensure_dir(DIR_FIG, "01_data")
TABLE_DIR <- ensure_dir(DIR_TAB, "01_data")

#### 1. LOAD DATA ####

wvs <- readRDS(file.path(DIR_DATA, "imputed_labeled.rds"))

stopifnot(is.data.frame(wvs), nrow(wvs) > 2000)
cat("Dataset loaded:", nrow(wvs), "rows x", ncol(wvs), "columns\n\n")


#### 2. INVENTORY: NUMBER OF QUESTIONS BY NUMBER OF CATEGORIES ####

# Label each variable with its block
assign_block <- function(v) {
  if (exists("Demographics") && v %in% Demographics) return("Sociodemographic")
  if (exists("Political")    && v %in% Political)    return("Political")
  "Value"
}

var_inventory <- tibble(
  variable = names(wvs),
  class    = map_chr(wvs, ~ class(.x)[1]),
  n_levels = map_int(wvs, ~ if (is.factor(.x)) nlevels(.x) else NA_integer_),
  ordered  = map_lgl(wvs, is.ordered),
  block    = map_chr(names(wvs), assign_block)
)

# Count by number of categories (categorical variables only).
# mclust is detached because it masks purrr::map().
detach("package:mclust", unload = TRUE)
cat_counts <- var_inventory %>%
  filter(!is.na(n_levels)) %>%
  count(n_levels, name = "n_variables") %>%
  mutate(
    type = case_when(
      n_levels == 2         ~ "Binary",
      n_levels %in% 3:5     ~ "Short Likert",
      n_levels %in% 6:8     ~ "Medium Likert",
      n_levels >= 9         ~ "10-point / long scale",
      TRUE                  ~ "Other"
    ),
    pct = round(n_variables / sum(n_variables) * 100, 1)
  ) %>%
  arrange(n_levels)

cat("--- Variables by number of categories ---\n")
print(cat_counts, n = Inf)

# Same, by block
cat("\n--- Number of categories x block ---\n")
print(
  var_inventory %>%
    filter(!is.na(n_levels)) %>%
    count(block, n_levels) %>%
    pivot_wider(names_from = block, values_from = n, values_fill = 0),
  n = Inf
)

# Non-categorical (numeric) variables
cat("\n--- Numeric (non-categorical) variables ---\n")
print(var_inventory %>% filter(is.na(n_levels)))

# Ordered vs nominal (relevant for the Gower metric)
cat("\n--- Ordered vs unordered ---\n")
print(var_inventory %>% filter(!is.na(n_levels)) %>% count(block, ordered))


#### 3. JOINT DISTRIBUTION OF AGE x SEX (Figure 2.1b) ####

df_age <- wvs %>%
  mutate(
    Q262_num  = as.numeric(as.character(Q262)),
    age_group = cut(
      Q262_num,
      breaks         = c(0, 20, 30, 40, 50, 60, 70, Inf),
      labels         = c("0-20", "20-30", "30-40", "40-50",
                         "50-60", "60-70", "70+"),
      right          = FALSE,
      include.lowest = TRUE
    ),
    sex = Q260
  ) %>%
  filter(!is.na(age_group), !is.na(sex)) %>%
  count(age_group, sex, name = "n")

n_tot <- sum(df_age$n)
stopifnot(n_tot > 0)

# Sex share WITHIN each age group
df_age <- df_age %>%
  group_by(age_group) %>%
  mutate(
    n_age      = sum(n),
    pct_within = n / n_age * 100
  ) %>%
  ungroup()

# Share of each age group in the whole sample (label above the bars)
df_age_tot <- df_age %>%
  group_by(age_group) %>%
  summarise(
    n_age = first(n_age),
    y_lab = max(n),
    .groups = "drop"
  ) %>%
  mutate(pct_age = round(n_age / n_tot * 100))

p_age_sex <- ggplot(df_age, aes(x = age_group, y = n, fill = sex)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  # share of men / women within the age group
  geom_text(
    aes(label = sprintf("%.1f%%", pct_within)),
    position = position_dodge(width = 0.8),
    vjust = 1.5, size = 2, colour = "white", fontface = "bold"
  ) +
  # share of the age group in the whole sample
  geom_text(
    data = df_age_tot,
    aes(x = age_group, y = y_lab, label = sprintf("%.f%%", pct_age)),
    inherit.aes = FALSE,
    vjust = -0.9, size = 3.5, colour = "grey25"
  ) +
  scale_y_continuous(
    name     = "Number of respondents",
    expand   = expansion(mult = c(0, 0.14))
  ) +
  scale_fill_manual(values = c("Male" = "#4E79A7", "Female" = "#E15759")) +
  labs(
    x     = "Age group",
    fill  = NULL,
    title = "Joint distribution of age and sex",
    subtitle = "Labels inside bars: sex composition within each age group. "
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "top",
    panel.grid.major.x = element_blank()
  )

print(p_age_sex)


#### 4. ETHNIC GROUPS (Q290): SAMPLE vs U.S. POPULATION 2017 (Figure 2.1a) ####

df_eth <- wvs %>%
  filter(!is.na(Q290)) %>%
  count(Q290, name = "n") %>%
  mutate(pct_sample = n / sum(n) * 100) %>%
  arrange(desc(n))

cat("\n--- Ethnic groups in the sample ---\n")
print(df_eth)

# Source: NCES Digest of Education Statistics 2021, Table 101.20
# (Census Bureau resident-population estimates), 2017 row.
# Racial categories exclude Hispanics, who are reported separately.
us_2017 <- tibble::tibble(
  group      = c("White", "Black", "Hispanic", "Asian",
                 "Pacific Islander", "American Indian/Alaska Native",
                 "Two or more races"),
  n_thousand = c(197798, 40589, 58403, 18348, 582, 2398, 7004),
  pct_census = c(60.8, 12.5, 18.0, 5.6, 0.2, 0.7, 2.2)
)

# --- WVS -> Census categories ------------------------------------------------
# Q290 has no Hispanic category and separates South/East Asian, which the
# Census merges into "Asian". We aggregate accordingly.
map_wvs_to_census <- c(
  "White"                = "White",
  "Black"                = "Black",
  "South Asian"          = "Asian",
  "East Asian"           = "Asian",
  "Arabic/Central Asian" = "Not in Census scheme"
)

df_eth_sample <- wvs %>%
  dplyr::filter(!is.na(Q290)) %>%
  dplyr::mutate(group = dplyr::recode(as.character(Q290), !!!map_wvs_to_census)) %>%
  dplyr::count(group, name = "n") %>%
  dplyr::mutate(pct_sample = n / sum(n) * 100)

# --- Join and comparability flag ---------------------------------------------
df_eth_cmp <- dplyr::full_join(df_eth_sample, us_2017, by = "group") %>%
  dplyr::mutate(
    comparable = !is.na(pct_sample) & !is.na(pct_census) &
      group != "Not in Census scheme",
    diff       = pct_sample - pct_census
  ) %>%
  dplyr::arrange(dplyr::desc(dplyr::coalesce(pct_census, pct_sample)))

cat("\n--- Ethnic composition: sample vs U.S. population 2017 ---\n")
print(df_eth_cmp)

# --- Plot --------------------------------------------------------------------
df_eth_long <- df_eth_cmp %>%
  dplyr::select(group, Sample = pct_sample, `U.S. population` = pct_census) %>%
  tidyr::pivot_longer(-group, names_to = "source", values_to = "pct") %>%
  # complete the grid: every group has both series (possibly NA)
  tidyr::complete(group, source)

lev_order <- df_eth_cmp$group[order(
  dplyr::coalesce(df_eth_cmp$pct_census, df_eth_cmp$pct_sample)
)]
df_eth_long$group  <- factor(df_eth_long$group, levels = lev_order)
df_eth_long$source <- factor(df_eth_long$source,
                             levels = c("Sample", "U.S. population"))

pd <- position_dodge2(width = 0.8, padding = 0.1)

df_eth_long <- df_eth_long %>%
  dplyr::filter(group %in% c("White", "Black", "Asian"))
p_eth <- ggplot(df_eth_long, aes(x = group, y = pct, fill = source)) +
  geom_col(position = pd, width = 0.75, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(is.na(pct), "", sprintf("%.1f%%", pct))),
    position = pd,
    hjust = -0.2, size = 3, na.rm = TRUE
  ) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  scale_fill_manual(
    values = c("Sample" = "#4E79A7", "U.S. population" = "#B07AA1"),
    drop = FALSE
  ) +
  labs(
    x = NULL, y = "% of population", fill = NULL,
    title    = "Ethnic composition: WVS sample vs U.S. population (2017)",
    "WVS South/East Asian aggregated into 'Asian'."
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

print(p_eth)


#### 5. RELIGIOSITY (Q173 self-identification, Q165 belief in God, Q289) ####

# Q173: "would you say you are a religious person / not religious / atheist"
df_q173 <- wvs %>%
  filter(!is.na(Q173)) %>%
  count(Q173, name = "n") %>%
  mutate(pct = n / sum(n) * 100)

df_Q289 <- wvs %>%
  filter(!is.na(Q289)) %>%
  count(Q289, name = "n") %>%
  mutate(pct = n / sum(n) *100)
cat("\n--- Q173: religious self-identification ---\n")
print(df_q173)
print(df_Q289)

# Q165: belief in God (binary yes/no), a different measure from Q173
df_q165 <- wvs %>%
  filter(!is.na(Q165)) %>%
  count(Q165, name = "n") %>%
  mutate(pct = n / sum(n) * 100)

cat("\n--- Q165: belief in God ---\n")
print(df_q165)

# Cross-tabulation: the two measures do not coincide
cat("\n--- Q173 x Q165 (row percentages) ---\n")
print(round(prop.table(table(wvs$Q173, wvs$Q165), margin = 1) * 100, 1))

p_q173 <- ggplot(df_q173, aes(x = reorder(Q173, -pct), y = pct)) +
  geom_col(fill = "#59A14F", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 3.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    x = NULL, y = "% of respondents",
    title = "Religious self-identification (Q173)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank())

print(p_q173)


#### 6. DISTRIBUTION OF EACH QUESTION ####

dist_all <- wvs %>%
  select(where(is.factor)) %>%
  pivot_longer(
    everything(),
    names_to  = "variable",
    values_to = "level",
    values_transform = list(level = as.character)
  ) %>%
  count(variable, level, name = "n") %>%
  group_by(variable) %>%
  mutate(
    pct       = n / sum(n) * 100,
    n_levels  = n_distinct(level),
    modal_pct = max(pct)
  ) %>%
  ungroup() %>%
  left_join(var_inventory %>% select(variable, block), by = "variable")

# Convenience function to inspect a single question, e.g. show_dist("Q184")
show_dist <- function(v) {
  dist_all %>%
    filter(variable == v) %>%
    select(variable, level, n, pct) %>%
    arrange(desc(n)) %>%
    print(n = Inf)
}

# Most unbalanced questions (one category dominates)
cat("\n--- 15 most unbalanced questions (modal share) ---\n")
print(
  dist_all %>%
    distinct(variable, block, n_levels, modal_pct) %>%
    arrange(desc(modal_pct)) %>%
    slice_head(n = 15),
  n = Inf
)

# Most balanced questions (potentially the most informative)
cat("\n--- 15 most balanced questions ---\n")
print(
  dist_all %>%
    distinct(variable, block, n_levels, modal_pct) %>%
    arrange(modal_pct) %>%
    slice_head(n = 15),
  n = Inf
)


#### 7. SAVE ####

ggsave(file.path(PLOTS_DIR, "age_sex_distribution.png"),
       p_age_sex, width = 10, height = 6, dpi = 300, bg = "white")

ggsave(file.path(PLOTS_DIR, "ethnic_composition.png"),
       p_eth, width = 9, height = 5.5, dpi = 300, bg = "white")

ggsave(file.path(PLOTS_DIR, "religious_identification.png"),
       p_q173, width = 8, height = 5, dpi = 300, bg = "white")

write.csv(var_inventory,
          file.path(TABLE_DIR, "variable_inventory.csv"), row.names = FALSE)
write.csv(cat_counts,
          file.path(TABLE_DIR, "category_counts.csv"), row.names = FALSE)
write.csv(dist_all,
          file.path(TABLE_DIR, "distribution_all_questions.csv"), row.names = FALSE)
