# CLAUDE.md

Replication code for the Master's thesis "Measuring cultural polarization across
social groups: a machine learning classification approach" (WVS Wave 7, U.S.
sample, R). See README.md for the pipeline and the script-to-thesis map.

## Working with the author

- Talk to the author in Italian. Everything written to the repository (code,
  comments, README, commit messages) is in English.
- The author is new to git and the terminal: explain each step in plain
  language before doing it.

## Rules that must not be broken

- Never `ordered = TRUE` in `factor()` for variables entering LASSO or logit:
  it produces polynomial contrasts and silently breaks the models.
  (`fct_ord()` in R/labels.R creates ordered factors for the labelled dataset
  only; every modelling script converts them to unordered factors or numeric
  scores before glm/glmnet.)
- `tigris::shift_geometry()` only on plotting objects, never on the spatial
  join object.
- Map `model.matrix` columns to variables with the `assign` attribute, not with
  regular expressions on the column names.
- `results/objects/clusters/all_clusters_aligned.rds` (`PATH_CLUSTERS`) is the
  canonical source of cluster labels. When a file exists in both
  `results/objects/` and `results/objects/clusters/`, `clusters/` wins.
  Only `scripts/02_clustering/07_align_clusters.R` may write it. It has 37
  columns: k = 2..6 for KM-PC{2,5,15,75,198}, KM-UMAP and PAM,
  `hdb_umap2_mp25` and `D_INTERVIEW`. The previous 59-column version is kept in
  `archive/results_objects_original/`.
- Sequential `lapply`, no future / parallel. Documented exception:
  `scripts/05_recovery/01_recovery_k3.R` and `04_recovery_all_k.R` keep
  `future` (multisession, `future.seed = TRUE`) so that the thesis numbers are
  reproduced exactly. Do not add parallelism anywhere else.
- Scripts must stay runnable piece by piece in RStudio. Every script starts
  with `source(here::here("R", "00_setup.R"))`; paths come from the `DIR_*`
  constants defined there (never absolute paths, never `setwd()`).
- Behaviour-preserving changes only: do not change statistical logic, seeds,
  parameters or model specifications without asking. List suspected bugs
  instead of fixing them (see README, Known limitations).

## Data protection

- WVS data cannot be redistributed. Nothing respondent-level may be committed:
  no raw or processed data, no `.rds` objects, no cluster labels, embeddings or
  distance matrices, and no HTML widgets that plot one point per respondent.
  `data/`, `results/objects/` and `results/figures/` are git-ignored; `docs/`
  publishes only the aggregate codebook (`docs/index.html`).
- Before every commit, check `git ls-files` for `.rds`, `.RData`, `.csv` data
  files and HTML widgets.

## Git setup

- The public history lives on the local branch `public-main`, pushed only to
  the remote `public` (`remote.public.push = refs/heads/public-main:refs/heads/main`).
- The old local branch `main` (remote `origin`, private) contains the raw WVS
  file in its history. NEVER check it out again from this folder: git would
  rewrite the working tree back to the old layout. Never push it to `public`.
- Never run `git restore`, `git checkout -- <path>`, `git reset --hard`,
  `git clean` or `git stash`. Never force-push. Never change repository
  visibility (the author does it).
- `archive/` holds superseded scripts and outputs (git-ignored, never delete).
