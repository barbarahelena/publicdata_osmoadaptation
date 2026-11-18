#!/usr/bin/env Rscript
# Select curatedMetagenomicData samples for catalogue-based osmoadaptation profiling.
#
# Builds three sample sets, prints the flow at every step, and writes the ENA/SRA
# accessions per set. Tune the thresholds at the top and rerun; the flow table is
# what goes in the supplement.
#
#   Rscript select_cmd_samples.R [outdir]

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(dplyr)
  library(tidyr)
  library(readr)
})

OUT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "results/cmd_selection"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- thresholds ------------------------------------------------------------
MIN_READS        <- 5e6     # raw reads; below this low-abundance families go to zero
MIN_READ_LEN     <- 75      # median read length, bp; 97% identity over 45 bp needs length
TARGET_DEPTH     <- 5e6     # read pairs to subsample to, downstream
MIN_PER_STUDY    <- 20      # samples per study after filtering
MIN_STUDIES      <- 2       # studies per stratum
PLATFORMS        <- c("IlluminaHiSeq", "IlluminaMiSeq", "IlluminaNovaSeq", "Illumina")

# ---- flow bookkeeping ------------------------------------------------------
flow <- tibble(set = character(), step = character(), n_samples = integer(),
               n_subjects = integer(), n_studies = integer())

note <- function(d, set, step) {
  force(d)  # evaluate upstream pipe steps first, or their flow rows get overwritten
  flow <<- bind_rows(flow, tibble(
    set = set, step = step,
    n_samples = nrow(d),
    n_subjects = dplyr::n_distinct(d$subject_id),
    n_studies = dplyr::n_distinct(d$study_name)))
  d
}

md <- sampleMetadata %>% as_tibble()
message("curatedMetagenomicData: ", nrow(md), " samples, ",
        dplyr::n_distinct(md$study_name), " studies")

# ---- shared technical filter ----------------------------------------------
base <- md %>%
  note("base", "all samples") %>%
  filter(body_site == "stool") %>%
  note("base", "stool only") %>%
  filter(!is.na(number_reads), number_reads >= MIN_READS) %>%
  note("base", sprintf("reads >= %.0e", MIN_READS)) %>%
  filter(!is.na(median_read_length), median_read_length >= MIN_READ_LEN) %>%
  note("base", sprintf("median read length >= %d bp", MIN_READ_LEN)) %>%
  filter(is.na(sequencing_platform) |
           grepl("Illumina", sequencing_platform, ignore.case = TRUE)) %>%
  note("base", "Illumina") %>%
  filter(!is.na(NCBI_accession)) %>%
  note("base", "has accession")

# study-level floor, applied after every stratum filter
apply_study_floor <- function(d, set) {
  d %>%
    group_by(study_name) %>% filter(dplyr::n() >= MIN_PER_STUDY) %>% ungroup() %>%
    note(set, sprintf("studies with >= %d samples", MIN_PER_STUDY))
}

# one sample per subject: keep the deepest
one_per_subject <- function(d, set) {
  d %>%
    group_by(subject_id) %>% slice_max(number_reads, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    note(set, "one sample per subject")
}

healthy <- function(d, set) {
  d %>%
    filter(is.na(study_condition) | study_condition == "control") %>%
    note(set, "controls only") %>%
    filter(is.na(antibiotics_current_use) | antibiotics_current_use == "no") %>%
    note(set, "no current antibiotics") %>%
    filter(is.na(pregnant) | pregnant == "no", is.na(lactating) | lactating == "no") %>%
    note(set, "not pregnant or lactating")
}

# ---- set 1: Bristol stool score -------------------------------------------
bristol <- base %>%
  filter(!is.na(bristol_score)) %>%
  note("bristol", "has bristol_score") %>%
  healthy("bristol") %>%
  one_per_subject("bristol") %>%
  apply_study_floor("bristol")

# ---- set 2: life course, longitudinal infants kept intact -----------------
# studies with repeated sampling of the same subject
longitudinal_studies <- base %>%
  group_by(study_name, subject_id) %>%
  summarise(n_visits = dplyr::n(), .groups = "drop") %>%
  group_by(study_name) %>%
  summarise(subjects_with_repeats = sum(n_visits > 1), .groups = "drop") %>%
  filter(subjects_with_repeats >= 10) %>%
  pull(study_name)

lifecourse <- base %>%
  note("lifecourse", "all stool, technical pass") %>%
  filter(is.na(antibiotics_current_use) | antibiotics_current_use == "no") %>%
  note("lifecourse", "no current antibiotics") %>%
  filter(is.na(study_condition) | study_condition == "control") %>%
  note("lifecourse", "controls only") %>%
  apply_study_floor("lifecourse")

lifecourse_long <- lifecourse %>%
  filter(study_name %in% longitudinal_studies) %>%
  note("lifecourse_long", "longitudinal studies") %>%
  group_by(subject_id) %>% filter(dplyr::n() > 1) %>% ungroup() %>%
  note("lifecourse_long", "subjects with >1 sample")

# ---- set 3: population comparison, matched on technical covariates --------
# match non-westernized and westernized within read-length and platform bins
population <- base %>%
  filter(!is.na(non_westernized)) %>%
  note("population", "has non_westernized") %>%
  healthy("population") %>%
  one_per_subject("population") %>%
  mutate(len_bin = cut(median_read_length,
                       breaks = c(0, 99, 125, 151, Inf),
                       labels = c("75-99", "100-125", "126-151", ">151"))) %>%
  group_by(len_bin) %>%
  filter(dplyr::n_distinct(non_westernized) == 2) %>%
  ungroup() %>%
  note("population", "read-length bins holding both groups") %>%
  apply_study_floor("population")

# ---- stratum coverage checks ----------------------------------------------
coverage <- function(d, var, set) {
  if (!var %in% names(d)) return(invisible(NULL))
  d %>%
    filter(!is.na(.data[[var]])) %>%
    count(.data[[var]], name = "n_samples") %>%
    left_join(d %>% filter(!is.na(.data[[var]])) %>%
                group_by(.data[[var]]) %>%
                summarise(n_studies = dplyr::n_distinct(study_name), .groups = "drop"),
              by = var) %>%
    mutate(set = set, variable = var) %>%
    rename(level = !!var) %>%
    mutate(level = as.character(level),
           usable = n_studies >= MIN_STUDIES)
}

checks <- bind_rows(
  coverage(bristol, "bristol_score", "bristol"),
  coverage(bristol, "gender", "bristol"),
  coverage(lifecourse, "age_category", "lifecourse"),
  coverage(lifecourse, "feeding_practice", "lifecourse"),
  coverage(lifecourse_long, "age_category", "lifecourse_long"),
  coverage(population, "non_westernized", "population"),
  coverage(population, "country", "population")
)

# BMI and age are continuous: report completeness instead
continuous <- bind_rows(
  tibble(set = "bristol", variable = "BMI",
         complete = sum(!is.na(bristol$BMI)), total = nrow(bristol)),
  tibble(set = "bristol", variable = "age",
         complete = sum(!is.na(bristol$age)), total = nrow(bristol)),
  tibble(set = "population", variable = "BMI",
         complete = sum(!is.na(population$BMI)), total = nrow(population)),
  tibble(set = "population", variable = "age",
         complete = sum(!is.na(population$age)), total = nrow(population))
) %>% mutate(fraction = round(complete / total, 3))

# ---- download volume ------------------------------------------------------
sets <- list(bristol = bristol, lifecourse_long = lifecourse_long,
             population = population)

volume <- do.call(bind_rows, lapply(names(sets), function(nm) {
  d <- sets[[nm]]
  tibble(set = nm, n_samples = nrow(d),
         # ~1 byte per base gzipped at roughly 3x compression of 2 bits/base
         est_download_GB = round(sum(pmin(d$number_reads, TARGET_DEPTH * 2) *
                                       d$median_read_length / 3) / 1e9, 1))
}))

# ---- write ----------------------------------------------------------------
for (nm in names(sets)) {
  d <- sets[[nm]]
  write_tsv(d, file.path(OUT, sprintf("samples_%s.tsv", nm)))
  d %>%
    select(study_name, sample_id, subject_id, NCBI_accession) %>%
    separate_rows(NCBI_accession, sep = ";") %>%
    mutate(NCBI_accession = trimws(NCBI_accession)) %>%
    filter(NCBI_accession != "") %>%
    write_tsv(file.path(OUT, sprintf("accessions_%s.tsv", nm)))
}
write_tsv(flow, file.path(OUT, "sample_flow.tsv"))
write_tsv(checks, file.path(OUT, "stratum_coverage.tsv"))
write_tsv(continuous, file.path(OUT, "continuous_completeness.tsv"))
write_tsv(volume, file.path(OUT, "download_volume.tsv"))

cat("\n== sample flow ==\n");            print(as.data.frame(flow), row.names = FALSE)
cat("\n== stratum coverage ==\n");       print(as.data.frame(checks), row.names = FALSE)
cat("\n== continuous completeness ==\n");print(as.data.frame(continuous), row.names = FALSE)
cat("\n== estimated download ==\n");     print(as.data.frame(volume), row.names = FALSE)
cat(sprintf("\nwritten to %s/\n", OUT))
cat("Strata with usable = FALSE come from fewer than", MIN_STUDIES,
    "studies: report them as single-cohort observations, not pooled estimates.\n")
