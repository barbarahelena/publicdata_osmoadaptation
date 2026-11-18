#!/usr/bin/env Rscript
# Step 1 of 2: select the samples to download (high vs low income comparison).
#
# Selection: stool, Illumina, >= MIN_READS reads, median read length >= MIN_READ_LEN,
# adults with sex and BMI available, controls only, antibiotics explicitly "no",
# not pregnant or lactating, one sample per person (the deepest).
#
# Writes (default dir results/cmd_download/):
#   samples_selected.tsv   one row per sample, with income group
#   runs.tsv               one row per run accession (a sample can have several)
#   selection_flow.tsv     sample counts after every filter
#
#   Rscript select_download_set.R [outdir]
# Next: Rscript ena_batches.R [outdir]

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(dplyr)
  library(tidyr)
  library(readr)
})

OUT <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "results/cmd_download"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- thresholds ------------------------------------------------------------
MIN_READS     <- 5e6
MIN_READ_LEN  <- 75
MIN_AGE       <- 18
MIN_PER_STUDY <- 1      # raise (e.g. 20) to drop small studies

# World Bank income group by country (FY2026, written from memory: check before
# publishing). Income is a country-level label, not a per-person measure.
INCOME <- c(
  # high
  GBR = "High", ISR = "High", NLD = "High", USA = "High", DNK = "High", IRL = "High",
  ITA = "High", DEU = "High", AUT = "High", FRA = "High", SWE = "High", CAN = "High",
  ESP = "High", LUX = "High", ISL = "High", NOR = "High", SVK = "High", HUN = "High",
  EST = "High", FIN = "High",
  # upper-middle
  CHN = "Upper-middle", KAZ = "Upper-middle", PER = "Upper-middle",
  # lower-middle
  IND = "Lower-middle", CMR = "Lower-middle",
  # low
  MDG = "Low")
INCOME_LEVELS <- c("High", "Upper-middle", "Lower-middle", "Low")

# ---- flow bookkeeping ------------------------------------------------------
flow <- tibble(step = character(), n_samples = integer(),
               n_subjects = integer(), n_studies = integer())
note <- function(d, step) {
  force(d)  # evaluate upstream pipe steps first, or their flow rows get overwritten
  flow <<- bind_rows(flow, tibble(
    step = step, n_samples = nrow(d),
    n_subjects = n_distinct(d$subject_id), n_studies = n_distinct(d$study_name)))
  d
}

sel <- sampleMetadata %>% as_tibble() %>%
  note("all samples") %>%
  filter(body_site == "stool") %>% note("stool only") %>%
  filter(!is.na(number_reads), number_reads >= MIN_READS) %>%
  note(sprintf("reads >= %.0e", MIN_READS)) %>%
  filter(!is.na(median_read_length), median_read_length >= MIN_READ_LEN) %>%
  note(sprintf("median read length >= %d bp", MIN_READ_LEN)) %>%
  filter(is.na(sequencing_platform) |
           grepl("Illumina", sequencing_platform, ignore.case = TRUE)) %>%
  note("Illumina") %>%
  filter(!is.na(NCBI_accession)) %>% note("has accession") %>%
  filter(!is.na(age), age >= MIN_AGE) %>% note(sprintf("age >= %d", MIN_AGE)) %>%
  filter(!is.na(gender)) %>% note("sex available") %>%
  filter(!is.na(BMI)) %>% note("BMI available") %>%
  filter(is.na(study_condition) | study_condition == "control") %>%
  note("controls only") %>%
  filter(antibiotics_current_use %in% "no") %>%
  note("antibiotics explicitly 'no'") %>%
  filter(is.na(pregnant) | pregnant == "no", is.na(lactating) | lactating == "no") %>%
  note("not pregnant or lactating") %>%
  group_by(subject_id) %>% slice_max(number_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>% note("one sample per subject") %>%
  # the same sequencing runs can be filed under two studies with different subject
  # ids (e.g. DhakanDB_2019 and GuptaA_2019): keep one copy of each
  group_by(NCBI_accession) %>% slice_max(number_reads, n = 1, with_ties = FALSE) %>%
  ungroup() %>% note("no duplicate runs across studies") %>%
  group_by(study_name) %>% filter(n() >= MIN_PER_STUDY) %>% ungroup() %>%
  note(sprintf("studies with >= %d samples", MIN_PER_STUDY)) %>%
  mutate(income = factor(unname(INCOME[country]), levels = INCOME_LEVELS))

unmapped <- sort(unique(sel$country[is.na(sel$income)]))
if (length(unmapped))
  stop("Countries missing from INCOME mapping: ", paste(unmapped, collapse = ", "))

runs <- sel %>%
  select(sample_id, subject_id, study_name, country, income, run_accession = NCBI_accession) %>%
  separate_rows(run_accession, sep = ";") %>%
  mutate(run_accession = trimws(run_accession)) %>%
  filter(run_accession != "") %>%
  distinct()

shared <- runs %>% group_by(run_accession) %>% filter(n_distinct(sample_id) > 1)
if (nrow(shared))
  stop(n_distinct(shared$run_accession), " run accessions still shared by several samples")

write_tsv(sel, file.path(OUT, "samples_selected.tsv"))
write_tsv(runs, file.path(OUT, "runs.tsv"))
write_tsv(flow, file.path(OUT, "selection_flow.tsv"))

cat("\n== selection flow ==\n"); print(as.data.frame(flow), row.names = FALSE)
cat("\n== by income group ==\n")
print(as.data.frame(sel %>% count(income, name = "samples") %>%
        left_join(sel %>% group_by(income) %>%
                    summarise(studies = n_distinct(study_name),
                              countries = n_distinct(country), .groups = "drop"),
                  by = "income")), row.names = FALSE)
cat(sprintf("\n%d samples, %d run accessions -> %s/\n", nrow(sel), nrow(runs), OUT))
cat("Next: Rscript scripts/ena_batches.R", OUT, "\n")
