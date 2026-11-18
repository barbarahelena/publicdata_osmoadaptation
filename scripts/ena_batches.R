#!/usr/bin/env Rscript
# Step 2 of 2: look up the real FASTQ files and sizes in ENA, split the samples
# into download batches of about BATCH_BYTES, and write a manifest, accession
# list and download script per batch.
#
#   Rscript ena_batches.R [outdir]      (outdir from select_download_set.R)
#
# Batching: samples are ordered so that every income group is spread evenly over
# the whole order, then cut into equal-sized batches. Each batch is therefore a
# representative mix, and the first batches are already usable for analysis.
# All runs of a sample stay in the same batch.
#
# Writes in outdir:
#   ena_files.tsv          cached ENA lookup, one row per FASTQ file
#   batches_summary.tsv    size and composition per batch
#   batch_NN/accessions.txt  run accessions
#   batch_NN/manifest.tsv    one row per file: url, md5, bytes, sample, income
#   batch_NN/download.sh     resumable, checksum-verified download

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

OUT         <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "results/cmd_download"
BATCH_BYTES <- 1e12     # target batch size (1 TB, decimal)
CHUNK       <- 100      # run accessions per ENA request
SEED        <- 1
set.seed(SEED)

runs <- read_tsv(file.path(OUT, "runs.tsv"), show_col_types = FALSE)
cache <- file.path(OUT, "ena_files.tsv")

# ---- ENA lookup ------------------------------------------------------------
ena_query <- function(accs) {
  q <- paste0("run_accession=", accs, collapse = " OR ")
  for (attempt in 1:4) {
    res <- suppressWarnings(system2("curl", c(
      "-sS", "-f", "-m", "180", "-X", "POST",
      "--data-urlencode", shQuote(paste0("query=", q)),
      "-d", "result=read_run",
      "-d", shQuote("fields=run_accession,library_layout,fastq_ftp,fastq_bytes,fastq_md5"),
      "-d", "format=tsv", "-d", "limit=0",
      "https://www.ebi.ac.uk/ena/portal/api/search"), stdout = TRUE, stderr = FALSE))
    if (length(res) >= 1 && grepl("^run_accession", res[1]))
      return(read_tsv(I(paste(res, collapse = "\n")),
                      col_types = cols(.default = "c"), na = ""))
    Sys.sleep(5 * attempt)
  }
  stop("ENA request failed for: ", paste(head(accs, 3), collapse = ","), " ...")
}

if (file.exists(cache)) {
  message("using cached ", cache, " (delete it to refresh)")
  ena <- read_tsv(cache, show_col_types = FALSE,
                  col_types = cols(.default = "c", bytes = "d"))
} else {
  accs <- unique(runs$run_accession)
  chunks <- split(accs, ceiling(seq_along(accs) / CHUNK))
  message("querying ENA for ", length(accs), " runs in ", length(chunks), " requests")
  raw <- bind_rows(lapply(seq_along(chunks), function(i) {
    if (i %% 10 == 1) message("  request ", i, "/", length(chunks))
    ena_query(chunks[[i]])
  }))
  # one row per file (a paired run has two, sometimes a third unpaired file)
  ena <- raw %>%
    filter(!is.na(fastq_ftp)) %>%
    separate_rows(fastq_ftp, fastq_bytes, fastq_md5, sep = ";") %>%
    mutate(fastq_bytes = as.numeric(fastq_bytes)) %>%
    rename(url = fastq_ftp, md5 = fastq_md5, bytes = fastq_bytes)
  write_tsv(ena, cache)
}

missing <- setdiff(runs$run_accession, ena$run_accession)
if (length(missing)) {
  write_lines(missing, file.path(OUT, "runs_without_fastq_in_ena.txt"))
  warning(length(missing), " runs have no FASTQ in ENA (listed in runs_without_fastq_in_ena.txt)")
}

# ---- samples with all their files ------------------------------------------
files <- runs %>% inner_join(ena, by = "run_accession", relationship = "many-to-many")

samples <- files %>%
  group_by(sample_id, subject_id, study_name, country, income) %>%
  summarise(bytes = sum(bytes), n_files = n(), .groups = "drop")

# a sample whose runs are only partly in ENA would be an incomplete sample: drop it
incomplete <- runs %>% group_by(sample_id) %>%
  summarise(any_missing = any(run_accession %in% missing), .groups = "drop") %>%
  filter(any_missing) %>% pull(sample_id)
if (length(incomplete)) {
  message("dropping ", length(incomplete), " samples with runs missing in ENA")
  samples <- filter(samples, !sample_id %in% incomplete)
  files   <- filter(files,   !sample_id %in% incomplete)
}

# ---- batches ---------------------------------------------------------------
# position of each sample inside its income group, scaled to (0,1]; sorting on it
# interleaves the groups so that each stretch of the order has the same mix
samples <- samples %>%
  group_by(income) %>%
  mutate(pos = (rank(runif(n()), ties.method = "first") - 0.5) / n()) %>%
  ungroup() %>%
  arrange(pos, income)

total <- sum(samples$bytes)
n_batches <- ceiling(total / BATCH_BYTES)
target <- total / n_batches
samples <- samples %>%
  mutate(cum_mid = cumsum(bytes) - bytes / 2,
         batch = pmin(n_batches, floor(cum_mid / target) + 1L))

batch_of <- samples %>% select(sample_id, batch)
files <- files %>% inner_join(batch_of, by = "sample_id")

# ---- write per batch -------------------------------------------------------
template <- file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "download_batch.sh")

for (b in seq_len(n_batches)) {
  dir <- file.path(OUT, sprintf("batch_%02d", b))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  fb <- files %>% filter(batch == b) %>%
    arrange(study_name, sample_id, run_accession, url) %>%
    transmute(run_accession, sample_id, subject_id, study_name, income,
              url, md5, bytes = format(bytes, scientific = FALSE, trim = TRUE))
  write_tsv(fb, file.path(dir, "manifest.tsv"))
  write_lines(unique(fb$run_accession), file.path(dir, "accessions.txt"))
  if (file.exists(template)) {
    file.copy(template, file.path(dir, "download.sh"), overwrite = TRUE)
    Sys.chmod(file.path(dir, "download.sh"), "755")
  } else warning("download_batch.sh not found next to this script: ", template)
}

summ <- samples %>%
  group_by(batch) %>%
  summarise(samples = n(), size_GB = round(sum(bytes) / 1e9),
            studies = n_distinct(study_name),
            High = sum(income == "High"), `Upper-middle` = sum(income == "Upper-middle"),
            `Lower-middle` = sum(income == "Lower-middle"), Low = sum(income == "Low"),
            .groups = "drop") %>%
  left_join(files %>% count(batch, name = "files"), by = "batch")
write_tsv(summ, file.path(OUT, "batches_summary.tsv"))
write_tsv(samples %>% select(sample_id, subject_id, study_name, country, income, bytes, batch),
          file.path(OUT, "sample_batches.tsv"))

cat("\n== batches ==\n"); print(as.data.frame(summ), row.names = FALSE)
cat(sprintf("\ntotal %.2f TB in %d batches of ~%.2f TB, %d samples, %d files\n",
            total / 1e12, n_batches, target / 1e12, nrow(samples), nrow(files)))
cat("Download a batch with:  ", file.path(OUT, "batch_01/download.sh"), " [target_dir] [parallel_jobs]\n")
