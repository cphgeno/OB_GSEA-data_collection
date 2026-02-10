suppressPackageStartupMessages({
    library(fs)
    library(optparse)
    library(tidyverse)
    library(R.utils)
})

option_list <- list(
  make_option("--output_dir",
    dest="output_dir", type="character",
    help="output directory where files will be saved"),
  make_option("--name",
    type = "character",
    help = "name of the analysis from omnibenchmark"),
  make_option("--input_dir",
    dest="input_dir", type="character",
    default = "OB_GSEA-data_input",
    help="output directory where files will be saved"),
  make_option("--genesets",
    type = "character",
    help = "name of the genesets for analysis"),
  make_option("--wreference",
    type = "character", default = 'None',
    help = "Provided dataset includes reference samples")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)


# get data directory
get_script_dir <- function() {
  commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = value, into = c("key", "value"), sep = "=", fill = "right"
    ) %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value) %>%
    dirname(.)
}


if (opts$input_dir == "OB_GSEA-data_input") {
  input_dir <- (file.path(get_script_dir(), opts$input_dir))
  gcompressed = TRUE
  gzip_ext = ".gz"
} else {
  input_dir <- opts$input_dir
  gcompressed = FALSE
  gzip_ext = ""
}


# identify the three required inputs for analysis
counts <- dir_ls(input_dir,
                  glob = paste0("*", opts$name, "-counts.tsv", gzip_ext),
                  recurse = FALSE)
metadata <- dir_ls(input_dir,
                  glob = paste0("*", opts$name, "-", opts$genesets, "-metadata.tsv", gzip_ext),
                  recurse = FALSE)
genesets <- dir_ls(input_dir,
                  glob = paste0("*", opts$genesets, ".gmt"),
                  recurse = FALSE)

# sanity check
if (length(counts) == 0) stop(paste("No *-counts.tsv file found in", input_dir))
if (length(metadata) == 0) stop(paste("No *-metadata.tsv file found in", input_dir))
if (length(genesets) == 0) stop(paste("No *.gmt file found in", input_dir))


if (!dir.exists(opts$output_dir)) {
  dir.create(opts$output_dir, recursive = TRUE)
}

# symlink genesets.gmt file
file.symlink(genesets, path(opts$output_dir, paste0(opts$name, '-GenesetsOI.gmt')))

if (gcompressed) {
  # if data from github repo, gunzip before symlinking
  R.utils::gunzip(
    filename = counts,
    overwrite = TRUE,
    remove = FALSE
  )
  counts = sub("\\.gz$", "", counts)
  R.utils::gunzip(
    filename = metadata,
    overwrite = TRUE,
    remove = FALSE
  )
  metadata = sub("\\.gz$", "", metadata)
}

# symlink counts, metadata and reference data based on type of analysis run
if (opts$wreference == 'VSdf') { # use input df as reference = Same-Cohort Reference
  file.symlink(counts, path(opts$output_dir, paste0(opts$name, '-counts.tsv')))
  file.symlink(metadata, path(opts$output_dir, paste0(opts$name, '-metadata.tsv')))
  file.symlink(counts, path(opts$output_dir, paste0(opts$name, '-REFERENCE-counts.tsv'))) # symlink counts as reference as well
} else if (opts$wreference == 'WT'){ # if metadata contains samples that have wildtype/control samples that should instead be used as reference
  meta <- read.table(metadata, sep = '\t', header = T)
  meta_ref <- meta %>% filter(annotation == 'WT')
  counts_df <- column_to_rownames(read.table(counts, sep = '\t', header = T, check.names = F), var = 'Geneid')
  counts_df_ref <- counts_df %>% select(meta_ref$filename)
  write.table(rownames_to_column(counts_df_ref, var = 'Geneid'), path(opts$output_dir, paste0(opts$name, '-REFERENCE-counts.tsv')), sep = '\t', quote = F, row.names = F)
  counts_df <- counts_df %>% select(-c(meta_ref$filename))
  write.table(rownames_to_column(counts_df, var = 'Geneid'), path(opts$output_dir, paste0(opts$name, '-counts.tsv')), sep = '\t', quote = F, row.names = F)
  write.table(meta %>% filter(!annotation == 'WT'), path(opts$output_dir, paste0(opts$name, '-metadata.tsv')), sep = '\t', quote = F, row.names = F)
} else { # Universal Reference
  counts_ref <- dir_ls(input_dir,
                      regexp = "UniversalReference-counts\\.tsv\\.gz$",
                      recurse = FALSE)
  if (length(counts_ref) == 0) stop(paste("No UniversalReference-counts.tsv.gz file found in", input_dir))
  if (gcompressed) {
    R.utils::gunzip(
      filename = counts_ref,
      overwrite = TRUE,
      remove = FALSE
    )
  }
  counts_ref = sub("\\.gz$", "", counts_ref)
  file.symlink(counts_ref, path(opts$output_dir, paste0(opts$name, '-REFERENCE-counts.tsv')))
  file.symlink(counts, path(opts$output_dir, paste0(opts$name, '-counts.tsv')))
  file.symlink(metadata, path(opts$output_dir, paste0(opts$name, '-metadata.tsv')))
}


message("All required input data symlinked to working dir")
