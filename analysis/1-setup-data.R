## DOWNLOAD PROCESSED DATA FROM ZENODO ----------------------------------------
#
# Downloads processed data files needed to run the analysis from the Zenodo
# record (https://doi.org/10.5281/zenodo.19139338, v1.1.0) into the data/
# directory. The data/ directory is gitignored except for the bundled example
# data (data/example-data-16s.RData) and the map data (data/gshhg.RData).
#
# Usage:
#   source('analysis/setup-data.R')               # core files only (~16 MB)
#   source('analysis/setup-data.R'); setup_data(partitions = TRUE)  # + partitions (~5.5 GB)
#
# The partition archives are only needed to reproduce the stability selection
# and nested validation steps. Model fitting and all vignettes work without them.

library(httr2)
library(fs)

ZENODO_BASE <- "https://zenodo.org/records/19139338/files"  # v1.1.0

# MD5 checksums from the Zenodo record
CHECKSUMS <- list(
  "ncog16s.RData"           = "f44227438d8d87eb1a630eb2baca3c27",
  "ncog18sv4.RData"         = "080b304353ca6bda357518ed88a5b34e",
  "ncog18sv9.RData"         = "c2231ff6c90df2b680e5e487fdab88a8",
  "density-estimates.RData" = "8bcd23bb6bc207784b18a711b043b659",
  "sightings.RData"         = "66c23555557fc5cddb85a3df8a48c98c",
  "sample_table.Rds"        = "55d087e492f5de7f8cbc17a5d7fd3a7f",
  "region_table.Rds"        = "3ec509e04ee47701e3a46703c929e25a",
  "obs_table.Rds"           = "9bf6c4da12d354c172cad05af718815c"
)

PARTITION_FILES <- list(
  "_combined-partitions.zip" = "a43ba0a093d49878fe5607954d1a07cb",
  "_partitions.zip"          = "6fa78986e43e41ca640d5208f88a47e3"
)

#' Download and optionally verify a single file from Zenodo
download_file <- function(filename, dest_dir, checksum = NULL) {
  dest <- path(dest_dir, filename)
  if (file_exists(dest)) {
    message(filename, " already exists, skipping")
    return(invisible(dest))
  }
  url <- paste0(ZENODO_BASE, "/", utils::URLencode(filename), "?download=1")
  message("Downloading ", filename, " ...")
  request(url) |> req_perform(path = dest)
  if (!is.null(checksum)) {
    actual <- tools::md5sum(dest)
    if (actual != checksum) {
      file_delete(dest)
      stop("MD5 mismatch for ", filename, ": expected ", checksum,
           ", got ", actual, ". File deleted.")
    }
    message("  checksum OK")
  }
  invisible(dest)
}

#' Download all processed data files from Zenodo
#'
#' @param partitions Logical. Also download partition archives (~5.5 GB)?
#'   Only needed to reproduce stability selection and nested validation.
#' @param dest_dir Destination directory. Defaults to 'data/'.
setup_data <- function(partitions = FALSE, dest_dir = "data") {
  dir_create(dest_dir)

  # Core processed data files
  for (nm in names(CHECKSUMS)) {
    download_file(nm, dest_dir, checksum = CHECKSUMS[[nm]])
  }

  # Partition archives (large, optional)
  if (partitions) {
    for (nm in names(PARTITION_FILES)) {
      download_file(nm, dest_dir, checksum = PARTITION_FILES[[nm]])
      zip_path <- path(dest_dir, nm)
      message("Unzipping ", nm, " ...")
      utils::unzip(zip_path, exdir = dest_dir)
      message("  done")
    }
  }

  message("\nData setup complete. Files written to: ", dest_dir)
}

# Run automatically when sourced
setup_data()
