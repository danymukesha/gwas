# Title: Extract and Harmonize AD GWAS Summary Statistics from GWAS Catalog Associations
# Description: This script reads the GWAS Catalog associations file, filters for Alzheimer's disease studies,
# downloads harmonized summary statistics, and standardizes them for meta-analysis.

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
if (!requireNamespace("furrr", quietly = TRUE)) install.packages("furrr")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("progressr", quietly = TRUE)) install.packages("progressr")

library(data.table)
library(dplyr)
library(httr)
library(future)
library(furrr)
library(tibble)
library(foreach)
library(progressr)

# Function to filter AD studies from GWAS Catalog associations
filter_studies <- function(gwas_data, trait = "Alzheimer") {
    ad_studies <- gwas_data |>
        filter(grepl(trait, `DISEASE/TRAIT`, ignore.case = TRUE)) |>
        dplyr::select(
            `STUDY ACCESSION`, `DISEASE/TRAIT`, PUBMEDID, `FIRST AUTHOR`,
            DATE, `INITIAL SAMPLE SIZE`, `REPLICATION SAMPLE SIZE`
        ) |>
        distinct()

    ad_studies
}


download_with_retry <- function(url, output_file, max_retries = 10) {
    attempt <- 1
    success <- FALSE
    while (attempt <= max_retries & !success) {
        try(
            {
                req <- httr2::request(url) |>
                    httr2::req_timeout(600) |>
                    httr2::req_progress()

                response <- httr2::req_perform(req)

                if (response$status_code %in% 200:299) {
                    writeBin(response$body, output_file)
                    success <- TRUE
                } else {
                    warning("Failed to download (status code: ", response$status_code, ")")
                }
            },
            silent = TRUE
        )

        if (!success) {
            message("Retrying... Attempt ", attempt)
            attempt <- attempt + 1
            Sys.sleep(5)
        }
    }
    if (!success) {
        stop("Failed to download after ", max_retries, " attempts.")
    }
}

download_summary_stats <- function(study_accession, output_dir = "gwas_data",
                                   harmonized_list_file = "harmonised_list.txt") {
    print(study_accession)
    with_progress({
        p <- progressor(steps = 3) # Three steps: harmonized list, raw download, harmonized download

        if (!dir.exists(output_dir)) dir.create(output_dir)

        if (!file.exists(harmonized_list_file)) {
            download.file("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt", harmonized_list_file)
            p(message = "Downloaded harmonized list.")
        }

        harmonized_list <- data.table::fread(harmonized_list_file, header = FALSE) |>
            tibble()

        path <- harmonized_list |> filter(grepl(study_accession, V1))

        if (nrow(path) == 0) {
            raw_url <- paste0("ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/", study_accession, "_raw.txt.gz")
            raw_file <- file.path(output_dir, paste0(study_accession, "_raw.txt.gz"))

            future({
                p(message = paste("Downloaded raw file:", raw_file))
                req <- httr2::request(raw_url) |>
                    httr2::req_timeout(300)
                response <- httr2::req_perform(req)

                if (response$status_code %in% 200:299) {
                    writeBin(response$body, raw_file)
                    message("Downloaded harmonized file: ", raw_file)
                } else {
                    message("Failed to download the file. Status code: ", response$status_code)
                }
            })

            return(raw_file)
        }

        harmonized_file <- file.path(output_dir, paste0(study_accession, "_harmonized.txt.gz"))

        if (!file.exists(harmonized_file)) {
            harmonized_url <- paste0(
                "ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
                sub("^./", "", path$V1)
            )

            download_with_retry(harmonized_url, harmonized_file)
        }


        rename_map <- c(
            "hm_rsid|rsid" = "SNP",
            "hm_chrom|chromosome" = "CHR",
            "hm_pos|base_pair_location" = "BP",
            "hm_effect_allele|effect_allele" = "A1",
            "hm_other_allele|other_allele" = "A2",
            "hm_beta|beta" = "BETA",
            "standard_error" = "SE",
            "p_value" = "P"
        )

        dt <- fread(harmonized_file, header = TRUE, sep = "\t")
        for (old_col in names(rename_map)) {
            possible_cols <- strsplit(old_col, "\\|")[[1]]
            existing_col <- possible_cols[possible_cols %in% colnames(dt)]
            if (length(existing_col) > 0) {
                dt <- dt |> dplyr::rename(!!rename_map[old_col] := existing_col[1])
            }
        }

        dt <- dt |>
            mutate(
                A1 = toupper(A1), A2 = toupper(A2),
                CHR = as.integer(gsub("chr", "", CHR)),
                P = as.numeric(P)
            ) |>
            filter(!is.na(P) & P > 0 & P <= 1)

        harmonized_file_out <- file.path(output_dir, paste0(study_accession, "_harmonized.txt"))
        fwrite(dt, harmonized_file_out, sep = "\t")

        p(message = paste("Saved harmonized data to:", harmonized_file_out))

        return(harmonized_file_out)
    })
}

main <- function(gwas_file = "gwas_associations.tsv") {
    gwas_data <- fread(gwas_file, sep = "\t", header = TRUE, quote = "")

    ad_studies <- filter_studies(gwas_data)
    if (nrow(ad_studies) == 0) {
        stop("No Alzheimer's disease studies found in GWAS Catalog associations")
    }

    ad_studies |> tibble()

    harmonized_files <- c()
    for (accession in ad_studies$`STUDY ACCESSION`) {
        file <- download_summary_stats(accession)
        if (!is.null(file)) {
            harmonized_files <- c(harmonized_files, file)
        }
    }

    message("Downloaded and harmonized ", length(harmonized_files), " studies")
    return(harmonized_files)
}

harmonized_files <- main()

cat("Harmonized files for meta-analysis:\n", paste(harmonized_files, collapse = "\n"), "\n")

# Downloaded and harmonized 235 studies
# Only 22 studies were correctly and completely harmonized.
