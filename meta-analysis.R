# Title: GWAS data processing and meta-analysis
# Description: In this part, we processes multiple GWAS (Genome-Wide Association Study) files,
# checks for missing columns, performs data cleaning, and executes an inverse-variance weighted
# meta-analysis to combine results across different studies. The script utilizes parallel
# processing and progress tracking to efficiently handle large data-sets. Missing column details
# are reported in a separate text file for troubleshooting. The final meta-analysis results,
# including beta coefficients, standard errors, and p-values, are saved to a text file for
# further interpretation.
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("furrr", quietly = TRUE)) install.packages("furrr")
if (!requireNamespace("foreach", quietly = TRUE)) install.packages("foreach")
if (!requireNamespace("progressr", quietly = TRUE)) install.packages("progressr")

library(data.table)
library(dplyr)
library(furrr)
library(foreach)
library(progressr)
library(foreach)

unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

load(file = "harmonized_files.Rdata")
file_paths <- harmonized_files[grep("harmonized", harmonized_files)]

missing_columns_report <- list()

read_gwas <- function(file, study_name) {
    dt <- fread(file, header = TRUE)
    required_cols <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
    missing_cols <- setdiff(required_cols, colnames(dt))
    if (length(missing_cols) > 0) {
        missing_columns_report <<-
            append(
                missing_columns_report,
                list(list(
                    file = file,
                    missing_columns = missing_cols
                ))
            )
        cat("Warning: Missing columns in file:", file, "\n")
        cat("Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    }

    dt[, A1 := toupper(A1)]
    dt[, A2 := toupper(A2)]
    dt$study_name <- study_name
    return(dt)
}

plan(multisession, workers = 5)

with_progress({
    p <- progressor(along = file_paths)
    gwas_list <- future_map(file_paths, function(file) {
        p()
        study_name <- sub(".*harmonized_(.*)\\.txt", "\\1", file)
        read_gwas(file, study_name)
    })
})

if (length(gwas_list) != length(file_paths)) {
    stop("The length of gwas_list and file_paths do not match!")
}

gwas_list_named <- Map(function(file_path, gwas_data) {
    data <- gwas_data
    data$file_path <- file_path
    return(data)
}, file_paths, gwas_list)

names(gwas_list_named) <- gsub("gwas_data/(.*)_harmonized.txt", "\\1", file_paths)

saveRDS(gwas_list, file = "gwas_list.rds", compress = "xz")

missing_report_df <-
    do.call(rbind, lapply(missing_columns_report, function(x) {
        data.frame(file = x$file, missing_columns = paste(x$missing_columns, collapse = ", "))
    }))

write.table(missing_report_df,
    file = "missing_columns_report.txt", sep = "\t",
    quote = FALSE, row.names = FALSE
)
cat("Missing columns report saved to 'missing_columns_report.txt'\n")
print(missing_report_df)

# Perform inverse-variance weighted meta-analysis.

meta_analysis <- function(gwas_list) {
    merged <- Reduce(function(x, y) {
        merge(x, y,
            by = c("SNP", "CHR", "BP", "A1", "A2"),
            all = TRUE
        )
    }, gwas_list)
    meta_results <- data.table(
        SNP = merged$SNP, CHR = merged$CHR, BP = merged$BP,
        A1 = merged$A1, A2 = merged$A2,
        BETA_meta = NA_real_, SE_meta = NA_real_, P_meta = NA_real_
    )

    for (i in 1:nrow(merged)) {
        betas <- unlist(merged[i, grep("BETA", colnames(merged)), with = FALSE])
        ses <- unlist(merged[i, grep("SE", colnames(merged)), with = FALSE])
        valid <- !is.na(betas) & !is.na(ses)
        betas <- betas[valid]
        ses <- ses[valid]

        if (length(betas) > 0) {
            weights <- 1 / (ses^2)
            BETA_meta <- sum(betas * weights) / sum(weights)
            SE_meta <- sqrt(1 / sum(weights))
            Z_meta <- BETA_meta / SE_meta
            P_meta <- 2 * pnorm(-abs(Z_meta))

            meta_results[i, BETA_meta := BETA_meta]
            meta_results[i, SE_meta := SE_meta]
            meta_results[i, P_meta := P_meta]
        }
    }

    return(meta_results)
}

meta_results <- meta_analysis(gwas_list)
fwrite(meta_results, "meta_analysis_results.txt", sep = "\t")
