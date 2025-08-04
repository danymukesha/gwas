# Title: Cross-Study Replication & Pleiotropy Analysis
# Description: Here we combine cross-study replication analysis with cross-trait pleiotropy detection
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)

library(parallel)
library(foreach)
library(doParallel)
library(Matrix)

unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

## Data preparation & filtering ====

gwas <- fread("../gwas_associations.tsv", sep = "\t", header = TRUE, quote = "")

# AD and related traits
ad_traits <- c(
    "Alzheimer disease",
    "Alzheimer's disease",
    "Late-onset Alzheimer's disease",
    "Early-onset Alzheimer's disease",
    "Dementia",
    "Memory decline",
    "Cognitive decline",
    "Memory performance",
    "Cognitive function",
    "Mild cognitive impairment"
)

ad_gwas <- gwas[
    grepl("Alzheimer|alzheimer|dementia|Dementia|memory|Memory|cognitive|Cognitive",
        `DISEASE/TRAIT`,
        ignore.case = TRUE
    ) |
        grepl("Alzheimer|alzheimer|dementia|Dementia|memory|Memory|cognitive|Cognitive",
            STUDY,
            ignore.case = TRUE
        )
]

cat("Total AD-related associations:", nrow(ad_gwas), "\n")
cat("Unique studies:", length(unique(ad_gwas$PUBMEDID)), "\n")
cat("Unique SNPs:", length(unique(ad_gwas$SNPS)), "\n")

## Cross-study replication analysis ====

perform_replication_analysis <- function(data) {
    data_clean <- data[
        !is.na(SNPS) &
            !is.na(`P-VALUE`) &
            `P-VALUE` != "" &
            !grepl("NR|nr", `P-VALUE`)
    ]
    data_clean$P_VALUE_NUM <- as.numeric(data_clean$`P-VALUE`)
    data_clean <- data_clean[!is.na(P_VALUE_NUM)]

    snp_replication <- data_clean[, .(
        n_studies = length(unique(PUBMEDID)),
        n_traits = length(unique(`DISEASE/TRAIT`)),
        studies = paste(unique(PUBMEDID), collapse = ","),
        traits = paste(unique(`DISEASE/TRAIT`), collapse = " | "),
        min_pvalue = min(P_VALUE_NUM, na.rm = TRUE),
        median_pvalue = median(P_VALUE_NUM, na.rm = TRUE),
        first_year = min(as.numeric(substr(DATE, 1, 4)), na.rm = TRUE),
        last_year = max(as.numeric(substr(DATE, 1, 4)), na.rm = TRUE),
        reported_genes = paste(unique(na.omit(`REPORTED GENE(S)`)), collapse = ", ")
    ), by = SNPS]

    snp_replication[, replication_status := fcase(
        n_studies == 1, "Single_study",
        n_studies >= 2 & n_studies <= 3, "Limited_replication",
        n_studies >= 4 & n_studies <= 6, "Moderate_replication",
        n_studies > 6, "High_replication"
    )]

    return(snp_replication)
}

replication_results <- perform_replication_analysis(ad_gwas)

replication_summary <- replication_results[, .N, by = replication_status][order(-N)]
print(replication_summary)

## Temporal discovery analysis ====

perform_temporal_analysis <- function(data) {
    data$year <- as.numeric(substr(data$DATE, 1, 4))
    data_clean <- data[!is.na(year) & year >= 2020 & year <= 2025]
    temporal_data <- data_clean[, .(
        total_associations = .N,
        unique_snps = length(unique(SNPS)),
        unique_studies = length(unique(PUBMEDID))
    ), by = year]
    temporal_data <- temporal_data[order(year)]
    temporal_data[, cumulative_snps := cumsum(unique_snps)]

    snp_first_year <- data_clean[, .(first_year = min(year)), by = SNPS]
    novel_per_year <- snp_first_year[, .(novel_snps = .N), by = first_year]

    temporal_data <- merge(temporal_data, novel_per_year,
        by.x = "year", by.y = "first_year", all.x = TRUE
    )
    temporal_data[is.na(novel_snps), novel_snps := 0]

    return(temporal_data)
}

temporal_results <- perform_temporal_analysis(ad_gwas)
print(temporal_results)

## Cross-trait pleiotropy analysis ====

perform_pleiotropy_analysis <- function(gwas_data, ad_data) {
    ad_snps <- unique(ad_data$SNPS)
    other_traits <- gwas_data[
        SNPS %in% ad_snps &
            !grepl("Alzheimer|alzheimer|dementia|Dementia|memory|Memory|cognitive|Cognitive",
                `DISEASE/TRAIT`,
                ignore.case = TRUE
            )
    ]

    if (nrow(other_traits) == 0) {
        cat("No pleiotropic associations found\n")
        return(NULL)
    }

    other_traits[, trait_category := fcase(
        grepl("diabetes|glucose|insulin|metabolic", `DISEASE/TRAIT`, ignore.case = TRUE), "Metabolic",
        grepl("cardiovascular|heart|coronary|blood pressure|hypertension", `DISEASE/TRAIT`, ignore.case = TRUE), "Cardiovascular",
        grepl("psychiatric|depression|schizophrenia|bipolar|anxiety", `DISEASE/TRAIT`, ignore.case = TRUE), "Psychiatric",
        grepl("autoimmune|inflammatory|arthritis|lupus", `DISEASE/TRAIT`, ignore.case = TRUE), "Autoimmune",
        grepl("cancer|tumor|carcinoma", `DISEASE/TRAIT`, ignore.case = TRUE), "Cancer",
        grepl("height|weight|BMI|obesity|anthropometric", `DISEASE/TRAIT`, ignore.case = TRUE), "Anthropometric",
        default = "Other"
    )]

    pleiotropy_summary <- other_traits[, .(
        n_other_traits = length(unique(`DISEASE/TRAIT`)),
        trait_categories = paste(unique(trait_category), collapse = ", "),
        other_traits = paste(unique(`DISEASE/TRAIT`), collapse = " | "),
        min_pvalue_other = min(as.numeric(`P-VALUE`), na.rm = TRUE)
    ), by = SNPS]

    ad_snp_summary <- ad_data[, .(
        ad_traits = paste(unique(`DISEASE/TRAIT`), collapse = " | "),
        ad_studies = length(unique(PUBMEDID)),
        min_pvalue_ad = min(`P-VALUE`, na.rm = TRUE),
        genes = paste(unique(na.omit(`REPORTED GENE(S)`)), collapse = ", ")
    ), by = SNPS]

    pleiotropy_results <- merge(ad_snp_summary, pleiotropy_summary, by = "SNPS")

    return(pleiotropy_results)
}

pleiotropy_results <- perform_pleiotropy_analysis(gwas, ad_gwas)
if (!is.null(pleiotropy_results)) {
    print(head(pleiotropy_results, 10))
}

## Population diversity analysis ====

perform_diversity_analysis <- function(data) {
    data$sample_description <- data$`INITIAL SAMPLE SIZE`

    data[, ancestry := fcase(
        grepl("European|Caucasian|White", sample_description, ignore.case = TRUE), "European",
        grepl("Asian|Chinese|Japanese|Korean", sample_description, ignore.case = TRUE), "East Asian",
        grepl("African|Black", sample_description, ignore.case = TRUE), "African",
        grepl("Hispanic|Latino", sample_description, ignore.case = TRUE), "Hispanic/Latino",
        grepl("admixed|mixed|multi", sample_description, ignore.case = TRUE), "Admixed",
        default = "Unspecified"
    )]

    data$year <- as.numeric(substr(data$DATE, 1, 4))

    diversity_trends <- data[, .(
        n_studies = length(unique(PUBMEDID)),
        n_associations = .N
    ), by = .(year, ancestry)]

    ancestry_summary <- data[, .(
        n_studies = length(unique(PUBMEDID)),
        n_associations = .N,
        unique_snps = length(unique(SNPS)),
        year_range = paste(min(year, na.rm = TRUE), max(year, na.rm = TRUE), sep = "-")
    ), by = ancestry][order(-n_studies)]

    return(list(trends = diversity_trends, summary = ancestry_summary))
}

diversity_results <- perform_diversity_analysis(ad_gwas)
print(diversity_results$summary)

## Create visualizations ====

create_visualizations <- function(replication_results, temporal_results, diversity_results) {
    p1 <- ggplot(
        replication_results[, .N, by = replication_status],
        aes(x = "", y = N, fill = replication_status)
    ) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        theme_void() +
        labs(
            title = "SNP Replication Status in AD GWAS",
            fill = "Replication Level"
        ) +
        scale_fill_brewer(palette = "Set3")

    p2 <- ggplot(temporal_results, aes(x = year)) +
        geom_line(aes(y = novel_snps, color = "Novel SNPs"), size = 1) +
        geom_line(aes(y = unique_studies * 10, color = "Studies (x10)"), size = 1) +
        scale_y_continuous(
            name = "Novel SNPs per Year",
            sec.axis = sec_axis(~ . / 10, name = "Number of Studies")
        ) +
        labs(
            title = "Temporal Trends in AD GWAS Discovery",
            x = "Year", color = "Metric"
        ) +
        theme_minimal()

    p3 <- ggplot(
        diversity_results$trends[ancestry != "Unspecified"],
        aes(x = year, y = n_studies, fill = ancestry)
    ) +
        geom_area(alpha = 0.7) +
        labs(
            title = "Population Diversity in AD GWAS Over Time",
            x = "Year", y = "Number of Studies", fill = "Ancestry"
        ) +
        theme_minimal()

    return(list(replication = p1, temporal = p2, diversity = p3))
}

plots <- create_visualizations(replication_results, temporal_results, diversity_results)

print(plots$replication)
print(plots$temporal)
print(plots$diversity)

# Summary tables ====

generate_summary_tables <- function(replication_results, temporal_results,
                                    pleiotropy_results, diversity_results) {
    total_studies <- length(unique(ad_gwas$PUBMEDID))
    total_snps <- length(unique(ad_gwas$SNPS))
    total_associations <- nrow(ad_gwas)
    year_range <- paste(min(ad_gwas$year, na.rm = TRUE),
        max(ad_gwas$year, na.rm = TRUE),
        sep = "-"
    )

    top_replicated <- replication_results[
        n_studies >= 3 & min_pvalue < 5e-8
    ][order(-n_studies, min_pvalue)][1:20]

    if (!is.null(pleiotropy_results)) {
        top_pleiotropic <- pleiotropy_results[
            n_other_traits >= 2
        ][order(-n_other_traits)][1:15]
    }

    population_table <- diversity_results$summary

    key_findings <- list(
        total_studies = total_studies,
        total_unique_snps = total_snps,
        highly_replicated = nrow(replication_results[n_studies >= 4]),
        single_study_only = nrow(replication_results[n_studies == 1]),
        pleiotropic_snps = if (!is.null(pleiotropy_results)) nrow(pleiotropy_results) else 0,
        european_dominance = round(100 * diversity_results$summary[ancestry == "European", n_studies] / total_studies, 1),
        peak_discovery_year = temporal_results[which.max(novel_snps), year]
    )

    return(list(
        overview = key_findings,
        top_replicated = top_replicated,
        top_pleiotropic = if (!is.null(pleiotropy_results)) top_pleiotropic else NULL,
        population_summary = population_table
    ))
}

summary_tables <- generate_summary_tables(
    replication_results, temporal_results,
    pleiotropy_results, diversity_results
)

cat("=== KEY FINDINGS FOR PUBLICATION ===\n")
cat("Total AD GWAS studies analyzed:", summary_tables$overview$total_studies, "\n")
cat("Total unique SNPs:", summary_tables$overview$total_unique_snps, "\n")
cat("Highly replicated SNPs (≥4 studies):", summary_tables$overview$highly_replicated, "\n")
cat("Single-study-only SNPs:", summary_tables$overview$single_study_only, "\n")
if (summary_tables$overview$pleiotropic_snps > 0) {
    cat("Pleiotropic SNPs:", summary_tables$overview$pleiotropic_snps, "\n")
}
cat("European ancestry dominance:", summary_tables$overview$european_dominance, "%\n")
cat("Peak discovery year:", summary_tables$overview$peak_discovery_year, "\n")
