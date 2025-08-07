library(data.table)
library(ggplot2)
library(stringr)
library(dplyr)
library(grid)

# wget -O gwas_associations.tsv "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"

gwas <- fread("gwas_associations.tsv", sep = "\t", header = TRUE, quote = "")
gwas
ad_associations <- gwas[
    grepl("Alzheimer", MAPPED_TRAIT, ignore.case = TRUE) & !is.na(PVALUE_MLOG) &
        PVALUE_MLOG > 7.3, .(SNPS, TRAIT = MAPPED_TRAIT, PVALUE_MLOG)
]

ad_associations_unique <- ad_associations[, .SD[which.max(PVALUE_MLOG)], by = SNPS]
significant_associations <- gwas[
    SNPS %in% ad_associations_unique$SNPS & !is.na(PVALUE_MLOG) &
        PVALUE_MLOG > 5, .(SNPS, TRAIT = MAPPED_TRAIT, PVALUE_MLOG)
]

significant_associations_unique <- significant_associations[, .SD[which.max(PVALUE_MLOG)], by = .(SNPS, TRAIT)]
significant_associations_unique[, is_AD := grepl("Alzheimer", TRAIT, ignore.case = TRUE)]
non_ad_associations <- significant_associations_unique[
    !significant_associations_unique$is_AD,
    .(SNPS, OTHER_TRAIT = TRAIT, OTHER_PVALUE_MLOG = PVALUE_MLOG)
]

ad_table <- ad_associations_unique[, .(SNPS, AD_TRAIT = TRAIT, AD_PVALUE_MLOG = PVALUE_MLOG)]
pleiotropy_table <- merge(ad_table, non_ad_associations, by = "SNPS", allow.cartesian = TRUE)


pleiotropy_table$AD_PVALUE_MLOG <- as.numeric(pleiotropy_table$AD_PVALUE_MLOG)
pleiotropy_table$OTHER_PVALUE_MLOG <- as.numeric(pleiotropy_table$OTHER_PVALUE_MLOG)

pleiotropy_table$highlight <- ifelse(
    pleiotropy_table$AD_PVALUE_MLOG > 10 & pleiotropy_table$OTHER_PVALUE_MLOG > 10,
    "Highly Significant", "Significant"
)

highest_ad_point <- pleiotropy_table |>
    filter(AD_PVALUE_MLOG == max(AD_PVALUE_MLOG))

highest_ad_point

highest_nonad_point <- pleiotropy_table |>
    filter(OTHER_PVALUE_MLOG == max(OTHER_PVALUE_MLOG))

highest_nonad_point

top_snp_points <- bind_rows(
    highest_ad_point |> dplyr::slice(1),
    highest_nonad_point
) |>
    mutate(SNP_label = case_when(
        SNPS == "rs814573" ~ "Top AD SNP (rs814573)",
        SNPS == "rs7412" ~ "Top non-AD SNP (rs7412)"
    ))

p1 <- ggplot(pleiotropy_table, aes(x = OTHER_PVALUE_MLOG, y = AD_PVALUE_MLOG)) +
    geom_point(aes(color = highlight), alpha = 0.7, size = 2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    geom_point(
        data = pleiotropy_table |> filter(AD_PVALUE_MLOG > 40 & OTHER_PVALUE_MLOG > 20),
        aes(x = OTHER_PVALUE_MLOG, y = AD_PVALUE_MLOG), color = "darkred", size = 2
    ) +
    geom_point(data = top_snp_points, aes(shape = SNP_label), color = "darkred", size = 4) +
    geom_text(data = top_snp_points, aes(label = SNPS), vjust = 1.8, hjust = 0.4, size = 3) +
    labs(
        x = "Non-AD Trait [-log(p-value)]", y = "AD Trait [-log(p-value)]",
        title = "A. Pleiotropy between AD and Non-AD Traits",
        subtitle = "Top SNPs highlighted with distinct shapes",
        caption = "Data from GWAS studies", color = "Significance Level",
        shape = "Key SNPs"
    ) +
    scale_color_manual(values = c("Highly Significant" = "darkred", "Significant" = "blue")) +
    scale_shape_manual(values = c("Top AD SNP (rs814573)" = 17, "Top non-AD SNP (rs7412)" = 18)) +
    theme_bw() +
    theme(
        panel.border = element_rect(colour = "gray", fill = NA, linewidth = 0.8),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12), axis.text = element_text(size = 10),
        axis.title = element_text(size = 12), legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "cm")
    )

ggsave(
    filename = paste0("figures/pleiotropy", ".btw_ADnNonAD.png"), plot = p1, width = 9, height = 6.5, dpi = 300
)

p2 <- highest_ad_point |>
    ggplot(aes(x = reorder(OTHER_TRAIT, OTHER_PVALUE_MLOG), y = OTHER_PVALUE_MLOG, fill = highlight)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
        title = "B. Pleiotropic Signals for Top AD SNP - rs814573",
        subtitle = "Ordered by –log10(p-value) of association",
        x = "Non-AD Trait", y = "-log10(p-value)", fill = "Significance"
    ) +
    scale_fill_manual(values = c("Highly Significant" = "darkred", "Significant" = "blue")) +
    theme_bw() +
    theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10), legend.position = "top"
    )

ggsave(
    filename = paste0("figures/highest_ad_point", ".btw_ADnNonAD.png"), plot = p2,
    width = 9, height = 5, dpi = 300
)

long_data <- data.table(
    SNP = c(pleiotropy_table$SNPS, pleiotropy_table$SNPS),
    Trait = c(pleiotropy_table$AD_TRAIT, pleiotropy_table$OTHER_TRAIT),
    PValue_MLOG = c(pleiotropy_table$AD_PVALUE_MLOG, pleiotropy_table$OTHER_PVALUE_MLOG),
    Type = c(rep("AD", nrow(pleiotropy_table)), rep("Non-AD", nrow(pleiotropy_table)))
)

long_data

pleiotropy_summary <- pleiotropy_table |>
    group_by(SNPS, AD_TRAIT, OTHER_TRAIT) |>
    summarize(
        AD_PVALUE_MLOG = mean(AD_PVALUE_MLOG),
        OTHER_PVALUE_MLOG = mean(OTHER_PVALUE_MLOG), .groups = "drop"
    )
pleiotropy_summary
