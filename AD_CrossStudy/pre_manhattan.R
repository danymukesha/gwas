library(data.table)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)

gwas <- fread("gwas_associations.tsv", sep = "\t", header = TRUE, quote = "")
valid_chromosomes <- as.character(c(1:22, "X", "Y", "MT"))
gwas[, CHR_ID_clean := ifelse(CHR_ID %in% valid_chromosomes, CHR_ID, NA_character_)]
gwas_cleaned <- gwas[!is.na(CHR_ID_clean)]
unique(gwas_cleaned$CHR_ID_clean)

read_rds(file = "good_ad_studies.rds")

alz <- gwas_cleaned[grepl("Alzheimer", `DISEASE/TRAIT`, ignore.case = TRUE)] |>
    as_tibble() #|>
# dplyr::filter(PUBMEDID %in% good_studies)

alz_gwas <- alz |>
    dplyr::mutate(
        SNP = paste0("rs", SNP_ID_CURRENT),
        CHR = as.factor(CHR_ID),
        BP = as.integer(CHR_POS),
        P = as.numeric(gsub(",", "", `P-VALUE`))
    ) |>
    filter(!is.na(CHR) & !is.na(BP) & !is.na(P)) |>
    dplyr::select(SNP, CHR, BP, P, `REPORTED GENE(S)`)

fwrite(alz_gwas, "alz_gwas_summary.tab", sep = "\t")

hits <- alz_gwas |>
    filter(P < 5e-8) |>
    mutate(STATUS = ifelse(P < 5e-9, "likely", "possible")) |>
    separate_rows(`REPORTED GENE(S)`, sep = ",") |>
    mutate(GENE = stringr::str_trim(`REPORTED GENE(S)`)) |>
    dplyr::select(SNP, STATUS, GENE)

fwrite(hits, "alz_hits.tab", sep = "\t") # read with fread()

cutoff <- 50
alz_gwas$log10Praw <- -log10(alz_gwas$P)
alz_gwas$log10P <- ifelse(alz_gwas$log10Praw > cutoff, cutoff, alz_gwas$log10Praw)
alz_gwas$Plevel <- NA
alz_gwas$Plevel[alz_gwas$P < 5E-08] <- "possible"
alz_gwas$Plevel[alz_gwas$P < 5E-09] <- "likely"

gwasFiltered <- subset(alz_gwas, log10P > 3.1)
snpsOfInterest <- hits$SNP

gwasToPlot <- merge(gwasFiltered, hits, by = "SNP", all.x = TRUE)
gwasToPlot <- gwasToPlot[order(gwasToPlot$CHR, gwasToPlot$BP), ]

gap <- 2e7
plotting <- gwasToPlot |>
    mutate(
        CHR = factor(CHR, levels = c(as.character(1:22), "X"))
    ) |>
    arrange(CHR, BP) |>
    group_by(CHR) |>
    summarise(chr_len = max(BP)) |>
    mutate(
        chr_len = as.numeric(chr_len),
        tot = cumsum(chr_len + gap) - (chr_len + gap)
    ) |>
    dplyr::select(-chr_len) |>
    left_join(gwasToPlot, ., by = c("CHR" = "CHR")) |>
    mutate(
        BPcum = BP + tot,
        is_highlight = ifelse(SNP %in% snpsOfInterest, "yes", "no"),
        is_annotate = ifelse(log10P > 7.3, "yes", "no")
    ) |>
    mutate(
        CHR = factor(CHR, levels = c(as.character(1:22), "X"))
    ) |>
    mutate(log10P = ifelse(is.infinite(log10P),
        max(log10P[is.finite(log10P)]) + 2, log10P
    ))

axisdf <- plotting |>
    group_by(CHR) |>
    summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop") |>
    mutate(CHR = factor(CHR, levels = c(as.character(1:22), "X"))) |>
    arrange(CHR)

thisManhattan <- ggplot(plotting, aes(x = BPcum, y = log10P)) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("grey", "black"), 23)) +
    scale_x_continuous(
        labels = axisdf$CHR, breaks = axisdf$center,
        expand = c(0.01, 0)
    ) +
    scale_y_continuous(
        expand = c(0, 0), limits = c(5, cutoff + 2)
    ) +
    geom_point(
        data = subset(plotting, is_highlight == "yes" & Plevel == "likely"),
        color = "red", size = 1.5
    ) +
    geom_point(
        data = subset(plotting, is_highlight == "yes" & Plevel == "possible"),
        color = "orange", size = 1.5
    ) +
    geom_label_repel(
        data = subset(plotting, is_annotate == "yes"),
        aes(label = GENE, fill = factor(STATUS)),
        alpha = 0.5, size = 2
    ) +
    scale_fill_manual(values = c("aquamarine", "cyan")) +
    # scale_y_discrete(guide = guide_axis(n.dodge = 2)) +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(
            # face = "bold", color = "#993333",
            # size = 8, vjust = 1.5
        )
    ) + # ylim(2,(plotting$log10P |> max() + 4)) +
    xlab("Chromosome") +
    ylab("-log10(P-value)")

ggsave(
    filename = paste0("figures/alz_gwas", ".manhattanPlot.png"),
    plot = thisManhattan +
        annotate("point",
            x = max(plotting$BPcum) * 0.90, y = cutoff - 1,
            color = "red", size = 2
        ) +
        annotate("text",
            x = max(plotting$BPcum) * 0.91, y = cutoff - 1,
            label = "P < 5e-9 (likely)", hjust = 0, size = 2.8
        ) +
        annotate("point",
            x = max(plotting$BPcum) * 0.90, y = cutoff - 2.2,
            color = "orange", size = 2
        ) +
        annotate("text",
            x = max(plotting$BPcum) * 0.91, y = cutoff - 2.2,
            label = "P < 5e-8 (possible)", hjust = 0, size = 2.8
        ),
    width = 10, height = 5, dpi = 300
)
