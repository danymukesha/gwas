# Usage: Rscript manhattanPlotter.R alz_gwas_summary.tab alz_hits.tab alzheimers_plot

args <- commandArgs(trailingOnly = TRUE)
gwasFile <- args[1]
hitsFile <- args[2]
output <- args[3]

library(data.table)
library(tidyverse)
library(ggrepel)

gwas <- fread(gwasFile)
hits <- fread(hitsFile)

gwas$log10Praw <- -log10(gwas$P)
gwas$log10P <- ifelse(gwas$log10Praw > 40, 40, gwas$log10Praw)
gwas$Plevel <- NA
gwas$Plevel[gwas$P < 5E-08] <- "possible"
gwas$Plevel[gwas$P < 5E-09] <- "likely"

gwasFiltered <- subset(gwas, log10P > 3.1)
snpsOfInterest <- hits$SNP

gwasToPlot <- merge(gwasFiltered, hits, by = "SNP", all.x = TRUE)
gwasToPlot <- gwasToPlot[order(gwasToPlot$CHR, gwasToPlot$BP), ]

plotting <- gwasToPlot |>
    group_by(CHR) |>
    summarise(chr_len = max(BP)) |>
    mutate(tot = cumsum(chr_len) - chr_len) |>
    select(-chr_len) |>
    left_join(gwasToPlot, ., by = c("CHR" = "CHR")) |>
    arrange(CHR, BP) |>
    mutate(
        BPcum = BP + tot,
        is_highlight = ifelse(SNP %in% snpsOfInterest, "yes", "no"),
        is_annotate = ifelse(log10P > 7.3, "yes", "no")
    )

axisdf <- plotting |>
    group_by(CHR) |>
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

thisManhattan <- ggplot(plotting, aes(x = BPcum, y = log10P)) +
    geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("grey", "black"), 22)) +
    scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) +
    geom_point(data = subset(plotting, is_highlight == "yes" & Plevel == "likely"), color = "red", size = 2) +
    geom_point(data = subset(plotting, is_highlight == "yes" & Plevel == "possible"), color = "orange", size = 2) +
    geom_label_repel(data = subset(plotting, is_annotate == "yes"), aes(label = GENE, fill = factor(STATUS)), alpha = 0.5, size = 2) +
    scale_fill_manual(values = c("aquamarine", "cyan")) +
    theme_bw() +
    theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) +
    xlab("Chromosome") +
    ylab("-log10(P-value)")

ggsave(filename = paste0(output, ".manhattanPlot.png"), plot = thisManhattan, width = 12, height = 5, dpi = 300)

# Rscript manhattanPlotter.R alz_gwas_summary.tab alz_hits.tab alzheimers_plot
