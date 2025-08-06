# Load data.table for fast processing
library(data.table)

# Read matched GWAS summary stats
trait1 <- gwas[
    grepl("Alzheimer", MAPPED_TRAIT, ignore.case = TRUE) & !is.na(PVALUE_MLOG) &
        PVALUE_MLOG > 7.3, .(SNPS, TRAIT = MAPPED_TRAIT, PVALUE_MLOG)
] |> dplyr::select(SNPS, ) # columns: SNP, Z, P, CHR, BP
trait2 <- fread("gwas_trait2.tsv") # same SNPs and columns

# Merge datasets on SNP
gwas_merged <- merge(
    trait1[, .(SNP, CHR, BP, Z1 = Z, P1 = P)],
    trait2[, .(SNP, Z2 = Z, P2 = P)],
    by = "SNP"
)

# Estimate variance and correlation under null
Z_mat <- as.matrix(gwas_merged[, .(Z1, Z2)])
P_mat <- as.matrix(gwas_merged[, .(P1, P2)])
varZ <- var.placo(Z_mat, P_mat)
corZ <- cor.pearson(Z_mat, P_mat, returnMatrix = FALSE)

# Apply PLACO+ to all SNPs
placo_results <- gwas_merged[,
    {
        res <- placo.plus(Z = c(Z1, Z2), VarZ = varZ, CorZ = corZ)
        .(CHR = CHR, BP = BP, T_PLACO = res$T.placo.plus, P_PLACO = res$p.placo.plus)
    },
    by = SNP
]


# Load ggplot2 for visualization
library(ggplot2)

# Create -log10(p) and cumulative position for Manhattan plot
placo_results[, logp := -log10(P_PLACO)]

# Order chromosomes
placo_results[, CHR := as.numeric(as.character(CHR))]
placo_results <- placo_results[order(CHR, BP)]

# Calculate cumulative base pair positions
placo_results[, pos_cum := BP]
chr_offsets <- placo_results[, .(chr_len = max(BP)), by = CHR][, cum_offset := cumsum(chr_len) - chr_len]
placo_results <- merge(placo_results, chr_offsets[, .(CHR, cum_offset)], by = "CHR")
placo_results[, pos_cum := BP + cum_offset]

# Label top hits
placo_results[, is_genomewide := P_PLACO < 5e-8]

# Manhattan plot
gg_manhattan <- ggplot(placo_results, aes(x = pos_cum, y = logp, color = factor(CHR %% 2))) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("#1b9e77", "#d95f02")) +
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed") +
    labs(
        title = "PLACO+ Manhattan Plot for Pleiotropic Association",
        x = "Genomic Position",
        y = expression(-log[10](italic(P))),
        color = "Chromosome"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

print(gg_manhattan)
