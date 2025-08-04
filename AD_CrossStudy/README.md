**Systematic Cross-Study Replication and Cross-Trait Pleiotropy Analysis
of Alzheimer’s Disease Genetic Associations in the GWAS Catalog**

# Abstract

**Background:** Genome-wide association studies (GWAS) have
significantly advanced our understanding of Alzheimer's disease (AD)
genetics, yet systematic evaluation of replication patterns and
cross-trait pleiotropy across the entire GWAS Catalog remains limited.
Understanding the reproducibility of genetic associations and their
broader phenotypic effects is crucial for prioritizing robust risk loci
and identifying shared biological pathways.

**Methods:** We performed a comprehensive analysis of 33,581 AD-related
associations from 268 studies and 25,279 unique SNPs in the GWAS Catalog
database through August 2025. Cross-study replication was assessed by
counting independent studies reporting each SNP, classifying variants as
single-study, limited (2-3 studies), moderate (4-6 studies), or highly
replicated (≥7 studies). Cross-trait pleiotropy was evaluated by
identifying AD-associated SNPs that also showed associations with non-AD
traits. Temporal discovery trends and population diversity were analyzed
to assess the maturity and representativeness of AD genetics
research[5][6].

**Results:** Of 25,279 unique AD-associated SNPs, 95.6% (24,159) were
reported in only single studies, while just 0.6% (150) showed high
replication across multiple independent studies. Cross-trait pleiotropy
analysis identified 13,950 AD-associated SNPs (55.2%) that also
influenced other phenotypes, revealing extensive shared genetic
architecture beyond neurodegeneration. Temporal analysis showed
accelerating discovery rates, with 2024 representing the peak year for
novel SNP identification (5,619 variants). Population ancestry analysis
revealed substantial European bias, with 70.9% of studies focusing on
European populations despite growing recognition of the need for diverse
genetic research.

**Conclusions:** This systematic analysis reveals that the vast majority
of reported AD genetic associations lack independent replication,
highlighting critical gaps in validation efforts. The extensive
pleiotropy observed suggests that AD shares genetic pathways with
numerous other traits, providing opportunities for therapeutic
re-purposing and mechanistic insights. The continued European ancestry
dominance and rapid discovery rates indicate both the field's maturity
in well-studied populations and the urgent need for more diverse genetic
studies. These findings provide essential guidance for prioritizing
robust AD risk loci and designing future replication studies to
strengthen the genetic foundation for precision medicine approaches in
AD.

## Methods and Materials

## Data Source

We utilized the NHGRI-EBI GWAS Catalog comprehensive summary statistics
release (gwas_associations.tsv; downloaded August 2025), which compiles
published associations from genome-wide association studies (GWAS)
across human traits and diseases. The dataset includes fields detailing
SNP identifiers, associated traits, p-values, sample sizes, population
backgrounds, publication metadata, and reported genes.

## Identification of AD Associations

Associations relevant to AD and related neurodegenerative or cognitive
traits were extracted using systematic string matching. The following
trait keywords were used for inclusion: “Alzheimer disease”,
“Alzheimer's disease”, “Late-onset Alzheimer's disease”, “Early-onset
Alzheimer's disease”, “dementia”, “memory decline”, “memory
performance”, “cognitive function”, and “mild cognitive impairment”.
Records were included if these patterns appeared in the reported
disease/trait or study description fields.

## Replication Analysis Across Studies

To assess replication, we grouped all AD-associated variants (SNPs) by
their unique identifiers. For each SNP, we enumerated the number of
independent studies (distinct PUBMEDIDs) reporting an association with
AD-related traits. SNPs were categorized based on the extent of
replication:

-   **Single-study**: Reported in only one GWAS.

-    **Limited replication**: Reported in 2–3 studies.

-    **Moderate replication**: Reported in 4–6 studies.

-    **High replication**: Reported in ≥7 independent studies.

## Temporal Analysis of SNP Discovery

We extracted the year of publication for each AD-associated SNP. Novel
discoveries and cumulative association counts were tabulated by
publication year to assess trends in SNP identification and the
saturation of genetic discoveries in AD GWAS.

## Cross-Trait Pleiotropy Assessment

Cross-trait pleiotropy was evaluated by identifying AD-associated SNPs
that were also significantly associated (in the GWAS Catalog) with
non-AD traits. For each such pleiotropic SNP, we summarized the number
and categories of additional traits, classified by phenotype domain
(e.g., metabolic, cardiovascular, psychiatric, autoimmune, cancer,
anthropometric, other).

## Population Diversity Characterization

Population ancestry was inferred from the initial sample size
descriptions using keyword matching (e.g., “European”, “East Asian”,
“African”, “Hispanic/Latino”, “admixed”, or “unspecified”). The number
of studies, associations, and unique SNPs were summarized by ancestry
group and by year, quantifying the representation of diverse populations
in AD GWAS.

> All data processing and analyses were conducted in R (version 4.x)
> using the `data.table`, `dplyr`, `tidyr`, `ggplot2`, and related
> packages. Summary statistics and tables were generated for
> publication, and figures illustrating replication distributions,
> discovery trends, pleiotropy, and population diversity were created
> for visualization. The complete analysis pipeline, including scripts
> and summary outputs, is available upon request.

source codee: `analysis_AD.R`
