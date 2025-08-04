# GWAS pleiotropy analysis: Study setup

The goal of this analysis is to identify pleiotropic SNPs (variants associated with multiple traits) using the GWAS Catalog summary statistics, 
focusing on how traits might be shared for each variant. Initially, I planned to investigate pleiotropy across multiple traits, but I decided 
to dig deeper into a single trait (e.g., Alzheimer's Disease) to provide a more focused study.

The source code used for the analysis is included in the code file. Please refer to it for specific details on how each analysis step is performed.

## Step 1: Download the GWAS Catalog data

### Official GWAS Catalog associations file

To start, you will need to download the most complete and up-to-date GWAS Catalog associations file. You can get it from the official source here:

* **URL**: [GWAS Catalog associations file](https://www.ebi.ac.uk/gwas/api/search/downloads/alternative)
* **File**: `gwas_catalog_v1.0.2-associations_e105_r2024-05-10.tsv` (or the latest version)

You can download the file to your working directory with the following command:

```bash
wget -O gwas_associations.tsv "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
```

### Harmonized GWAS Data

For further analysis, you will often need harmonized GWAS data. This data is updated daily, and you can download the latest list of harmonized studies
from the following link:

[Download Harmonized Studies List](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt)

You can download it manually by clicking the link above, or use the terminal with either of the following commands:

* **Using `curl`**:

  ```bash
  curl -O https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt
  ```

* **Using `wget`**:

  ```bash
  wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/harmonised_list.txt
  ```

> Formatted files are harmonized summary statistics from GWAS studies. These files are cleaned and standardized to ensure consistency across datasets,
> making them easier to work with for your analysis. The data includes key information such as SNP identifiers and p-values that have been harmonized
> across studies, which simplifies further steps like meta-analysis or functional annotation.

---

## Step 2: Running the GWAS meta-analysis

Once the GWAS Catalog data and harmonized data have been downloaded, the analysis will proceed with reading the GWAS files, checking for missing columns,
and performing an inverse-variance weighted meta-analysis. The analysis includes:

1. **Loading GWAS files**: by reading each file using `fread` and checking for the presence of required columns.
3. **Meta-analysis**: we perfom an inverse-variance weighted meta-analysis to combine results from multiple studies.
4. **Missing columns report**: any missing required columns will be logged in a report file for troubleshooting.

---

## Step 3: GWAS pleiotropy analysis - Focus on Alzheimer's disease

The next objective of the analysis is to investigate the pleiotropy of SNPs; i.e., identify variants that are associated with multiple traits. 
While initially the intention was to explore pleiotropy across multiple traits, we focus here on **Alzheimer's Disease (AD)** as a more targeted study.

In this analysis, we will examine how different traits might share associations with common SNPs. This will help identify potential pleiotropic variants 
that could influence the development of AD as well as other traits.

The following analysis steps will be performed:

* **Data processing**: we load and clean GWAS summary statistics.
* **Pleiotropy identification**: we search for SNPs that overlap between AD and other traits.
* **Final reporting**: in the end, we aggregate data across multiple studies to identify consistent findings.

---

## Step 4: Further exploration and functional annotation

Once the meta-analysis is complete, the next steps will include functional annotation of the identified SNPs and further investigation of the biological 
significance of pleiotropic variants. Functional annotation could include checking for overlaps with known regulatory regions, pathways, or protein-coding genes.

