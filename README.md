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

