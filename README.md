# gwas2crispr

## GWAS-to-CRISPR Data Pipeline for High-Throughput SNP Target Extraction

`gwas2crispr` is an R package designed to provide a reproducible, streamlined, and fully automated pipeline for extracting, processing, and exporting significant GWAS SNPs associated with human traits or diseases (identified by their EFO ID). It prepares data for downstream CRISPR screening and guide RNA design.

---

## Features

* Fetches significant GWAS SNPs from the GWAS Catalog using EFO trait identifiers.
* Annotates SNPs with relevant genomic contexts (e.g., associated genes).
* Exports comprehensive metadata tables (CSV).
* Provides genomic intersection files (BED format).
* Generates FASTA files with user-specified flanking regions for each SNP.
* Designed to integrate seamlessly with downstream functional genomics analyses, including QTL, cCREs, and CRISPR on/off-target tools.

---

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Biostrings", "BSgenome.Hsapiens.UCSC.hg38"))

install.packages("devtools")
devtools::install_github("yourusername/gwas2crispr")
```

Replace `yourusername` with your GitHub username.

---

## Usage

### From R

```r
library(gwas2crispr)

# Prostate cancer example
run_gwas2crispr("EFO_0001663", p_cut = 5e-8, flank_bp = 200, out_prefix = "prostate")

# Lung cancer example
run_gwas2crispr("EFO_0000707", p_cut = 1e-6, flank_bp = 300, out_prefix = "lung")
```

### Command Line

```bash
Rscript path/to/gwas2crispr.R -e EFO_0001663 -p 5e-8 -f 200 -o prostate
```

Replace `path/to/gwas2crispr.R` with the actual path to your installed script.

---

## Example Output

After running, the package generates three output files:

* **CSV**: Comprehensive metadata for SNPs
* **BED**: Genomic intervals suitable for downstream intersection analyses
* **FASTA**: Sequences flanking each SNP (user-defined length)

Output example:

```
prostate_snps_full.csv
prostate_snps_hg38.bed
prostate_snps_flank200.fa
```

---

## Dependencies

The package relies on several CRAN and Bioconductor packages:

* gwasrapidd
* httr
* dplyr
* purrr
* tibble
* tidyr
* readr
* Biostrings
* BSgenome.Hsapiens.UCSC.hg38
* optparse

All dependencies are automatically handled during installation.

---

## License

This package is released under the MIT License.

---

## Author

* **Othman S. I. Mohammed**
  [Leopard.ly](https://leopard.ly)
  Email: [admin@leopard.ly](mailto:admin@leopard.ly)

---

## Citation

If you use `gwas2crispr` in your research, please cite:

```
Mohammed OSI. gwas2crispr: GWAS-to-CRISPR Data Pipeline for High-Throughput SNP Target Extraction. 2025. GitHub repository. https://github.com/yourusername/gwas2crispr
```

---

For questions or support, please open an issue in this repository.
