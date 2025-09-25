---
title: "gwas2crispr: from GWAS associations to CRISPR-ready targets"
tags:
  - R
  - GWAS
  - CRISPR
  - Bioinformatics
  - Reproducibility
authors:
  - name: Othman Mohammed
    affiliation: 1
    orcid: 0009-0000-5548-790X
affiliations:
  - name: LEOPARD.LY LTD, Birmingham, United Kingdom
    index: 1
date: 2025-09-25
---

## Summary
Genome‑wide association studies (GWAS) have catalogued hundreds of thousands of variants associated with human traits and diseases, but moving from lists of significant single‑nucleotide polymorphisms (SNPs) to functional hypotheses still requires bespoke data wrangling. Researchers typically query the NHGRI–EBI GWAS Catalog, filter by significance, harmonise coordinates to a common genome build, annotate nearby genes and then prepare inputs for guide‑RNA design—steps that are often repeated for each project. **gwas2crispr** is an open‑source R package that automates this entire workflow. Given a trait identifier from the Experimental Factor Ontology (EFO) and a p‑value threshold, it queries the GWAS Catalog via the gwasrapidd client, retrieves all genome‑wide significant associations for that trait, standardises variant coordinates to GRCh38 and parses metadata such as risk allele, effect size and gene context. Results are returned as a data frame and may be exported as three complementary files: (1) a **CSV** summary containing SNP identifiers and annotations; (2) a **BED** file of one‑base SNP loci for genomic visualisation and intersection; and (3) a **FASTA** file of user‑defined flanking sequences (default 200 bp either side) ready for CRISPR guide‑design tools. A command‑line interface wraps the R functionality to allow non‑interactive use on any system. In a prostate‑cancer case study (EFO_0001663) with a 5×10⁻⁸ significance threshold and 200 bp flanks, gwas2crispr retrieved 1 309 unique SNPs mapped to 805 genes, generated a comprehensive target set and highlighted clusters on chromosomes 8, 2 and 6. By removing manual steps and embracing reproducibility, **gwas2crispr** enables rapid transition from statistical association to experimental design.

## Statement of need
The NHGRI–EBI GWAS Catalog and national biobanks offer programmatic access to tens of thousands of trait‑associated variants, and numerous web servers exist for designing CRISPR guides. Yet there is a conspicuous gap between these resources: no open, community‑maintained tool automatically converts a set of significant GWAS hits into formatted inputs for genome editing. Most researchers assemble ad‑hoc scripts to download associations, filter by significance, lift over coordinates, map SNPs to genes and extract flanking sequences before using guide‑design tools. These custom pipelines are error‑prone, irreproducible and difficult to share. **gwas2crispr** addresses this unmet need by providing a CRAN‑compliant, tested and archived R package that encapsulates the hand‑off from population genetics to CRISPR. It abstracts away data wrangling, ensures results are anchored to the GRCh38 genome build and generates standard output formats that integrate seamlessly with genome browsers and guide‑design platforms. By condensing many manual steps into a single, scriptable command, the package promotes reproducible research and lowers the barrier for functional follow‑up of GWAS discoveries.

## Functionality
* **Inputs:** an Experimental Factor Ontology (EFO) trait identifier and a significance threshold (default 5 × 10⁻⁸).
* **Core process:** queries the GWAS Catalog via gwasrapidd, filters associations by p‑value, harmonises variant coordinates to GRCh38 and aggregates metadata (risk allele, effect size, gene mapping, study identifiers).
* **Outputs:** a **CSV** file summarising all significant SNPs and their annotations; a **BED** file containing 1‑bp SNP coordinates on GRCh38; and a **FASTA** file of reference sequences flanking each SNP (default 200 bp, user configurable). When required packages are unavailable, FASTA generation is skipped with a message.
* **Interfaces:** R API (functions exported by the gwas2crispr package) and a command‑line script located in `inst/scripts/gwas2crispr.R`, enabling the pipeline to run non‑interactively: `Rscript gwas2crispr.R -e EFO_0001663 -p 5e-8 -f 200 -o prostate`.
* **Installation:** distributed under the MIT license and available on CRAN. Install with `install.packages("gwas2crispr")`.
* **Repository:** https://github.com/leopard0ly/gwas2crispr
* **Archive DOI:** 10.5281/zenodo.16878244

## Acknowledgements
The author thanks the NHGRI–EBI GWAS Catalog team for maintaining the open infrastructure that enables gwas2crispr and acknowledges colleagues who provided feedback during development. No funding was received for this work, and the author declares no competing interests.

## References
Buniello, A., et al. (2019). The NHGRI–EBI GWAS Catalog of published genome‑wide association studies, targeted arrays and summary statistics 2019. *Nucleic Acids Research*, **47**(D1), D1005–D1012. https://doi.org/10.1093/nar/gky1120  
Jinek, M., et al. (2012). A programmable dual‑RNA‑guided DNA endonuclease in adaptive bacterial immunity. *Science*, **337**(6096), 816–821. https://doi.org/10.1126/science.1225829  
Labun, K., et al. (2019). CHOPCHOP v3: expanding the CRISPR web toolbox beyond genome editing. *Nucleic Acids Research*, **47**(W1), W171–W174. https://doi.org/10.1093/nar/gkz199  
Magno, R., & Maia, A.T. (2020). gwasrapidd: an R package to query, download and wrangle GWAS Catalog data. *Bioinformatics*, **36**(2), 649–650. https://doi.org/10.1093/bioinformatics/btz669  
Malone, J., et al. (2010). Modelling sample variables with an Experimental Factor Ontology. *Bioinformatics*, **26**(8), 1112–1118. https://doi.org/10.1093/bioinformatics/btq099
