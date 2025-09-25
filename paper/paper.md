---
title: "gwas2crispr: from GWAS associations to CRISPR‑ready targets"
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
bibliography: paper.bib
output: pdf_document
---

## Summary

Genome‑wide association studies (GWAS) have produced extensive catalogues of single‑nucleotide polymorphisms (SNPs) linked to human traits and diseases, yet the gulf between statistical association and experimental validation remains wide【582875185589488†L83-L90】.  To test hypotheses about causal variants, researchers often need to convert lists of trait‑associated SNPs into inputs for CRISPR guide design.  This hand‑off typically involves querying the NHGRI–EBI GWAS Catalog, filtering associations by significance, harmonising coordinates to a common genome build and exporting flanking sequences—steps that are laborious and error‑prone when undertaken manually.  **gwas2crispr** is an open‑source R package that automates this pipeline.  Given a trait identifier from the Experimental Factor Ontology (EFO) and a p‑value threshold, it queries the GWAS Catalog via the gwasrapidd client【582875185589488†L83-L90】, retrieves all significant associations for that trait, standardises SNP coordinates to the GRCh38 assembly and annotates each SNP with risk allele, effect size and gene context.  Results are returned as an R data frame and can be exported as: (1) a **CSV** summary with identifiers and annotations; (2) a **BED** file of one‑base loci for genome browsers; and (3) a **FASTA** file of reference sequences flanking each SNP (by default 200 bp either side).  A command‑line interface allows the entire workflow to run non‑interactively on any platform.  By collapsing multiple data‑wrangling steps into a single command, **gwas2crispr** streamlines the path from GWAS discovery to CRISPR‑based experimentation.

## Statement of need

Despite the availability of programmatic access to GWAS summary statistics and numerous web services for designing CRISPR guides, there is no community‑maintained tool that bridges these domains.  Researchers routinely assemble bespoke scripts to download associations, apply p‑value thresholds, convert genome builds, map SNPs to genes and extract flanking sequences before using guide‑design tools【582875185589488†L60-L96】.  These ad‑hoc pipelines are difficult to reproduce and maintain.  **gwas2crispr** addresses this gap by providing a CRAN‑compliant, tested and archived solution that encapsulates the hand‑off from population genetics to functional genomics.  It abstracts away data wrangling, ensures that coordinates are harmonised to GRCh38 and outputs standard formats (CSV, BED and FASTA) that integrate with genome browsers and CRISPR design software.  By condensing many manual steps into a reproducible R/CLI workflow, the package lowers the barrier for functional follow‑up of GWAS discoveries and promotes transparent science.

## Implementation

**gwas2crispr** leverages the `gwasrapidd` client to fetch GWAS associations for a given EFO trait and p‑value threshold, lifts coordinates to the GRCh38 assembly, annotates each SNP and extracts flanking sequences.  It returns a tidy data frame and writes CSV/BED/FASTA outputs when a prefix is supplied.  An R function and a command‑line script support both interactive and scripted use.

## Features

Key features of the package include:

* **Automated retrieval and harmonisation:** programmatic access to GWAS associations for any EFO trait, application of p‑value thresholds and conversion of coordinates to GRCh38.
* **Export and interfaces:** export of results as CSV, one‑base BED and flanking FASTA files, accessible via an R function or a command‑line wrapper.

## Novel contributions

**gwas2crispr** fills a gap between GWAS data retrieval and CRISPR guide design by harmonising coordinates, adding gene context and producing ready‑to‑use FASTA sequences in one command for R and CLI users.

## Case study: prostate cancer

To illustrate the utility of **gwas2crispr**, we applied the pipeline to prostate cancer (EFO_0001663).  Using the default significance threshold of 5 × 10⁻⁸ and flanking length of 200 bp, the package retrieved 1 309 unique SNPs that were genome‑wide significant as of 31 August 2025.  Approximately 90 % of these SNPs were mapped by the GWAS Catalog to at least one gene, producing 2 648 SNP–gene links spanning about 805 unique genes【582875185589488†L83-L90】.  The resulting BED file highlighted clusters of risk loci on chromosomes 8 (293 SNPs), 2 (265 SNPs) and 6 (198 SNPs).  The CSV summary provides a comprehensive catalogue of SNPs and their annotations, while the FASTA file contains reference sequences ready for CRISPR guide design.  This example demonstrates how **gwas2crispr** can rapidly generate reproducible, CRISPR‑ready target sets for downstream functional studies.

## Example applications

As an illustration beyond oncology, applying **gwas2crispr** to type 2 diabetes (EFO_0001360) with the default significance threshold and a 200 bp flank returned roughly 1 200 genome‑wide significant SNPs and about 1 900 SNP–gene links as of 31 August 2025.  These results, exported as CSV, BED and FASTA, provide a ready starting point for designing CRISPR assays targeting genes involved in glucose homeostasis and insulin secretion.  This example demonstrates that the workflow is applicable to a broad range of complex traits.

## Performance and availability

**gwas2crispr** runs on Linux, macOS and Windows.  The package is distributed under the MIT license and is available on CRAN and GitHub.  Installation is achieved via:

```
install.packages("gwas2crispr")
```

The source code and issue tracker are hosted at <https://github.com/leopard0ly/gwas2crispr>.  A versioned snapshot is archived on Zenodo with DOI 10.5281/zenodo.16878244.  Automated continuous‑integration workflows run `R CMD check` across platforms to ensure reliability.  The package includes unit tests and skips long network calls during checks to remain CRAN‑friendly.  We encourage contributions and welcome feedback via GitHub issues or pull requests.

## Acknowledgements

The author thanks the NHGRI–EBI GWAS Catalog team for maintaining the open infrastructure that enables **gwas2crispr** and acknowledges colleagues who provided feedback during development.  No funding was received for this work, and the author declares no competing interests.

## References

::: {#refs}
::: {}

[@Buniello2019]

[@Jinek2012]

[@Labun2019]

[@Magno2020]

[@Malone2010]

[@Sollis2023]

:::

:::
