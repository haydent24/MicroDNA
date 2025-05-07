# MicroDNA Circle Detection

This project provides a pipeline for detecting microDNA circles in a human genome using soft-clip analysis from BAM files. MicroDNA circles are small, extrachromosomal DNA molecules that can be identified by analyzing soft-clipped reads in sequencing data. The pipeline includes tools for identifying breakpoints, generating consensus sequences, and aligning sequences to validate microDNA circles.

---

## Features

- **Soft-Clip Analysis**: Detects soft-clipped reads in BAM files to identify potential breakpoints of microDNA circles.
- **Consensus Sequence Generation**: Clusters soft-clipped reads and generates consensus sequences for breakpoints.
- **Alignment Validation**: Uses Smith-Waterman alignment to validate microDNA circles by aligning consensus sequences to the reference genome.
- **Customizable Parameters**: Allows users to adjust thresholds for clustering, alignment scoring, and circle size.

---

## Requirements

- Python 3.8+
- Required Python libraries:
  - `pysam`
  - `argparse`
  - `numpy`
  - `matplotlib`

Install the required libraries using:
```bash
pip install pysam numpy matplotlib