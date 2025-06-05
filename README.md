# CNV_identification
# 🧬 Gene-to-Bin Annotation Matcher

This script annotates genomic bins with gene names and Arabidopsis homologs by checking positional overlaps between two datasets. It supports both **exact** and **±1000 bp** matching windows for more flexible matching.

---

## 📂 Input Files

1. **`updated_merged_chromosomes_bins.csv`**
   - Contains chromosome bin data with median coverage values.
   - Expected columns:
     ```
     chrom, start, end, median.coverage.window_1, ..., median.coverage.window_18
     ```

2. **`annotated_genes_plot.txt`**
   - Tab-delimited file with gene annotations.
   - Expected columns:
     ```
     chrom, start, end, gene_name, arabidopsis_gene
     ```

---

## ✅ Features

- Matches gene entries **exactly** to bin positions based on `chrom`, `start`, and `end`.
- If no exact match is found, the script looks for a bin within **±1000 bp** range.
- Genes that cannot be matched (even within ±1000 bp) are saved for manual inspection.
- Adds the `gene_name` and `arabidopsis_gene` to the corresponding bins in a new output CSV.
- Handles common formatting issues (e.g., whitespaces, mixed data types).

---

## 🛠 Usage

Run the script in your Python environment:

```bash
python 2_overlap.py
