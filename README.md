# Genome Analysis Pipeline

Analyze a personal VCF file across five open genomic databases and produce a single JSON report.

## What this is

A command-line tool that takes a standard VCF file (e.g., from 23andMe, AncestryDNA, or clinical sequencing) and runs it through five analysis layers: polygenic risk scores, pharmacogenetics, ClinVar pathogenic variants, ancestry inference, and GWAS trait associations. The output is a unified JSON file designed for downstream consumption by LLMs or custom dashboards.

## Quick start

**Prerequisites:** Python 3.8+

```bash
pip install requests pgscatalog-core
```

`pgscatalog-core` is optional but required for downloading PGS scoring files from source via `download_pgs.py`.

**Option A: Download pre-built databases from GitHub Releases (fastest)**

```bash
# Download from Releases page
gh release download v1.0.0 --repo eyev0/genome-analysis

# Extract PGS scoring files
tar xzf pgs_scoring_files_grch37.tar.gz

# Move ClinVar VCF to cache
mkdir -p cache
mv clinvar_grch37.vcf.gz cache/

# Run
python genome_analysis.py --vcf your_file.vcf --pgs-dir ./pgs_scoring_files
```

Or download manually from [Releases](https://github.com/eyev0/genome-analysis/releases/tag/v1.0.0) and extract to the same locations.

**Option B: Download from original sources**

```bash
python download_pgs.py                                          # Download 44 PGS files from PGS Catalog (~1.3GB)
python genome_analysis.py --vcf your_file.vcf --pgs-dir ./pgs_scoring_files
# ClinVar VCF (~190MB) auto-downloads to cache/ on first run
```

**Output:** `results/genome_report.json`

**Useful flags:**

- `--layer N` -- run only specific layer(s), e.g. `--layer 2 --layer 3`
- `--output-dir DIR` -- output directory (default: `./results`)
- `--cache-dir DIR` -- cache directory for auto-downloaded databases like ClinVar (default: `./cache`)
- `-v` -- verbose logging

## Layers

| Layer | Name | What it analyzes | Source |
|-------|------|-----------------|--------|
| 1 | Polygenic Risk Scores | 44 conditions across 6 health domains | PGS Catalog |
| 2 | Pharmacogenetics | 13 genes, drug-gene interactions | CPIC Guidelines |
| 3 | ClinVar Variants | 65 curated + full ClinVar VCF cross-reference | NCBI ClinVar |
| 4 | Ancestry | Population proportions from 40 AIMs | Published AIMs panels |
| 5 | GWAS Associations | 60 curated trait associations | GWAS Catalog |

## Data sources

- **PGS Catalog** -- <https://www.pgscatalog.org/> -- polygenic score repository; 44 harmonized GRCh37 scoring files
- **ClinVar** -- <https://www.ncbi.nlm.nih.gov/clinvar/> -- variant-disease database; auto-downloaded on first run (~190 MB)
- **CPIC** -- <https://cpicpgx.org/> -- clinical pharmacogenomics guidelines
- **Pre-downloaded PGS files** -- see Releases

## Output

The pipeline writes `genome_report.json` with the following top-level keys:

- `meta` -- sample ID, variant counts, reference build, run date
- `layer1_prs` -- 44 polygenic risk scores with category tags (cardiometabolic, neurodegeneration, mental health, metabolism, oncology, immunology)
- `layer2_pharmacogenetics` -- 13 genes with diplotypes, phenotypes, and drug recommendations
- `layer3_clinvar` -- carrier variants, ClinVar VCF hits, APOE haplotype, FTD analysis
- `layer4_ancestry` -- population proportions (EUR, AFR, EAS, AMR, SAS)
- `layer5_gwas` -- trait associations grouped by category

## Limitations

- **Educational/informational only** -- not clinical-grade. Do not make medical decisions based on this output.
- Genotyping chip coverage is ~630K out of millions of possible variants; many scores will have low overlap.
- APOE-defining variants (rs429358, rs7412) may be missing from some VCF files due to chip design.
- C9orf72 repeat expansion (relevant to FTD/ALS) cannot be detected from VCF data.
- PRS scores are raw sums with no population percentile normalization.
