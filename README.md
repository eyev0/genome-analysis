# Genome Analysis Pipeline

Analyze a personal VCF file across five open genomic databases and produce a single JSON report.

## What this is

A command-line tool that takes a standard VCF file (e.g., from 23andMe, AncestryDNA, or clinical sequencing) and runs it through five analysis layers: polygenic risk scores, pharmacogenetics, ClinVar pathogenic variants, ancestry inference, and GWAS trait associations. The output is a unified JSON file designed for downstream consumption by LLMs or custom dashboards.

## Quick start

**Prerequisites:** Python 3.8+, Homebrew (macOS) for Layer 4 external tools.

### One-shot install

```bash
make install         # Python deps + brew tools + haplogrep JAR
make check           # verify binaries
```

What `make install` does:
- `pip install -r requirements.txt` (requests, numpy, scipy, pgscatalog-core, yhaplo)
- `brew install plink2 bcftools admixture eigensoft openjdk` (Layer 4 binaries)
- Downloads `haplogrep3` JAR to `~/.local/bin/` and writes a wrapper script

Use `make install-python` alone if you only need Layers 1–3, 5 (no Layer 4 ADMIXTURE/haplogroups).

### Get the databases

**Option A: Download pre-built databases from GitHub Releases (fastest)**

```bash
gh release download v1.0.0 --repo eyev0/genome-analysis
tar xzf pgs_large_grch37.tar.gz
tar xzf pgs_small_grch37.tar.gz
mkdir -p cache && mv clinvar_grch37.vcf.gz cache/
```

Or download manually from [Releases](https://github.com/eyev0/genome-analysis/releases/tag/v1.0.0).

**Option B: Download from original sources**

```bash
make pgs             # PGS Catalog scoring files (~1.3 GB)
make references      # 1000G + HGDP + SGDP (+ AADR) reference panel for ADMIXTURE (slow, GBs)
                     # ClinVar VCF (~190 MB) auto-downloads on first Layer 3 run
```

`make references` is only required for Layer 4 ADMIXTURE sub-analysis. Other sub-analyses (AIMs, K-calculators, haplogroups) work without it.

### Run

```bash
make analyze VCF=path/to/your.vcf
# or directly:
python genome_analysis.py --vcf your.vcf --pgs-dir ./pgs_scoring_files \
    --output-dir ./reports --cache-dir ./cache
```

**Output:** `reports/genome_report.json`

**Useful flags:**

- `--layer N` — run only specific layer(s), e.g. `--layer 2 --layer 3`
- `--output-dir DIR` — output directory (default: `./reports`)
- `--cache-dir DIR` — cache directory (default: `./cache`)
- `-v` — verbose logging

## Layers

| Layer | Name | What it analyzes | Source |
|-------|------|-----------------|--------|
| 1 | Polygenic Risk Scores | 44 conditions across 6 health domains | PGS Catalog |
| 2 | Pharmacogenetics | 13 genes, drug-gene interactions | CPIC Guidelines |
| 3 | ClinVar Variants | 65 curated + full ClinVar VCF cross-reference | NCBI ClinVar |
| 4 | Ancestry (4 sub-analyses) | Continental + regional + micro-ethnicity + haplogroups | 1000G, HGDP, SGDP, AADR, ISOGG, PhyloTree |
| 5 | GWAS Associations | 60 curated trait associations | GWAS Catalog |

### Layer 4 sub-analyses

| Sub-analysis | Method | Resolution | Requires |
|--------------|--------|------------|----------|
| **AIMs** | 40 high-Fst markers, log-likelihood | Continent (EUR/AFR/EAS/SAS/AMR) | Python only |
| **ADMIXTURE** | Supervised ADMIXTURE vs 1000G+HGDP+SGDP(+AADR) | Region/country | `plink2`, `admixture`, reference panel built |
| **K-calculators** | GEDmatch-style Eurogenes K13 / MDLP K23 (Python solver) | Micro-ethnicity components | `scipy` (offline K=7 fallback if `.par` files unreachable) |
| **Haplogroups** | Y-DNA (yhaplo) + mtDNA (haplogrep) | Direct paternal/maternal lines | `yhaplo`, `haplogrep` |

Each sub-analysis gracefully skips with a clear reason if its dependency is missing — Layer 4 always produces partial output.

## Data sources

- **PGS Catalog** — <https://www.pgscatalog.org/> — polygenic score repository; 44 harmonized GRCh37 scoring files
- **ClinVar** — <https://www.ncbi.nlm.nih.gov/clinvar/> — variant-disease database; auto-downloaded on first run (~190 MB)
- **CPIC** — <https://cpicpgx.org/> — clinical pharmacogenomics guidelines
- **1000 Genomes Phase 3** — <http://ftp.1000genomes.ebi.ac.uk/> — 2,504 samples, 26 populations
- **HGDP** — Sanger CGP HGDP — 929 samples, 54 populations
- **SGDP** — Simons Genome Diversity Project — 300 samples, 142 populations
- **AADR** — <https://reich.hms.harvard.edu/> — ancient DNA resource (optional; requires `convertf` from EIGENSOFT)
- **yhaplo** — <https://github.com/23andMe/yhaplo> — Y-DNA haplogroup caller (ISOGG tree)
- **haplogrep3** — <https://github.com/genepi/haplogrep3> — mtDNA haplogroup caller (PhyloTree)
- **Pre-downloaded PGS files** — see Releases

## Output

The pipeline writes `genome_report.json` with the following top-level keys:

- `meta` — sample ID, variant counts, reference build, run date
- `layer1_prs` — 44 polygenic risk scores with category tags
- `layer2_pharmacogenetics` — 13 genes with diplotypes, phenotypes, drug recommendations
- `layer3_clinvar` — carrier variants, ClinVar VCF hits, APOE haplotype, FTD analysis
- `layer4_ancestry` — `{aims, admixture, k_calculators, haplogroups}` (each may be `{skipped: true, reason: ...}` if its tools are missing)
- `layer5_gwas` — trait associations grouped by category

## Makefile targets

| Target | Action |
|--------|--------|
| `make install` | Full setup: Python deps + brew tools + haplogrep |
| `make install-python` | Python deps only |
| `make install-tools` | Homebrew binaries only |
| `make install-haplogrep` | haplogrep3 JAR + wrapper |
| `make check` | Verify all binaries on PATH |
| `make pgs` | Download PGS scoring files |
| `make references` | Build reference panel for ADMIXTURE (slow) |
| `make analyze VCF=path` | Run full pipeline on a VCF |
| `make clean-cache` | Remove cache directory |

## Limitations

- **Educational/informational only** — not clinical-grade. Do not make medical decisions based on this output.
- Genotyping chip coverage is ~630K out of millions of possible variants; many scores will have low overlap.
- APOE-defining variants (rs429358, rs7412) may be missing from some VCF files due to chip design.
- C9orf72 repeat expansion (relevant to FTD/ALS) cannot be detected from VCF data.
- PRS scores are raw sums with no population percentile normalization.
- **Layer 4 K-calculators** use a built-in K=7 fallback when canonical Eurogenes K13 / MDLP K23 `.par` files are unreachable. The fallback is not a full GEDmatch substitute — it demonstrates the math but uses ~40 AIMs vs the thousands of SNPs in real calculators.
- **Layer 4 haplogroups** illuminate only two of your thousands of ancestral lines (direct paternal Y, direct maternal mt). They say nothing about autosomal ancestry.
- **No relative matching.** GEDmatch-style one-to-many matching requires a participant database (open alternatives don't exist at scale). This pipeline computes admixture and haplogroups only.
