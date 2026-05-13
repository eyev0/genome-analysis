# CLAUDE.md

## Structure
- `genome_analysis.py` — CLI entry point, report generation, orchestration (~189 lines)
- `core.py` — `Variant`/`GenomeData` dataclasses, VCF parser, download utilities
- `layers/prs.py` — Layer 1: Polygenic Risk Scores (PGS Catalog)
- `layers/pharmacogenetics.py` — Layer 2: Pharmacogenetics (CPIC)
- `layers/clinvar.py` — Layer 3: ClinVar pathogenic variants + APOE + FTD analysis
- `layers/ancestry.py` — Layer 4: Ancestry — AIMs, supervised ADMIXTURE, K13/K23 calculators, Y-DNA + mtDNA haplogroups
- `layers/gwas.py` — Layer 5: GWAS trait associations
- `download_pgs.py` — Standalone PGS scoring file downloader

## Key types
- `GenomeData` (core.py): Parsed VCF. Has `.variants` (dict), `.rsid_index` (rsid→key), `.pos_index` (chrom:pos→key)
- `Variant` (core.py): Single VCF row. Has `.genotype_alleles`, `.effect_allele_count`

## Layer pattern
Each layer module exports `run_layerN_*(genome: GenomeData, ..., logger) -> dict`.
Layer functions return dicts, main assembles them into genome_report.json.

## Commands
```
pip install requests pgscatalog-core  # optional but recommended
python download_pgs.py                # download 44 PGS scoring files
python genome_analysis.py --vcf FILE --pgs-dir ./pgs_scoring_files --output-dir ./genome_report --cache-dir ./cache
python genome_analysis.py --vcf FILE --layer 3  # run single layer
```

### External tools for Layer 4 (optional)
- `brew install plink2 bcftools admixture` (or: `conda install -c bioconda admixture plink2 bcftools`)
- `pip install yhaplo scipy`
- haplogrep: download from https://github.com/genepi/haplogrep3/releases (JAR) and put on PATH as `haplogrep`
- `python download_references.py` — to build the reference panel used by ADMIXTURE

## Data
- VCF files, pgs_scoring_files/, cache/, genome_report/ are gitignored
- ClinVar VCF (~190MB) auto-downloads to cache/ on first Layer 3 run
- Reference build: GRCh37/hg19 throughout
- cache/reference_panel/, cache/k_calculators/, cache/ancestry_user/ — gitignored Layer 4 working directories

## Layer 4 sub-analyses
- AIMs: quick continental sanity check (continental level)
- ADMIXTURE: regional/sub-continental decomposition vs. 1000G+HGDP+SGDP(+AADR) reference
- K-calculators: GEDmatch-style Eurogenes K13 / MDLP K23 micro-ethnicity components
- Haplogroups: Y-DNA + mtDNA direct paternal/maternal line classification
