# Genome Report Summary Analysis

You are analyzing a personal genome report in JSON format. The report was generated from a VCF file (genotyping chip, ~630K variants, GRCh37/hg19 build) and contains 5 analysis layers:

- **layer1_prs**: Polygenic risk scores for 44 conditions. Raw scores (no population percentile). `category` tags: cardio, metabolism, neuro, mental, oncology, immune. Higher coverage % = more reliable.
- **layer2_pharmacogenetics**: 13 genes with diplotypes, phenotypes, and drug recommendations. Anything other than "Normal Metabolizer" is actionable.
- **layer3_clinvar**: Curated pathogenic variants + full ClinVar VCF cross-reference (41K ACMG variants). `apoe` = Alzheimer risk haplotype. `ftd_analysis` = frontotemporal dementia genes.
- **layer4_ancestry**: Population proportions from ancestry-informative markers.
- **layer5_gwas**: Trait associations from GWAS catalog. Only carrier variants shown (dosage > 0).

## Task

Analyze this genome report and provide a **concise actionable summary** covering:

1. **Clinically significant findings** — anything pathogenic, drug-response variants that change dosing, carrier status for serious conditions
2. **Top pharmacogenetic results** — non-normal metabolizer phenotypes and what drugs are affected
3. **Notable risk factors** — APOE status, FTD genes (especially if family history), highest-impact GWAS hits
4. **Interesting traits** — caffeine metabolism, bitter taste, lactase persistence, etc.

Keep it under 500 words. Lead with what matters most. Flag anything that warrants follow-up with a doctor.

## Report

```json
<PASTE genome_report.json HERE>
```
