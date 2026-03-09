# Genome Report Deep Dive Analysis

You are a genetics-literate analyst reviewing a personal genome report. The report is JSON from a genotyping chip (~630K variants, GRCh37/hg19). It has 5 layers — PRS, pharmacogenetics, ClinVar variants, ancestry, and GWAS associations.

## Key context

- PRS scores are **raw** (not population-normalized). They are useful for relative comparison between conditions but cannot be directly interpreted as "high" or "low" risk without reference population data. Note the coverage % — scores with <10% coverage are less reliable.
- Pharmacogenetics diplotypes follow CPIC star-allele nomenclature. Non-normal phenotypes have direct clinical implications for drug prescribing.
- ClinVar `clinvar_vcf_hits` = pathogenic variants found by scanning the full ClinVar database (41K ACMG variants). Empty = no known pathogenic variants in actionable genes.
- APOE genotype is critical for Alzheimer risk assessment. If "undetermined", the SNPs (rs429358, rs7412) were not on the genotyping chip.
- GWAS `or_beta` > 1.0 = increased risk for that trait; < 1.0 = decreased risk. `pvalue` indicates statistical confidence.

## Task

Provide a comprehensive analysis organized as:

### 1. Clinical Priority Findings
Anything that could change medical decisions: pathogenic variants, actionable pharmacogenetics, carrier status for serious Mendelian conditions. For each finding, explain the clinical implication and recommended action.

### 2. Pharmacogenetics Profile
For every gene analyzed: what the diplotype means, which drugs are affected, and what dosing adjustments are recommended. Pay special attention to non-normal metabolizers. Create a simple reference table.

### 3. Disease Risk Landscape
Synthesize PRS scores and GWAS hits by health domain (cardiovascular, metabolic, neurological, mental health, oncology, immune). Where PRS and GWAS overlap for the same condition, note whether they point in the same direction.

### 4. Neurodegeneration Focus
Detailed analysis of Alzheimer/FTD/Parkinson risk: APOE status, FTD gene panel (MAPT, GRN, LRRK2, GBA1), C9orf72 limitations, relevant PRS scores (PGS003992 for Alzheimer, PGS000903 for Parkinson). Note what the VCF can and cannot tell us.

### 5. Biohacking & Optimization
Actionable lifestyle/supplement insights: COMT status (dopamine metabolism), MTHFR if present, vitamin D/B12/iron PRS scores, caffeine metabolism, lactose tolerance, bitter taste perception. Practical recommendations.

### 6. Ancestry Context
How ancestry composition affects interpretation of PRS scores (most are trained on European populations). Note any limitations.

### 7. Gaps & Recommended Follow-up
What this report cannot tell: APOE if undetermined, C9orf72 repeat expansion, CNVs, structural variants, full pharmacogenomic panel. Recommend specific clinical tests if findings warrant it.

Be thorough but structured. Use tables where appropriate. Cite specific rsIDs and gene names.

## Report

```json
<PASTE genome_report.json HERE>
```
