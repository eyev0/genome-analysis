# Clinically Validated Monogenic Variants Checkable in Consumer Genomics VCF Files

## Comprehensive Reference for VCF-Based Genetic Screening

Last updated: 2026-03-09

---

## TABLE OF CONTENTS

1. [ACMG Secondary Findings (SF v3.2) Gene List](#1-acmg-secondary-findings-sf-v32)
2. [High-Penetrance Variants with Known rsIDs](#2-high-penetrance-variants-with-known-rsids)
3. [Pharmacogenomics (PGx) Variants](#3-pharmacogenomics-pgx-variants)
4. [Common Well-Studied SNPs](#4-common-well-studied-snps)
5. [FTD/Neurodegeneration Variants](#5-ftd-and-neurodegeneration-variants)
6. [Databases and Resources](#6-databases-and-resources)
7. [Technical Limitations](#7-technical-limitations-of-vcf-based-screening)

---

## 1. ACMG SECONDARY FINDINGS (SF v3.2)

The ACMG SF v3.2 list (approved Feb 2023) contains **81 genes** across 4 phenotype categories.
These are genes where incidental findings of pathogenic/likely pathogenic variants are considered
medically actionable and should be reported back.

**IMPORTANT NOTE FOR VCF-BASED SCREENING:** The ACMG SF list does NOT specify individual rsIDs
because pathogenic variants in these genes are spread across the entire coding sequence. Most are
rare, private mutations found only in individual families. A VCF-based screen should cross-reference
against the full ClinVar database for pathogenic/likely pathogenic variants in these genes.

### 1A. Cancer Predisposition Genes (28 genes)

| Gene | Disease/Syndrome | Inheritance | Penetrance |
|------|-----------------|-------------|------------|
| APC | Familial adenomatous polyposis | AD | ~100% colon cancer by age 40 |
| BRCA1 | Hereditary breast/ovarian cancer | AD | 60-80% breast, 40-60% ovarian |
| BRCA2 | Hereditary breast/ovarian cancer | AD | 45-70% breast, 15-30% ovarian |
| PALB2 | Hereditary breast cancer | AD | 35-60% breast cancer |
| TP53 | Li-Fraumeni syndrome | AD | ~90% any cancer by age 60 |
| RET | MEN2A/2B, familial medullary thyroid | AD | ~90% medullary thyroid cancer |
| RB1 | Retinoblastoma | AD | ~90% |
| PTEN | PTEN hamartoma tumor syndrome | AD | 85% breast, 35% thyroid |
| STK11 | Peutz-Jeghers syndrome | AD | 85% any cancer by age 70 |
| MLH1 | Lynch syndrome (HNPCC) | AD | 40-80% colorectal cancer |
| MSH2 | Lynch syndrome (HNPCC) | AD | 40-80% colorectal cancer |
| MSH6 | Lynch syndrome | AD | 10-30% colorectal cancer |
| PMS2 | Lynch syndrome | AD | 15-20% colorectal cancer |
| MUTYH | MYH-associated polyposis | AR (biallelic) | ~50-80% colorectal (biallelic) |
| MEN1 | Multiple endocrine neoplasia type 1 | AD | >95% |
| VHL | Von Hippel-Lindau syndrome | AD | ~90% |
| NF2 | Neurofibromatosis type 2 | AD | ~90% |
| SDHB | Hereditary paraganglioma-pheochromocytoma | AD | 30-70% |
| SDHD | Hereditary paraganglioma-pheochromocytoma | AD | 50-80% |
| SDHC | Paraganglioma syndrome 3 | AD | Variable |
| SDHAF2 | Paraganglioma syndrome 2 | AD | Variable |
| MAX | Hereditary pheochromocytoma | AD | Variable |
| TMEM127 | Hereditary pheochromocytoma | AD | Variable |
| TSC1 | Tuberous sclerosis 1 | AD | ~85% |
| TSC2 | Tuberous sclerosis 2 | AD | ~85% |
| WT1 | Wilms tumor | AD | Variable |
| BMPR1A | Juvenile polyposis syndrome | AD | 40-70% |
| SMAD4 | Juvenile polyposis/HHT | AD | 40-70% |

### 1B. Cardiovascular Genes (40 genes)

#### Cardiomyopathy Genes
| Gene | Disease | Inheritance |
|------|---------|-------------|
| MYH7 | Hypertrophic/dilated cardiomyopathy | AD |
| MYBPC3 | Hypertrophic cardiomyopathy | AD |
| TNNT2 | Hypertrophic/dilated cardiomyopathy | AD |
| TPM1 | Hypertrophic cardiomyopathy | AD |
| TNNI3 | Hypertrophic cardiomyopathy | AD |
| MYL3 | Hypertrophic cardiomyopathy | AD |
| MYL2 | Hypertrophic cardiomyopathy | AD |
| ACTC1 | Hypertrophic cardiomyopathy | AD |
| PRKAG2 | Hypertrophic cardiomyopathy (glycogen storage) | AD |
| LMNA | Dilated cardiomyopathy 1A | AD |
| DES | Dilated cardiomyopathy/myofibrillar myopathy | AD |
| TNNC1 | Dilated cardiomyopathy | AD |
| RBM20 | Dilated cardiomyopathy | AD |
| BAG3 | Dilated cardiomyopathy | AD |
| TTN | Dilated cardiomyopathy (truncating variants only) | AD |
| FLNC | Dilated cardiomyopathy/myofibrillar myopathy | AD |
| SCN5A | Brugada syndrome / Long QT 3 / DCM | AD |

#### Arrhythmogenic Cardiomyopathy Genes
| Gene | Disease | Inheritance |
|------|---------|-------------|
| DSP | ARVC 8 | AD |
| PKP2 | ARVC 9 | AD |
| DSG2 | ARVC 10 | AD |
| DSC2 | ARVC 11 | AD |
| TMEM43 | ARVC 5 | AD |

#### Arrhythmia / Channelopathy Genes
| Gene | Disease | Inheritance |
|------|---------|-------------|
| KCNQ1 | Long QT syndrome 1 | AD |
| KCNH2 | Long QT syndrome 2 | AD |
| RYR2 | Catecholaminergic polymorphic VT | AD |
| CASQ2 | Catecholaminergic polymorphic VT 2 | AR |
| TRDN | CPVT 5 / Long QT | AR |
| CALM1 | CPVT 4 / Long QT 14 | AD |
| CALM2 | Long QT 15 | AD |
| CALM3 | Long QT 16 | AD |

#### Aortopathy / Vascular Genes
| Gene | Disease | Inheritance |
|------|---------|-------------|
| FBN1 | Marfan syndrome | AD |
| TGFBR1 | Loeys-Dietz syndrome 1A | AD |
| TGFBR2 | Loeys-Dietz syndrome 1B | AD |
| SMAD3 | Loeys-Dietz syndrome 3 | AD |
| COL3A1 | Vascular Ehlers-Danlos (type IV) | AD |
| MYH11 | Familial thoracic aortic aneurysm | AD |
| ACTA2 | Familial thoracic aortic aneurysm | AD |

#### Vascular Malformation Genes
| Gene | Disease | Inheritance |
|------|---------|-------------|
| ACVRL1 | Hereditary hemorrhagic telangiectasia 2 | AD |
| ENG | Hereditary hemorrhagic telangiectasia 1 | AD |

#### Malignant Hyperthermia
| Gene | Disease | Inheritance |
|------|---------|-------------|
| RYR1 | Malignant hyperthermia susceptibility | AD |
| CACNA1S | Malignant hyperthermia susceptibility | AD |

### 1C. Inborn Errors of Metabolism (5 genes)
| Gene | Disease | Inheritance |
|------|---------|-------------|
| GAA | Pompe disease | AR |
| GLA | Fabry disease | XL |
| OTC | Ornithine transcarbamylase deficiency | XL |
| BTD | Biotinidase deficiency | AR |
| ATP7B | Wilson disease | AR |

### 1D. Metabolic / Miscellaneous (8 genes)
| Gene | Disease | Inheritance |
|------|---------|-------------|
| LDLR | Familial hypercholesterolemia 1 | AD (semi-dominant) |
| APOB | Familial hypercholesterolemia 2 | AD |
| PCSK9 | Familial hypercholesterolemia 3 | AD |
| HFE | Hereditary hemochromatosis (C282Y homozygotes only) | AR |
| TTR | Hereditary transthyretin amyloidosis | AD |
| HNF1A | MODY (maturity-onset diabetes of the young) | AD |
| RPE65 | RPE65-related retinal dystrophy | AR |

---

## 2. HIGH-PENETRANCE VARIANTS WITH KNOWN rsIDs

These are specific, well-characterized variants with established rsIDs that can be directly
looked up in a VCF file by rsID or genomic position.

### 2A. BRCA1/BRCA2 Ashkenazi Jewish Founder Mutations

| rsID | Gene | Variant (HGVS) | Legacy Name | Condition | Penetrance | Actionable |
|------|------|----------------|-------------|-----------|------------|------------|
| rs80357914 | BRCA1 | c.68_69del | 185delAG | Breast/ovarian cancer | 60-80% breast | YES - enhanced screening, prophylactic surgery |
| rs80357906 | BRCA1 | c.5266dup | 5382insC | Breast/ovarian cancer | 60-80% breast | YES |
| rs80359550 | BRCA2 | c.5946del | 6174delT | Breast/ovarian/prostate cancer | 45-70% breast | YES |

**Notes:** ~1 in 40 Ashkenazi Jews carry one of these three mutations. Together they account for
>90% of BRCA-related cancers in Ashkenazi populations. These are INDELS, not SNPs -- they require
VCF files that include indel calls (WGS and most modern arrays can detect these). 23andMe now
tests for 44 BRCA1/2 variants with FDA clearance (expanded from original 3 in 2023).

### 2B. Hemochromatosis (HFE)

| rsID | Gene | Variant | Amino Acid | Condition | Penetrance | Actionable |
|------|------|---------|------------|-----------|------------|------------|
| rs1800562 | HFE | c.845G>A | C282Y | Hereditary hemochromatosis | Males: 24-43% iron overload; Females: 1-14% | YES - serum ferritin monitoring, phlebotomy |
| rs1799945 | HFE | c.187C>G | H63D | Hereditary hemochromatosis | Very low (~1-2% compound het) | MODERATE - monitor ferritin if compound het |

**Genotype-risk breakdown:**
- C282Y/C282Y (homozygous): Main clinical genotype. ~28% of males develop clinical disease
- C282Y/H63D (compound het): ~1-2% penetrance for iron overload disease
- H63D/H63D: ~6.7% documented iron overload, very rarely clinical disease
- ACMG only reports C282Y homozygotes

### 2C. Factor V Leiden & Prothrombin

| rsID | Gene | Variant | Condition | Risk (Het) | Risk (Hom) | Actionable |
|------|------|---------|-----------|------------|------------|------------|
| rs6025 | F5 | c.1601G>A (R506Q) | Factor V Leiden thrombophilia | OR 3-5x VTE | OR 18-80x VTE | YES - avoid OCPs, anticoagulation if events |
| rs1799963 | F2 | c.*97G>A | Prothrombin G20210A | OR 2-3x VTE | OR 5-10x VTE | YES - avoid OCPs, anticoagulation if events |

**Notes:** ~10% of FVL carriers develop abnormal clots in their lifetime. Double heterozygotes
(FVL + PT G20210A) have OR ~5.24 for VTE, approximating FVL homozygosity.

### 2D. Familial Hypercholesterolemia - Specific Known Variants

| rsID | Gene | Variant | Condition | Penetrance | Actionable |
|------|------|---------|-----------|------------|------------|
| rs5742904 | APOB | c.10580G>A (R3527Q) | FH type 2 | High (variable) | YES - statins, PCSK9i |
| rs28362286 | PCSK9 | c.2037C>A (C679X) | PROTECTIVE - Low LDL | Loss-of-function | INFORMATIONAL - 40% lower LDL |
| rs11591147 | PCSK9 | c.137G>T (R46L) | PROTECTIVE - Low LDL | Loss-of-function | INFORMATIONAL - 15-28% lower LDL |

**Notes on LDLR:** LDLR has >2,100 known pathogenic variants. No single rsID predominates.
FH screening is best done via full ClinVar cross-reference against all known LDLR pathogenic variants.

### 2E. Alpha-1 Antitrypsin Deficiency (SERPINA1)

| rsID | Gene | Variant | Allele | Condition | Risk | Actionable |
|------|------|---------|--------|-----------|------|------------|
| rs28929474 | SERPINA1 | c.1096G>A (E342K) | Pi*Z | A1AT deficiency (COPD, liver) | ZZ homozygotes: high risk emphysema | YES - avoid smoking, augmentation therapy |
| rs17580 | SERPINA1 | c.863A>T (E264V) | Pi*S | A1AT deficiency (mild) | SS: mild; SZ: moderate risk | MODERATE - smoking cessation critical |

**Notes:** ZZ homozygotes have only 10-15% of normal AAT levels. Smoking is the major environmental
modifier -- nonsmoking ZZ individuals may have normal lifespan. On ACMG list via SERPINA1.

---

## 3. PHARMACOGENOMICS (PGx) VARIANTS

These are the most practically impactful variants for a VCF-based screen. Unlike disease-causing
variants, PGx variants are COMMON and directly affect drug dosing.

### 3A. CYP2D6 (Codeine, Tamoxifen, SSRIs, Tramadol, Ondansetron, TCAs)

CYP2D6 is complex -- it has >100 star alleles, gene deletions, and duplications. Standard
VCF analysis CANNOT detect copy number variants (gene deletion = *5, duplications = xN).
However, the key SNP-based alleles are:

| rsID | Star Allele | Function | Key Info |
|------|-------------|----------|----------|
| rs3892097 | *4 | No function | Most common null allele in Europeans (~20% allele freq). Splicing defect |
| rs35742686 | *3 | No function | Frameshift deletion. ~1-2% in Europeans |
| rs5030655 | *6 | No function | Deletion, 0-1% frequency |
| rs1065852 | *10 | Decreased function | Most common in East Asians (~40% allele freq). Also part of *4 haplotype |
| rs28371706 | *17 | Decreased function | Most common in Africans (~20% allele freq) |
| rs28371725 | *41 | Decreased function | Splicing defect. ~8% in Europeans |
| rs16947 | *2 | Normal function (may be decreased) | Recent evidence suggests 2-fold reduced expression |
| rs1135840 | *2/*10 defining | Context-dependent | Part of multiple haplotypes |
| rs5030656 | *9 | Decreased function | Deletion of amino acid |

**CPIC Tier 1 (minimum recommended panel):** *2, *3, *4, *5 (deletion), *6, *9, *10, *17, *29, *41
**Limitation:** *5 (whole gene deletion) and gene duplications are NOT detectable from SNP-based VCF.

### 3B. CYP2C19 (Clopidogrel, PPIs, Voriconazole, SSRIs, TCAs)

| rsID | Star Allele | Function | Clinical Impact |
|------|-------------|----------|-----------------|
| rs4244285 | *2 | No function | Most common LoF allele. ~15% Europeans, ~30% East Asians |
| rs4986893 | *3 | No function | Premature stop codon. Rare in Europeans, ~5% East Asians |
| rs28399504 | *4 | No function | Rare |
| rs56337013 | *5 | No function | Rare |
| rs72552267 | *6 | No function | Rare |
| rs72558186 | *7 | No function | Rare |
| rs41291556 | *8 | No function | Rare |
| rs12248560 | *17 | Increased function | Ultrarapid metabolism. ~21% Europeans |

**Key clinical scenario:** CYP2C19 *2/*2 = poor metabolizer = clopidogrel INEFFECTIVE (use
prasugrel/ticagrelor instead). This is one of the most impactful PGx findings.

### 3C. CYP2C9 (Warfarin, Phenytoin, NSAIDs, Sulfonylureas)

| rsID | Star Allele | Function | Clinical Impact |
|------|-------------|----------|-----------------|
| rs1799853 | *2 | Decreased function | ~30% reduced metabolism. ~13% Europeans |
| rs1057910 | *3 | Decreased function | ~80% reduced metabolism. ~7% Europeans |
| rs56165452 | *5 | No function | Rare. Common in African Americans |
| rs28371686 | *6 | No function | Rare |
| rs9332131 | *8 | No function | Frameshift. ~0.5% in Africans |
| rs7900194 | *11 | Decreased function | Rare |

**Key clinical scenario:** CYP2C9 *2/*3 or *3/*3 carriers need significantly lower warfarin doses.
Combined with VKORC1 genotype for optimal warfarin dosing algorithms.

### 3D. VKORC1 (Warfarin)

| rsID | Gene | Variant | Effect | Clinical Impact |
|------|------|---------|--------|-----------------|
| rs9923231 | VKORC1 | c.-1639G>A | Reduced expression | AA genotype: ~50% lower warfarin dose needed. ~37% Europeans carry A allele |

**This is the single most important pharmacogenomic variant for warfarin dosing.**

### 3E. CYP3A5 (Tacrolimus)

| rsID | Gene | Star Allele | Function | Clinical Impact |
|------|------|-------------|----------|-----------------|
| rs776746 | CYP3A5 | *3 | No function (splicing defect) | ~85-95% of Europeans are *3/*3 (non-expressors). Expressors need higher tacrolimus doses |

### 3F. DPYD (5-Fluorouracil, Capecitabine)

| rsID | Star Allele | Function | Clinical Impact |
|------|-------------|----------|-----------------|
| rs3918290 | *2A (c.1905+1G>A) | No function | LIFE-THREATENING toxicity if treated with 5-FU. ~1-2% Europeans |
| rs55886062 | *13 (c.1679T>G) | No function | Severe toxicity risk |
| rs67376798 | c.2846A>T (D949V) | Decreased function | Moderate toxicity risk |
| rs75017182 | c.1129-5923C>G (HapB3) | Decreased function | Moderate toxicity risk. ~4% Europeans |

**CRITICAL:** DPYD testing before fluoropyrimidine chemotherapy is now MANDATED in the EU and
recommended by CPIC/NCCN. A single dose of 5-FU can be lethal in DPYD-deficient patients.

### 3G. TPMT (Azathioprine, 6-Mercaptopurine, Thioguanine)

| rsID | Star Allele | Function | Clinical Impact |
|------|-------------|----------|-----------------|
| rs1800462 | *2 | No function | Rare |
| rs1800460 | *3B | No function | Part of *3A haplotype |
| rs1142345 | *3C | No function | Part of *3A haplotype. Most common LoF allele |

**Notes:** ~10% of population is heterozygous (intermediate metabolizer). ~0.3% are homozygous
poor metabolizers -- standard thiopurine doses can cause FATAL myelosuppression.

### 3H. NUDT15 (Azathioprine, 6-Mercaptopurine)

| rsID | Gene | Variant | Function | Clinical Impact |
|------|------|---------|----------|-----------------|
| rs116855232 | NUDT15 | c.415C>T (R139C) | No function (*3) | 7.86x risk of leukopenia with thiopurines. ~10% in East Asians, rare in Europeans |

### 3I. UGT1A1 (Irinotecan, Atazanavir)

| rsID | Gene | Variant | Function | Clinical Impact |
|------|------|---------|----------|-----------------|
| rs8175347 | UGT1A1 | TA repeat in promoter | *28 = 7 repeats (decreased) | Reduced clearance of irinotecan. ~10% Europeans are *28/*28 |

**Note:** rs8175347 is a microsatellite (TA repeat), not a simple SNP. Many VCF files may not
capture this correctly. Some arrays report it, WGS typically can call it.

### 3J. SLCO1B1 (Statins, especially Simvastatin)

| rsID | Gene | Variant | Star Allele | Clinical Impact |
|------|------|---------|-------------|-----------------|
| rs4149056 | SLCO1B1 | c.521T>C (V174A) | *5 | CC genotype: 17x risk simvastatin myopathy. ~15% Europeans carry C allele |

**Highly actionable:** CC carriers should avoid simvastatin or use lowest dose. Consider
alternative statins (rosuvastatin, pravastatin).

### 3K. IFNL3/IL28B (Historical - Hepatitis C)

| rsID | Gene | Variant | Clinical Impact |
|------|------|---------|-----------------|
| rs12979860 | IFNL3 (IL28B) | C>T | CC genotype: 2-3x better response to interferon-based HCV treatment. Less relevant in DAA era |

**Note:** With direct-acting antivirals (DAAs), this variant is less clinically impactful than
it was historically, but still informative for HCV clearance probability.

### Summary of CPIC Level A Gene-Drug Pairs (Key Examples)

| Gene | Drug(s) | Level | Clinical Action |
|------|---------|-------|-----------------|
| CYP2D6 | Codeine, tramadol | A | Avoid in PM and UM |
| CYP2D6 | Tamoxifen | A | Avoid in PM, use aromatase inhibitor |
| CYP2D6 | SSRIs (paroxetine, fluvoxamine) | A | Dose adjustment or alternative |
| CYP2D6 | Ondansetron | A | Alternative in UM |
| CYP2C19 | Clopidogrel | A | Avoid in PM (use prasugrel/ticagrelor) |
| CYP2C19 | PPIs | A | Dose adjustment |
| CYP2C19 | Voriconazole | A | Dose adjustment or alternative |
| CYP2C19 | SSRIs (sertraline, escitalopram) | A | Dose adjustment |
| CYP2C9 + VKORC1 | Warfarin | A | Dose calculation algorithm |
| CYP3A5 | Tacrolimus | A | Dose adjustment |
| DPYD | 5-FU, capecitabine | A | Dose reduction or contraindicated in PM |
| TPMT + NUDT15 | Azathioprine, 6-MP | A | Dose reduction or contraindicated in PM |
| UGT1A1 | Irinotecan | A | Dose reduction in *28/*28 |
| SLCO1B1 | Simvastatin | A | Dose limit or alternative statin |
| CYP2B6 | Efavirenz | A | Dose adjustment |
| CFTR | Ivacaftor | A | Indicated based on specific variants |
| G6PD | Rasburicase | A | Contraindicated in deficiency |

---

## 4. COMMON WELL-STUDIED SNPs WITH STRONG CLINICAL EVIDENCE

### 4A. APOE - Alzheimer's Disease Risk

| rsID | Gene | Allele | Risk |
|------|------|--------|------|
| rs429358 | APOE | T>C defines epsilon-4 | See haplotype table below |
| rs7412 | APOE | C>T defines epsilon-2 | See haplotype table below |

**APOE Haplotype Determination (from two SNPs):**

| rs429358 | rs7412 | APOE Allele |
|----------|--------|-------------|
| T | T | epsilon-2 |
| T | C | epsilon-3 (most common, reference) |
| C | C | epsilon-4 (risk allele) |

**Risk by Genotype (relative to e3/e3):**

| Genotype | Alzheimer's OR | Lifetime Risk (approx.) | Actionable? |
|----------|---------------|------------------------|-------------|
| e2/e2 | 0.6x | ~5% | Protective |
| e2/e3 | 0.6x | ~5% | Protective |
| e3/e3 | 1.0x (reference) | ~10-12% | Reference |
| e3/e4 | 3.2-3.7x | ~25-30% | MODERATE - lifestyle, monitoring |
| e4/e4 | 12-15x | ~50-60% | YES - aggressive risk reduction |
| e2/e4 | 2.6x | ~20% | Moderate |

**Notes:** APOE e4/e4 is now considered by some researchers as "semi-deterministic" rather than
merely a risk factor. Recent evidence suggests it may function closer to a monogenic cause in
homozygotes. APOE e4 also increases cardiovascular disease risk.

### 4B. MTHFR - Folate Metabolism

| rsID | Gene | Variant | Effect | Clinical Significance |
|------|------|---------|--------|----------------------|
| rs1801133 | MTHFR | c.677C>T (A222V) | TT: 30% enzyme activity | Elevated homocysteine, NTD risk |
| rs1801131 | MTHFR | c.1298A>C (E429A) | Reduced activity | Weaker effect, LD with C677T |

**Clinical reality:** CDC states that MTHFR genotyping is NOT recommended for clinical decision-
making. All MTHFR genotypes respond to standard folic acid supplementation (400 mcg/day). The
TT genotype results in only ~16% lower blood folate levels. However, knowing the genotype can
motivate individuals to ensure adequate folate intake.

**Penetrance for hyperhomocysteinemia:** TT genotype + low folate = elevated homocysteine
(~25% higher). With adequate folate, homocysteine is typically normal.

### 4C. Thrombophilia Panel (expanded)

| rsID | Gene | Variant | Condition | OR (Het) | OR (Hom) | Actionable |
|------|------|---------|-----------|----------|----------|------------|
| rs6025 | F5 | R506Q | Factor V Leiden | 3-5x VTE | 18-80x VTE | YES |
| rs1799963 | F2 | G20210A | Prothrombin mutation | 2-3x VTE | 5-10x VTE | YES |
| rs1800562 | HFE | C282Y | (see hemochromatosis) | -- | -- | YES |
| rs1799945 | HFE | H63D | (see hemochromatosis) | -- | -- | MODERATE |

### 4D. Hereditary Hemochromatosis (expanded)

See Section 2B above for full details.

---

## 5. FTD AND NEURODEGENERATION VARIANTS

### 5A. LRRK2 - Parkinson's Disease

| rsID | Gene | Variant | Condition | Penetrance | Actionable |
|------|------|---------|-----------|------------|------------|
| rs34637584 | LRRK2 | c.6055G>A (G2019S) | Parkinson's disease | 25-42.5% by age 80 | EMERGING - clinical trials available (LRRK2 inhibitors) |

**Notes:** Most common genetic cause of Parkinson's. ~1% of PD cases in Europeans, much higher
in Ashkenazi Jews (~15-20% of PD) and North African Berbers (~30-40% of PD). Age-dependent
penetrance. Polygenic risk score modifies penetrance.

### 5B. GBA1 - Parkinson's Disease / Gaucher Disease

| rsID | Gene | Variant | Condition | OR for PD | Actionable |
|------|------|---------|-----------|-----------|------------|
| rs76763715 | GBA1 | c.1226A>G (N370S/N409S) | PD risk / Gaucher carrier | OR 3.3-5.6 | EMERGING - GCase-targeted therapies in trials |
| rs421016 | GBA1 | c.1448T>C (L444P/L483P) | PD risk / Gaucher carrier | OR 4.9-9.7 | EMERGING |
| rs2230288 | GBA1 | c.1093G>A (E326K/E365K) | PD risk (not Gaucher) | OR ~2.0 | INFORMATIONAL |

**Penetrance for PD:** Age-dependent: ~2% by age 60-65, ~8-11% by age 80-85, ~30% by age 80
(varies by study). N370S and L444P carriers have similar penetrance.

**Note on Gaucher disease:** N370S and L444P are the two most common Gaucher disease mutations.
HOMOZYGOUS N370S causes type 1 (non-neuronopathic) Gaucher. Being a carrier (heterozygous)
increases PD risk without causing Gaucher.

### 5C. MAPT - Frontotemporal Dementia

| rsID | Gene | Variant | Condition | Penetrance | VCF Detectable |
|------|------|---------|-----------|------------|----------------|
| rs63751273 | MAPT | c.902C>T (P301L) | FTD | ~100% (AD) | YES - if in VCF |
| rs63750424 | MAPT | c.1216C>T (R406W) | FTD / AD-like | ~100% (AD) | YES - if in VCF |
| rs143624519 | MAPT | c.454G>A (A152T) | FTD/AD risk factor | OR 2.3-3.0 (NOT fully penetrant) | YES |

**Notes:** MAPT P301L is a rare, highly penetrant FTD mutation (mean onset ~51 years). A152T is
a risk factor, not deterministic. These are rare variants -- they will only be present in WGS
data, not consumer SNP arrays.

### 5D. GRN - Frontotemporal Dementia

| rsID | Gene | Variant | Condition | Penetrance | VCF Detectable |
|------|------|---------|-----------|------------|----------------|
| Various (>70 known) | GRN | Mostly LOF (nonsense, frameshift, splice) | FTD | ~90% by age 70 | Depends on variant |

**Notes:** GRN mutations cause FTD via haploinsufficiency. Most are private (family-specific)
loss-of-function mutations. No single common rsID. Best approach: cross-reference VCF against
ClinVar for all GRN pathogenic variants.

### 5E. C9orf72 - FTD/ALS

**NOT RELIABLY DETECTABLE FROM STANDARD VCF FILES.**

The C9orf72 pathogenic variant is a GGGGCC hexanucleotide repeat expansion (typically >30 repeats,
often hundreds to thousands). This CANNOT be detected by:
- SNP arrays (23andMe, Ancestry, etc.)
- Standard short-read WGS variant calling
- Most VCF files

**Detection requires:** Repeat-primed PCR, Southern blot, or targeted long-read sequencing.
Some specialized bioinformatics tools (ExpansionHunter, STRipy) can estimate repeat length
from short-read WGS BAM files, but this is NOT in standard VCF output.

---

## 6. DATABASES AND RESOURCES

### 6A. ClinVar (Primary Resource)

- **URL:** https://www.ncbi.nlm.nih.gov/clinvar/
- **Download VCF:** `wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`
- **GRCh37:** `wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz`
- **Contents:** ~2.5M variant submissions. Includes clinical significance (Pathogenic, Likely
  Pathogenic, VUS, Benign, etc.), disease associations, review status (star rating)
- **Update frequency:** Weekly (Mondays)
- **Key fields in VCF:** CLNSIG (significance), CLNDN (disease name), CLNREVSTAT (review status)
- **Filtering strategy:** Filter for CLNSIG = "Pathogenic" or "Likely_pathogenic" AND
  CLNREVSTAT includes "reviewed_by_expert_panel" or "criteria_provided,_multiple_submitters"
- **LIMITATION:** ClinVar VCF only contains variants with dbSNP rsIDs. Not all ClinVar
  variants have rsIDs. The XML or TSV downloads are more complete.

### 6B. PharmGKB / CPIC

- **URL:** https://www.pharmgkb.org/ and https://cpicpgx.org/
- **Downloads:** https://www.pharmgkb.org/downloads (Creative Commons license)
- **PharmCAT:** https://github.com/PharmGKB/PharmCAT - Tool that takes a VCF and produces
  a pharmacogenomics report with CPIC-level recommendations
- **Contents:** Clinical annotations linking variants to drug response, dosing guidelines,
  drug labels
- **Key data:** Gene-drug pair spreadsheet at https://files.cpicpgx.org/data/report/current/pair/cpic_gene-drug_pairs.xlsx
- **Allele definitions:** https://files.cpicpgx.org/data/report/current/allele_function_reference/

### 6C. SNPedia

- **URL:** https://www.snpedia.com/
- **Contents:** ~111,700 SNPs with functional annotations from peer-reviewed literature
- **Key features:** "Magnitude" score (0-10 significance), "Repute" (Good/Bad)
- **Integration:** Used by Promethease for consumer genome interpretation
- **API:** Available for bulk queries
- **Best for:** Quick lookup of individual rsIDs with plain-language summaries

### 6D. Additional Resources

| Resource | URL | Best For |
|----------|-----|---------|
| dbSNP | https://www.ncbi.nlm.nih.gov/snp/ | rsID lookup, allele frequencies |
| gnomAD | https://gnomad.broadinstitute.org/ | Population allele frequencies |
| ClinGen | https://clinicalgenome.org/ | Gene-disease validity, variant curation |
| OMIM | https://omim.org/ | Gene-disease relationships |
| InSiGHT (Lynch) | https://www.insight-group.org/ | MMR gene variant database |
| LOVD | https://www.lovd.nl/ | Locus-specific variant databases |
| PharmVar | https://www.pharmvar.org/ | Star allele definitions for PGx genes |
| HGMD (commercial) | https://www.hgmd.cf.ac.uk/ | Comprehensive mutation database (requires license) |

---

## 7. TECHNICAL LIMITATIONS OF VCF-BASED SCREENING

### What CAN be detected from a standard VCF:
1. **Single nucleotide variants (SNVs)** - all rsID-based lookups
2. **Small indels** (typically <50bp) - includes BRCA founder mutations
3. **Multi-nucleotide variants** if called by the variant caller
4. **Most pharmacogenomics SNPs** (except structural variants)
5. **Common risk alleles** with known rsIDs

### What CANNOT be reliably detected:
1. **Repeat expansions** - C9orf72, HTT (Huntington), DMPK (myotonic dystrophy), FMR1 (Fragile X)
2. **Copy number variants** - CYP2D6 gene deletion (*5) and duplications, EPCAM deletions,
   large BRCA1/2 rearrangements, SMN1 deletions (SMA)
3. **Structural variants** - inversions, translocations
4. **Methylation changes** - MLH1 promoter methylation (Lynch)
5. **Mosaic variants** at low allele fraction
6. **Phasing limitations** - Cannot determine if two variants are on the same or different
   chromosomes (cis vs trans) without phased VCF or family data

### Consumer Genomics Platform Limitations:
- **23andMe/Ancestry (SNP arrays):** ~650K-700K SNPs. Good for common PGx variants, APOE,
  Factor V Leiden, HFE. Limited for rare pathogenic variants in ACMG genes.
  23andMe tests 44 BRCA variants with FDA clearance.
- **Nebula Genomics:** Low-pass WGS (~0.4x). Good coverage but imputation-based. May miss
  rare variants. Better than arrays for structural assessment.
- **Whole Genome Sequencing (30x+):** Best for comprehensive screening. Can detect all SNVs,
  small indels. With specialized tools can estimate some repeat expansions.

### Recommended Approach for Python Script:

1. **Layer 1 - Direct rsID lookup:** Check all variants in this document by rsID
   (fast, works with any VCF)
2. **Layer 2 - ClinVar cross-reference:** Download ClinVar VCF, intersect with personal VCF
   for all P/LP variants in ACMG genes
3. **Layer 3 - PharmCAT:** Run PharmCAT for comprehensive pharmacogenomics
4. **Layer 4 - Position-based lookup:** For variants without rsIDs, check by chr:pos:ref:alt

---

## APPENDIX A: QUICK-REFERENCE rsID TABLE FOR DIRECT VCF LOOKUP

This is the master list of specific rsIDs suitable for direct lookup in any VCF file.

### Tier 1: Highest Clinical Impact (always check these)

| rsID | Gene | Condition | Category |
|------|------|-----------|----------|
| rs429358 | APOE | Alzheimer's risk (e4 allele) | Risk/Neuro |
| rs7412 | APOE | Alzheimer's risk (e2 allele) | Risk/Neuro |
| rs6025 | F5 | Factor V Leiden | Thrombophilia |
| rs1799963 | F2 | Prothrombin G20210A | Thrombophilia |
| rs1800562 | HFE | Hemochromatosis C282Y | Metabolism |
| rs1799945 | HFE | Hemochromatosis H63D | Metabolism |
| rs80357914 | BRCA1 | 185delAG (Ashkenazi) | Cancer |
| rs80357906 | BRCA1 | 5382insC (Ashkenazi) | Cancer |
| rs80359550 | BRCA2 | 6174delT (Ashkenazi) | Cancer |
| rs3918290 | DPYD | *2A - 5-FU toxicity | PGx CRITICAL |
| rs55886062 | DPYD | *13 - 5-FU toxicity | PGx CRITICAL |
| rs67376798 | DPYD | D949V - 5-FU toxicity | PGx |
| rs75017182 | DPYD | HapB3 - 5-FU toxicity | PGx |
| rs4149056 | SLCO1B1 | Statin myopathy | PGx |
| rs9923231 | VKORC1 | Warfarin sensitivity | PGx |
| rs4244285 | CYP2C19 | *2 - clopidogrel resistance | PGx |
| rs4986893 | CYP2C19 | *3 - clopidogrel resistance | PGx |
| rs12248560 | CYP2C19 | *17 - ultrarapid metabolism | PGx |
| rs3892097 | CYP2D6 | *4 - poor metabolizer | PGx |
| rs1799853 | CYP2C9 | *2 - warfarin sensitivity | PGx |
| rs1057910 | CYP2C9 | *3 - warfarin sensitivity | PGx |
| rs28929474 | SERPINA1 | Pi*Z - A1AT deficiency | Metabolism |

### Tier 2: Important Clinical Variants

| rsID | Gene | Condition | Category |
|------|------|-----------|----------|
| rs17580 | SERPINA1 | Pi*S - A1AT deficiency | Metabolism |
| rs35742686 | CYP2D6 | *3 - no function | PGx |
| rs5030655 | CYP2D6 | *6 - no function | PGx |
| rs1065852 | CYP2D6 | *10 - decreased function | PGx |
| rs28371706 | CYP2D6 | *17 - decreased function | PGx |
| rs28371725 | CYP2D6 | *41 - decreased function | PGx |
| rs16947 | CYP2D6 | *2 - normal/decreased | PGx |
| rs5030656 | CYP2D6 | *9 - decreased function | PGx |
| rs28399504 | CYP2C19 | *4 - no function | PGx |
| rs56337013 | CYP2C19 | *5 - no function | PGx |
| rs72552267 | CYP2C19 | *6 - no function | PGx |
| rs56165452 | CYP2C9 | *5 - no function | PGx |
| rs28371686 | CYP2C9 | *6 - no function | PGx |
| rs776746 | CYP3A5 | *3 - non-expressor | PGx |
| rs1800462 | TPMT | *2 - no function | PGx |
| rs1800460 | TPMT | *3B - no function | PGx |
| rs1142345 | TPMT | *3C - no function | PGx |
| rs116855232 | NUDT15 | *3 - thiopurine toxicity | PGx |
| rs8175347 | UGT1A1 | *28 - irinotecan toxicity | PGx |
| rs34637584 | LRRK2 | G2019S - Parkinson's risk | Neuro |
| rs76763715 | GBA1 | N370S - Parkinson's/Gaucher | Neuro |
| rs421016 | GBA1 | L444P - Parkinson's/Gaucher | Neuro |
| rs2230288 | GBA1 | E326K - Parkinson's risk | Neuro |
| rs5742904 | APOB | R3527Q - FH | Cardiovascular |
| rs11591147 | PCSK9 | R46L - protective (low LDL) | Cardiovascular |
| rs28362286 | PCSK9 | C679X - protective (low LDL) | Cardiovascular |
| rs12979860 | IFNL3 | IL28B - HCV response | PGx |
| rs1801133 | MTHFR | C677T - folate metabolism | Metabolism |
| rs1801131 | MTHFR | A1298C - folate metabolism | Metabolism |

### Tier 3: FTD/Neurodegeneration (rare, mainly WGS only)

| rsID | Gene | Variant | Condition | Notes |
|------|------|---------|-----------|-------|
| rs63751273 | MAPT | P301L | FTD | Rare, highly penetrant |
| rs63750424 | MAPT | R406W | FTD | Rare, highly penetrant |
| rs143624519 | MAPT | A152T | FTD/AD risk | OR 2.3-3.0, risk factor |

---

## APPENDIX B: STRATEGY FOR COMPREHENSIVE VCF SCREENING SCRIPT

### Recommended Architecture:

```
1. PARSE VCF → Build index by rsID AND by chr:pos:ref:alt

2. DIRECT rsID CHECK (this document)
   ├── Tier 1 variants (always report)
   ├── Tier 2 variants (report with context)
   └── Tier 3 variants (report if found)

3. CLINVAR CROSS-REFERENCE
   ├── Download: clinvar.vcf.gz (GRCh37 or GRCh38 to match your VCF)
   ├── Filter: CLNSIG in (Pathogenic, Likely_pathogenic)
   ├── Filter: Gene in ACMG_SF_v3.2_genes
   ├── Intersect with personal VCF
   └── Report matches with disease context

4. PHARMACOGENOMICS
   ├── Check all PGx rsIDs from Section 3
   ├── Determine star alleles (diplotype calling)
   ├── Map to metabolizer phenotype
   └── Report drug implications per CPIC

5. APOE HAPLOTYPING
   ├── Read rs429358 and rs7412 genotypes
   ├── Determine e2/e3/e4 diplotype
   └── Report risk stratification

6. REPORT GENERATION
   ├── Pathogenic findings (ClinVar)
   ├── PGx drug-gene interactions
   ├── Risk alleles (APOE, GBA, LRRK2)
   └── Carrier status (BRCA, HFE, GBA, etc.)
```

### Key Implementation Notes:
- Always normalize VCF chromosomes (strip "chr" prefix or add it to match reference)
- Handle multi-allelic sites (split or iterate)
- Consider VCF quality filters (GQ, DP) before reporting
- For PGx, default to *1 (wild-type) for alleles not tested
- ClinVar matching should be done by position+alleles, not just rsID (rsIDs can merge/change)
- Include genome build (GRCh37 vs GRCh38) awareness throughout

---

## Sources

- [ACMG SF v3.2 - Genetics in Medicine](https://www.gimjournal.org/article/S1098-3600(23)00879-1/fulltext)
- [ClinVar ACMG Gene Table](https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/)
- [ClinGen ACMG SF Genes](https://search.clinicalgenome.org/kb/genes/acmgsf)
- [CPIC Genes-Drugs](https://cpicpgx.org/genes-drugs/)
- [PharmGKB](https://www.pharmgkb.org/)
- [PharmCAT - GitHub](https://github.com/PharmGKB/PharmCAT)
- [CPIC CYP2D6 Recommendations](https://cpicpgx.org/wp-content/uploads/2021/09/CPIC_CYP2D6_v2.pdf)
- [PharmVar CYP2C19](https://pmc.ncbi.nlm.nih.gov/articles/PMC7769975/)
- [CYP2D6 Overview - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK574601/)
- [Pharmacogenomics Overview - StatPearls](https://www.ncbi.nlm.nih.gov/books/NBK617055/)
- [BRCA Ashkenazi Founder Mutations](https://pmc.ncbi.nlm.nih.gov/articles/PMC9184654/)
- [23andMe BRCA FDA Clearance](https://investors.23andme.com/news-releases/news-release-details/23andme-granted-new-fda-clearance-report-additional-brca)
- [Factor V Leiden - GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1368/)
- [FVL/PT G20210A Thrombosis Risk](https://ashpublications.org/blood/article/143/23/2425/515331)
- [HFE Hemochromatosis Penetrance](https://pmc.ncbi.nlm.nih.gov/articles/PMC1880832/)
- [APOE and Alzheimer's](https://pmc.ncbi.nlm.nih.gov/articles/PMC8096522/)
- [LRRK2 G2019S Penetrance](https://pmc.ncbi.nlm.nih.gov/articles/PMC8975556/)
- [GBA and Parkinson's](https://pmc.ncbi.nlm.nih.gov/articles/PMC10548077/)
- [GBA E326K and PD Risk](https://pmc.ncbi.nlm.nih.gov/articles/PMC5901859/)
- [MAPT P301L - Alzforum](https://www.alzforum.org/mutations/mapt-p301l)
- [MAPT A152T Risk Factor](https://pmc.ncbi.nlm.nih.gov/articles/PMC3657164/)
- [Alpha-1 Antitrypsin Deficiency - GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1519/)
- [SERPINA1 Variants](https://pmc.ncbi.nlm.nih.gov/articles/PMC7997584/)
- [C9orf72 Repeat Expansion Detection Limitations](https://molecularneurodegeneration.biomedcentral.com/articles/10.1186/s13024-018-0274-4)
- [Lynch Syndrome - GeneReviews](https://www.ncbi.nlm.nih.gov/books/NBK1211/)
- [FH Variants LDLR/APOB/PCSK9](https://www.nature.com/articles/gim2017151)
- [MTHFR - CDC](https://www.cdc.gov/folic-acid/data-research/mthfr/index.html)
- [DPYD Genotyping Recommendations](https://www.jmdjournal.org/article/S1525-1578(24)00154-5/fulltext)
- [ClinVar Data Downloads](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
- [SNPedia](https://www.snpedia.com/)
