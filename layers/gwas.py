import logging
from collections import defaultdict
from typing import Dict, List

from core import GenomeData

# ============================================================================
# SECTION 8: LAYER 5 — GWAS TRAIT ASSOCIATIONS
# ============================================================================

# Curated GWAS trait associations — top validated SNPs from NHGRI-EBI GWAS Catalog
# Selected: genome-wide significant (p < 5e-8), well-replicated, on common microarrays
# Format: rsID → list of {trait, category, risk_allele, or_beta, pvalue, gene, source}
GWAS_CURATED = {
    # === Cardiovascular ===
    "rs10455872": [{"trait": "Coronary artery disease", "cat": "Cardiovascular",
                    "risk": "G", "or": 1.63, "p": 3e-15, "gene": "LPA"}],
    "rs4977574": [{"trait": "Coronary artery disease", "cat": "Cardiovascular",
                   "risk": "G", "or": 1.33, "p": 2e-58, "gene": "CDKN2A/B (9p21)"}],
    "rs1333049": [{"trait": "Coronary artery disease", "cat": "Cardiovascular",
                   "risk": "C", "or": 1.36, "p": 1e-44, "gene": "9p21.3"}],
    "rs3184504": [{"trait": "Coronary artery disease / Blood pressure", "cat": "Cardiovascular",
                   "risk": "T", "or": 1.13, "p": 1e-13, "gene": "SH2B3"}],
    "rs1746048": [{"trait": "Coronary artery disease", "cat": "Cardiovascular",
                   "risk": "C", "or": 1.09, "p": 5e-9, "gene": "CXCL12"}],
    "rs17249754": [{"trait": "Hypertension / Blood pressure", "cat": "Cardiovascular",
                    "risk": "A", "or": 1.0, "p": 7e-24, "gene": "ATP2B1",
                    "note": "Beta = -0.9 mmHg per allele"}],
    "rs1799945": [{"trait": "Iron levels / Hemochromatosis", "cat": "Hematological",
                   "risk": "G", "or": 1.65, "p": 2e-75, "gene": "HFE (H63D)"}],
    "rs2200733": [{"trait": "Atrial fibrillation", "cat": "Cardiovascular",
                   "risk": "T", "or": 1.72, "p": 4e-29, "gene": "PITX2"}],
    "rs13376333": [{"trait": "Atrial fibrillation", "cat": "Cardiovascular",
                    "risk": "T", "or": 1.18, "p": 1e-9, "gene": "KCNN3"}],
    # === Metabolic ===
    "rs7903146": [{"trait": "Type 2 diabetes", "cat": "Metabolic",
                   "risk": "T", "or": 1.40, "p": 1e-150, "gene": "TCF7L2"}],
    "rs1801282": [{"trait": "Type 2 diabetes (protective)", "cat": "Metabolic",
                   "risk": "C", "or": 0.87, "p": 3e-17, "gene": "PPARG (Pro12Ala)"}],
    "rs5219": [{"trait": "Type 2 diabetes", "cat": "Metabolic",
                "risk": "T", "or": 1.15, "p": 7e-17, "gene": "KCNJ11 (E23K)"}],
    "rs1558902": [{"trait": "Body mass index", "cat": "Metabolic",
                   "risk": "A", "or": 1.0, "p": 5e-120, "gene": "FTO",
                   "note": "Beta = +0.37 kg/m² per allele"}],
    "rs6567160": [{"trait": "Body mass index", "cat": "Metabolic",
                   "risk": "C", "or": 1.0, "p": 2e-48, "gene": "MC4R",
                   "note": "Beta = +0.30 kg/m² per allele"}],
    "rs8050136": [{"trait": "Obesity", "cat": "Metabolic",
                   "risk": "A", "or": 1.22, "p": 3e-20, "gene": "FTO"}],
    "rs13266634": [{"trait": "Type 2 diabetes", "cat": "Metabolic",
                    "risk": "C", "or": 0.88, "p": 5e-14, "gene": "SLC30A8"}],
    "rs4988235": [{"trait": "Lactose intolerance", "cat": "Metabolic",
                   "risk": "G", "or": 8.0, "p": 3e-180, "gene": "LCT/MCM6",
                   "note": "GG = lactose intolerant in EUR"}],
    # === Cancer ===
    "rs10993994": [{"trait": "Prostate cancer", "cat": "Cancer",
                    "risk": "T", "or": 1.25, "p": 8e-29, "gene": "MSMB"}],
    "rs6983267": [{"trait": "Colorectal cancer / Prostate cancer", "cat": "Cancer",
                   "risk": "G", "or": 1.21, "p": 1e-18, "gene": "8q24"}],
    "rs1447295": [{"trait": "Prostate cancer", "cat": "Cancer",
                   "risk": "A", "or": 1.28, "p": 7e-16, "gene": "8q24"}],
    "rs10490924": [{"trait": "Age-related macular degeneration", "cat": "Ophthalmological",
                    "risk": "T", "or": 2.69, "p": 2e-89, "gene": "ARMS2 (A69S)"}],
    "rs1061170": [{"trait": "Age-related macular degeneration", "cat": "Ophthalmological",
                   "risk": "C", "or": 2.45, "p": 4e-74, "gene": "CFH (Y402H)"}],
    "rs11571833": [{"trait": "Breast cancer", "cat": "Cancer",
                    "risk": "T", "or": 1.26, "p": 2e-8, "gene": "BRCA2 (K3326X)"}],
    # === Immune / Autoimmune ===
    "rs2187668": [{"trait": "Celiac disease", "cat": "Immune/Autoimmune",
                   "risk": "T", "or": 7.0, "p": 1e-300, "gene": "HLA-DQ2.5"}],
    "rs6457620": [{"trait": "Rheumatoid arthritis", "cat": "Immune/Autoimmune",
                   "risk": "T", "or": 2.07, "p": 1e-120, "gene": "HLA-DRB1"}],
    "rs2476601": [{"trait": "Type 1 diabetes / Rheumatoid arthritis", "cat": "Immune/Autoimmune",
                   "risk": "A", "or": 1.76, "p": 8e-35, "gene": "PTPN22 (R620W)"}],
    "rs3135388": [{"trait": "Multiple sclerosis", "cat": "Immune/Autoimmune",
                   "risk": "A", "or": 2.97, "p": 1e-100, "gene": "HLA-DRB1*15:01"}],
    "rs2395029": [{"trait": "Psoriasis", "cat": "Immune/Autoimmune",
                   "risk": "G", "or": 4.1, "p": 3e-36, "gene": "HLA-C"}],
    "rs11209026": [{"trait": "Crohn disease / IBD", "cat": "Immune/Autoimmune",
                    "risk": "A", "or": 0.43, "p": 2e-30, "gene": "IL23R (R381Q)",
                    "note": "Protective variant"}],
    "rs2066847": [{"trait": "Crohn disease", "cat": "Immune/Autoimmune",
                   "risk": "C", "or": 3.1, "p": 2e-40, "gene": "NOD2 (L1007fs)"}],
    # === Neurological ===
    "rs429358": [{"trait": "Alzheimer disease", "cat": "Neurological",
                  "risk": "C", "or": 3.68, "p": 1e-500, "gene": "APOE (e4)"}],
    "rs7412": [{"trait": "Alzheimer disease (protective)", "cat": "Neurological",
                "risk": "T", "or": 0.62, "p": 5e-50, "gene": "APOE (e2)"}],
    "rs6656401": [{"trait": "Alzheimer disease", "cat": "Neurological",
                   "risk": "A", "or": 1.18, "p": 3e-14, "gene": "CR1"}],
    "rs6733839": [{"trait": "Alzheimer disease", "cat": "Neurological",
                   "risk": "T", "or": 1.20, "p": 2e-27, "gene": "BIN1"}],
    "rs2514218": [{"trait": "Schizophrenia", "cat": "Neurological",
                   "risk": "T", "or": 1.08, "p": 4e-12, "gene": "DRD2"}],
    # === Anthropometric ===
    "rs1805007": [{"trait": "Red hair / fair skin", "cat": "Anthropometric",
                   "risk": "T", "or": 5.7, "p": 3e-60, "gene": "MC1R (R151C)"}],
    "rs1805008": [{"trait": "Red hair / fair skin", "cat": "Anthropometric",
                   "risk": "T", "or": 3.2, "p": 1e-30, "gene": "MC1R (R160W)"}],
    "rs12913832": [{"trait": "Eye color (blue vs brown)", "cat": "Anthropometric",
                    "risk": "G", "or": 30.0, "p": 1e-300, "gene": "HERC2/OCA2",
                    "note": "GG = blue eyes with ~80% probability"}],
    "rs16891982": [{"trait": "Skin pigmentation", "cat": "Anthropometric",
                    "risk": "G", "or": 7.0, "p": 3e-50, "gene": "SLC45A2 (L374F)"}],
    "rs1426654": [{"trait": "Skin pigmentation", "cat": "Anthropometric",
                   "risk": "A", "or": 30.0, "p": 1e-100, "gene": "SLC24A5 (A111T)",
                   "note": "Fixed in Europeans, absent in Africans"}],
    "rs4778138": [{"trait": "Eye color", "cat": "Anthropometric",
                   "risk": "A", "or": 2.0, "p": 5e-15, "gene": "OCA2"}],
    "rs12203592": [{"trait": "Freckles / hair color", "cat": "Anthropometric",
                    "risk": "T", "or": 2.9, "p": 8e-40, "gene": "IRF4"}],
    "rs1800407": [{"trait": "Eye color (blue-green)", "cat": "Anthropometric",
                   "risk": "T", "or": 2.5, "p": 3e-20, "gene": "OCA2 (R419Q)"}],
    "rs17822931": [{"trait": "Earwax type (wet/dry)", "cat": "Anthropometric",
                    "risk": "T", "or": 50.0, "p": 1e-100, "gene": "ABCC11",
                    "note": "CC = dry earwax (common in East Asians)"}],
    "rs4481887": [{"trait": "Asparagus anosmia", "cat": "Anthropometric",
                   "risk": "A", "or": 1.4, "p": 2e-10, "gene": "OR2M7"}],
    "rs713598": [{"trait": "Bitter taste perception (PTC)", "cat": "Anthropometric",
                  "risk": "G", "or": 5.0, "p": 1e-30, "gene": "TAS2R38"}],
    # === Hematological ===
    "rs855791": [{"trait": "Iron levels / Hemoglobin", "cat": "Hematological",
                  "risk": "A", "or": 1.0, "p": 2e-65, "gene": "TMPRSS6",
                  "note": "Beta = -0.15 SD iron per allele"}],
    "rs1800562": [{"trait": "Hemochromatosis / High iron", "cat": "Hematological",
                   "risk": "A", "or": 4.5, "p": 1e-100, "gene": "HFE (C282Y)"}],
    "rs28929474": [{"trait": "Alpha-1 antitrypsin deficiency", "cat": "Respiratory",
                    "risk": "T", "or": 5.0, "p": 1e-40, "gene": "SERPINA1 (Z allele)"}],
    # === Renal ===
    "rs2231142": [{"trait": "Gout / Uric acid levels", "cat": "Renal",
                   "risk": "T", "or": 1.73, "p": 2e-90, "gene": "ABCG2 (Q141K)"}],
    "rs12498742": [{"trait": "Uric acid / Gout", "cat": "Renal",
                    "risk": "A", "or": 1.0, "p": 5e-200, "gene": "SLC2A9",
                    "note": "Strongest genetic determinant of serum urate"}],
    # === Pharmacogenomic ===
    "rs9923231": [{"trait": "Warfarin dose requirement", "cat": "Pharmacogenomic",
                   "risk": "T", "or": 1.0, "p": 1e-80, "gene": "VKORC1",
                   "note": "T allele = lower warfarin dose"}],
    "rs4149056": [{"trait": "Statin myopathy risk", "cat": "Pharmacogenomic",
                   "risk": "C", "or": 4.5, "p": 4e-9, "gene": "SLCO1B1",
                   "note": "CC genotype = highest simvastatin myopathy risk"}],
    "rs12248560": [{"trait": "Clopidogrel response", "cat": "Pharmacogenomic",
                    "risk": "T", "or": 1.0, "p": 1e-10, "gene": "CYP2C19*17",
                    "note": "Ultrarapid metabolizer"}],
    "rs4244285": [{"trait": "Clopidogrel poor response", "cat": "Pharmacogenomic",
                   "risk": "A", "or": 1.0, "p": 1e-15, "gene": "CYP2C19*2",
                   "note": "Reduced efficacy of clopidogrel"}],
    "rs762551": [{"trait": "Caffeine metabolism", "cat": "Pharmacogenomic",
                  "risk": "C", "or": 1.0, "p": 5e-30, "gene": "CYP1A2*1F",
                  "note": "CC = slow caffeine metabolism"}],
    "rs671": [{"trait": "Alcohol flush reaction", "cat": "Pharmacogenomic",
               "risk": "A", "or": 75.0, "p": 1e-200, "gene": "ALDH2*2",
               "note": "Nearly exclusive to East Asian populations"}],
    # === Hepatic ===
    "rs738409": [{"trait": "Non-alcoholic fatty liver disease", "cat": "Hepatic",
                  "risk": "G", "or": 3.26, "p": 5e-21, "gene": "PNPLA3 (I148M)"}],
    "rs4149056": [{"trait": "Bilirubin levels", "cat": "Hepatic",
                   "risk": "C", "or": 1.0, "p": 1e-20, "gene": "SLCO1B1"}],
    # === Sleep ===
    "rs1801260": [{"trait": "Chronotype (morning/evening preference)", "cat": "Behavioral",
                   "risk": "C", "or": 1.3, "p": 2e-10, "gene": "CLOCK (3111T>C)"}],
    # === Longevity ===
    "rs2802292": [{"trait": "Longevity / Lifespan", "cat": "Longevity",
                   "risk": "G", "or": 1.17, "p": 2e-10, "gene": "FOXO3"}],
}


def match_gwas_curated(genome: GenomeData, logger: logging.Logger) -> List[Dict]:
    """Match user variants against curated GWAS database."""
    results = []

    for rsid, associations in GWAS_CURATED.items():
        if rsid not in genome.rsid_index:
            continue

        key = genome.rsid_index[rsid]
        variant = genome.variants[key]
        a1, a2 = variant.genotype_alleles

        for assoc in associations:
            risk_allele = assoc["risk"]
            dosage = 0
            if a1.upper() == risk_allele.upper():
                dosage += 1
            if a2.upper() == risk_allele.upper():
                dosage += 1

            results.append({
                "rsid": rsid,
                "gene": assoc["gene"],
                "trait": assoc["trait"],
                "category": assoc["cat"],
                "risk_allele": risk_allele,
                "or_beta": assoc["or"],
                "pvalue": assoc["p"],
                "dosage": dosage,
                "genotype": f"{a1}/{a2}",
                "note": assoc.get("note", ""),
            })

    logger.info(f"  Matched {len(results)} GWAS associations from "
                f"{len(GWAS_CURATED)} curated variants")
    return results


def run_layer5_gwas(genome: GenomeData, cache_dir: str,
                    logger: logging.Logger) -> Dict:
    """Run GWAS analysis using curated database. Returns results dict."""
    logger.info("=" * 60)
    logger.info("LAYER 5: GWAS Trait Associations (Curated)")
    logger.info("=" * 60)

    results = match_gwas_curated(genome, logger)

    # Separate carriers vs non-carriers
    carriers = [r for r in results if r["dosage"] > 0]
    non_carriers = [r for r in results if r["dosage"] == 0]

    # Group carriers by category
    categories = defaultdict(list)
    for r in carriers:
        categories[r["category"]].append(r)

    logger.info(f"Layer 5 complete: {len(carriers)} carrier associations across {len(categories)} categories")
    return {"layer": 5, "total": len(carriers),
            "categories": {k: len(v) for k, v in categories.items()},
            "results": carriers}
