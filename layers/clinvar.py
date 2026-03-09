import gzip
import logging
import os
from typing import Dict, List, Optional, Set

from core import GenomeData, download_file, normalize_chrom

# ============================================================================
# SECTION 6: LAYER 3 — CLINVAR PATHOGENIC VARIANTS
# ============================================================================

# Curated database of clinically significant variants commonly found on
# consumer genotyping microarrays. Sources: ClinVar, OMIM, ACMG SF v3.2
# Format: rsID → {gene, significance, conditions, risk_allele, zygosity_note}
CLINVAR_CURATED = {
    # === ACMG Actionable Genes (Secondary Findings v3.2) ===
    # Cardiac
    "rs121913626": {"gene": "MYBPC3", "sig": "Pathogenic", "cond": "Hypertrophic cardiomyopathy",
                    "risk": "T", "note": "Common HCM mutation in European populations"},
    "rs3218713": {"gene": "MYBPC3", "sig": "Pathogenic", "cond": "Hypertrophic cardiomyopathy",
                  "risk": "A", "note": "HCM founder mutation"},
    "rs397516083": {"gene": "SCN5A", "sig": "Pathogenic", "cond": "Brugada syndrome / Long QT",
                    "risk": "A", "note": "Cardiac arrhythmia risk"},
    "rs199473282": {"gene": "KCNQ1", "sig": "Pathogenic", "cond": "Long QT syndrome type 1",
                    "risk": "T", "note": "Cardiac arrhythmia risk"},
    "rs199473656": {"gene": "KCNH2", "sig": "Pathogenic", "cond": "Long QT syndrome type 2",
                    "risk": "A", "note": "Cardiac arrhythmia risk"},
    # Thrombophilia
    "rs6025": {"gene": "F5", "sig": "Risk factor", "cond": "Factor V Leiden thrombophilia",
               "risk": "T", "note": "VCF alleles C>T at chr1:169519049. ~5% EUR carrier freq"},
    "rs1799963": {"gene": "F2", "sig": "Risk factor", "cond": "Prothrombin G20210A thrombophilia",
                  "risk": "A", "note": "~2% EUR carrier freq. Venous thrombosis risk"},
    # Hemochromatosis
    "rs1800562": {"gene": "HFE", "sig": "Pathogenic", "cond": "Hereditary hemochromatosis",
                  "risk": "A", "note": "C282Y mutation. ~10% EUR carrier freq. Homozygous = disease"},
    "rs1799945": {"gene": "HFE", "sig": "Risk factor", "cond": "Hereditary hemochromatosis",
                  "risk": "G", "note": "H63D variant. Compound het with C282Y increases risk"},
    # BRCA
    "rs80357906": {"gene": "BRCA1", "sig": "Pathogenic", "cond": "Hereditary breast/ovarian cancer",
                   "risk": "T", "note": "5382insC Ashkenazi founder mutation"},
    "rs28897743": {"gene": "BRCA2", "sig": "Pathogenic", "cond": "Hereditary breast/ovarian cancer",
                   "risk": "T", "note": "6174delT Ashkenazi founder mutation"},
    "rs80357713": {"gene": "BRCA1", "sig": "Pathogenic", "cond": "Hereditary breast/ovarian cancer",
                   "risk": "A", "note": "185delAG Ashkenazi founder mutation"},
    # Cancer predisposition
    "rs121913039": {"gene": "TP53", "sig": "Pathogenic", "cond": "Li-Fraumeni syndrome",
                    "risk": "A", "note": "R248W hotspot. Very rare"},
    "rs121913343": {"gene": "TP53", "sig": "Pathogenic", "cond": "Li-Fraumeni syndrome",
                    "risk": "T", "note": "R175H hotspot. Very rare"},
    "rs11571833": {"gene": "BRCA2", "sig": "Pathogenic", "cond": "Breast cancer susceptibility",
                   "risk": "T", "note": "K3326X truncation. Moderate risk increase"},
    "rs1800056": {"gene": "ATM", "sig": "Likely pathogenic", "cond": "Breast cancer susceptibility",
                  "risk": "C", "note": "ATM missense variant. Risk allele is ALT=C"},
    # Lynch syndrome
    "rs63750447": {"gene": "MLH1", "sig": "Pathogenic", "cond": "Lynch syndrome (colorectal cancer)",
                   "risk": "A", "note": "Common MLH1 pathogenic variant"},
    "rs267607918": {"gene": "MSH2", "sig": "Pathogenic", "cond": "Lynch syndrome",
                    "risk": "T", "note": "MSH2 pathogenic variant"},

    # === Metabolic / Common health ===
    "rs4880": {"gene": "SOD2", "sig": "Drug response", "cond": "Oxidative stress susceptibility",
               "risk": "T", "note": "A16V (Ala/Val). Val allele = mitochondrial targeting change"},
    "rs4680": {"gene": "COMT", "sig": "Drug response", "cond": "Catechol-O-methyltransferase activity",
               "risk": "A", "note": "Val158Met. Met/Met = low COMT, higher dopamine"},
    "rs1801133": {"gene": "MTHFR", "sig": "Drug response", "cond": "MTHFR deficiency / folate metabolism",
                  "risk": "A", "note": "C677T. TT genotype = 70% reduced enzyme activity"},
    "rs1801131": {"gene": "MTHFR", "sig": "Drug response", "cond": "MTHFR deficiency",
                  "risk": "G", "note": "A1298C. Compound het with C677T relevant"},
    "rs1805007": {"gene": "MC1R", "sig": "Risk factor", "cond": "Red hair / skin cancer susceptibility",
                  "risk": "T", "note": "R151C. Red hair variant, increased melanoma risk"},
    "rs1805008": {"gene": "MC1R", "sig": "Risk factor", "cond": "Red hair / skin cancer susceptibility",
                  "risk": "T", "note": "R160W. Increased melanoma risk"},
    "rs4988235": {"gene": "MCM6/LCT", "sig": "Association", "cond": "Lactase persistence/non-persistence",
                  "risk": "G", "note": "C/T(-13910). CC = lactose intolerant (non-persistence)"},
    "rs182549": {"gene": "MCM6/LCT", "sig": "Association", "cond": "Lactase persistence",
                 "risk": "T", "note": "Secondary lactase persistence marker"},
    # Celiac
    "rs2187668": {"gene": "HLA-DQ2.5", "sig": "Risk factor", "cond": "Celiac disease",
                  "risk": "T", "note": "HLA-DQ2.5 tag SNP. Required but not sufficient"},
    "rs7454108": {"gene": "HLA-DQ8", "sig": "Risk factor", "cond": "Celiac disease",
                  "risk": "C", "note": "HLA-DQ8 tag SNP"},
    # Alpha-1 antitrypsin
    "rs28929474": {"gene": "SERPINA1", "sig": "Pathogenic", "cond": "Alpha-1 antitrypsin deficiency",
                   "risk": "T", "note": "Z allele (Pi*Z). 1 in 25 EUR carriers. COPD/liver risk"},
    "rs17580": {"gene": "SERPINA1", "sig": "Pathogenic", "cond": "Alpha-1 antitrypsin deficiency",
                "risk": "A", "note": "S allele (Pi*S). Milder than Z. Compound ZS = moderate risk"},
    # Glucose-6-phosphate dehydrogenase
    "rs1050828": {"gene": "G6PD", "sig": "Pathogenic", "cond": "G6PD deficiency (favism)",
                  "risk": "T", "note": "Val68Met (G6PD A-). Hemolytic anemia with oxidative stress"},
    # Cystic fibrosis
    "rs75039782": {"gene": "CFTR", "sig": "Pathogenic", "cond": "Cystic fibrosis",
                   "risk": "A", "note": "G551D mutation. Ivacaftor-responsive"},
    "rs113993960": {"gene": "CFTR", "sig": "Pathogenic", "cond": "Cystic fibrosis",
                    "risk": "ATCT", "note": "F508del. Most common CF mutation (70% of alleles)"},
    # Pharmacogenetics (complement Layer 2)
    "rs9923231": {"gene": "VKORC1", "sig": "Drug response", "cond": "Warfarin dose requirement",
                  "risk": "T", "note": "-1639G>A. T allele = lower warfarin dose needed"},
    "rs12248560": {"gene": "CYP2C19", "sig": "Drug response", "cond": "CYP2C19 ultrarapid metabolizer",
                   "risk": "T", "note": "*17 allele. Increased clopidogrel activation"},
    "rs4244285": {"gene": "CYP2C19", "sig": "Drug response", "cond": "CYP2C19 poor metabolizer",
                  "risk": "A", "note": "*2 allele. Reduced clopidogrel efficacy"},
    "rs3892097": {"gene": "CYP2D6", "sig": "Drug response", "cond": "CYP2D6 poor metabolizer",
                  "risk": "A", "note": "*4 allele. Reduced codeine/tramadol activation"},
    # Hearing loss
    "rs80338939": {"gene": "GJB2", "sig": "Pathogenic", "cond": "Autosomal recessive deafness",
                   "risk": "T", "note": "35delG. Most common hereditary hearing loss mutation in EUR"},
    # Wilson disease
    "rs76151636": {"gene": "ATP7B", "sig": "Pathogenic", "cond": "Wilson disease",
                   "risk": "A", "note": "H1069Q. Common EUR mutation for copper metabolism disorder"},
    # Familial Mediterranean Fever
    "rs61752717": {"gene": "MEFV", "sig": "Pathogenic", "cond": "Familial Mediterranean fever",
                   "risk": "A", "note": "M694V mutation. Common in Mediterranean populations"},
    # Sickle cell
    "rs334": {"gene": "HBB", "sig": "Pathogenic", "cond": "Sickle cell disease / trait",
              "risk": "T", "note": "HbS (E6V). Carrier = sickle cell trait. Homozygous = disease"},
    # APOE (Alzheimer)
    "rs429358": {"gene": "APOE", "sig": "Risk factor", "cond": "Alzheimer disease susceptibility",
                 "risk": "C", "note": "APOE e4 determinant. CC = e4/e4 (highest risk)"},
    "rs7412": {"gene": "APOE", "sig": "Protective", "cond": "Alzheimer disease",
               "risk": "T", "note": "APOE e2 determinant. e2 allele is protective"},
    # === FTD / Neurodegeneration (family history relevant) ===
    # MAPT — Frontotemporal Dementia
    "rs63751273": {"gene": "MAPT", "sig": "Pathogenic", "cond": "Frontotemporal dementia (FTDP-17)",
                   "risk": "T", "note": "P301L. ~100% penetrant. Mean onset ~51 years. Rare, WGS only"},
    "rs63750424": {"gene": "MAPT", "sig": "Pathogenic", "cond": "Frontotemporal dementia / AD-like",
                   "risk": "T", "note": "R406W. ~100% penetrant. Later onset than P301L"},
    "rs143624519": {"gene": "MAPT", "sig": "Risk factor", "cond": "FTD / Alzheimer risk",
                    "risk": "A", "note": "A152T. OR 2.3-3.0. Risk factor, not deterministic"},
    # LRRK2 — Parkinson's
    "rs34637584": {"gene": "LRRK2", "sig": "Pathogenic", "cond": "Parkinson disease",
                   "risk": "A", "note": "G2019S. 25-42% penetrance by age 80. Most common genetic cause of PD. ~1% EUR PD, ~15-20% Ashkenazi PD"},
    # GBA1 — Parkinson's / Gaucher carrier
    "rs76763715": {"gene": "GBA1", "sig": "Risk factor", "cond": "Parkinson disease / Gaucher carrier",
                   "risk": "G", "note": "N370S (N409S). OR 3.3-5.6 for PD. Homozygous = type 1 Gaucher disease"},
    "rs421016": {"gene": "GBA1", "sig": "Risk factor", "cond": "Parkinson disease / Gaucher carrier",
                 "risk": "C", "note": "L444P (L483P). OR 4.9-9.7 for PD. Homozygous = Gaucher disease"},
    "rs2230288": {"gene": "GBA1", "sig": "Risk factor", "cond": "Parkinson disease",
                  "risk": "A", "note": "E326K (E365K). OR ~2.0 for PD. Does NOT cause Gaucher"},
    # === Protective variants ===
    "rs11591147": {"gene": "PCSK9", "sig": "Protective", "cond": "Low LDL cholesterol",
                   "risk": "T", "note": "R46L loss-of-function. 15-28% lower LDL. Reduced cardiovascular risk"},
    "rs28362286": {"gene": "PCSK9", "sig": "Protective", "cond": "Low LDL cholesterol",
                   "risk": "A", "note": "C679X loss-of-function. ~40% lower LDL"},
    # Familial Hypercholesterolemia
    "rs5742904": {"gene": "APOB", "sig": "Pathogenic", "cond": "Familial hypercholesterolemia type 2",
                  "risk": "A", "note": "R3527Q. High penetrance for elevated LDL. Statin-responsive"},
    # === Additional critical PGx (complement Layer 2) ===
    "rs3918290": {"gene": "DPYD", "sig": "Drug response", "cond": "5-FU/capecitabine toxicity",
                  "risk": "G", "note": "*2A. LIFE-THREATENING toxicity. EU-mandated test before fluoropyrimidines"},
    "rs55886062": {"gene": "DPYD", "sig": "Drug response", "cond": "5-FU/capecitabine toxicity",
                   "risk": "C", "note": "*13. Severe toxicity risk"},
    "rs4149056": {"gene": "SLCO1B1", "sig": "Drug response", "cond": "Statin myopathy",
                  "risk": "C", "note": "*5. CC genotype: 17x risk of simvastatin myopathy"},
    "rs116855232": {"gene": "NUDT15", "sig": "Drug response", "cond": "Thiopurine toxicity",
                    "risk": "T", "note": "*3. 7.86x risk of leukopenia with azathioprine/6-MP. ~10% East Asian"},
    # Age-related macular degeneration
    "rs10490924": {"gene": "ARMS2", "sig": "Risk factor", "cond": "Age-related macular degeneration",
                   "risk": "T", "note": "A69S. Strong AMD risk factor, OR ~2.7"},
    "rs1061170": {"gene": "CFH", "sig": "Risk factor", "cond": "Age-related macular degeneration",
                  "risk": "C", "note": "Y402H. Complement factor H variant, OR ~2.5"},
    # Alcohol flush
    "rs671": {"gene": "ALDH2", "sig": "Drug response", "cond": "Alcohol sensitivity / flush reaction",
              "risk": "A", "note": "ALDH2*2. Very common in East Asian populations"},
    # Caffeine metabolism
    "rs762551": {"gene": "CYP1A2", "sig": "Drug response", "cond": "Caffeine metabolism",
                 "risk": "C", "note": "CYP1A2*1F. CC = slow caffeine metabolizer"},
    # Bitter taste
    "rs713598": {"gene": "TAS2R38", "sig": "Association", "cond": "Bitter taste perception (PTC)",
                 "risk": "G", "note": "PAV haplotype. G allele = taster"},
    "rs1726866": {"gene": "TAS2R38", "sig": "Association", "cond": "Bitter taste perception",
                  "risk": "G", "note": "PAV haplotype component"},
    # Asparagus anosmia
    "rs4481887": {"gene": "OR2M7", "sig": "Association", "cond": "Asparagus anosmia (smell)",
                  "risk": "A", "note": "Ability to smell asparagus metabolites in urine"},
    # ADD1 hypertension
    "rs4961": {"gene": "ADD1", "sig": "Risk factor", "cond": "Salt-sensitive hypertension",
               "risk": "T", "note": "G460W. T allele = increased salt sensitivity"},
}

# ACMG Secondary Findings v3.2 — 81+ genes for ClinVar VCF filtering
# Plus FTD-relevant genes (MAPT, GRN, LRRK2, GBA1) for family history analysis
ACMG_SF_GENES = {
    # Cancer predisposition (28)
    "APC", "BRCA1", "BRCA2", "PALB2", "TP53", "RET", "RB1", "PTEN", "STK11",
    "MLH1", "MSH2", "MSH6", "PMS2", "MUTYH", "MEN1", "VHL", "NF2",
    "SDHB", "SDHD", "SDHC", "SDHAF2", "MAX", "TMEM127",
    "TSC1", "TSC2", "WT1", "BMPR1A", "SMAD4",
    # Cardiovascular (40)
    "MYH7", "MYBPC3", "TNNT2", "TPM1", "TNNI3", "MYL3", "MYL2", "ACTC1",
    "PRKAG2", "LMNA", "DES", "TNNC1", "RBM20", "BAG3", "TTN", "FLNC", "SCN5A",
    "DSP", "PKP2", "DSG2", "DSC2", "TMEM43",
    "KCNQ1", "KCNH2", "RYR2", "CASQ2", "TRDN", "CALM1", "CALM2", "CALM3",
    "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "COL3A1", "MYH11", "ACTA2",
    "ACVRL1", "ENG", "RYR1", "CACNA1S",
    # Inborn errors of metabolism (5)
    "GAA", "GLA", "OTC", "BTD", "ATP7B",
    # Other metabolic/misc (8)
    "LDLR", "APOB", "PCSK9", "HFE", "TTR", "HNF1A", "RPE65", "SERPINA1",
    # FTD-relevant (not in ACMG but clinically important)
    "MAPT", "GRN", "LRRK2", "GBA1",
}

CLINVAR_VCF_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz"


def download_clinvar_vcf(cache_dir: str, logger) -> str:
    """Download ClinVar VCF (GRCh37) to cache. Returns path."""
    dest = os.path.join(cache_dir, "clinvar_grch37.vcf.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        logger.info(f"ClinVar VCF cached: {dest} ({os.path.getsize(dest) // 1_000_000}MB)")
        return dest
    logger.info("Downloading ClinVar VCF (GRCh37, ~190MB)...")
    download_file(CLINVAR_VCF_URL, dest, logger)
    return dest


def parse_and_match_clinvar_vcf(clinvar_path: str, genome: GenomeData,
                                 logger) -> List[Dict]:
    """Parse ClinVar VCF, find pathogenic/likely_pathogenic variants in ACMG genes matching user genome."""
    hits = []
    acmg_pathogenic_checked = 0
    seen_keys = set()  # avoid duplicates

    opener = gzip.open if clinvar_path.endswith(".gz") else open
    with opener(clinvar_path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.split("\t", 8)
            if len(parts) < 8:
                continue

            chrom_raw, pos_str, vid, ref, alt, qual, filt, info_str = parts[:8]

            # Quick filter on INFO field before full parse
            if "Pathogenic" not in info_str:
                continue

            # Parse INFO
            info = {}
            for field in info_str.split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    info[k] = v

            # Filter: Pathogenic or Likely_pathogenic
            clnsig = info.get("CLNSIG", "")
            if "Pathogenic" not in clnsig:
                continue

            # Filter: in ACMG genes
            geneinfo = info.get("GENEINFO", "")
            gene_names = [g.split(":")[0] for g in geneinfo.split("|") if g]
            matching_genes = [g for g in gene_names if g in ACMG_SF_GENES]
            if not matching_genes:
                continue

            acmg_pathogenic_checked += 1

            # Normalize chromosome
            norm_chrom = normalize_chrom(chrom_raw)
            if not norm_chrom:
                continue
            pos = int(pos_str)
            gene = matching_genes[0]

            # Match against user genome by position
            for user_alt in alt.split(","):
                pos_keys = genome.pos_index.get((norm_chrom, pos), [])
                for vkey in pos_keys:
                    user_var = genome.variants[vkey]
                    if user_var.ref == ref and user_var.alt == user_alt:
                        dosage = user_var.effect_allele_count
                        if dosage > 0:
                            dedup_key = f"{norm_chrom}:{pos}:{ref}:{user_alt}"
                            if dedup_key not in seen_keys:
                                seen_keys.add(dedup_key)
                                hits.append({
                                    "chrom": norm_chrom,
                                    "pos": pos,
                                    "rsid": vid if vid.startswith("rs") else None,
                                    "ref": ref,
                                    "alt": user_alt,
                                    "gene": gene,
                                    "clnsig": clnsig.replace("_", " "),
                                    "condition": info.get("CLNDN", "not_provided").replace("_", " "),
                                    "review_status": info.get("CLNREVSTAT", "").replace("_", " "),
                                    "dosage": dosage,
                                    "zygosity": "Homozygous" if dosage == 2 else "Heterozygous",
                                    "genotype": user_var.gt,
                                })

            # Also try rsID match as fallback
            if vid.startswith("rs"):
                vkey = genome.rsid_index.get(vid)
                if vkey:
                    user_var = genome.variants[vkey]
                    dosage = user_var.effect_allele_count
                    if dosage > 0:
                        dedup_key = f"rsid:{vid}"
                        if dedup_key not in seen_keys:
                            seen_keys.add(dedup_key)
                            hits.append({
                                "chrom": norm_chrom,
                                "pos": pos,
                                "rsid": vid,
                                "ref": ref,
                                "alt": alt.split(",")[0],
                                "gene": gene,
                                "clnsig": clnsig.replace("_", " "),
                                "condition": info.get("CLNDN", "not_provided").replace("_", " "),
                                "review_status": info.get("CLNREVSTAT", "").replace("_", " "),
                                "dosage": dosage,
                                "zygosity": "Homozygous" if dosage == 2 else "Heterozygous",
                                "genotype": user_var.gt,
                            })

    logger.info(f"ClinVar VCF: {acmg_pathogenic_checked} pathogenic ACMG variants checked, "
                f"{len(hits)} found in user genome")
    return hits


def match_clinvar_curated(genome: GenomeData, logger: logging.Logger) -> List[Dict]:
    """Match user variants against curated ClinVar database."""
    results = []

    for rsid, info in CLINVAR_CURATED.items():
        if rsid not in genome.rsid_index:
            continue

        key = genome.rsid_index[rsid]
        variant = genome.variants[key]
        a1, a2 = variant.genotype_alleles
        risk_allele = info["risk"]

        # Count risk allele dosage
        dosage = 0
        if len(risk_allele) == 1:  # SNP
            if a1.upper() == risk_allele.upper():
                dosage += 1
            if a2.upper() == risk_allele.upper():
                dosage += 1
        else:
            # For indels, check ALT match
            if variant.alt.upper() == risk_allele.upper():
                dosage = variant.effect_allele_count

        zygosity = "Homozygous" if dosage == 2 else "Heterozygous" if dosage == 1 else "Reference"

        results.append({
            "rsid": rsid,
            "gene": info["gene"],
            "significance": info["sig"],
            "condition": info["cond"],
            "risk_allele": risk_allele,
            "dosage": dosage,
            "zygosity": zygosity,
            "genotype": f"{a1}/{a2}",
            "note": info["note"],
            "found": True,
        })

    logger.info(f"  Matched {len(results)} of {len(CLINVAR_CURATED)} curated ClinVar variants")
    return results


def run_layer3_clinvar(genome: GenomeData, cache_dir: str,
                       logger: logging.Logger) -> Dict:
    """Run ClinVar analysis using curated database. Returns results dict."""
    logger.info("=" * 60)
    logger.info("LAYER 3: ClinVar Pathogenic Variants (Curated)")
    logger.info("=" * 60)

    results = match_clinvar_curated(genome, logger)

    # Separate by significance and carrier status
    carriers = [r for r in results if r["dosage"] > 0]
    non_carriers = [r for r in results if r["dosage"] == 0]

    pathogenic = [r for r in carriers if "athogenic" in r["significance"]]
    drug_resp = [r for r in carriers if r["significance"] == "Drug response"]
    risk_factor = [r for r in carriers if r["significance"] in ("Risk factor", "Association")]
    protective = [r for r in carriers if r["significance"] == "Protective"]

    # APOE haplotyping
    apoe_e4 = next((r for r in results if r["rsid"] == "rs429358"), None)
    apoe_e2 = next((r for r in results if r["rsid"] == "rs7412"), None)
    if apoe_e4 and apoe_e2:
        e4_dose = apoe_e4["dosage"]
        e2_dose = apoe_e2["dosage"]
        if e4_dose == 0 and e2_dose == 0:
            apoe_type = "e3/e3 (most common, average risk)"
        elif e4_dose == 1 and e2_dose == 0:
            apoe_type = "e3/e4 (~3x Alzheimer risk vs e3/e3)"
        elif e4_dose == 2 and e2_dose == 0:
            apoe_type = "e4/e4 (~12x Alzheimer risk vs e3/e3)"
        elif e4_dose == 0 and e2_dose == 1:
            apoe_type = "e2/e3 (reduced Alzheimer risk, slightly higher lipid issues)"
        elif e4_dose == 0 and e2_dose == 2:
            apoe_type = "e2/e2 (lowest Alzheimer risk, type III hyperlipoproteinemia risk)"
        elif e4_dose == 1 and e2_dose == 1:
            apoe_type = "e2/e4 (mixed — e4 risk partially offset by e2)"
        else:
            apoe_type = f"Unusual combination (e4 dosage={e4_dose}, e2 dosage={e2_dose})"
        apoe_dict = {
            "genotype": apoe_type,
            "rs429358": {"genotype": apoe_e4["genotype"], "e4_dosage": e4_dose},
            "rs7412": {"genotype": apoe_e2["genotype"], "e2_dosage": e2_dose},
        }
    else:
        missing = []
        if not apoe_e4:
            missing.append("rs429358")
        if not apoe_e2:
            missing.append("rs7412")
        apoe_dict = {
            "genotype": "undetermined",
            "note": f"APOE haplotype cannot be determined — {', '.join(missing)} not found in VCF. "
                    "Consider clinical APOE genotyping (blood test) for definitive results.",
        }

    # ClinVar VCF cross-reference (full database)
    clinvar_vcf_hits = []
    try:
        clinvar_path = download_clinvar_vcf(cache_dir, logger)
        clinvar_vcf_hits = parse_and_match_clinvar_vcf(clinvar_path, genome, logger)
    except Exception as e:
        logger.warning(f"ClinVar VCF cross-reference failed: {e}")

    # FTD-specific analysis
    ftd_analysis = {
        "mapt_variants": [r for r in carriers if r.get("gene") == "MAPT"],
        "grn_variants": [h for h in clinvar_vcf_hits if h.get("gene") == "GRN"],
        "lrrk2_variants": [r for r in carriers if r.get("gene") == "LRRK2"],
        "gba1_variants": [r for r in carriers if r.get("gene") == "GBA1"],
        "c9orf72_note": "C9orf72 hexanucleotide repeat expansion CANNOT be detected from VCF files. "
                        "It is the most common genetic cause of FTD (~25-40% of familial cases) and ALS. "
                        "Detection requires repeat-primed PCR or long-read sequencing.",
    }

    logger.info(f"Layer 3 complete: {len(carriers)} carrier variants, {len(pathogenic)} pathogenic")
    return {
        "layer": 3,
        "curated_results": carriers,
        "non_carriers": non_carriers,
        "apoe": apoe_dict,
        "pathogenic_count": len(pathogenic),
        "drug_response_count": len(drug_resp),
        "risk_factor_count": len(risk_factor),
        "protective_count": len(protective),
        "total_screened": len(CLINVAR_CURATED),
        "total_found": len(results),
        "clinvar_vcf_hits": clinvar_vcf_hits,
        "ftd_analysis": ftd_analysis,
    }
