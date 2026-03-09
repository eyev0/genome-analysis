import logging
from typing import Dict, List

from core import GenomeData

# ============================================================================
# SECTION 5: LAYER 2 — PHARMACOGENETICS (CPIC)
# ============================================================================

# Curated pharmacogene database: gene → defining SNPs → star alleles → phenotypes → drugs
# Based on CPIC guidelines (cpicpgx.org) and PharmGKB
PHARMACO_DB = {
    "CYP2D6": {
        "defining_snps": {
            "rs3892097": {"*4": "A"},   # splicing defect, loss of function
            "rs1065852": {"*10": "A"},  # reduced function (P34S)
            "rs5030655": {"*6": "del"}, # frameshift, null
            "rs16947":   {"*2": "A"},   # normal function
            "rs1135840":  {"*2": "G"},  # normal function
            "rs35742686": {"*3": "del"},  # frameshift, no function. ~1-2% EUR
            "rs5030656":  {"*9": "del"},  # deletion, decreased function
            "rs28371706": {"*17": "T"},   # decreased function. ~20% AFR
            "rs28371725": {"*41": "T"},   # splicing defect, decreased. ~8% EUR
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Normal Metabolizer",
            "*1/*4": "Intermediate Metabolizer",
            "*1/*10": "Intermediate Metabolizer",
            "*4/*4": "Poor Metabolizer",
            "*4/*10": "Poor Metabolizer",
            "*10/*10": "Intermediate Metabolizer",
            "*1/*41": "Intermediate Metabolizer",
            "*1/*17": "Intermediate Metabolizer",
            "*4/*41": "Poor Metabolizer",
            "*41/*41": "Intermediate Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*3/*4": "Poor Metabolizer",
            "*1/*9": "Intermediate Metabolizer",
        },
        "drugs": [
            {"drug": "Codeine", "action_pm": "Avoid; use alternative analgesic", "action_im": "Consider reduced dose or alternative", "action_nm": "Standard dosing"},
            {"drug": "Tramadol", "action_pm": "Avoid; use alternative", "action_im": "Consider alternative", "action_nm": "Standard dosing"},
            {"drug": "Tamoxifen", "action_pm": "Consider aromatase inhibitor", "action_im": "Consider higher dose or alternative", "action_nm": "Standard dosing"},
            {"drug": "Ondansetron", "action_pm": "Standard dosing", "action_im": "Standard dosing", "action_nm": "Standard dosing"},
            {"drug": "Aripiprazole", "action_pm": "Reduce dose by 50%", "action_im": "Reduce dose by 25%", "action_nm": "Standard dosing"},
        ],
    },
    "CYP2C19": {
        "defining_snps": {
            "rs4244285": {"*2": "A"},   # splicing defect, null
            "rs4986893": {"*3": "A"},   # stop codon, null
            "rs12248560": {"*17": "T"}, # ultra-rapid
            "rs28399504": {"*4": "G"},  # no function. Rare
            "rs56337013": {"*5": "T"},  # no function. Rare
            "rs72552267": {"*6": "A"},  # no function. Rare
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*1/*17": "Rapid Metabolizer",
            "*17/*17": "Ultra-rapid Metabolizer",
            "*2/*2": "Poor Metabolizer",
            "*2/*3": "Poor Metabolizer",
            "*2/*17": "Intermediate Metabolizer",
            "*1/*4": "Intermediate Metabolizer",
            "*1/*5": "Intermediate Metabolizer",
            "*1/*6": "Intermediate Metabolizer",
            "*4/*4": "Poor Metabolizer",
        },
        "drugs": [
            {"drug": "Clopidogrel", "action_pm": "Use alternative (prasugrel/ticagrelor)", "action_im": "Use alternative", "action_nm": "Standard dosing"},
            {"drug": "Omeprazole", "action_pm": "Standard dosing (increased efficacy)", "action_im": "Standard dosing", "action_nm": "Standard dosing"},
            {"drug": "Escitalopram", "action_pm": "Reduce dose by 50%", "action_im": "Reduce dose by 25%", "action_nm": "Standard dosing"},
            {"drug": "Sertraline", "action_pm": "Consider alternative", "action_im": "Standard dosing", "action_nm": "Standard dosing"},
            {"drug": "Voriconazole", "action_pm": "Use alternative antifungal", "action_im": "Standard dosing", "action_nm": "Standard dosing"},
        ],
    },
    "CYP2C9": {
        "defining_snps": {
            "rs1799853": {"*2": "T"},   # R144C, reduced function
            "rs1057910": {"*3": "C"},   # I359L, reduced function
            "rs56165452": {"*5": "T"},   # no function. ~1% AFR
            "rs28371686": {"*6": "A"},   # no function. Rare
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*2/*2": "Poor Metabolizer",
            "*2/*3": "Poor Metabolizer",
            "*3/*3": "Poor Metabolizer",
            "*1/*5": "Intermediate Metabolizer",
            "*1/*6": "Intermediate Metabolizer",
            "*5/*5": "Poor Metabolizer",
        },
        "drugs": [
            {"drug": "Warfarin", "action_pm": "Reduce dose significantly (≤2mg/day)", "action_im": "Reduce initial dose by 25-50%", "action_nm": "Standard dosing algorithm"},
            {"drug": "Phenytoin", "action_pm": "Reduce dose by 50%", "action_im": "Reduce dose by 25%", "action_nm": "Standard dosing"},
            {"drug": "Celecoxib", "action_pm": "Reduce dose by 50%", "action_im": "Consider lower dose", "action_nm": "Standard dosing"},
            {"drug": "Flurbiprofen", "action_pm": "Reduce dose by 50%", "action_im": "Reduce dose by 25%", "action_nm": "Standard dosing"},
        ],
    },
    "CYP3A5": {
        "defining_snps": {
            "rs776746": {"*3": "G"},  # splicing defect (most common non-expressor)
        },
        "phenotypes": {
            "*1/*1": "Expressor (extensive)",
            "*1/*3": "Intermediate Expressor",
            "*3/*3": "Non-expressor",
        },
        "drugs": [
            {"drug": "Tacrolimus", "action_pm": "Standard dosing", "action_im": "Increase dose by 1.5-2x", "action_nm": "Increase dose by 1.5-2x"},
        ],
    },
    "DPYD": {
        "defining_snps": {
            "rs3918290": {"*2A": "G"},  # splice, null (IVS14+1G>A) — ALT=G in VCF
            "rs55886062": {"*13": "C"}, # I560S, null — ALT=C in VCF
            "rs67376798": {"HapB3": "A"}, # D949V, reduced — ALT=A in VCF
            "rs75017182": {"c.1129-5923C>G": "C"}, # reduced — ALT=C in VCF
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2A": "Intermediate Metabolizer (50% DPD activity)",
            "*2A/*2A": "Poor Metabolizer (DPD deficient)",
        },
        "drugs": [
            {"drug": "5-Fluorouracil", "action_pm": "CONTRAINDICATED — fatal toxicity risk", "action_im": "Reduce dose by 50%", "action_nm": "Standard dosing"},
            {"drug": "Capecitabine", "action_pm": "CONTRAINDICATED — fatal toxicity risk", "action_im": "Reduce dose by 50%", "action_nm": "Standard dosing"},
        ],
    },
    "TPMT": {
        "defining_snps": {
            "rs1800462": {"*2": "G"},   # A80P
            "rs1800460": {"*3B": "A"},  # A154T
            "rs1142345": {"*3C": "G"},  # Y240C
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*3C": "Intermediate Metabolizer",
            "*3B/*3C": "Poor Metabolizer",
            "*3C/*3C": "Poor Metabolizer",
        },
        "drugs": [
            {"drug": "Azathioprine", "action_pm": "Reduce dose by 90% or use alternative", "action_im": "Reduce dose by 30-70%", "action_nm": "Standard dosing"},
            {"drug": "6-Mercaptopurine", "action_pm": "Reduce dose by 90% or use alternative", "action_im": "Reduce dose by 30-70%", "action_nm": "Standard dosing"},
            {"drug": "Thioguanine", "action_pm": "Reduce dose by 90% or use alternative", "action_im": "Reduce dose by 30-70%", "action_nm": "Standard dosing"},
        ],
    },
    "NUDT15": {
        "defining_snps": {
            "rs116855232": {"*3": "T"},  # R139C, null
            "rs186364861": {"*2": "A"},  # V18I, reduced
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*3/*3": "Poor Metabolizer",
        },
        "drugs": [
            {"drug": "Azathioprine", "action_pm": "Reduce dose by 90% or use alternative", "action_im": "Reduce dose by 50%", "action_nm": "Standard dosing"},
            {"drug": "6-Mercaptopurine", "action_pm": "Reduce dose by 90% or use alternative", "action_im": "Reduce dose by 50%", "action_nm": "Standard dosing"},
        ],
    },
    "SLCO1B1": {
        "defining_snps": {
            "rs4149056": {"*5": "C"},  # V174A, reduced transport (T>C in VCF)
        },
        "phenotypes": {
            "*1/*1": "Normal Function",
            "*1/*5": "Decreased Function",
            "*5/*5": "Poor Function",
        },
        "drugs": [
            {"drug": "Simvastatin", "action_pm": "Use alternative statin or low dose (≤20mg)", "action_im": "Use ≤20mg or alternative statin", "action_nm": "Standard dosing"},
            {"drug": "Atorvastatin", "action_pm": "Consider lower dose", "action_im": "Monitor for myopathy", "action_nm": "Standard dosing"},
        ],
    },
    "UGT1A1": {
        "defining_snps": {
            "rs8175347": {"*28": "T"},   # TA repeat proxy. Decreased expression
            "rs887829": {"*80": "T"},    # better proxy for *28 on arrays. High LD
        },
        "phenotypes": {
            "*1/*1": "Normal Metabolizer",
            "*1/*28": "Intermediate Metabolizer",
            "*28/*28": "Poor Metabolizer",
            "*1/*80": "Intermediate Metabolizer",
            "*80/*80": "Poor Metabolizer",
        },
        "drugs": [
            {"drug": "Irinotecan", "action_pm": "Reduce dose by 30%", "action_im": "Consider reduced dose", "action_nm": "Standard dosing"},
            {"drug": "Atazanavir", "action_pm": "Monitor for jaundice", "action_im": "Monitor for jaundice", "action_nm": "Standard dosing"},
        ],
    },
    "VKORC1": {
        "defining_snps": {
            "rs9923231": {"Low-dose": "T"},    # -1639G>A → VCF shows C>T on opposite strand
            "rs9934438": {"Low-dose": "A"},    # 1173C>T → VCF shows G>A
            "rs8050894": {"variable": "G"},    # VCF shows C>G
        },
        "phenotypes": {
            # VCF alleles are C/T for rs9923231 (= G/A on gene strand)
            # CC = GG (normal), CT = GA (moderate), TT = AA (high sensitivity)
            "CC": "Normal Warfarin Sensitivity",
            "CT": "Moderate Warfarin Sensitivity (lower dose needed)",
            "TC": "Moderate Warfarin Sensitivity (lower dose needed)",
            "TT": "High Warfarin Sensitivity (significantly lower dose needed)",
        },
        "drugs": [
            {"drug": "Warfarin", "action_pm": "Start at 50-75% of standard dose", "action_im": "Start at 75% of standard dose", "action_nm": "Standard dosing algorithm"},
            {"drug": "Acenocoumarol", "action_pm": "Reduce initial dose", "action_im": "Consider reduced dose", "action_nm": "Standard dosing"},
        ],
    },
    "COMT": {
        "defining_snps": {
            "rs4680": {"Met": "A"},  # Val158Met
        },
        "phenotypes": {
            "Val/Val": "High COMT Activity (rapid dopamine clearance)",
            "Val/Met": "Intermediate COMT Activity",
            "Met/Met": "Low COMT Activity (slow dopamine clearance)",
        },
        "drugs": [
            {"drug": "Catechol-O-methyltransferase substrates", "action_pm": "Enhanced response to dopaminergic drugs", "action_im": "Variable response", "action_nm": "Standard response"},
        ],
    },
    "IFNL3": {
        "defining_snps": {
            "rs12979860": {"CC": "C"},  # favorable response to IFN-α therapy
        },
        "phenotypes": {
            "CC": "Favorable HCV treatment response",
            "CT": "Intermediate HCV treatment response",
            "TT": "Unfavorable HCV treatment response",
        },
        "drugs": [
            {"drug": "Peginterferon alfa-2a/2b", "action_pm": "Standard duration may not be effective", "action_im": "Consider extended course", "action_nm": "Standard treatment course"},
        ],
    },
    "HLA-B": {
        "defining_snps": {
            "rs2395029": {"*5701": "G"},  # proxy for HLA-B*57:01
        },
        "phenotypes": {
            "carrier": "HLA-B*57:01 carrier (abacavir hypersensitivity risk)",
            "non-carrier": "Non-carrier",
        },
        "drugs": [
            {"drug": "Abacavir", "action_pm": "CONTRAINDICATED (hypersensitivity)", "action_im": "CONTRAINDICATED (hypersensitivity)", "action_nm": "Standard dosing"},
        ],
    },
}

def call_simple_diplotype(gene_name: str, gene_data: Dict,
                          genome: GenomeData) -> Dict:
    """
    Simplified diplotype calling based on defining SNPs.
    Returns {allele1, allele2, phenotype, snps_found}.
    """
    defining = gene_data.get("defining_snps", {})
    found_variants = {}

    for rsid, allele_map in defining.items():
        if rsid in genome.rsid_index:
            key = genome.rsid_index[rsid]
            variant = genome.variants[key]
            a1, a2 = variant.genotype_alleles
            found_variants[rsid] = {"a1": a1, "a2": a2, "gt": variant.gt}

    # Special handling for VKORC1 and COMT (not star-allele based)
    if gene_name == "VKORC1":
        rs = "rs9923231"
        if rs in found_variants:
            a1, a2 = found_variants[rs]["a1"], found_variants[rs]["a2"]
            gt_str = f"{a1}{a2}"
            pheno_key = gt_str if gt_str in gene_data["phenotypes"] else None
            # Try reversed
            if not pheno_key:
                pheno_key = f"{a2}{a1}" if f"{a2}{a1}" in gene_data["phenotypes"] else None
            phenotype = gene_data["phenotypes"].get(pheno_key, f"Genotype: {gt_str}")
            return {"allele1": a1, "allele2": a2, "phenotype": phenotype,
                    "snps_found": found_variants, "gene": gene_name}
        return {"allele1": "?", "allele2": "?", "phenotype": "No data",
                "snps_found": {}, "gene": gene_name}

    if gene_name == "COMT":
        rs = "rs4680"
        if rs in found_variants:
            a1, a2 = found_variants[rs]["a1"], found_variants[rs]["a2"]
            # A = Met, G = Val
            label1 = "Met" if a1 == "A" else "Val"
            label2 = "Met" if a2 == "A" else "Val"
            gt_label = f"{label1}/{label2}"
            phenotype = gene_data["phenotypes"].get(gt_label, f"Genotype: {gt_label}")
            return {"allele1": label1, "allele2": label2, "phenotype": phenotype,
                    "snps_found": found_variants, "gene": gene_name}
        return {"allele1": "?", "allele2": "?", "phenotype": "No data",
                "snps_found": {}, "gene": gene_name}

    if gene_name in ("IFNL3", "HLA-B"):
        # Simple genotype-based
        for rsid in defining:
            if rsid in found_variants:
                a1, a2 = found_variants[rsid]["a1"], found_variants[rsid]["a2"]
                gt_str = f"{a1}{a2}"
                phenotypes = gene_data["phenotypes"]
                phenotype = "Unknown"
                for k, v in phenotypes.items():
                    if k in gt_str or gt_str in k:
                        phenotype = v
                        break
                # HLA-B special: check carrier status
                if gene_name == "HLA-B":
                    allele_info = list(defining[rsid].values())[0]
                    if a1 == allele_info or a2 == allele_info:
                        phenotype = phenotypes.get("carrier", "Carrier")
                    else:
                        phenotype = phenotypes.get("non-carrier", "Non-carrier")
                return {"allele1": a1, "allele2": a2, "phenotype": phenotype,
                        "snps_found": found_variants, "gene": gene_name}
        return {"allele1": "?", "allele2": "?", "phenotype": "No data",
                "snps_found": {}, "gene": gene_name}

    # Standard star-allele based genes
    # Determine which non-*1 alleles are present
    detected_alleles_hap1 = ["*1"]
    detected_alleles_hap2 = ["*1"]

    for rsid, allele_map in defining.items():
        if rsid not in found_variants:
            continue
        a1, a2 = found_variants[rsid]["a1"], found_variants[rsid]["a2"]
        for star_name, defining_allele in allele_map.items():
            if a1 == defining_allele:
                detected_alleles_hap1.append(star_name)
            if a2 == defining_allele:
                detected_alleles_hap2.append(star_name)

    # Pick highest-priority non-*1 allele for each haplotype
    # (simplified: just take first non-*1)
    allele1 = next((a for a in detected_alleles_hap1 if a != "*1"), "*1")
    allele2 = next((a for a in detected_alleles_hap2 if a != "*1"), "*1")

    # Sort for consistent lookup
    pair = "/".join(sorted([allele1, allele2]))
    phenotype = gene_data["phenotypes"].get(
        f"{allele1}/{allele2}",
        gene_data["phenotypes"].get(
            f"{allele2}/{allele1}",
            gene_data["phenotypes"].get(pair, f"Diplotype: {allele1}/{allele2}")
        )
    )

    return {"allele1": allele1, "allele2": allele2, "phenotype": phenotype,
            "snps_found": found_variants, "gene": gene_name}

def get_drug_action(drug_entry: Dict, phenotype: str) -> str:
    """Pick the right drug action based on phenotype category."""
    pheno_lower = phenotype.lower()
    if any(k in pheno_lower for k in ["poor", "deficient", "non-express",
                                       "low comt", "high warfarin", "carrier",
                                       "unfavorable", "contraindicated"]):
        return drug_entry.get("action_pm", "See CPIC guidelines")
    elif any(k in pheno_lower for k in ["intermediate", "decreased",
                                         "moderate warfarin", "reduced"]):
        return drug_entry.get("action_im", "See CPIC guidelines")
    else:
        return drug_entry.get("action_nm", "Standard dosing")

def run_layer2_pharmacogenetics(genome: GenomeData,
                                 logger: logging.Logger) -> Dict:
    """Run pharmacogenetic analysis. Returns results dict."""
    logger.info("=" * 60)
    logger.info("LAYER 2: Pharmacogenetics (CPIC)")
    logger.info("=" * 60)

    results = []
    for gene_name, gene_data in PHARMACO_DB.items():
        logger.info(f"  Analyzing {gene_name}...")
        diplotype = call_simple_diplotype(gene_name, gene_data, genome)

        drug_recs = []
        for drug_entry in gene_data.get("drugs", []):
            action = get_drug_action(drug_entry, diplotype["phenotype"])
            drug_recs.append({
                "drug": drug_entry["drug"],
                "action": action,
            })

        results.append({
            "gene": gene_name,
            "allele1": diplotype["allele1"],
            "allele2": diplotype["allele2"],
            "phenotype": diplotype["phenotype"],
            "snps_found": len(diplotype["snps_found"]),
            "snps_total": len(gene_data.get("defining_snps", {})),
            "drugs": drug_recs,
        })
        logger.info(f"    {diplotype['allele1']}/{diplotype['allele2']} → {diplotype['phenotype']}")

    logger.info(f"Layer 2 complete: {len(results)} genes analyzed")
    return {"layer": 2, "results": results}
