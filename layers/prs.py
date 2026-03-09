import gzip
import logging
import os
from typing import Dict, List, Tuple

from core import GenomeData

# ============================================================================
# SECTION 4: LAYER 1 — POLYGENIC RISK SCORES (PRS)
# ============================================================================

# Curated SNP-level PRS data from published GWAS meta-analyses
# Each condition has top SNPs with validated effect alleles and log(OR) weights
# Sources: CARDIoGRAMplusC4D, DIAGRAM, GIANT, AFGen, IGAP, PGC, GLGC, ICBP
PRS_DATABASE = {
    "Coronary Artery Disease": {
        "source": "CARDIoGRAMplusC4D 2015 + Khera 2018",
        "snps": {
            "rs10455872": {"effect": "G", "weight": 0.49},   # LPA
            "rs4977574": {"effect": "G", "weight": 0.29},    # 9p21 CDKN2A/B
            "rs12526453": {"effect": "C", "weight": 0.12},   # PHACTR1
            "rs6725887": {"effect": "C", "weight": 0.17},    # WDR12
            "rs2505083": {"effect": "C", "weight": 0.10},    # KIAA1462
            "rs9349379": {"effect": "G", "weight": 0.13},    # PHACTR1
            "rs12190287": {"effect": "C", "weight": 0.10},   # TCF21
            "rs1333049": {"effect": "C", "weight": 0.28},    # 9p21
            "rs11206510": {"effect": "T", "weight": 0.15},   # PCSK9
            "rs515135": {"effect": "G", "weight": 0.12},     # APOB
            "rs3184504": {"effect": "T", "weight": 0.13},    # SH2B3
            "rs2895811": {"effect": "C", "weight": 0.07},    # HHIPL1
            "rs12936587": {"effect": "G", "weight": 0.07},   # RASD1/SMCR3
            "rs7412": {"effect": "T", "weight": -0.19},      # APOE (protective)
            "rs429358": {"effect": "C", "weight": 0.10},     # APOE e4
            "rs1746048": {"effect": "C", "weight": 0.09},    # CXCL12
            "rs17114036": {"effect": "A", "weight": -0.14},  # PPAP2B (protective)
            "rs264": {"effect": "G", "weight": 0.06},        # LPL
            "rs3798220": {"effect": "C", "weight": 0.51},    # LPA
        },
    },
    "Type 2 Diabetes": {
        "source": "DIAGRAM 2018 + Mahajan et al",
        "snps": {
            "rs7903146": {"effect": "T", "weight": 0.35},    # TCF7L2 (strongest)
            "rs1801282": {"effect": "C", "weight": -0.14},   # PPARG
            "rs5219": {"effect": "T", "weight": 0.14},       # KCNJ11
            "rs13266634": {"effect": "C", "weight": -0.12},  # SLC30A8
            "rs4402960": {"effect": "T", "weight": 0.14},    # IGF2BP2
            "rs10811661": {"effect": "T", "weight": 0.20},   # CDKN2A/B
            "rs7756992": {"effect": "G", "weight": 0.17},    # CDKAL1
            "rs1111875": {"effect": "C", "weight": 0.13},    # HHEX
            "rs10946398": {"effect": "C", "weight": 0.14},   # CDKAL1
            "rs8050136": {"effect": "A", "weight": 0.15},    # FTO
            "rs2237892": {"effect": "C", "weight": 0.26},    # KCNQ1
            "rs231362": {"effect": "G", "weight": 0.08},     # KCNQ1
            "rs7961581": {"effect": "C", "weight": 0.09},    # TSPAN8/LGR5
            "rs12779790": {"effect": "G", "weight": 0.12},   # CDC123/CAMK1D
            "rs7578597": {"effect": "T", "weight": -0.10},   # THADA
            "rs4607103": {"effect": "C", "weight": 0.08},    # ADAMTS9
            "rs10923931": {"effect": "T", "weight": 0.14},   # NOTCH2
        },
    },
    "Atrial Fibrillation": {
        "source": "AFGen 2018 + Nielsen et al",
        "snps": {
            "rs2200733": {"effect": "T", "weight": 0.51},    # PITX2 (strongest)
            "rs10033464": {"effect": "T", "weight": 0.30},   # PITX2
            "rs6843082": {"effect": "G", "weight": 0.14},    # ZFHX3
            "rs2106261": {"effect": "T", "weight": 0.19},    # ZFHX3
            "rs13376333": {"effect": "T", "weight": 0.17},   # KCNN3
            "rs3807989": {"effect": "A", "weight": 0.12},    # CAV1
            "rs7164883": {"effect": "G", "weight": 0.10},    # KCNN2
            "rs2040862": {"effect": "T", "weight": 0.08},    # NEURL
            "rs1152591": {"effect": "A", "weight": 0.07},    # SYNE2
            "rs10821415": {"effect": "C", "weight": 0.09},   # C9orf3
            "rs6838973": {"effect": "T", "weight": 0.08},    # PRRX1
            "rs7508": {"effect": "G", "weight": 0.06},       # HCN4
        },
    },
    "Alzheimer Disease": {
        "source": "IGAP 2019 + Kunkle et al",
        "snps": {
            "rs429358": {"effect": "C", "weight": 1.15},     # APOE e4 (very strong)
            "rs7412": {"effect": "T", "weight": -0.47},      # APOE e2 (protective)
            "rs6656401": {"effect": "A", "weight": 0.17},    # CR1
            "rs6733839": {"effect": "T", "weight": 0.18},    # BIN1
            "rs35349669": {"effect": "T", "weight": 0.10},   # INPP5D
            "rs190982": {"effect": "G", "weight": -0.08},    # MEF2C
            "rs2718058": {"effect": "G", "weight": -0.07},   # NME8
            "rs1476679": {"effect": "C", "weight": -0.09},   # ZCWPW1
            "rs10948363": {"effect": "G", "weight": 0.11},   # CD2AP
            "rs11771145": {"effect": "A", "weight": -0.10},  # EPHA1
            "rs9271192": {"effect": "C", "weight": 0.11},    # HLA-DRB5
            "rs28834970": {"effect": "C", "weight": 0.11},   # PTK2B
            "rs11218343": {"effect": "C", "weight": -0.17},  # SORL1
            "rs10498633": {"effect": "T", "weight": -0.10},  # SLC24A4
            "rs8093731": {"effect": "T", "weight": -0.36},   # DSG2
        },
    },
    "Body Mass Index": {
        "source": "GIANT 2018 + Yengo et al",
        "snps": {
            "rs1558902": {"effect": "A", "weight": 0.37},    # FTO (strongest)
            "rs6567160": {"effect": "C", "weight": 0.30},    # MC4R
            "rs13021737": {"effect": "G", "weight": 0.22},   # TMEM18
            "rs10938397": {"effect": "G", "weight": 0.18},   # GNPDA2
            "rs543874": {"effect": "G", "weight": 0.19},     # SEC16B
            "rs2867125": {"effect": "C", "weight": 0.16},    # TMEM18
            "rs571312": {"effect": "A", "weight": 0.23},     # MC4R
            "rs10182181": {"effect": "G", "weight": 0.15},   # ADCY3
            "rs12446632": {"effect": "G", "weight": 0.16},   # GPRC5B
            "rs2815752": {"effect": "A", "weight": 0.13},    # NEGR1
            "rs7359397": {"effect": "T", "weight": 0.13},    # SH2B1
            "rs987237": {"effect": "G", "weight": 0.14},     # TFAP2B
            "rs9816226": {"effect": "T", "weight": 0.10},    # ETV5
            "rs3101336": {"effect": "C", "weight": 0.12},    # NEGR1
            "rs7138803": {"effect": "A", "weight": 0.10},    # BCDIN3D
        },
    },
    "LDL Cholesterol": {
        "source": "GLGC 2013 + Willer et al",
        "snps": {
            "rs6511720": {"effect": "T", "weight": -0.22},   # LDLR
            "rs515135": {"effect": "G", "weight": 0.16},     # APOB
            "rs2228671": {"effect": "T", "weight": -0.15},   # LDLR
            "rs629301": {"effect": "T", "weight": -0.18},    # SORT1/CELSR2
            "rs11206510": {"effect": "T", "weight": 0.12},   # PCSK9
            "rs4420638": {"effect": "G", "weight": 0.19},    # APOC1/APOE
            "rs429358": {"effect": "C", "weight": 0.11},     # APOE
            "rs7412": {"effect": "T", "weight": -0.37},      # APOE e2
            "rs1367117": {"effect": "A", "weight": 0.10},    # APOB
            "rs6544713": {"effect": "T", "weight": 0.07},    # ABCG8
            "rs3846663": {"effect": "T", "weight": -0.06},   # HMGCR
            "rs2650000": {"effect": "A", "weight": 0.05},    # HNF1A
        },
    },
    "Systolic Blood Pressure": {
        "source": "ICBP 2018 + Evangelou et al",
        "snps": {
            "rs13082711": {"effect": "T", "weight": -0.50},  # SLC4A7
            "rs3184504": {"effect": "T", "weight": 0.62},    # SH2B3
            "rs1458038": {"effect": "T", "weight": 0.56},    # FGF5
            "rs17249754": {"effect": "G", "weight": -0.90},  # ATP2B1
            "rs11191548": {"effect": "T", "weight": 0.80},   # CYP17A1
            "rs12946454": {"effect": "T", "weight": 0.45},   # PLCD3
            "rs2681472": {"effect": "A", "weight": 0.52},    # ATP2B1
            "rs381815": {"effect": "T", "weight": 0.40},     # PLEKHA7
            "rs1530440": {"effect": "T", "weight": 0.42},    # c10orf107
            "rs4373814": {"effect": "G", "weight": -0.30},   # CACNB2
            "rs13107325": {"effect": "T", "weight": 0.60},   # SLC39A8
            "rs1799945": {"effect": "G", "weight": 0.38},    # HFE
        },
    },
    "Schizophrenia": {
        "source": "PGC3 2022 + Trubetskoy et al",
        "snps": {
            "rs2514218": {"effect": "T", "weight": 0.12},    # DRD2
            "rs1233578": {"effect": "G", "weight": 0.09},    # GRIA1
            "rs4391122": {"effect": "G", "weight": 0.11},    # FURIN/FES
            "rs2007044": {"effect": "G", "weight": 0.08},    # CACNA1C
            "rs115329265": {"effect": "A", "weight": 0.13},  # MHC region
            "rs140505938": {"effect": "C", "weight": 0.09},  # GRIN2A
            "rs2851447": {"effect": "G", "weight": 0.10},    # MIR137
            "rs11191580": {"effect": "T", "weight": 0.08},   # AS3MT
            "rs12704290": {"effect": "G", "weight": 0.09},   # GRM3
            "rs7432375": {"effect": "G", "weight": 0.07},    # SNAP91
        },
    },
    "Major Depressive Disorder": {
        "source": "PGC MDD 2019 + Howard et al",
        "snps": {
            "rs1432639": {"effect": "G", "weight": 0.03},    # VRK2
            "rs10514299": {"effect": "T", "weight": 0.05},   # TMEM161B
            "rs1354115": {"effect": "A", "weight": 0.03},    # SORCS3
            "rs2568958": {"effect": "A", "weight": 0.03},    # NEGR1
            "rs301806": {"effect": "A", "weight": 0.04},     # RERE
            "rs1475120": {"effect": "T", "weight": 0.03},    # VRK2
            "rs10786831": {"effect": "G", "weight": 0.03},   # LHFPL2
            "rs2422321": {"effect": "A", "weight": 0.03},    # OLFM4
            "rs61902811": {"effect": "A", "weight": 0.03},   # ZKSCAN4/NKAPL
        },
    },
    "Prostate Cancer": {
        "source": "PRACTICAL 2018 + Schumacher et al",
        "snps": {
            "rs10993994": {"effect": "T", "weight": 0.20},   # MSMB
            "rs1447295": {"effect": "A", "weight": 0.25},    # 8q24
            "rs6983267": {"effect": "G", "weight": 0.18},    # 8q24
            "rs16901979": {"effect": "A", "weight": 0.44},   # 8q24
            "rs4242382": {"effect": "A", "weight": 0.22},    # 8q24
            "rs620861": {"effect": "G", "weight": -0.11},    # LMTK2
            "rs1859962": {"effect": "G", "weight": 0.18},    # 17q24.3
            "rs17632542": {"effect": "T", "weight": -0.20},  # KLK3
            "rs2735839": {"effect": "G", "weight": -0.19},   # KLK3
            "rs7931342": {"effect": "G", "weight": 0.15},    # 11q13
            "rs10896449": {"effect": "G", "weight": 0.12},   # 11q13
        },
    },
}

# Population mean/SD approximations for European reference (from published studies)
PRS_POPULATION_REF = {
    "Coronary Artery Disease": {"note": "Higher = more risk. Typical EUR range: -0.5 to 2.5"},
    "Type 2 Diabetes": {"note": "Higher = more risk. TCF7L2 rs7903146 TT confers ~2x risk"},
    "Atrial Fibrillation": {"note": "Higher = more risk. PITX2 variants are strongest drivers"},
    "Alzheimer Disease": {"note": "APOE e4 dominates. One e4 allele ~3x risk, two e4 ~12x risk"},
    "Body Mass Index": {"note": "Higher = tendency toward higher BMI. FTO is strongest single locus"},
    "LDL Cholesterol": {"note": "Higher = genetically higher LDL. APOE e2 strongly lowers LDL"},
    "Systolic Blood Pressure": {"note": "Score in mmHg effect. Higher = tendency to higher BP"},
    "Schizophrenia": {"note": "Higher = more genetic loading. Effect sizes are very small per SNP"},
    "Major Depressive Disorder": {"note": "Effect sizes very small per SNP. Highly polygenic"},
    "Prostate Cancer": {"note": "Higher = more risk. 8q24 region contributes multiple signals"},
}

# Map PGS IDs to health domains for JSON report grouping
PGS_CATEGORIES = {
    # Cardio
    "PGS000018": "cardio", "PGS000013": "cardio", "PGS000296": "cardio",
    "PGS000016": "cardio", "PGS000039": "cardio", "PGS002997": "cardio",
    "PGS000034": "cardio",
    # Metabolic
    "PGS000014": "metabolism", "PGS000807": "metabolism", "PGS000033": "metabolism",
    "PGS000027": "metabolism", "PGS003037": "metabolism", "PGS000115": "metabolism",
    "PGS003147": "metabolism", "PGS002937": "metabolism", "PGS005276": "metabolism",
    "PGS002356": "metabolism", "PGS002813": "metabolism", "PGS000882": "metabolism",
    "PGS003154": "metabolism", "PGS003426": "metabolism", "PGS004332": "metabolism",
    "PGS002768": "metabolism", "PGS005222": "metabolism",
    # Neuro
    "PGS003992": "neuro", "PGS000812": "neuro", "PGS000903": "neuro",
    "PGS002760": "neuro", "PGS000137": "neuro", "PGS001797": "neuro",
    # Mental health
    "PGS003497": "mental", "PGS002786": "mental", "PGS004451": "mental",
    "PGS002746": "mental", "PGS002790": "mental", "PGS000135": "mental",
    "PGS005211": "mental", "PGS002231": "mental",
    # Oncology
    "PGS004834": "oncology",
    # Biometric / other
    "PGS001229": "metabolism", "PGS000037": "immune",
    # Immune
    "PGS000017": "immune", "PGS002344": "immune", "PGS000728": "immune",
}

# Map curated PRS condition names to categories
CURATED_CATEGORIES = {
    "Coronary Artery Disease": "cardio",
    "Type 2 Diabetes": "metabolism",
    "Atrial Fibrillation": "cardio",
    "Alzheimer Disease": "neuro",
    "Body Mass Index (BMI)": "metabolism",
    "LDL Cholesterol": "metabolism",
    "Systolic Blood Pressure": "cardio",
    "Schizophrenia": "mental",
    "Major Depressive Disorder": "mental",
    "Prostate Cancer": "oncology",
}


def score_prs_curated(genome: GenomeData, snps: Dict,
                      logger: logging.Logger) -> Tuple[float, int, int, List[Dict]]:
    """
    Calculate PRS from curated SNP database.
    Returns (score, matched, total, details_list).
    """
    total_score = 0.0
    matched = 0
    details = []

    for rsid, info in snps.items():
        if rsid not in genome.rsid_index:
            continue

        key = genome.rsid_index[rsid]
        variant = genome.variants[key]
        effect_allele = info["effect"]
        weight = info["weight"]

        a1, a2 = variant.genotype_alleles
        dosage = 0
        if a1.upper() == effect_allele.upper():
            dosage += 1
        if a2.upper() == effect_allele.upper():
            dosage += 1

        contribution = dosage * weight
        total_score += contribution
        matched += 1

        if dosage > 0:
            details.append({
                "rsid": rsid,
                "effect_allele": effect_allele,
                "weight": weight,
                "dosage": dosage,
                "contribution": contribution,
                "genotype": f"{a1}/{a2}",
            })

    details.sort(key=lambda d: abs(d["contribution"]), reverse=True)
    return (total_score, matched, len(snps), details)


# ---------------------------------------------------------------------------
# PGS CATALOG SCORING FILE PARSER
# ---------------------------------------------------------------------------
# Reads harmonized scoring files from PGS Catalog (format: gzipped TSV)
# Header lines start with '#', data columns include:
#   rsID, chr_name, chr_position, effect_allele, other_allele, effect_weight
# Harmonized files (GRCh37) have columns:
#   hm_chr, hm_pos, hm_inferOtherAllele, hm_match_chr, hm_match_pos
# Download URL pattern:
#   https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS00XXXX/
#     ScoringFiles/Harmonized/PGS00XXXX_hmPOS_GRCh37.txt.gz

def parse_pgs_scoring_file(filepath: str,
                           logger: logging.Logger) -> Dict:
    """
    Parse a PGS Catalog scoring file (.txt.gz or .txt).
    Returns dict with metadata and list of variant weights.
    """
    metadata = {}
    variants = []
    header_cols = None

    open_fn = gzip.open if filepath.endswith(".gz") else open

    with open_fn(filepath, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue

            # Header metadata lines
            if line.startswith("#"):
                if "=" in line:
                    key, _, val = line[1:].partition("=")
                    metadata[key.strip().lower()] = val.strip()
                continue

            # Column header line (first non-# non-empty line)
            if header_cols is None:
                header_cols = line.split("\t")
                # Normalize column names to lowercase
                header_cols = [c.strip().lower() for c in header_cols]
                continue

            # Data rows
            fields = line.split("\t")
            if len(fields) != len(header_cols):
                continue

            row = dict(zip(header_cols, fields))
            variant_entry = {}

            # rsID
            rsid = row.get("rsid", row.get("snpid", "")).strip()
            if rsid and rsid.startswith("rs"):
                variant_entry["rsid"] = rsid

            # Chromosome and position — prefer harmonized columns
            chrom = row.get("hm_chr", row.get("chr_name", "")).strip()
            pos_str = row.get("hm_pos", row.get("chr_position", "")).strip()
            if chrom:
                variant_entry["chr"] = chrom.replace("chr", "")
            if pos_str:
                try:
                    variant_entry["pos"] = int(pos_str)
                except ValueError:
                    pass

            # Alleles
            ea = row.get("effect_allele", "").strip().upper()
            oa = row.get("other_allele",
                         row.get("reference_allele",
                                 row.get("hm_inferotherallele", ""))).strip().upper()
            if ea:
                variant_entry["effect_allele"] = ea
            if oa:
                variant_entry["other_allele"] = oa

            # Weight
            weight_str = row.get("effect_weight", "").strip()
            if not weight_str:
                continue
            try:
                variant_entry["weight"] = float(weight_str)
            except ValueError:
                continue

            # Must have either rsID or chr+pos
            if "rsid" not in variant_entry and \
               ("chr" not in variant_entry or "pos" not in variant_entry):
                continue

            variants.append(variant_entry)

    # Extract trait name from metadata
    trait = metadata.get("trait_mapped",
                         metadata.get("pgs_name",
                                      metadata.get("reported_trait",
                                                    os.path.basename(filepath))))
    pgs_id = metadata.get("pgs_id", "unknown")
    n_variants = metadata.get("variants_number",
                              metadata.get("number_of_variants", str(len(variants))))

    logger.info(f"  PGS file parsed: {pgs_id} | {trait} | "
                f"{len(variants)} weights loaded (header says {n_variants})")

    return {
        "pgs_id": pgs_id,
        "trait": trait,
        "metadata": metadata,
        "variants": variants,
    }


def score_pgs_catalog(genome: GenomeData, pgs_data: Dict,
                      logger: logging.Logger) -> Dict:
    """
    Score a genome against a PGS Catalog scoring file.
    Matches by rsID first, falls back to chr:pos.
    Returns dict with score, matched count, details.
    """
    total_score = 0.0
    matched = 0
    total = len(pgs_data["variants"])
    top_contributors = []  # keep track of largest absolute contributions

    for v in pgs_data["variants"]:
        variant = None
        rsid = v.get("rsid")

        # Try rsID match first
        if rsid and rsid in genome.rsid_index:
            key = genome.rsid_index[rsid]
            variant = genome.variants[key]
        # Fall back to chr:pos match
        elif "chr" in v and "pos" in v:
            chrom = v["chr"]
            pos = v["pos"]
            keys = genome.pos_index.get((chrom, pos), [])
            if keys:
                variant = genome.variants[keys[0]]

        if variant is None:
            continue

        effect_allele = v["effect_allele"]
        weight = v["weight"]

        a1, a2 = variant.genotype_alleles
        dosage = 0
        if a1.upper() == effect_allele:
            dosage += 1
        if a2.upper() == effect_allele:
            dosage += 1

        contribution = dosage * weight
        total_score += contribution
        matched += 1

        # Track top contributors (keep top 50 by absolute contribution)
        if dosage > 0:
            entry = {
                "rsid": rsid or f"{variant.chrom}:{variant.pos}",
                "effect_allele": effect_allele,
                "weight": weight,
                "dosage": dosage,
                "contribution": contribution,
                "genotype": f"{a1}/{a2}",
            }
            top_contributors.append(entry)

    # Sort and keep top 20 contributors for display
    top_contributors.sort(key=lambda d: abs(d["contribution"]), reverse=True)
    top_contributors = top_contributors[:20]

    coverage_pct = matched / total * 100 if total else 0

    logger.info(f"    Score: {total_score:.4f} | Matched: {matched:,}/{total:,} "
                f"({coverage_pct:.1f}%)")

    return {
        "condition": pgs_data["trait"],
        "source": f"PGS Catalog {pgs_data['pgs_id']}",
        "pgs_id": pgs_data["pgs_id"],
        "score": total_score,
        "matched": matched,
        "total": total,
        "coverage": coverage_pct,
        "details": top_contributors,
        "note": f"Full PGS with {total:,} variants. "
                f"Matched {matched:,} in your VCF ({coverage_pct:.1f}%).",
    }


def discover_pgs_files(pgs_dir: str, logger: logging.Logger) -> List[str]:
    """Find all PGS scoring files in directory (*.txt.gz or *.txt)."""
    if not pgs_dir or not os.path.isdir(pgs_dir):
        return []
    files = []
    for f in sorted(os.listdir(pgs_dir)):
        if f.endswith(".txt.gz") or (f.endswith(".txt") and "PGS" in f.upper()):
            files.append(os.path.join(pgs_dir, f))
    logger.info(f"  Found {len(files)} PGS scoring file(s) in {pgs_dir}")
    return files


def run_layer1_prs(genome: GenomeData, cache_dir: str,
                   logger: logging.Logger, pgs_dir: str = None) -> Dict:
    """
    Run PRS analysis.
    If pgs_dir contains PGS Catalog scoring files, use those (full PRS).
    Otherwise fall back to curated top-SNP database.
    """
    logger.info("=" * 60)
    logger.info("LAYER 1: Polygenic Risk Scores (PRS)")
    logger.info("=" * 60)

    results = []
    pgs_files = discover_pgs_files(pgs_dir, logger) if pgs_dir else []
    used_pgs_catalog = len(pgs_files) > 0

    if used_pgs_catalog:
        # ---- PGS Catalog mode: full scoring files ----
        logger.info(f"  Using PGS Catalog scoring files from: {pgs_dir}")
        for fpath in pgs_files:
            logger.info(f"  Processing: {os.path.basename(fpath)}")
            try:
                pgs_data = parse_pgs_scoring_file(fpath, logger)
                result = score_pgs_catalog(genome, pgs_data, logger)
                result["category"] = PGS_CATEGORIES.get(result.get("pgs_id", ""), "other")
                results.append(result)
            except Exception as e:
                logger.error(f"  Failed to process {fpath}: {e}")
        method_label = "PGS Catalog harmonized scoring files (GRCh37)"
    else:
        # ---- Fallback: curated top-SNP PRS ----
        logger.info("  No PGS Catalog files found — using curated top-SNP database")
        for condition, cdata in PRS_DATABASE.items():
            snps = cdata["snps"]
            source = cdata["source"]
            logger.info(f"  {condition} ({len(snps)} SNPs)")

            score, matched, total, details = score_prs_curated(genome, snps, logger)
            coverage = matched / total * 100 if total else 0
            logger.info(f"    Score: {score:.4f} | Matched: {matched}/{total} ({coverage:.0f}%)")

            results.append({
                "condition": condition, "source": source,
                "score": score, "matched": matched, "total": total,
                "coverage": coverage, "details": details,
                "note": PRS_POPULATION_REF.get(condition, {}).get("note", ""),
                "category": CURATED_CATEGORIES.get(condition, "other"),
            })
        method_label = "Curated top-SNP PRS from published GWAS meta-analyses"

    logger.info(f"Layer 1 complete: {len(results)} conditions scored")
    return {"layer": 1, "results": results}
