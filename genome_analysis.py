#!/usr/bin/env python3
"""
Personal Genome Analysis Pipeline
===================================
Comprehensive analysis of a VCF file against 5 open databases:
  Layer 1: Polygenic Risk Scores (PGS Catalog)
  Layer 2: Pharmacogenetics (CPIC guidelines)
  Layer 3: ClinVar pathogenic variants
  Layer 4: Ancestry inference (AIMs panel)
  Layer 5: GWAS trait associations

Usage:
    python genome_analysis.py --vcf input.vcf --output-dir ./results
    python genome_analysis.py --vcf input.vcf --layer 3   # run only ClinVar
"""

import argparse
import csv
import gzip
import io
import json
import logging
import math
import os
import shutil
import sys
import time
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from urllib.request import urlopen, Request
from urllib.error import URLError

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

# ============================================================================
# SECTION 1: CONFIGURATION & DATA STRUCTURES
# ============================================================================

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(message)s"

@dataclass
class Variant:
    """Single variant from VCF."""
    chrom: str       # normalized: "1"-"22", "X", "Y", "MT"
    pos: int
    ref: str
    alt: str
    rsid: Optional[str] = None
    gt: Optional[str] = None   # e.g. "0|1", "1|1"

    @property
    def key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @property
    def genotype_alleles(self) -> Tuple[str, str]:
        """Return actual allele pair from GT field."""
        if not self.gt:
            return (self.ref, self.ref)
        sep = "|" if "|" in self.gt else "/"
        parts = self.gt.split(sep)
        alleles_list = [self.ref] + self.alt.split(",")
        try:
            a = alleles_list[int(parts[0])]
            b = alleles_list[int(parts[1])]
            return (a, b)
        except (ValueError, IndexError):
            return (self.ref, self.ref)

    @property
    def effect_allele_count(self) -> int:
        """Count of ALT alleles (0, 1, or 2)."""
        if not self.gt:
            return 0
        sep = "|" if "|" in self.gt else "/"
        parts = self.gt.split(sep)
        try:
            return int(parts[0]) + int(parts[1])
        except ValueError:
            return 0

    @property
    def is_snp(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1 and self.alt != "."

@dataclass
class GenomeData:
    """In-memory VCF representation with multiple indices."""
    variants: Dict[str, Variant] = field(default_factory=dict)       # key → Variant
    rsid_index: Dict[str, str] = field(default_factory=dict)         # rsID → key
    pos_index: Dict[Tuple[str, int], List[str]] = field(default_factory=lambda: defaultdict(list))
    total: int = 0
    sample_id: str = "SAMPLE"


# ============================================================================
# SECTION 2: VCF PARSER
# ============================================================================

def normalize_chrom(chrom: str) -> Optional[str]:
    """Normalize chromosome name: chr1 → 1, chrX → X, chrM → MT."""
    c = chrom.replace("chr", "").upper()
    if c in [str(i) for i in range(1, 23)]:
        return c
    if c == "X": return "X"
    if c == "Y": return "Y"
    if c in ("M", "MT"): return "MT"
    return None

def parse_vcf(vcf_path: str, logger: logging.Logger) -> GenomeData:
    """
    Parse VCF file into GenomeData with lookup indices.
    Handles .vcf and .vcf.gz. Reads once, caches in memory.
    """
    genome = GenomeData()
    open_fn = gzip.open if vcf_path.endswith(".gz") else open

    logger.info(f"Parsing VCF: {vcf_path}")
    t0 = time.time()

    with open_fn(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                if len(header) > 9:
                    genome.sample_id = header[9]
                continue
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            chrom_raw, pos_str, rsid_raw, ref, alt = fields[0:5]
            fmt_col, sample_col = fields[8], fields[9]

            chrom = normalize_chrom(chrom_raw)
            if chrom is None:
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            # Extract GT
            fmt_keys = fmt_col.split(":")
            gt = None
            if "GT" in fmt_keys:
                gt_idx = fmt_keys.index("GT")
                sample_parts = sample_col.split(":")
                if gt_idx < len(sample_parts):
                    gt = sample_parts[gt_idx]

            rsid = rsid_raw if rsid_raw.startswith("rs") else None

            # Handle multi-allelic: only take first ALT for now
            first_alt = alt.split(",")[0] if "," in alt else alt

            variant = Variant(
                chrom=chrom, pos=pos, ref=ref, alt=first_alt,
                rsid=rsid, gt=gt
            )

            key = variant.key
            genome.variants[key] = variant
            if rsid:
                genome.rsid_index[rsid] = key
            genome.pos_index[(chrom, pos)].append(key)
            genome.total += 1

            if genome.total % 200_000 == 0:
                logger.info(f"  ...parsed {genome.total:,} variants")

    elapsed = time.time() - t0
    logger.info(f"VCF parsed: {genome.total:,} variants, {len(genome.rsid_index):,} with rsIDs ({elapsed:.1f}s)")
    return genome


# ============================================================================
# SECTION 3: DOWNLOAD MANAGER
# ============================================================================

def download_file(url: str, dest: str, logger: logging.Logger,
                  chunk_size: int = 1024*1024, max_retries: int = 5) -> str:
    """Download file with progress, resume support, and retries."""
    if os.path.exists(dest):
        logger.info(f"  Cached: {os.path.basename(dest)}")
        return dest

    logger.info(f"  Downloading: {url}")
    tmp = dest + ".tmp"

    for attempt in range(1, max_retries + 1):
        try:
            # Resume from partial download if .tmp exists
            existing_size = os.path.getsize(tmp) if os.path.exists(tmp) else 0
            headers = {}
            if existing_size > 0:
                headers["Range"] = f"bytes={existing_size}-"
                logger.info(f"  Resuming from {existing_size / 1e6:.1f} MB (attempt {attempt}/{max_retries})")

            if HAS_REQUESTS and url.startswith("http"):
                r = requests.get(url, stream=True, timeout=(30, 300),
                                 headers=headers)
                r.raise_for_status()

                if r.status_code == 206:  # Partial content (resume)
                    content_range = r.headers.get("content-range", "")
                    total = int(content_range.split("/")[-1]) if "/" in content_range else 0
                    downloaded = existing_size
                else:
                    total = int(r.headers.get("content-length", 0))
                    downloaded = 0
                    existing_size = 0  # Server didn't support range, start fresh

                mode = "ab" if existing_size > 0 else "wb"
                with open(tmp, mode) as f:
                    for chunk in r.iter_content(chunk_size=chunk_size):
                        f.write(chunk)
                        downloaded += len(chunk)
                        if total and downloaded % (10 * chunk_size) == 0:
                            pct = downloaded / total * 100
                            logger.info(f"    {downloaded / 1e6:.1f} MB / {total / 1e6:.1f} MB ({pct:.0f}%)")
            else:
                req = Request(url)
                req.add_header("User-Agent", "GenomeAnalysisPipeline/1.0")
                if existing_size > 0:
                    req.add_header("Range", f"bytes={existing_size}-")
                resp = urlopen(req, timeout=300)
                mode = "ab" if existing_size > 0 else "wb"
                with open(tmp, mode) as f:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)

            shutil.move(tmp, dest)
            size_mb = os.path.getsize(dest) / 1e6
            logger.info(f"  Saved: {os.path.basename(dest)} ({size_mb:.1f} MB)")
            return dest

        except Exception as e:
            # Keep .tmp for resume on next attempt
            logger.warning(f"  Download attempt {attempt}/{max_retries} failed: {e}")
            if attempt == max_retries:
                logger.error(f"  Download failed after {max_retries} attempts: {e}")
                raise

def ensure_cache_dir(cache_dir: str) -> str:
    os.makedirs(cache_dir, exist_ok=True)
    return cache_dir


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


# ============================================================================
# SECTION 7: LAYER 4 — ANCESTRY INFERENCE (AIMs)
# ============================================================================

# Curated Ancestry-Informative Markers with allele frequencies from 1000 Genomes Phase 3
# Format: rsid → {alt_freq_EUR, alt_freq_EAS, alt_freq_AFR, alt_freq_SAS, alt_freq_AMR}
# Selected for high Fst (population differentiation) across superpopulations
AIMS_PANEL = {
    # Each entry: "allele" = the nucleotide whose frequency is given for each population
    # Frequencies are from 1000 Genomes Phase 3 superpopulations
    # Skin/hair pigmentation-related AIMs (very high Fst)
    "rs1426654": {"allele": "A", "EUR": 0.999, "EAS": 0.007, "AFR": 0.044, "SAS": 0.726, "AMR": 0.613},  # SLC24A5
    "rs16891982": {"allele": "G", "EUR": 0.959, "EAS": 0.003, "AFR": 0.014, "SAS": 0.077, "AMR": 0.393},  # SLC45A2
    "rs1545397": {"allele": "A", "EUR": 0.149, "EAS": 0.733, "AFR": 0.030, "SAS": 0.274, "AMR": 0.233},  # OCA2
    "rs12913832": {"allele": "A", "EUR": 0.789, "EAS": 0.007, "AFR": 0.015, "SAS": 0.088, "AMR": 0.256},  # HERC2
    "rs1805007": {"allele": "T", "EUR": 0.102, "EAS": 0.000, "AFR": 0.002, "SAS": 0.007, "AMR": 0.033},  # MC1R
    # Duffy antigen
    "rs2814778": {"allele": "T", "EUR": 0.993, "EAS": 1.000, "AFR": 0.003, "SAS": 0.995, "AMR": 0.777},  # DARC
    # EDAR — East Asian hair thickness
    "rs3827760": {"allele": "G", "EUR": 0.005, "EAS": 0.850, "AFR": 0.000, "SAS": 0.008, "AMR": 0.445},
    # ADH1B — alcohol metabolism
    "rs1229984": {"allele": "T", "EUR": 0.028, "EAS": 0.739, "AFR": 0.006, "SAS": 0.107, "AMR": 0.138},
    # ABCC11 — earwax type
    "rs17822931": {"allele": "T", "EUR": 0.082, "EAS": 0.918, "AFR": 0.015, "SAS": 0.079, "AMR": 0.354},
    # LCT — lactase persistence (G = REF = non-persistence allele)
    "rs4988235": {"allele": "G", "EUR": 0.735, "EAS": 0.003, "AFR": 0.067, "SAS": 0.253, "AMR": 0.281},
    # Various high-Fst AIMs from published panels
    "rs2065160": {"allele": "T", "EUR": 0.088, "EAS": 0.019, "AFR": 0.805, "SAS": 0.121, "AMR": 0.286},
    "rs1834640": {"allele": "A", "EUR": 0.530, "EAS": 0.093, "AFR": 0.969, "SAS": 0.317, "AMR": 0.647},
    "rs1871534": {"allele": "G", "EUR": 0.161, "EAS": 0.003, "AFR": 0.784, "SAS": 0.068, "AMR": 0.298},
    "rs3916235": {"allele": "T", "EUR": 0.127, "EAS": 0.000, "AFR": 0.743, "SAS": 0.045, "AMR": 0.212},
    "rs730570": {"allele": "A", "EUR": 0.834, "EAS": 0.111, "AFR": 0.053, "SAS": 0.327, "AMR": 0.414},
    "rs4411548": {"allele": "T", "EUR": 0.203, "EAS": 0.810, "AFR": 0.063, "SAS": 0.385, "AMR": 0.292},
    "rs260690": {"allele": "A", "EUR": 0.761, "EAS": 0.293, "AFR": 0.131, "SAS": 0.453, "AMR": 0.470},
    "rs2250072": {"allele": "A", "EUR": 0.284, "EAS": 0.840, "AFR": 0.087, "SAS": 0.446, "AMR": 0.377},
    "rs310644": {"allele": "T", "EUR": 0.677, "EAS": 0.199, "AFR": 0.098, "SAS": 0.324, "AMR": 0.395},
    "rs12498138": {"allele": "A", "EUR": 0.042, "EAS": 0.005, "AFR": 0.664, "SAS": 0.060, "AMR": 0.150},
    "rs174570": {"allele": "C", "EUR": 0.302, "EAS": 0.610, "AFR": 0.132, "SAS": 0.340, "AMR": 0.280},
    "rs1800404": {"allele": "C", "EUR": 0.285, "EAS": 0.847, "AFR": 0.118, "SAS": 0.386, "AMR": 0.341},
    "rs3811801": {"allele": "G", "EUR": 0.785, "EAS": 0.215, "AFR": 0.113, "SAS": 0.446, "AMR": 0.470},
    "rs6497268": {"allele": "A", "EUR": 0.073, "EAS": 0.004, "AFR": 0.622, "SAS": 0.048, "AMR": 0.149},
    "rs10497191": {"allele": "C", "EUR": 0.116, "EAS": 0.579, "AFR": 0.023, "SAS": 0.223, "AMR": 0.196},
    "rs6451722": {"allele": "G", "EUR": 0.600, "EAS": 0.112, "AFR": 0.876, "SAS": 0.388, "AMR": 0.563},
    "rs7657799": {"allele": "G", "EUR": 0.143, "EAS": 0.704, "AFR": 0.030, "SAS": 0.290, "AMR": 0.249},
    "rs2789823": {"allele": "A", "EUR": 0.686, "EAS": 0.175, "AFR": 0.110, "SAS": 0.368, "AMR": 0.382},
    "rs1079597": {"allele": "C", "EUR": 0.160, "EAS": 0.023, "AFR": 0.752, "SAS": 0.181, "AMR": 0.261},
    "rs6003": {"allele": "C", "EUR": 0.096, "EAS": 0.001, "AFR": 0.518, "SAS": 0.043, "AMR": 0.159},
    "rs1800414": {"allele": "A", "EUR": 0.005, "EAS": 0.534, "AFR": 0.001, "SAS": 0.012, "AMR": 0.074},
    "rs7495174": {"allele": "A", "EUR": 0.835, "EAS": 0.246, "AFR": 0.025, "SAS": 0.324, "AMR": 0.382},
    "rs1393350": {"allele": "A", "EUR": 0.219, "EAS": 0.015, "AFR": 0.024, "SAS": 0.059, "AMR": 0.098},
    "rs12821256": {"allele": "C", "EUR": 0.118, "EAS": 0.001, "AFR": 0.002, "SAS": 0.007, "AMR": 0.036},
    "rs4959270": {"allele": "A", "EUR": 0.410, "EAS": 0.792, "AFR": 0.157, "SAS": 0.473, "AMR": 0.377},
    "rs1408799": {"allele": "C", "EUR": 0.314, "EAS": 0.028, "AFR": 0.830, "SAS": 0.242, "AMR": 0.418},
    "rs2402130": {"allele": "A", "EUR": 0.310, "EAS": 0.028, "AFR": 0.026, "SAS": 0.101, "AMR": 0.135},
    "rs12203592": {"allele": "T", "EUR": 0.152, "EAS": 0.000, "AFR": 0.005, "SAS": 0.019, "AMR": 0.048},
    "rs1042602": {"allele": "A", "EUR": 0.372, "EAS": 0.005, "AFR": 0.061, "SAS": 0.167, "AMR": 0.209},
    "rs6119471": {"allele": "T", "EUR": 0.015, "EAS": 0.000, "AFR": 0.549, "SAS": 0.006, "AMR": 0.083},
}

POPULATION_NAMES = {
    "EUR": "European",
    "EAS": "East Asian",
    "AFR": "African",
    "SAS": "South Asian",
    "AMR": "Americas (Admixed)",
}

def calc_ancestry_likelihood(genome: GenomeData, aims: Dict,
                             logger: logging.Logger) -> Dict:
    """
    Estimate ancestry proportions using log-likelihood method.
    For each AIM, compute P(genotype | population) and sum log-likelihoods.
    Each AIM entry has an "allele" field specifying which nucleotide the
    frequency values correspond to. We count that allele in the genotype.
    """
    populations = ["EUR", "EAS", "AFR", "SAS", "AMR"]
    log_likes = {pop: 0.0 for pop in populations}
    used_markers = 0
    eps = 0.001  # avoid log(0)

    for rsid, pop_freqs in aims.items():
        if rsid not in genome.rsid_index:
            continue

        key = genome.rsid_index[rsid]
        variant = genome.variants[key]
        a1, a2 = variant.genotype_alleles

        # Count copies of the specified frequency allele
        freq_allele = pop_freqs.get("allele", variant.alt).upper()
        dosage = (1 if a1.upper() == freq_allele else 0) + \
                 (1 if a2.upper() == freq_allele else 0)

        used_markers += 1
        for pop in populations:
            p = pop_freqs.get(pop, 0.5)
            p = max(eps, min(1 - eps, p))  # clamp

            # P(genotype | pop) assuming HWE: p^2, 2pq, q^2
            if dosage == 0:
                prob = (1 - p) ** 2
            elif dosage == 1:
                prob = 2 * p * (1 - p)
            else:
                prob = p ** 2

            log_likes[pop] += math.log(max(prob, eps))

    if used_markers == 0:
        logger.warning("  No AIM markers found in VCF!")
        return {pop: 1.0 / len(populations) for pop in populations}

    # Convert log-likelihoods to proportions via softmax
    max_ll = max(log_likes.values())
    exp_likes = {pop: math.exp(ll - max_ll) for pop, ll in log_likes.items()}
    total = sum(exp_likes.values())
    proportions = {pop: exp_likes[pop] / total for pop in populations}

    logger.info(f"  Used {used_markers}/{len(aims)} AIMs for ancestry estimation")
    return proportions

def run_layer4_ancestry(genome: GenomeData,
                        logger: logging.Logger) -> Dict:
    """Run ancestry inference. Returns results dict."""
    logger.info("=" * 60)
    logger.info("LAYER 4: Ancestry Inference (AIMs)")
    logger.info("=" * 60)

    proportions = calc_ancestry_likelihood(genome, AIMS_PANEL, logger)

    # Count how many AIMs were found
    found = sum(1 for rs in AIMS_PANEL if rs in genome.rsid_index)

    # Build AIM details list
    aim_details = []
    for rsid, pop_freqs in sorted(AIMS_PANEL.items()):
        if rsid not in genome.rsid_index:
            continue
        key = genome.rsid_index[rsid]
        variant = genome.variants[key]
        a1, a2 = variant.genotype_alleles
        aim_details.append({
            "rsid": rsid,
            "genotype": f"{a1}/{a2}",
            "freq_allele": pop_freqs.get("allele", ""),
            "EUR": pop_freqs["EUR"],
            "EAS": pop_freqs["EAS"],
            "AFR": pop_freqs["AFR"],
            "SAS": pop_freqs["SAS"],
            "AMR": pop_freqs["AMR"],
        })

    logger.info(f"Layer 4 complete: {found} AIMs found")
    return {"layer": 4, "proportions": proportions, "aims_found": found,
            "aim_details": aim_details}


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


# ============================================================================
# SECTION 9: JSON REPORT GENERATOR
# ============================================================================

def generate_json_report(layer_results: Dict, output_dir: str, genome: GenomeData,
                         logger: logging.Logger) -> str:
    """Generate unified genome_report.json from all layer results."""
    logger.info("=" * 60)
    logger.info("GENERATING JSON REPORT")
    logger.info("=" * 60)

    report = {
        "meta": {
            "date": datetime.now().isoformat(),
            "sample_id": genome.sample_id,
            "total_variants": genome.total,
            "variants_with_rsid": len(genome.rsid_index),
            "reference_build": "GRCh37/hg19",
            "tool": "genome_analysis.py",
            "disclaimer": "For informational/educational purposes only. Not a substitute for clinical genetic testing or medical advice.",
        },
        "layer1_prs": layer_results.get(1, {}),
        "layer2_pharmacogenetics": layer_results.get(2, {}),
        "layer3_clinvar": layer_results.get(3, {}),
        "layer4_ancestry": layer_results.get(4, {}),
        "layer5_gwas": layer_results.get(5, {}),
    }

    output_file = os.path.join(output_dir, "genome_report.json")
    with open(output_file, "w") as f:
        json.dump(report, f, indent=2, default=str, ensure_ascii=False)

    logger.info(f"JSON report: {output_file}")
    return output_file


# ============================================================================
# SECTION 10: MAIN ENTRY POINT & CLI
# ============================================================================

def setup_logging(verbose: bool = False) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT, stream=sys.stderr)
    return logging.getLogger("genome_analysis")

def main():
    parser = argparse.ArgumentParser(
        description="Personal Genome Analysis Pipeline — 5-layer analysis from VCF",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Layers:
  1  Polygenic Risk Scores (PGS Catalog)
  2  Pharmacogenetics (CPIC)
  3  ClinVar Pathogenic Variants
  4  Ancestry Inference (AIMs)
  5  GWAS Trait Associations

Output: genome_report.json (single unified JSON file)

Examples:
  python genome_analysis.py --vcf my.vcf --output-dir ./results
  python genome_analysis.py --vcf my.vcf --layer 2 --layer 3
        """
    )
    parser.add_argument("--vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--output-dir", default="./results", help="Output directory for genome_report.json")
    parser.add_argument("--cache-dir", default="./cache", help="Cache directory for downloaded databases")
    parser.add_argument("--pgs-dir", default=None,
                        help="Directory with PGS Catalog scoring files (.txt.gz). "
                             "If provided, Layer 1 uses full PGS instead of curated SNPs.")
    parser.add_argument("--layer", type=int, action="append",
                        help="Run specific layer(s) only (1-5). Omit to run all.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()
    logger = setup_logging(args.verbose)

    if not os.path.exists(args.vcf):
        logger.error(f"VCF file not found: {args.vcf}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    ensure_cache_dir(args.cache_dir)

    layers_to_run = set(args.layer) if args.layer else {1, 2, 3, 4, 5}

    logger.info("=" * 60)
    logger.info("PERSONAL GENOME ANALYSIS PIPELINE")
    logger.info(f"VCF: {args.vcf}")
    logger.info(f"Layers: {sorted(layers_to_run)}")
    logger.info(f"Output: {args.output_dir}")
    logger.info(f"Cache: {args.cache_dir}")
    logger.info("=" * 60)

    # Parse VCF
    genome = parse_vcf(args.vcf, logger)

    layer_results = {}
    t_start = time.time()

    # Layer 1: PRS
    if 1 in layers_to_run:
        try:
            layer_results[1] = run_layer1_prs(genome, args.cache_dir, logger,
                                              pgs_dir=args.pgs_dir)
        except Exception as e:
            logger.error(f"Layer 1 failed: {e}")
            layer_results[1] = {"layer": 1, "error": str(e)}

    # Layer 2: Pharmacogenetics
    if 2 in layers_to_run:
        try:
            layer_results[2] = run_layer2_pharmacogenetics(genome, logger)
        except Exception as e:
            logger.error(f"Layer 2 failed: {e}")
            layer_results[2] = {"layer": 2, "error": str(e)}

    # Layer 3: ClinVar
    if 3 in layers_to_run:
        try:
            layer_results[3] = run_layer3_clinvar(genome, args.cache_dir, logger)
        except Exception as e:
            logger.error(f"Layer 3 failed: {e}")
            layer_results[3] = {"layer": 3, "error": str(e)}

    # Layer 4: Ancestry
    if 4 in layers_to_run:
        try:
            layer_results[4] = run_layer4_ancestry(genome, logger)
        except Exception as e:
            logger.error(f"Layer 4 failed: {e}")
            layer_results[4] = {"layer": 4, "error": str(e)}

    # Layer 5: GWAS
    if 5 in layers_to_run:
        try:
            layer_results[5] = run_layer5_gwas(genome, args.cache_dir, logger)
        except Exception as e:
            logger.error(f"Layer 5 failed: {e}")
            layer_results[5] = {"layer": 5, "error": str(e)}

    # Always generate JSON report
    generate_json_report(layer_results, args.output_dir, genome, logger)

    elapsed = time.time() - t_start
    logger.info("=" * 60)
    logger.info(f"PIPELINE COMPLETE in {elapsed:.0f}s")
    logger.info(f"Output files in: {args.output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
