#!/usr/bin/env python3
"""
Core data structures, VCF parser, and download utilities.
"""

import gzip
import logging
import os
import shutil
import time
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
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
