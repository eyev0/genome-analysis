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
import json
import logging
import os
import sys
import time
from datetime import datetime
from typing import Dict

from core import (
    LOG_FORMAT, Variant, GenomeData,
    normalize_chrom, parse_vcf, download_file, ensure_cache_dir,
)
from layers.prs import run_layer1_prs
from layers.pharmacogenetics import run_layer2_pharmacogenetics
from layers.clinvar import run_layer3_clinvar
from layers.ancestry import run_layer4_ancestry
from layers.gwas import run_layer5_gwas


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
    parser.add_argument("--output-dir", default="./reports", help="Output directory for genome_report.json")
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
