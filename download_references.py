#!/usr/bin/env python3
"""
Download and merge GRCh37/hg19 reference panels for supervised ADMIXTURE.

Final outputs:
    cache/reference_panel/merged.bed
    cache/reference_panel/merged.bim
    cache/reference_panel/merged.fam
    cache/reference_panel/merged.pop

Prerequisites:
    macOS: brew install bcftools plink2

AADR is optional. It is distributed as EIGENSTRAT and conversion is best done
with EIGENSOFT convertf. Ancient AADR genotypes are often low coverage and can
increase missingness/batch effects, so use --skip-aadr for a present-day-only
reference panel. The downstream user sample should append "-" to merged.pop.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from core import LOG_FORMAT, download_file, ensure_cache_dir


AUTOSOMES = tuple(str(i) for i in range(1, 23))
PopMap = Dict[Tuple[str, str], str]

ONEKG_BASE = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
ONEKG_PANEL = "integrated_call_samples_v3.20130502.ALL.panel"
ONEKG_VCF = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

HGDP_BASE = "ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516"
HGDP_VCF = "hgdp_wgs.20190516.full.chr{chrom}.vcf.gz"
HGDP_METADATA = (
    f"{HGDP_BASE}/hgdp_wgs.20190516.metadata.txt",
    f"{HGDP_BASE}/hgdp_wgs.20190516.sample_info.txt",
    f"{HGDP_BASE}/sample_info.txt",
)

SGDP_TAR = (
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
    "vcf_variants/vcfs.variants.public_samples.279samples.tar"
)
SGDP_METADATA = (
    "https://sf-web-assets-prod.s3.amazonaws.com/wp-content/uploads/"
    "2023/05/05_02_2023_SGDP_metainformation_update.txt",
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/"
    "SGDP_metadata.279public.21signedLetter.44Fan.samples.txt",
)

AADR_BASE = (
    "https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/"
    "curated_releases/V54/V54.1.p1/SHARE/public.dir/dataverse"
)
AADR = {
    "geno": "aadr_v54.1.p1_1240K_public.geno.gz",
    "snp": "aadr_v54.1.p1_1240K_public.snp",
    "ind": "aadr_v54.1.p1_1240K_public.ind",
    "anno": "aadr_v54.1.p1_1240K_public.anno",
}


@dataclass(frozen=True)
class SourceResult:
    name: str
    fid: str
    prefix: Path
    pop_map: PopMap


def logger() -> logging.Logger:
    logging.basicConfig(level=logging.INFO, format=LOG_FORMAT)
    return logging.getLogger("download_references")


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def bed_exists(prefix: Path) -> bool:
    return all(prefix.with_suffix(ext).exists() for ext in (".bed", ".bim", ".fam"))


def pgen_exists(prefix: Path) -> bool:
    return all(prefix.with_suffix(ext).exists() for ext in (".pgen", ".pvar", ".psam"))


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lines), encoding="utf-8")


def run(cmd: Sequence[object], log: logging.Logger) -> None:
    log.info("  $ %s", " ".join(str(x) for x in cmd))
    subprocess.run([str(x) for x in cmd], check=True)


def require_tools() -> None:
    missing = [tool for tool in ("bcftools", "plink2") if shutil.which(tool) is None]
    if missing:
        raise SystemExit(
            "Missing required command(s): "
            + ", ".join(missing)
            + "\n\nInstall on macOS:\n  brew install bcftools plink2\n"
        )


def dl(url: str, dest: Path, log: logging.Logger, retries: int = 5) -> Path:
    mkdir(dest.parent)
    return Path(download_file(url, str(dest), log, max_retries=retries))


def dl_first(urls: Sequence[str], dest: Path, log: logging.Logger, required: bool) -> Optional[Path]:
    if nonempty(dest):
        log.info("  Cached: %s", dest.name)
        return dest
    last: Optional[Exception] = None
    for url in urls:
        try:
            return dl(url, dest, log, retries=1)
        except Exception as exc:
            last = exc
            log.warning("  Candidate failed: %s", url)
    if required:
        raise RuntimeError(f"Could not download {dest.name}: {last}")
    log.warning("  Optional metadata not available: %s", dest.name)
    return None


def fields(line: str) -> List[str]:
    line = line.strip()
    if "\t" in line:
        return [x.strip() for x in line.split("\t")]
    if "," in line and line.count(",") >= 2:
        return [x.strip() for x in next(csv.reader([line]))]
    return re.split(r"\s+", line)


def key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def pick(headers: Sequence[str], candidates: Sequence[str]) -> Optional[int]:
    normalized = [key(h) for h in headers]
    for candidate in candidates:
        want = key(candidate)
        for idx, header in enumerate(normalized):
            if header == want:
                return idx
    for candidate in candidates:
        want = key(candidate)
        for idx, header in enumerate(normalized):
            if want in header:
                return idx
    return None


def parse_table_metadata(path: Optional[Path], fid: str, log: logging.Logger) -> PopMap:
    if path is None or not path.exists():
        return {}
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    sample_names = ("sample", "sample_id", "sample_name", "id", "genetic_id", "genetic id")
    pop_names = ("pop", "population", "population_id", "ethnicity", "ethnic", "group")
    for header_idx, line in enumerate(lines[:50]):
        header = fields(line)
        sample_col, pop_col = pick(header, sample_names), pick(header, pop_names)
        if sample_col is None or pop_col is None:
            continue
        out: PopMap = {}
        for row in lines[header_idx + 1 :]:
            row_fields = fields(row)
            if max(sample_col, pop_col) < len(row_fields):
                sample = row_fields[sample_col]
                pop = row_fields[pop_col].replace(" ", "_")
                if sample and pop:
                    out[(fid, sample)] = pop
        log.info("  Parsed %d labels from %s", len(out), path.name)
        return out
    log.warning("  Could not infer metadata columns in %s", path.name)
    return {}


def parse_1000g_panel(path: Path) -> PopMap:
    out: PopMap = {}
    with path.open("r", encoding="utf-8") as handle:
        header = [key(x) for x in fields(handle.readline())]
        sample_col = header.index("sample") if "sample" in header else 0
        pop_col = header.index("pop") if "pop" in header else 1
        for line in handle:
            row = fields(line)
            if max(sample_col, pop_col) < len(row):
                out[("1kg", row[sample_col])] = row[pop_col]
    return out


def parse_ind(path: Path, fid: str) -> PopMap:
    out: PopMap = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        row = fields(line)
        if len(row) >= 3:
            out[(fid, row[0])] = row[2].replace(" ", "_")
    return out


def fam_labels(prefix: Path, fid: str) -> PopMap:
    out: PopMap = {}
    fam = prefix.with_suffix(".fam")
    if not fam.exists():
        return out
    for line in fam.read_text(encoding="utf-8", errors="replace").splitlines():
        row = line.split()
        if len(row) >= 2:
            label = row[0] if row[0] != fid else fid
            out[(row[0], row[1])] = label
            out[(fid, row[1])] = label
    return out


def has_index(vcf: Path) -> bool:
    return Path(str(vcf) + ".tbi").exists() or Path(str(vcf) + ".csi").exists()


def vcf_to_pgen(vcf: Path, out: Path, fid: str, threads: int, log: logging.Logger) -> Path:
    if pgen_exists(out):
        log.info("  Cached pgen: %s", out.name)
        return out
    filtered = out.parent / f"{out.name}.biallelic-snps.vcf.gz"
    if not nonempty(filtered):
        run(
            ["bcftools", "view", "--threads", threads, "-m", "2", "-M", "2", "-v", "snps", "-Oz", "-o", filtered, vcf],
            log,
        )
        run(["bcftools", "index", "-f", filtered], log)
    run(
        [
            "plink2", "--vcf", filtered, "--const-fid", fid, "--vcf-half-call", "missing",
            "--snps-only", "just-acgt", "--max-alleles", "2", "--maf", "0.05",
            "--set-all-var-ids", "@:#:$r:$a", "--new-id-max-allele-len", "200", "truncate",
            "--make-pgen", "--out", out,
        ],
        log,
    )
    return out


def merge_pfiles_to_bed(pgens: Sequence[Path], out: Path, log: logging.Logger) -> Path:
    if bed_exists(out):
        log.info("  Cached bed: %s", out.name)
        return out
    if len(pgens) == 1:
        run(["plink2", "--pfile", pgens[0], "--make-bed", "--out", out], log)
        return out
    merge_list = out.parent / f"{out.name}.pmerge-list"
    write_lines(merge_list, (f"{p}\n" for p in pgens[1:]))
    run(["plink2", "--pfile", pgens[0], "--pmerge-list", merge_list, "--make-bed", "--out", out], log)
    return out


def build_vcf_source(
    name: str,
    fid: str,
    vcfs: Sequence[Path],
    work_dir: Path,
    threads: int,
    log: logging.Logger,
    pop_map: Optional[PopMap] = None,
) -> SourceResult:
    out = work_dir / name
    if not bed_exists(out):
        pgens = [vcf_to_pgen(vcf, work_dir / f"chr{i + 1}", fid, threads, log) for i, vcf in enumerate(vcfs)]
        merge_pfiles_to_bed(pgens, out, log)
    return SourceResult(name, fid, out, pop_map or {})


def build_1000g(root: Path, threads: int, log: logging.Logger) -> SourceResult:
    log.info("Preparing 1000 Genomes Phase 3")
    source_dir, work_dir = mkdir(root / "sources" / "1000g"), mkdir(root / "work" / "1000g")
    panel = dl(f"{ONEKG_BASE}/{ONEKG_PANEL}", source_dir / ONEKG_PANEL, log)
    vcfs = []
    for chrom in AUTOSOMES:
        name = ONEKG_VCF.format(chrom=chrom)
        vcfs.append(dl(f"{ONEKG_BASE}/{name}", source_dir / name, log))
        dl(f"{ONEKG_BASE}/{name}.tbi", source_dir / f"{name}.tbi", log)
    return build_vcf_source("1000g", "1kg", vcfs, work_dir, threads, log, parse_1000g_panel(panel))


def build_hgdp(root: Path, threads: int, log: logging.Logger) -> SourceResult:
    log.info("Preparing HGDP")
    source_dir, work_dir = mkdir(root / "sources" / "hgdp"), mkdir(root / "work" / "hgdp")
    metadata = dl_first(HGDP_METADATA, source_dir / "hgdp_metadata.txt", log, required=False)
    vcfs = []
    for chrom in AUTOSOMES:
        name = HGDP_VCF.format(chrom=chrom)
        vcfs.append(dl(f"{HGDP_BASE}/{name}", source_dir / name, log))
        dl(f"{HGDP_BASE}/{name}.tbi", source_dir / f"{name}.tbi", log)
    return build_vcf_source("hgdp", "hgdp", vcfs, work_dir, threads, log, parse_table_metadata(metadata, "hgdp", log))


def sgdp_vcf(root: Path, tar_path: Path, threads: int, log: logging.Logger) -> Path:
    work_dir = mkdir(root / "work" / "sgdp")
    merged = work_dir / "sgdp.merged.vcf.gz"
    if nonempty(merged):
        return merged

    extract_dir = mkdir(root / "sources" / "sgdp" / "vcfs")
    marker = extract_dir / ".extract-complete"
    if not marker.exists():
        log.info("  Extracting SGDP tarball")
        with tarfile.open(tar_path, "r:*") as archive:
            for member in archive.getmembers():
                if member.isfile() and re.search(r"\.vcf(\.gz|\.bgz)?$|\.tbi$", member.name):
                    archive.extract(member, extract_dir)
        marker.touch()

    vcfs = sorted(p for p in extract_dir.rglob("*") if p.is_file() and re.search(r"\.vcf(\.gz|\.bgz)$", p.name))
    plain = sorted(p for p in extract_dir.rglob("*.vcf") if p.is_file())
    if len(vcfs) == 1:
        return vcfs[0]
    if not vcfs and len(plain) == 1:
        return plain[0]
    if not vcfs:
        raise RuntimeError(f"No SGDP VCFs found under {extract_dir}")

    for vcf in vcfs:
        if not has_index(vcf):
            run(["bcftools", "index", "-f", vcf], log)
    run(["bcftools", "merge", "--threads", threads, "-m", "none", "-Oz", "-o", merged, *vcfs], log)
    run(["bcftools", "index", "-f", merged], log)
    return merged


def build_sgdp(root: Path, threads: int, log: logging.Logger) -> SourceResult:
    log.info("Preparing SGDP")
    source_dir, work_dir = mkdir(root / "sources" / "sgdp"), mkdir(root / "work" / "sgdp")
    metadata = dl_first(SGDP_METADATA, source_dir / "sgdp_metadata.txt", log, required=False)
    tar_path = dl(SGDP_TAR, source_dir / "vcfs.variants.public_samples.279samples.tar", log)
    out = work_dir / "sgdp"
    if not bed_exists(out):
        pgen = vcf_to_pgen(sgdp_vcf(root, tar_path, threads, log), work_dir / "sgdp_pgen", "sgdp", threads, log)
        merge_pfiles_to_bed([pgen], out, log)
    return SourceResult("sgdp", "sgdp", out, parse_table_metadata(metadata, "sgdp", log))


def build_aadr(root: Path, log: logging.Logger) -> Optional[SourceResult]:
    log.info("Preparing AADR v54.1.p1")
    source_dir, work_dir = mkdir(root / "sources" / "aadr"), mkdir(root / "work" / "aadr")
    for filename in AADR.values():
        dl(f"{AADR_BASE}/{filename}", source_dir / filename, log)
    ind = source_dir / AADR["ind"]
    out, raw = work_dir / "aadr", work_dir / "aadr_raw"
    if bed_exists(out):
        return SourceResult("aadr", "aadr", out, parse_ind(ind, "aadr"))

    convertf = shutil.which("convertf")
    if convertf is None:
        log.warning("  AADR downloaded but not converted: install EIGENSOFT convertf or rerun with --skip-aadr")
        return None

    par = work_dir / "convertf.par"
    write_lines(
        par,
        [
            f"genotypename: {source_dir / AADR['geno']}\n",
            f"snpname: {source_dir / AADR['snp']}\n",
            f"indivname: {ind}\n",
            "outputformat: PACKEDPED\n",
            f"genotypeoutname: {raw.with_suffix('.bed')}\n",
            f"snpoutname: {raw.with_suffix('.bim')}\n",
            f"indivoutname: {raw.with_suffix('.fam')}\n",
            "familynames: NO\n",
        ],
    )
    run([convertf, "-p", par], log)
    run(
        [
            "plink2", "--bfile", raw, "--autosome", "--snps-only", "just-acgt",
            "--max-alleles", "2", "--maf", "0.05", "--make-bed", "--out", out,
        ],
        log,
    )
    return SourceResult("aadr", "aadr", out, parse_ind(ind, "aadr"))


def merge_sources(results: Sequence[SourceResult], root: Path, log: logging.Logger) -> Path:
    log.info("Merging source panels")
    merge_dir, combined = mkdir(root / "work" / "merge"), root / "work" / "merge" / "combined"
    if bed_exists(combined):
        return combined
    pgens = []
    for result in results:
        pgen = merge_dir / result.name
        if not pgen_exists(pgen):
            run(["plink2", "--bfile", result.prefix, "--make-pgen", "--out", pgen], log)
        pgens.append(pgen)
    return merge_pfiles_to_bed(pgens, combined, log)


def final_filter_and_prune(combined: Path, final: Path, root: Path, log: logging.Logger) -> None:
    log.info("Filtering and LD-pruning merged panel")
    if bed_exists(final):
        return
    filtered, prune = root / "work" / "merged.filtered", root / "work" / "merged.ld"
    if not bed_exists(filtered):
        run(
            [
                "plink2", "--bfile", combined, "--autosome", "--snps-only", "just-acgt",
                "--max-alleles", "2", "--maf", "0.05", "--make-bed", "--out", filtered,
            ],
            log,
        )
    if not nonempty(prune.with_suffix(".prune.in")):
        run(["plink2", "--bfile", filtered, "--indep-pairwise", "50", "5", "0.1", "--out", prune], log)
    run(["plink2", "--bfile", filtered, "--extract", prune.with_suffix(".prune.in"), "--make-bed", "--out", final], log)


def write_pop(final: Path, results: Sequence[SourceResult], log: logging.Logger) -> None:
    fam, pop = final.with_suffix(".fam"), final.with_suffix(".pop")
    rows = [line.split() for line in fam.read_text(encoding="utf-8").splitlines() if line.strip()]
    if nonempty(pop) and len(pop.read_text(encoding="utf-8").splitlines()) == len(rows):
        log.info("  Cached pop labels: %s", pop)
        return

    labels: PopMap = {}
    for result in results:
        labels.update(fam_labels(result.prefix, result.fid))
        labels.update(result.pop_map)

    missing, output = 0, []
    for row in rows:
        fid, iid = row[0], row[1]
        label = labels.get((fid, iid)) or labels.get(("*", iid))
        if not label:
            label = fid if fid not in ("0", "") else "-"
            missing += 1
        output.append(f"{label}\n")
    write_lines(pop, output)
    log.info("  Wrote %s%s", pop, f" ({missing} fallback labels)" if missing else "")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a merged GRCh37/hg19 reference panel for ADMIXTURE.")
    parser.add_argument("--cache-dir", default="./cache", help="Cache root directory (default: ./cache)")
    parser.add_argument("--skip-aadr", action="store_true", help="Skip optional AADR EIGENSTRAT conversion")
    parser.add_argument(
        "--threads",
        type=int,
        default=max(1, min(8, (os.cpu_count() or 2) // 2)),
        help="Threads for bcftools/plink2 where supported (default: half CPUs, max 8)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = parse_args(argv)
    log = logger()
    require_tools()
    root = mkdir(Path(ensure_cache_dir(args.cache_dir)).resolve() / "reference_panel")
    final = root / "merged"
    if bed_exists(final) and nonempty(final.with_suffix(".pop")):
        log.info("Final merged panel already exists: %s", final)
        return

    results = [build_1000g(root, args.threads, log), build_hgdp(root, args.threads, log), build_sgdp(root, args.threads, log)]
    if args.skip_aadr:
        log.info("Skipping AADR by request")
    else:
        aadr = build_aadr(root, log)
        if aadr is not None:
            results.append(aadr)
    combined = merge_sources(results, root, log)
    final_filter_and_prune(combined, final, root, log)
    write_pop(final, results, log)


if __name__ == "__main__":
    main()
