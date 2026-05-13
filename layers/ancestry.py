import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple

from core import GenomeData, download_file

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

EUROGENES_K13_LABELS = "North_Atlantic Baltic West_Med West_Asian East_Med Red_Sea South_Asian East_Asian Siberian Amerindian Oceanian Northeast_African Sub_Saharan".split()
MDLP_K23_LABELS = "Amerindian Ancestral_Altaic South_Central_Asian Arctic South_Indian Australoid Austronesian Caucasian Archaic_African East_African East_Siberian European_Early_Farmers Khoisan Melano_Polynesian Archaic_Human North_African Paleo_Siberian African_Pygmy South_East_Asian Subsaharian Tungus_Altaic European_Hunters_Gatherers East_Asian".split()

K_CALCULATOR_LABELS = {"Eurogenes_K13": EUROGENES_K13_LABELS, "MDLP_K23": MDLP_K23_LABELS}
Y_HAPLOGROUP_ANNOTATIONS = {
    "R1A": {"geographic_origin": "Eastern Europe and Central/South Asia", "age_kya": 22, "description": "Common in Slavic, Indo-Iranian, and Central Asian paternal lineages."},
    "R1B": {"geographic_origin": "Western Europe", "age_kya": 20, "description": "The dominant paternal lineage in much of Atlantic and western Europe."},
    "I1": {"geographic_origin": "Northern Europe", "age_kya": 27, "description": "Frequent in Scandinavia and associated with northern European founder expansions."},
    "I2": {"geographic_origin": "Southeastern and Eastern Europe", "age_kya": 22, "description": "A European paternal lineage with strong Balkan and eastern European branches."},
    "J1": {"geographic_origin": "Near East and Arabian Peninsula", "age_kya": 20, "description": "Common in Arabian, Levantine, and some Caucasus paternal lineages."},
    "J2": {"geographic_origin": "Near East, Anatolia, and Mediterranean", "age_kya": 27, "description": "Linked to Near Eastern and Mediterranean Neolithic and Bronze Age expansions."},
    "G": {"geographic_origin": "Caucasus, Anatolia, and Europe", "age_kya": 26, "description": "Often associated with early Neolithic farmer ancestry in Europe."},
    "E1B1B": {"geographic_origin": "Northeast Africa and Mediterranean", "age_kya": 25, "description": "Found across North/East Africa, the Near East, and southern Europe."},
    "Q": {"geographic_origin": "North Asia and the Americas", "age_kya": 31, "description": "A Siberian-rooted lineage important in Native American paternal ancestry."},
    "N": {"geographic_origin": "Northern Eurasia", "age_kya": 20, "description": "Frequent in Uralic-speaking and northern Eurasian populations."},
}
MT_HAPLOGROUP_ANNOTATIONS = {
    "H": {"geographic_origin": "Europe and Near East", "age_kya": 20, "description": "The most common maternal lineage in Europe."},
    "U": {"geographic_origin": "West Eurasia", "age_kya": 45, "description": "An old West Eurasian maternal lineage with European and South Asian branches."},
    "J": {"geographic_origin": "Near East and Europe", "age_kya": 35, "description": "A West Eurasian maternal lineage expanded into Europe after the Ice Age."},
    "T": {"geographic_origin": "Near East and Europe", "age_kya": 25, "description": "A West Eurasian lineage common in Europe, Anatolia, and the Near East."},
    "K": {"geographic_origin": "West Eurasia", "age_kya": 16, "description": "A branch of U, common in Europe and the Near East."},
    "V": {"geographic_origin": "Western Europe", "age_kya": 15, "description": "A European maternal lineage enriched in some Atlantic and Saami populations."},
    "X": {"geographic_origin": "West Eurasia and Native America", "age_kya": 30, "description": "A low-frequency lineage found in West Eurasia and some Native American groups."},
    "L": {"geographic_origin": "Sub-Saharan Africa", "age_kya": 150, "description": "A deep African maternal lineage ancestral to all non-African mtDNA branches."},
    "M": {"geographic_origin": "Asia and Africa", "age_kya": 60, "description": "One of the two major non-African maternal macro-haplogroups."},
}


def _run_aims(genome: GenomeData, logger: logging.Logger) -> Dict:
    """Estimate ancestry proportions using the existing AIM log-likelihood method."""
    populations = ["EUR", "EAS", "AFR", "SAS", "AMR"]
    log_likes = {pop: 0.0 for pop in populations}
    used_markers = 0
    eps = 0.001  # avoid log(0)

    for rsid, pop_freqs in AIMS_PANEL.items():
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
        proportions = {pop: 1.0 / len(populations) for pop in populations}
    else:
        # Convert log-likelihoods to proportions via softmax
        max_ll = max(log_likes.values())
        exp_likes = {pop: math.exp(ll - max_ll) for pop, ll in log_likes.items()}
        total = sum(exp_likes.values())
        proportions = {pop: exp_likes[pop] / total for pop in populations}
        logger.info(f"  Used {used_markers}/{len(AIMS_PANEL)} AIMs for ancestry estimation")

    found = sum(1 for rs in AIMS_PANEL if rs in genome.rsid_index)

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

    logger.info(f"  AIMs complete: {found} markers found")
    return {"proportions": proportions, "aims_found": found, "aim_details": aim_details}


def _write_genome_vcf(genome: GenomeData, path: str,
                      chroms: Optional[Sequence[str]] = None,
                      allowed_rsids: Optional[set] = None) -> int:
    chrom_filter = set(chroms) if chroms else None
    count = 0
    with open(path, "w", encoding="utf-8") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=genome_analysis_layer4\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        f.write(f"{genome.sample_id}\n")
        for variant in genome.variants.values():
            if chrom_filter and variant.chrom not in chrom_filter:
                continue
            if allowed_rsids is not None and variant.rsid not in allowed_rsids:
                continue
            if not variant.is_snp:
                continue
            gt = variant.gt if variant.gt and "." not in variant.gt else "0/0"
            rsid = variant.rsid or "."
            f.write(
                f"{variant.chrom}\t{variant.pos}\t{rsid}\t{variant.ref}\t"
                f"{variant.alt}\t.\tPASS\t.\tGT\t{gt}\n"
            )
            count += 1
    return count


def _run_command(cmd: List[str], logger: logging.Logger,
                 cwd: Optional[str] = None) -> subprocess.CompletedProcess:
    logger.info(f"  Running: {' '.join(cmd)}")
    try:
        return subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        if len(detail) > 1200:
            detail = detail[-1200:]
        raise RuntimeError(f"{cmd[0]} failed: {detail or exc}") from exc


def _read_bim_rsids(bim_path: str) -> set:
    rsids = set()
    with open(bim_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2 and parts[1] != ".":
                rsids.add(parts[1])
    return rsids


def _count_bim_snps(bim_path: str) -> int:
    with open(bim_path, "r", encoding="utf-8") as f:
        return sum(1 for _ in f)


def _fam_samples(fam_path: str) -> List[Tuple[str, str]]:
    samples = []
    with open(fam_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                samples.append((parts[0], parts[1]))
    return samples


def _read_reference_populations(pop_path: str, fam_path: str) -> Tuple[Dict[str, str], List[str]]:
    fam = _fam_samples(fam_path)
    labels_by_sample = {}
    order = []
    with open(pop_path, "r", encoding="utf-8") as f:
        for idx, line in enumerate(f):
            parts = line.split()
            if not parts:
                continue
            if len(parts) == 1:
                if idx >= len(fam):
                    continue
                fid, iid = fam[idx]
                sample_keys = [iid, f"{fid}:{iid}"]
                label = parts[0]
            elif len(parts) == 2:
                sample_keys = [parts[0]]
                label = parts[1]
            else:
                sample_keys = [parts[1], f"{parts[0]}:{parts[1]}"]
                label = parts[-1]
            for sample_key in sample_keys:
                labels_by_sample[sample_key] = label
            if label != "-" and label not in order:
                order.append(label)
    return labels_by_sample, order


def _write_supervised_pop(ref_pop: str, ref_fam: str, merged_fam: str,
                          out_pop: str) -> Tuple[List[str], int]:
    labels_by_sample, order = _read_reference_populations(ref_pop, ref_fam)
    unknown = 0
    with open(out_pop, "w", encoding="utf-8") as out:
        for fid, iid in _fam_samples(merged_fam):
            label = labels_by_sample.get(iid) or labels_by_sample.get(f"{fid}:{iid}")
            if label is None:
                label = "-"
                unknown += 1
            out.write(f"{label}\n")
    return order, unknown


def _run_admixture(genome: GenomeData, cache_dir: str, logger: logging.Logger) -> Dict:
    ref_prefix = os.path.join(cache_dir, "reference_panel", "merged")
    required = [f"{ref_prefix}.{ext}" for ext in ("bed", "bim", "fam")]
    if not all(os.path.exists(path) for path in required):
        return {"skipped": True, "reason": "reference panel not built; run download_references.py"}

    ref_pop = f"{ref_prefix}.pop"
    if not os.path.exists(ref_pop):
        return {"skipped": True, "reason": "reference population labels missing; expected cache/reference_panel/merged.pop"}

    plink2 = shutil.which("plink2")
    if not plink2:
        return {"skipped": True, "reason": "plink2 not found; install with `brew install plink2`"}

    admixture = shutil.which("admixture")
    if not admixture:
        return {"skipped": True, "reason": "admixture not found; install with `brew install brewsci/bio/admixture` or from the ADMIXTURE project"}

    out_dir = os.path.join(cache_dir, "ancestry_user")
    os.makedirs(out_dir, exist_ok=True)

    reference_rsids = _read_bim_rsids(f"{ref_prefix}.bim")
    extract_path = os.path.join(out_dir, "reference_snps.txt")
    with open(extract_path, "w", encoding="utf-8") as f:
        for rsid in sorted(reference_rsids):
            f.write(f"{rsid}\n")

    user_vcf = os.path.join(out_dir, "user.vcf")
    written = _write_genome_vcf(genome, user_vcf, allowed_rsids=reference_rsids)
    if written == 0:
        return {"skipped": True, "reason": "no user SNPs overlap the reference panel"}

    user_prefix = os.path.join(out_dir, "user")
    _run_command([
        plink2, "--vcf", user_vcf, "--extract", extract_path,
        "--double-id", "--max-alleles", "2", "--make-bed", "--out", user_prefix,
    ], logger)

    user_bim = f"{user_prefix}.bim"
    if not os.path.exists(user_bim) or _count_bim_snps(user_bim) == 0:
        return {"skipped": True, "reason": "no SNPs remained after PLINK conversion"}

    merge_list = os.path.join(out_dir, "pmerge_list.txt")
    with open(merge_list, "w", encoding="utf-8") as f:
        f.write(f"{ref_prefix}\n")
        f.write(f"{user_prefix}\n")

    merged_prefix = os.path.join(out_dir, "merged")
    _run_command([
        plink2, "--pmerge-list", merge_list, "bfile",
        "--make-bed", "--out", merged_prefix,
    ], logger)

    pop_order, unknown = _write_supervised_pop(
        ref_pop, f"{ref_prefix}.fam", f"{merged_prefix}.fam", f"{merged_prefix}.pop"
    )
    if unknown == 0:
        logger.warning("  No unknown sample found in merged .fam; ADMIXTURE will still run")
    k = len(pop_order)
    if k < 2:
        return {"skipped": True, "reason": "reference panel has fewer than two population labels"}

    _run_command([admixture, "--supervised", "merged.bed", str(k)], logger, cwd=out_dir)

    q_path = os.path.join(out_dir, f"merged.{k}.Q")
    if not os.path.exists(q_path):
        raise FileNotFoundError(f"ADMIXTURE output not found: {q_path}")

    with open(q_path, "r", encoding="utf-8") as f:
        rows = [line.split() for line in f if line.strip()]
    if not rows:
        raise ValueError("ADMIXTURE .Q output is empty")

    q = [float(x) for x in rows[-1]]
    proportions = {
        pop_order[i]: q[i]
        for i in range(min(len(pop_order), len(q)))
    }
    return {"proportions": proportions, "K": k, "n_snps_used": _count_bim_snps(f"{merged_prefix}.bim")}


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


def _download_first(urls: Sequence[str], dest: str,
                    logger: logging.Logger) -> Optional[str]:
    if os.path.exists(dest):
        return dest
    for url in urls:
        try:
            return download_file(url, dest, logger, max_retries=1)
        except Exception as exc:
            logger.warning(f"  Could not download {url}: {exc}")
    return None
def _calculator_urls(name: str, ext: str) -> List[str]:
    if name == "Eurogenes_K13":
        steven_ext = "13.F" if ext == "par" else "alleles"
        return [f"https://raw.githubusercontent.com/stevenliuyi/admix/master/admix/data/K13.{steven_ext}", f"https://raw.githubusercontent.com/wegene-llc/admix/master/calculator/{name}.{ext}", f"https://raw.githubusercontent.com/stevemkopp/dna-tools/master/{name}.{ext}"]
    return [
        f"https://raw.githubusercontent.com/wegene-llc/admix/master/calculator/{name}.{ext}",
        f"https://raw.githubusercontent.com/stevemkopp/dna-tools/master/{name}.{ext}",
    ]
def _parse_alleles_file(path: str) -> List[Tuple[str, str]]:
    records = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"[\s,]+", line.strip())
            rsid = next((p for p in parts if p.startswith("rs")), None)
            allele = next((p.upper() for p in reversed(parts) if p.upper() in {"A", "C", "G", "T"}), None)
            if rsid and allele:
                records.append((rsid, allele))
    return records


def _parse_par_file(path: str, alleles: List[Tuple[str, str]],
                    default_labels: List[str]) -> Tuple[List[str], List[Tuple[str, str, List[float]]]]:
    allele_by_rsid = {rsid: allele for rsid, allele in alleles}
    labels = list(default_labels)
    rows = []
    allele_idx = 0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = re.split(r"[\s,]+", line.strip())
            rsid = next((p for p in parts if p.startswith("rs")), None)
            if rsid:
                rsid_index = parts.index(rsid)
                freqs = [float(p) for p in parts[rsid_index + 1:] if _is_float(p) and 0.0 <= float(p) <= 1.0]
                allele = allele_by_rsid.get(rsid)
            else:
                freqs = [float(p) for p in parts if _is_float(p) and 0.0 <= float(p) <= 1.0]
                if not freqs:
                    header = [p for p in parts if p.lower() not in {"snp", "rsid", "marker", "allele"}]
                    if len(header) > 1:
                        labels = header
                    continue
                if allele_idx >= len(alleles):
                    continue
                rsid, allele = alleles[allele_idx]
                allele_idx += 1
            if not allele or not freqs:
                continue
            if len(freqs) != len(labels):
                labels = [f"Component_{i + 1}" for i in range(len(freqs))]
            rows.append((rsid, allele, [max(0.0, min(1.0, freq)) for freq in freqs]))
    return labels, rows


def _fallback_k7_panel() -> Tuple[List[str], List[Tuple[str, str, List[float]]]]:
    """
    Offline fallback K=7 demonstration panel. Full Eurogenes/MDLP .par files
    lack a stable canonical URL, so this projects built-in 1000G AIM frequencies
    into broad 1000G/HGDP-like components. It is not a full GEDmatch substitute.
    """
    labels = "Northwest_European Mediterranean Caucasus_Near_East Baltic_Slavic East_Asian South_Asian Sub_Saharan_African".split()
    rows = []
    for rsid, freqs in AIMS_PANEL.items():
        eur, eas, afr = freqs["EUR"], freqs["EAS"], freqs["AFR"]
        sas = freqs["SAS"]
        component_freqs = [eur, 0.75 * eur + 0.25 * sas, 0.50 * eur + 0.50 * sas, 0.90 * eur + 0.05 * sas + 0.05 * eas, eas, sas, afr]
        rows.append((rsid, freqs["allele"], [max(0.0, min(1.0, x)) for x in component_freqs]))
    return labels, rows


def _load_calculator(name: str, calc_root: str,
                     logger: logging.Logger) -> Optional[Tuple[List[str], List[Tuple[str, str, List[float]]]]]:
    calc_dir = os.path.join(calc_root, name)
    os.makedirs(calc_dir, exist_ok=True)
    par_path = os.path.join(calc_dir, f"{name}.par")
    alleles_path = os.path.join(calc_dir, f"{name}.alleles")

    par = _download_first(_calculator_urls(name, "par"), par_path, logger)
    alleles = _download_first(_calculator_urls(name, "alleles"), alleles_path, logger)
    if not par or not alleles:
        return None

    allele_records = _parse_alleles_file(alleles)
    labels, rows = _parse_par_file(par, allele_records, K_CALCULATOR_LABELS[name])
    return (labels, rows) if rows else None


def _dosage_for_allele(genome: GenomeData, rsid: str, allele: str) -> Optional[int]:
    key = genome.rsid_index.get(rsid)
    if not key:
        return None
    variant = genome.variants[key]
    allele = allele.upper()
    observed = {variant.ref.upper(), variant.alt.upper()}
    if allele not in observed:
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}.get(allele)
        if complement in observed:
            allele = complement
        else:
            return None
    a1, a2 = variant.genotype_alleles
    return (1 if a1.upper() == allele else 0) + (1 if a2.upper() == allele else 0)


def _solve_components(np, scipy_minimize, f_matrix, dosage_vector,
                      logger: logging.Logger) -> List[float]:
    k = f_matrix.shape[1]
    x0 = np.full(k, 1.0 / k)
    if scipy_minimize is not None:
        result = scipy_minimize(
            lambda q: float(np.sum((dosage_vector - f_matrix.dot(q)) ** 2)), x0,
            method="SLSQP", bounds=[(0.0, 1.0)] * k,
            constraints={"type": "eq", "fun": lambda q: float(np.sum(q) - 1.0)}, options={"maxiter": 500, "ftol": 1e-10},
        )
        if result.success:
            q = np.maximum(result.x, 0.0)
        else:
            logger.warning(f"  scipy SLSQP failed ({result.message}); using clipped least squares")
            q = np.maximum(np.linalg.lstsq(f_matrix, dosage_vector, rcond=None)[0], 0.0)
    else:
        logger.warning("  scipy not available; using numpy clipped least-squares approximation")
        q = np.maximum(np.linalg.lstsq(f_matrix, dosage_vector, rcond=None)[0], 0.0)
    total = float(np.sum(q))
    if total <= 0:
        q = x0
    else:
        q = q / total
    return [float(x) for x in q]


def _score_calculator(labels: List[str], rows: List[Tuple[str, str, List[float]]],
                      genome: GenomeData, np, scipy_minimize,
                      logger: logging.Logger) -> Dict:
    matrix_rows = []
    dosages = []
    for rsid, allele, freqs in rows:
        dosage = _dosage_for_allele(genome, rsid, allele)
        if dosage is None:
            continue
        matrix_rows.append(freqs)
        dosages.append(dosage / 2.0)

    if not dosages:
        return {"skipped": True, "reason": "no overlapping calculator SNPs", "n_snps_used": 0}

    f_matrix = np.array(matrix_rows, dtype=float)
    dosage_vector = np.array(dosages, dtype=float)
    q = _solve_components(np, scipy_minimize, f_matrix, dosage_vector, logger)
    return {
        "components": {labels[i]: q[i] for i in range(min(len(labels), len(q)))},
        "n_snps_used": len(dosages),
    }


def _run_k_calculators(genome: GenomeData, cache_dir: str, logger: logging.Logger) -> Dict:
    try:
        import numpy as np
    except ImportError:
        return {"skipped": True, "reason": "numpy not installed; K calculators require numpy"}

    try:
        from scipy.optimize import minimize as scipy_minimize
    except ImportError:
        scipy_minimize = None

    calc_root = os.path.join(cache_dir, "k_calculators")
    os.makedirs(calc_root, exist_ok=True)

    calculators = {}
    for name in K_CALCULATOR_LABELS:
        loaded = _load_calculator(name, calc_root, logger)
        if not loaded:
            calculators[name] = {"skipped": True, "reason": "calculator .par/.alleles files unavailable"}
            continue
        labels, rows = loaded
        calculators[name] = _score_calculator(labels, rows, genome, np, scipy_minimize, logger)

    if all(value.get("skipped") for value in calculators.values()):
        labels, rows = _fallback_k7_panel()
        calculators["Fallback_K7_AIMs"] = _score_calculator(
            labels, rows, genome, np, scipy_minimize, logger
        )

    return {"calculators": calculators}


def _annotation_for_haplogroup(haplogroup: str, annotations: Dict[str, Dict]) -> Dict:
    if not haplogroup:
        return {}
    normalized = haplogroup.upper()
    for prefix in sorted(annotations, key=len, reverse=True):
        if normalized.startswith(prefix):
            return annotations[prefix]
    return {}


def _parse_yhaplo_output(out_dir: str) -> Optional[str]:
    pattern = re.compile(r"\b(E1b1b|R1a|R1b|[A-Z][0-9][A-Za-z0-9._-]*)\b", re.IGNORECASE)
    for root, _, files in os.walk(out_dir):
        for filename in files:
            path = os.path.join(root, filename)
            try:
                with open(path, "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if "haplogroup" not in line.lower() and "\t" not in line:
                            continue
                        matches = pattern.findall(line)
                        if matches:
                            return matches[-1]
            except OSError:
                continue
    return None


def _parse_haplogrep_output(path: str) -> Tuple[Optional[str], Optional[float]]:
    if not os.path.exists(path):
        return None, None
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f if line.strip()]
    if len(lines) < 2:
        return None, None
    header = re.split(r"[\t,;]", lines[0])
    row = re.split(r"[\t,;]", lines[1])
    lower = [h.lower() for h in header]
    hg_idx = next((i for i, h in enumerate(lower) if "haplogroup" in h or h == "hg"), None)
    q_idx = next((i for i, h in enumerate(lower) if "quality" in h or h == "quality_score"), None)
    haplogroup = row[hg_idx].strip() if hg_idx is not None and hg_idx < len(row) else None
    quality = None
    if q_idx is not None and q_idx < len(row):
        try:
            quality = float(row[q_idx])
        except ValueError:
            quality = None
    return haplogroup, quality


def _run_haplogroups(genome: GenomeData, cache_dir: str, logger: logging.Logger) -> Dict:
    del cache_dir
    results = {}
    with tempfile.TemporaryDirectory(prefix="genome_haplogroups_") as tmp:
        y_vcf = os.path.join(tmp, "y.vcf")
        y_count = _write_genome_vcf(genome, y_vcf, chroms=["Y"])
        if y_count == 0:
            results["Y"] = {"skipped": True, "reason": "no Y-chromosome data (likely female)"}
        else:
            yhaplo = shutil.which("yhaplo")
            if not yhaplo:
                results["Y"] = {"skipped": True, "reason": "yhaplo not installed"}
            else:
                y_out = os.path.join(tmp, "yhaplo")
                os.makedirs(y_out, exist_ok=True)
                _run_command([yhaplo, "-i", y_vcf, "-o", y_out], logger)
                haplogroup = _parse_yhaplo_output(y_out)
                results["Y"] = {
                    "haplogroup": haplogroup,
                    "annotation": _annotation_for_haplogroup(haplogroup or "", Y_HAPLOGROUP_ANNOTATIONS),
                    "n_variants_used": y_count,
                }

        mt_vcf = os.path.join(tmp, "mt.vcf")
        mt_count = _write_genome_vcf(genome, mt_vcf, chroms=["MT"])
        if mt_count == 0:
            results["mt"] = {"skipped": True, "reason": "no mitochondrial data"}
        else:
            haplogrep = shutil.which("haplogrep")
            if not haplogrep:
                results["mt"] = {"skipped": True, "reason": "haplogrep not installed"}
            else:
                mt_out = os.path.join(tmp, "haplogrep.txt")
                _run_command([
                    haplogrep, "classify", "--in", mt_vcf,
                    "--format", "vcf", "--out", mt_out,
                ], logger)
                haplogroup, quality = _parse_haplogrep_output(mt_out)
                results["mt"] = {
                    "haplogroup": haplogroup,
                    "quality": quality,
                    "annotation": _annotation_for_haplogroup(haplogroup or "", MT_HAPLOGROUP_ANNOTATIONS),
                    "n_variants_used": mt_count,
                }
    return results


def _run_subsection(label: str, func, logger: logging.Logger) -> Dict:
    logger.info("=" * 60)
    logger.info(label)
    logger.info("=" * 60)
    try:
        return func()
    except Exception as exc:
        logger.error(f"  {label} failed: {exc}")
        return {"error": str(exc)}


def run_layer4_ancestry(genome: GenomeData, cache_dir: str,
                        logger: logging.Logger) -> Dict:
    """Run ancestry inference. Returns results dict."""
    logger.info("=" * 60)
    logger.info("LAYER 4: Ancestry Inference")
    logger.info("=" * 60)

    results = {
        "layer": 4,
        "aims": _run_subsection("Layer 4 / AIMs", lambda: _run_aims(genome, logger), logger),
        "admixture": _run_subsection("Layer 4 / ADMIXTURE", lambda: _run_admixture(genome, cache_dir, logger), logger),
        "k_calculators": _run_subsection("Layer 4 / K Calculators", lambda: _run_k_calculators(genome, cache_dir, logger), logger),
        "haplogroups": _run_subsection("Layer 4 / Haplogroups", lambda: _run_haplogroups(genome, cache_dir, logger), logger),
    }

    aims_found = results.get("aims", {}).get("aims_found", 0)
    logger.info(f"Layer 4 complete: {aims_found} AIMs found")
    return results
