import logging
import math
from typing import Dict, List

from core import GenomeData

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
