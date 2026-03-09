#!/usr/bin/env python3
"""
download_pgs.py — Скачивает PGS Catalog scoring files через pgscatalog-core.

Зависимости:
    pip install pgscatalog-utils

Использование:
    python3 download_pgs.py                          # скачать все рекомендованные
    python3 download_pgs.py -o ./my_pgs_dir          # в свою папку
    python3 download_pgs.py -i PGS000018 PGS000014   # только конкретные
    python3 download_pgs.py --list-traits             # показать что будет скачано
    python3 download_pgs.py --overwrite               # перезаписать существующие

После скачивания запусти:
    python3 genome_analysis.py --vcf yyq9p5.vcf --pgs-dir ./pgs_scoring_files \\
        --output-dir ./genome_report --cache-dir ./cache
"""

import argparse
import os
import socket
import sys

from pgscatalog.core import GenomeBuild, ScoringFile

# --- Рекомендованные PGS ID ---
RECOMMENDED_PGS = [
    # ═══ КАРДИОМЕТАБОЛИЧЕСКИЕ ═══
    ("PGS000018", "Coronary Artery Disease (metaGRS, 1.7M var)"),
    ("PGS000013", "Coronary Artery Disease (Khera, 6.6M var)"),
    ("PGS000296", "Coronary Artery Disease (LDpred)"),
    ("PGS000016", "Atrial Fibrillation"),
    ("PGS000039", "Stroke"),
    ("PGS002997", "Hypertension (7.6M var)"),
    ("PGS000014", "Type 2 Diabetes"),
    ("PGS000807", "Type 2 Diabetes (newer)"),
    ("PGS000033", "T2D (insulin resistance SNPs)"),
    ("PGS000034", "Systolic Blood Pressure"),
    # ═══ НЕЙРОДЕГЕНЕРАЦИЯ (семейная история — болезнь Пика / FTD) ═══
    ("PGS003992", "Alzheimer's Disease (1.1M var, лучший)"),
    ("PGS000812", "Alzheimer's Disease (базовый)"),
    ("PGS000903", "Parkinson's Disease (1.8K var)"),
    ("PGS002760", "Epilepsy (836K var)"),
    ("PGS000137", "Glaucoma (2.7K var)"),
    ("PGS001797", "Primary Open-Angle Glaucoma (885K var)"),
    # ═══ МЕНТАЛЬНОЕ ЗДОРОВЬЕ ═══
    ("PGS003497", "Depression (980K var)"),
    ("PGS002786", "Bipolar Disorder (949K var)"),
    ("PGS004451", "Anxiety Disorders (1.06M var)"),
    ("PGS002746", "ADHD (514K var)"),
    ("PGS002790", "Autism Spectrum Disorder (917K var)"),
    ("PGS000135", "Schizophrenia"),
    ("PGS005211", "Alcohol Use Disorder (2.3M var)"),
    # ═══ БИОХАКИНГ / МЕТАБОЛИЗМ ═══
    ("PGS003037", "LDL Cholesterol (8.9M var, улучшенный)"),
    ("PGS000115", "LDL Cholesterol (базовый)"),
    ("PGS003147", "Triglycerides (8.9M var)"),
    ("PGS002937", "Glucose (7.6M var)"),
    ("PGS005276", "Fasting Insulin adj. BMI (1.1M var)"),
    ("PGS000027", "Body Mass Index (BMI)"),
    ("PGS002356", "Waist-Hip Ratio (1.1M var)"),
    ("PGS002813", "Visceral Adipose Tissue (1.1M var)"),
    ("PGS000882", "Vitamin D — 25(OH)D (1.09M var)"),
    ("PGS003154", "Vitamin B12 (315K var)"),
    ("PGS003426", "Ferritin / Iron (1.08M var)"),
    ("PGS004332", "Calcium (1.06M var)"),
    ("PGS002768", "Osteoporosis (1.09M var)"),
    ("PGS002231", "Educational Attainment / IQ proxy (951K var)"),
    ("PGS005222", "Frailty Index (1.7M var)"),
    # ═══ ОНКОЛОГИЯ ═══
    ("PGS001229", "Prostate Cancer"),
    ("PGS000037", "Breast Cancer"),
    ("PGS004834", "Skin Cancer / Melanoma (1.3M var)"),
    # ═══ ИММУНОЛОГИЯ / ВОСПАЛЕНИЕ ═══
    ("PGS000017", "Inflammatory Bowel Disease (6.9M var)"),
    ("PGS002344", "Psoriasis (1.1M var)"),
    ("PGS000728", "Chronic Kidney Disease (2.0M var)"),
]

BUILD = GenomeBuild.GRCh37
EBI_FTP_HOST = "ftp.ebi.ac.uk"


def check_ebi_connectivity() -> bool:
    """Быстрая проверка доступности EBI FTP сервера."""
    try:
        sock = socket.create_connection((EBI_FTP_HOST, 443), timeout=5)
        sock.close()
        return True
    except (socket.timeout, OSError):
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Download PGS Catalog scoring files for genome_analysis.py"
    )
    parser.add_argument(
        "-o", "--output-dir",
        default="./pgs_scoring_files",
        help="Output directory (default: ./pgs_scoring_files)",
    )
    parser.add_argument(
        "-i", "--ids",
        nargs="+",
        default=None,
        help="Specific PGS IDs to download (default: all recommended)",
    )
    parser.add_argument(
        "--list-traits",
        action="store_true",
        help="Show recommended PGS IDs and exit",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files (re-download if newer version available)",
    )
    parser.add_argument(
        "--no-check",
        action="store_true",
        help="Skip EBI FTP connectivity check (try downloading anyway)",
    )
    args = parser.parse_args()

    if args.list_traits:
        print("Recommended PGS scoring files:")
        print(f"{'PGS ID':<14} {'Description'}")
        print("-" * 65)
        for pgs_id, desc in RECOMMENDED_PGS:
            print(f"{pgs_id:<14} {desc}")
        print(f"\nTotal: {len(RECOMMENDED_PGS)} files")
        print(f"Download all: python3 {sys.argv[0]}")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    # Определяем что качать
    if args.ids:
        to_download = [(pid, "") for pid in args.ids]
    else:
        to_download = RECOMMENDED_PGS

    print("=" * 60)
    print("  PGS Catalog Scoring File Downloader")
    print(f"  Files to download: {len(to_download)}")
    print(f"  Output: {os.path.abspath(args.output_dir)}")
    print(f"  Build: {BUILD.value}")
    print("=" * 60)
    print()

    # Проверяем доступность EBI FTP до начала скачивания
    if not args.no_check:
        print("Checking EBI FTP connectivity...", end=" ", flush=True)
    if not args.no_check and not check_ebi_connectivity():
        print("UNREACHABLE")
        print()
        print("  ftp.ebi.ac.uk недоступен с этой сети.")
        print("  Все файлы PGS Catalog хостятся на этом сервере.")
        print()
        print("  Варианты решения:")
        print("    1. Используй VPN (подключись к серверу в EU/UK/US)")
        print("    2. Скачай файлы через браузер на другой сети:")
        for pgs_id, desc in to_download:
            expected = os.path.join(args.output_dir, f"{pgs_id}_hmPOS_GRCh37.txt.gz")
            if not os.path.exists(expected):
                print(f"       https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/ScoringFiles/Harmonized/{pgs_id}_hmPOS_GRCh37.txt.gz")
        print(f"    3. Положи скачанные .txt.gz файлы в {args.output_dir}/")
        print()
        print("  После скачивания запусти этот скрипт снова — он найдёт файлы и пропустит их.")
        sys.exit(1)
    elif not args.no_check:
        print("OK")
        print()

    ok, fail, skip = 0, 0, 0

    for pgs_id, desc in to_download:
        expected_file = os.path.join(
            args.output_dir, f"{pgs_id}_hmPOS_GRCh37.txt.gz"
        )

        # Пропускаем если уже скачан (и не --overwrite)
        if not args.overwrite and os.path.exists(expected_file) and os.path.getsize(expected_file) > 100:
            size_kb = os.path.getsize(expected_file) // 1024
            print(f"✓ {pgs_id} — already downloaded ({size_kb} KB)")
            skip += 1
            continue

        label = f"{pgs_id} — {desc}" if desc else pgs_id
        print(f"▶ {label}")

        try:
            sf = ScoringFile(pgs_id, target_build=BUILD)
            print(f"  Downloading from PGS Catalog...")
            sf.download(args.output_dir, overwrite=args.overwrite)
            size_kb = os.path.getsize(expected_file) // 1024
            print(f"  ✓ OK ({size_kb} KB)")
            ok += 1
        except Exception as e:
            print(f"  ✗ FAILED: {e}")
            fail += 1

        print()

    print("=" * 60)
    print(f"  Done! Downloaded: {ok} | Skipped: {skip} | Failed: {fail}")
    print(f"  Files: {os.path.abspath(args.output_dir)}/")
    print("=" * 60)

    if ok + skip > 0:
        print()
        print("Запуск анализа:")
        print(f"  python3 genome_analysis.py \\")
        print(f"    --vcf yyq9p5.vcf \\")
        print(f"    --pgs-dir {args.output_dir} \\")
        print(f"    --output-dir ./genome_report \\")
        print(f"    --cache-dir ./cache")

    if fail > 0:
        print()
        print("Для файлов которые не скачались:")
        print("  1. Открой https://www.pgscatalog.org/score/<PGS_ID>/")
        print("  2. Скачай Harmonized GRCh37 файл через браузер")
        print(f"  3. Положи в {args.output_dir}/")


if __name__ == "__main__":
    main()
