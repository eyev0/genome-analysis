# Инструкция по скачиванию PGS Catalog файлов

## Что это

PGS Catalog — открытая база полигенных скоров. Наш скрипт `genome_analysis.py` умеет читать
их гармонизированные scoring-файлы (формат GRCh37) и считать скор по твоему VCF.

Файлы весят от 100 КБ до ~200 МБ (для полногеномных скоров с миллионами вариантов).
Чем больше вариантов в файле — тем точнее скор, но из 632K SNP в твоём VCF совпадёт
типично 5-30% от полногеномного файла. Это нормально — направление скора останется информативным.

## Как скачать

### Шаблон URL

```
https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{PGS_ID}/ScoringFiles/Harmonized/{PGS_ID}_hmPOS_GRCh37.txt.gz
```

Замени `{PGS_ID}` на нужный идентификатор (например, `PGS000018`).

### Команда для скачивания одного файла

```bash
wget https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000018/ScoringFiles/Harmonized/PGS000018_hmPOS_GRCh37.txt.gz
```

### Скрипт для скачивания всех рекомендованных

```bash
#!/bin/bash
# download_pgs.sh — скачивает рекомендованные PGS scoring files

DEST_DIR="./pgs_scoring_files"
mkdir -p "$DEST_DIR"

PGS_IDS=(
  PGS000018   # CAD (metaGRS, ~1.7M variants)
  PGS000013   # CAD (Khera 2018, 6.6M variants)
  PGS000296   # CAD (LDpred)
  PGS000014   # Type 2 Diabetes
  PGS000807   # Type 2 Diabetes (newer)
  PGS000027   # BMI
  PGS000812   # Alzheimer's Disease
  PGS000016   # Atrial Fibrillation
  PGS000035   # Atrial Fibrillation (alternative)
  PGS000115   # LDL Cholesterol
  PGS002274   # LDL Cholesterol (newer)
  PGS000135   # Schizophrenia
  PGS000338   # Atrial Fibrillation (Kany 2023)
  PGS000037   # Breast Cancer
  PGS001229   # Prostate Cancer
  PGS000034   # Systolic Blood Pressure
  PGS000039   # Stroke
  PGS001327   # Diabetes (general)
  PGS000033   # T2D (insulin resistance SNPs)
)

BASE_URL="https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores"

for PGS_ID in "${PGS_IDS[@]}"; do
  # Skip comment lines
  [[ "$PGS_ID" == \#* ]] && continue

  URL="${BASE_URL}/${PGS_ID}/ScoringFiles/Harmonized/${PGS_ID}_hmPOS_GRCh37.txt.gz"
  OUTFILE="${DEST_DIR}/${PGS_ID}_hmPOS_GRCh37.txt.gz"

  if [ -f "$OUTFILE" ]; then
    echo "SKIP: $PGS_ID (already exists)"
  else
    echo "Downloading $PGS_ID..."
    wget -q --show-progress -O "$OUTFILE" "$URL"
    if [ $? -ne 0 ]; then
      echo "FAIL: $PGS_ID — URL may be incorrect, check pgscatalog.org"
      rm -f "$OUTFILE"
    fi
  fi
done

echo ""
echo "Done! Files saved to $DEST_DIR"
echo "Run: python3 genome_analysis.py --vcf your.vcf --pgs-dir $DEST_DIR --output-dir ./results --cache-dir ./cache"
```

## Рекомендованные PGS файлы

### Приоритет 1: Клинически важные (скачай эти в первую очередь)

| PGS ID | Условие | Комментарий |
|--------|---------|-------------|
| **PGS000018** | Ишемическая болезнь сердца (CAD) | metaGRS, ~1.7M вариантов, хорошо валидирован в EUR |
| **PGS000013** | CAD | Khera 2018, ~6.6M вариантов, один из самых цитируемых |
| **PGS000014** | Сахарный диабет 2 типа | Базовый T2D скор |
| **PGS000807** | Сахарный диабет 2 типа | Более новая версия |
| **PGS000027** | Индекс массы тела (BMI) | Полногеномный BMI скор |
| **PGS000812** | Болезнь Альцгеймера | Основной AD скор |
| **PGS000016** | Фибрилляция предсердий | Базовый AF скор |

### Приоритет 2: Дополнительные клинические

| PGS ID | Условие | Комментарий |
|--------|---------|-------------|
| **PGS000115** | LDL холестерин | Генетически предсказанный уровень LDL |
| **PGS002274** | LDL холестерин | Более новый, возможно больше вариантов |
| **PGS000135** | Шизофрения | PGC-based |
| **PGS000034** | Систолическое АД | Артериальное давление |
| **PGS001229** | Рак простаты | Для мужчин |
| **PGS000037** | Рак молочной железы | Для женщин в семье |

### Приоритет 3: Расширенные

| PGS ID | Условие | Комментарий |
|--------|---------|-------------|
| **PGS000296** | CAD (LDpred) | Альтернативный метод скоринга |
| **PGS000035** | Фибрилляция предсердий | Альтернативный AF скор |
| **PGS000338** | Фибрилляция предсердий | Kany 2023, новее |
| **PGS000039** | Инсульт | Ишемический инсульт |
| **PGS001327** | Диабет (общий) | Объединённый скор |
| **PGS000033** | T2D (инсулинорезистентность) | Подтип T2D по механизму |

## Как запустить после скачивания

```bash
# 1. Скачай файлы (или запусти download_pgs.sh)
# 2. Положи .txt.gz файлы в папку pgs_scoring_files/
# 3. Запусти:

python3 genome_analysis.py \
  --vcf yyq9p5.vcf \
  --pgs-dir ./pgs_scoring_files \
  --output-dir ./genome_report \
  --cache-dir ./cache

# Или только Layer 1:
python3 genome_analysis.py \
  --vcf yyq9p5.vcf \
  --pgs-dir ./pgs_scoring_files \
  --output-dir ./genome_report \
  --cache-dir ./cache \
  --layer 1
```

## Что ожидать

- Скрипт автоматически найдёт все `.txt.gz` файлы в `--pgs-dir`
- Для каждого файла посчитает скор и покажет coverage (% совпавших вариантов)
- При 632K SNP из Genotek ожидай ~5-30% coverage для полногеномных скоров
- Это нормально — даже при 10% coverage направление скора информативно
- Для полного coverage нужна импутация (Impute.me или Michigan Imputation Server)

## Если файл не скачивается

Некоторые PGS ID могут не иметь гармонизированных файлов для GRCh37.
В таком случае:

1. Открой `https://www.pgscatalog.org/score/{PGS_ID}/` в браузере
2. Проверь раздел "Downloads" → "Harmonized Scoring Files"
3. Если GRCh37 нет — попробуй GRCh38 (но наш VCF в hg19, могут не совпасть позиции)
4. Или поищи альтернативный PGS ID для того же заболевания на pgscatalog.org

## Дополнительно: поиск других скоров

На сайте https://www.pgscatalog.org/browse/traits/ можно найти PGS для любого заболевания.
Фильтруй по:
- Ancestry: European (для лучшей применимости)
- Number of variants: чем больше, тем точнее (но и дольше скачивать)
- Development method: LDpred2, PRS-CS, metaGRS — современные методы лучше простого clumping+thresholding
