.PHONY: help install install-python install-tools install-haplogrep check clean-cache pgs references analyze

CACHE_DIR ?= ./cache
PGS_DIR ?= ./pgs_scoring_files
OUTPUT_DIR ?= ./reports
VCF ?=
HAPLOGREP_VERSION ?= 3.2.2
HAPLOGREP_DIR ?= $(HOME)/.local/bin
HAPLOGREP_JAR ?= $(HAPLOGREP_DIR)/haplogrep.jar

help:
	@echo "Targets:"
	@echo "  install            Install Python deps + system tools + haplogrep"
	@echo "  install-python     pip install -r requirements.txt"
	@echo "  install-tools      brew install plink2 bcftools admixture eigensoft"
	@echo "  install-haplogrep  Download haplogrep3 JAR + wrapper to $(HAPLOGREP_DIR)"
	@echo "  check              Verify all binaries on PATH"
	@echo "  pgs                Download PGS scoring files"
	@echo "  references         Build reference panel for ADMIXTURE (slow, GBs)"
	@echo "  analyze VCF=path   Run full pipeline on a VCF"
	@echo "  clean-cache        Remove cache directory"

install: install-python install-tools install-haplogrep check

install-python:
	pip install -r requirements.txt

install-tools:
	@command -v brew >/dev/null 2>&1 || { echo "Homebrew required. Install from https://brew.sh"; exit 1; }
	brew install plink2 bcftools admixture eigensoft openjdk

install-haplogrep:
	@mkdir -p $(HAPLOGREP_DIR)
	@if [ ! -f $(HAPLOGREP_JAR) ]; then \
		echo "Downloading haplogrep $(HAPLOGREP_VERSION)..."; \
		curl -L -o $(HAPLOGREP_JAR) \
			https://github.com/genepi/haplogrep3/releases/download/v$(HAPLOGREP_VERSION)/haplogrep3-$(HAPLOGREP_VERSION).jar; \
	fi
	@printf '#!/bin/sh\nexec java -jar %s "$$@"\n' "$(HAPLOGREP_JAR)" > $(HAPLOGREP_DIR)/haplogrep
	@chmod +x $(HAPLOGREP_DIR)/haplogrep
	@echo "haplogrep installed at $(HAPLOGREP_DIR)/haplogrep"
	@echo "Ensure $(HAPLOGREP_DIR) is on PATH"

check:
	@echo "Python deps:"
	@python -c "import requests, numpy, scipy, yhaplo; print('  OK')" || echo "  MISSING (run: make install-python)"
	@echo "System tools:"
	@for t in plink2 bcftools admixture convertf haplogrep java; do \
		if command -v $$t >/dev/null 2>&1; then \
			echo "  $$t: $$(command -v $$t)"; \
		else \
			echo "  $$t: MISSING"; \
		fi; \
	done

pgs:
	python download_pgs.py

references:
	python download_references.py --cache-dir $(CACHE_DIR)

analyze:
	@test -n "$(VCF)" || { echo "Usage: make analyze VCF=path/to.vcf"; exit 1; }
	python genome_analysis.py --vcf $(VCF) --pgs-dir $(PGS_DIR) \
		--output-dir $(OUTPUT_DIR) --cache-dir $(CACHE_DIR)

clean-cache:
	rm -rf $(CACHE_DIR)
