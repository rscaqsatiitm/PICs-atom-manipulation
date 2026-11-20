# Run every simulation module in dependency order. MEEP is optional (grating FDTD).
PY ?= python3
MEEPPY ?= $(shell command -v meep-python 2>/dev/null)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

.PHONY: all test check check-numbers m1 m2 m3 m3g m4 m5 m7

all: m1 m2 m3g m3 m4 m5 m7

test:
	$(PY) -m pytest -q

# every number in report/sections must occur in docs, README or results (exit 1 otherwise)
check-numbers:
	$(PY) report/check_numbers.py

check: test check-numbers

m1:
	$(PY) scripts/m1_platform.py
m2:
	$(PY) scripts/m2_couplers_mmi.py
m3g:
	@if [ -n "$(MEEPPY)" ]; then PYTHONPATH=. $(MEEPPY) scripts/m3_grating_fdtd.py; \
	else echo "MEEP interpreter not set (MEEPPY=...); skipping grating FDTD"; fi
m3:
	$(PY) scripts/m3_opa.py
m4:
	$(PY) scripts/m4_switches.py
m5:
	$(PY) scripts/m5_heaters.py
m7:
	$(PY) scripts/m7_budget.py

.PHONY: report
report: check-numbers
	cd report && $(PY) build_docs_volume.py && $(PY) make_appendices.py && \
	latexmk -pdf -interaction=nonstopmode report.tex && latexmk -pdf -interaction=nonstopmode docs_volume.tex
