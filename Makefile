# kmerstash — laptop-portable DNA k-mer presence screening
#
#   make            build the release binary
#   make demo       gen a synthetic dataset + screen it (see the ALERT table)
#   make bench      print the three experiment CSVs (throughput/error/separation)
#   make notebook   regenerate kmerstash_demo.ipynb
#   make run-notebook  execute the notebook end-to-end (nbconvert)
#   make test       cargo unit tests
#   make clean

BIN  = target/release/kmerstash
SRC  = src/main.rs src/kmer.rs Cargo.toml

# conda python carries nbformat / pandas / matplotlib (same env pdf_search uses).
CONDA_PYTHON ?= /var/home/fraser/anaconda3/envs/py314/bin/python3
PYTHON3 := $(shell command -v $(CONDA_PYTHON) 2>/dev/null || command -v python3)

.PHONY: all build demo bench notebook run-notebook test clean

all: build

build: $(BIN)
$(BIN): $(SRC)
	cargo build --release

demo: build
	./$(BIN) gen --out data
	@echo
	./$(BIN) screen --ref data/panel.fa --sample data/sample.fa

bench: build
	@echo "### throughput";  ./$(BIN) bench --kind throughput
	@echo "### error";       ./$(BIN) bench --kind error
	@echo "### separation (first rows)"; ./$(BIN) bench --kind separation | head

notebook: build
	$(PYTHON3) build_notebook.py

run-notebook: notebook
	$(PYTHON3) -m jupyter nbconvert --to notebook --execute --inplace \
		--ExecutePreprocessor.timeout=600 kmerstash_demo.ipynb
	@echo "executed kmerstash_demo.ipynb"

test:
	cargo test --release

clean:
	cargo clean
	rm -rf data
