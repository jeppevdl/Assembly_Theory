# Makefile for KEGG data processing
JULIA:=julia --project=.

# Directories
DATA_DIR := data
ADJMAT_DIR := $(DATA_DIR)/adjmat
LOOKUP_DIR := $(DATA_DIR)/lookup
METADATA_DIR := $(DATA_DIR)/metadata

# Files
KEGG_SMALL := $(DATA_DIR)/kegg-small.tsv
KEGG_CPD_JSON := $(METADATA_DIR)/kegg_cpd.latest.json
KEGG_RN_JSON := $(METADATA_DIR)/kegg_rn.latest.json

# Phony targets
.PHONY: prepare clean setup

# Setup project
setup:
	$(JULIA) -e "using Pkg; Pkg.instantiate()"

# Default target
prepare: setup directories $(KEGG_SMALL) $(KEGG_CPD_JSON) $(KEGG_RN_JSON) reconstruct

# Create directory structure
directories:
	mkdir -p $(ADJMAT_DIR) $(LOOKUP_DIR) $(METADATA_DIR)

# Generate kegg-small.tsv
$(KEGG_SMALL):
	$(JULIA) generate_kegg-small.jl

# Get KEGG compound data
$(KEGG_CPD_JSON):
	$(JULIA) get_kegg_compound_data.jl

# Get KEGG reaction data
$(KEGG_RN_JSON):
	$(JULIA) get_kegg_reaction_data.jl

# Process AutoKEGGRec
reconstruct: $(KEGG_SMALL) $(KEGG_CPD_JSON) $(KEGG_RN_JSON)
	cut -f1 $(KEGG_SMALL) | tail -n+2 > $(DATA_DIR)/kegg-small.lst
	$(JULIA) AutoKEGGRec.jl -i $(DATA_DIR)/kegg-small.lst 

# Clean up generated files
clean:
	rm -rf $(DATA_DIR)
