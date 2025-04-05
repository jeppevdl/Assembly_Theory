#!/bin/bash

# Define the working directories and JSON files
MAP="pathway_no_MA"
WORK_DIR="$HOME/OneDrive/Documenten/Bioinformatics/Tweede master/Master Thesis/Assembly_Theory/GEMs"
ASSEMBLY_GO_DIR="$WORK_DIR/bin/assembly_go"
INPUT_FILE="$WORK_DIR/data/pathway_complexities/complexities_$MAP.json"
OUTPUT_FILE="$WORK_DIR/data/bash_MA_output/bash_output_$MAP.json"
TIMEOUT_DURATION=30  # Timeout in seconds

# Ensure jq is installed
if ! command -v jq &> /dev/null; then
    echo "Error: jq is not installed. Please install jq to proceed."
    exit 1
fi

# Read the input JSON
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found."
    exit 1
fi

# Extract compound IDs from JSON
COMPOUND_IDS=$(jq -r '.columns[0][]' "$INPUT_FILE" | tr -d '\r' | tail -n +5691)

# Change to the assembly directory
cd "$ASSEMBLY_GO_DIR" || exit 1

# Initialize results associative array
declare -A results

# Loop over each compound ID
echo "cpd\tma" > "$WORK_DIR/data/bash_MA_output/bash_MA_$MAP.tsv"
for compoundId in $COMPOUND_IDS; do
    echo "Processing compound ID: $compoundId"
    COMMAND="./assembly molfiles/$compoundId.mol"
    
    # Run the command with a timeout and send SIGINT if it times out
    timeout --foreground -s SIGINT $TIMEOUT_DURATION $COMMAND > "output/output_$compoundId.txt" 2> "output/error_$compoundId.txt"
    EXIT_STATUS=$?
    
    # Extract the first number found in the output file
    result=$(grep -o '[0-9]\+' "output/output_$compoundId.txt" | head -n 1)
    
    if [ -n "$result" ]; then
        echo -e "$compoundId\t$result" >> "$WORK_DIR/data/bash_MA_output/bash_MA_$MAP.tsv"
    else
        echo -e "$compoundId\tna" >> "$WORK_DIR/data/bash_MA_output/bash_MA_$MAP.tsv"
    fi
    
    # Check for errors
    if [ -s "output/error_$compoundId.txt" ]; then
        echo "Error output for compound ID $compoundId:"
        cat "output/error_$compoundId.txt"
    fi
done