#!/bin/bash

# --- Configuration ---
DATA_DIR="data"
# Replace the URL below with your actual Google Drive direct link or Zenodo URL
# Example: Zenodo direct download link
IS_REFERENCE_URL="https://zenodo.org/record/18996276/files/reference_genome.gb"
FQ1_URL="https://zenodo.org/record/18996276/files/reads_R1.fastq.gz"
FQ2_URL="https://zenodo.org/record/18996276/files/reads_R2.fastq.gz"
IS_ELEMENTS_URL="https://zenodo.org/record/18996276/files/is_elements.fasta"


# --- Setup ---
echo "🚀 Starting data download for ISdetector..."

# Create data directory if it doesn't exist
if [ ! -d "$DATA_DIR" ]; then
    echo "📁 Creating $DATA_DIR directory..."
    mkdir -p "$DATA_DIR"
fi

# --- Download Files ---
# Function to download if file doesn't exist
download_file() {
    local URL=$1
    local FILENAME=$2
    if [ -f "$DATA_DIR/$FILENAME" ]; then
        echo "✅ $FILENAME already exists, skipping."
    else
        echo "📥 Downloading $FILENAME..."
        curl -L $URL -o "$DATA_DIR/$FILENAME"
    fi
}

download_file "$IS_REFERENCE_URL" "reference_genome.gb"
download_file "$FQ1_URL" "reads_R1.fastq.gz"
download_file "$FQ2_URL" "reads_R2.fastq.gz"
download_file "$IS_ELEMENTS_URL" "is_elements.fasta"

echo "✨ All data files are ready in the /$DATA_DIR folder!"
