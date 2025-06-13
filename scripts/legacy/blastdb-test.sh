#!/bin/bash

# ====================================================================================
#
# build_dbs_for_multiple_species.sh
#
# Description:
# This script automates the creation of custom BLAST databases for multiple species
# by repeatedly calling the 'build_custom_blast_db.py' script.
#
# It is configured for four species: Human, Mouse, Zebrafish, and Yeast, using
# URLs from Ensembl release 112. You can easily add or modify the species list
# and corresponding URLs in the configuration section.
#
# For each species, the script will:
# 1. Create a dedicated output directory.
# 2. Call 'build_custom_blast_db.py' to download the required genome, cDNA, ncRNA,
#    and GTF files, and then process them into a single FASTA file.
# 3. Call 'makeblastdb' from the NCBI BLAST+ suite to create a nucleotide BLAST
#    database from the generated FASTA file.
#
# Requirements:
# - A Python 3 environment with the necessary libraries for the Python script
#   (requests, gffutils, pyfaidx, numpy).
# - NCBI BLAST+ command-line tools installed and in your system's PATH,
#   specifically the 'makeblastdb' executable.
# - The 'build_custom_blast_db.py' script must be in the same directory as this
#   shell script, or you must provide a full path to it.
#
# ====================================================================================

# --- Configuration ---

# Path to your Python interpreter and the build script.
# Modify these if they are not in the current directory or your system's PATH.
PYTHON_EXE="python3"
BUILD_SCRIPT_PATH="./build_custom_blast_db.py"

# Number of flanking nucleotides for intron extraction.
FLANKING_NT=40

# --- Species Data Configuration (Ensembl Release 112) ---

# Common FTP base URL for Ensembl Release 112
ENSEMBL_BASE_URL="https://ftp.ensembl.org/pub/release-112"

# Array of species names (used for directory naming and logging)
SPECIES_NAMES=(
    "homo_sapiens"
    "mus_musculus"
    "danio_rerio"
    "saccharomyces_cerevisiae"
)

# --- URLs for each species ---

# Human (Homo sapiens, GRCh38)
HOMO_SAPIENS_GENOME_URL="${ENSEMBL_BASE_URL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
HOMO_SAPIENS_CDNA_URL="${ENSEMBL_BASE_URL}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
HOMO_SAPIENS_NCRNA_URL="${ENSEMBL_BASE_URL}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
HOMO_SAPIENS_GTF_URL="${ENSEMBL_BASE_URL}/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"

# Mouse (Mus musculus, GRCm39)
MUS_MUSCULUS_GENOME_URL="${ENSEMBL_BASE_URL}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
MUS_MUSCULUS_CDNA_URL="${ENSEMBL_BASE_URL}/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
MUS_MUSCULUS_NCRNA_URL="${ENSEMBL_BASE_URL}/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz"
MUS_MUSCULUS_GTF_URL="${ENSEMBL_BASE_URL}/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz"

# Zebrafish (Danio rerio, GRCz11)
DANIO_RERIO_GENOME_URL="${ENSEMBL_BASE_URL}/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
DANIO_RERIO_CDNA_URL="${ENSEMBL_BASE_URL}/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz"
DANIO_RERIO_NCRNA_URL="${ENSEMBL_BASE_URL}/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz"
DANIO_RERIO_GTF_URL="${ENSEMBL_BASE_URL}/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz"

# Yeast (Saccharomyces cerevisiae, R64-1-1)
SACCHAROMYCES_CEREVISIAE_GENOME_URL="${ENSEMBL_BASE_URL}/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
SACCHAROMYCES_CEREVISIAE_CDNA_URL="${ENSEMBL_BASE_URL}/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
SACCHAROMYCES_CEREVISIAE_NCRNA_URL="${ENSEMBL_BASE_URL}/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"
SACCHAROMYCES_CEREVISIAE_GTF_URL="${ENSEMBL_BASE_URL}/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.112.gtf.gz"

# --- Main Processing Loop ---

# Check if the build script exists
if [ ! -f "$BUILD_SCRIPT_PATH" ]; then
    echo "Error: Build script not found at '$BUILD_SCRIPT_PATH'. Please check the path."
    exit 1
fi

# Check if makeblastdb is installed
if ! command -v makeblastdb &> /dev/null; then
    echo "Error: 'makeblastdb' command not found. Please install NCBI BLAST+ and ensure it's in your PATH."
    exit 1
fi

echo "Starting multi-species BLAST database build process..."
echo "======================================================"

for SPECIES in "${SPECIES_NAMES[@]}"; do
    echo -e "\n--- Processing species: $SPECIES ---\n"

    # Dynamically construct variable names for the URLs
    GENOME_URL_VAR="${SPECIES^^}_GENOME_URL"
    CDNA_URL_VAR="${SPECIES^^}_CDNA_URL"
    NCRNA_URL_VAR="${SPECIES^^}_NCRNA_URL"
    GTF_URL_VAR="${SPECIES^^}_GTF_URL"

    # Use indirect variable expansion to get the actual URLs
    GENOME_URL=${!GENOME_URL_VAR}
    CDNA_URL=${!CDNA_URL_VAR}
    NCRNA_URL=${!NCRNA_URL_VAR}
    GTF_URL=${!GTF_URL_VAR}

    # Define the output directory for the current species
    OUTPUT_DIR="./${SPECIES}_db_output"

    # Step 1: Run the Python script to download and prepare the FASTA file
    echo "Step 1: Running build_custom_blast_db.py for $SPECIES..."
    "$PYTHON_EXE" "$BUILD_SCRIPT_PATH" \
      --genome-fasta-url "$GENOME_URL" \
      --cdna-fasta-url "$CDNA_URL" \
      --ncrna-fasta-url "$NCRNA_URL" \
      --gtf-url "$GTF_URL" \
      --output-dir "$OUTPUT_DIR" \
      --flanking-nt "$FLANKING_NT"

    # Check if the Python script executed successfully
    if [ $? -ne 0 ]; then
        echo "Error: Python script failed for $SPECIES. Skipping to the next species."
        continue
    fi

    # Step 2: Run makeblastdb to create the BLAST database
    CONCATENATED_FASTA_PATH="${OUTPUT_DIR}/custom_rna_db.fa"
    BLAST_DB_PATH="${OUTPUT_DIR}/${SPECIES}_blast_db"

    if [ -s "$CONCATENATED_FASTA_PATH" ]; then
        echo -e "\nStep 2: Creating BLAST database for $SPECIES..."
        makeblastdb \
          -in "$CONCATENATED_FASTA_PATH" \
          -dbtype nucl \
          -parse_seqids \
          -out "$BLAST_DB_PATH" \
          -title "${SPECIES} Custom RNA Database"

        if [ $? -eq 0 ]; then
            echo "Successfully created BLAST database for $SPECIES at: $BLAST_DB_PATH"
        else
            echo "Error: makeblastdb failed for $SPECIES."
        fi
    else
        echo "Error: The combined FASTA file was not created or is empty. Cannot run makeblastdb for $SPECIES."
    fi

    echo "--- Finished processing $SPECIES ---"
done

echo -e "\n======================================================"
echo "All species have been processed. Process complete."





python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz \
  --gtf-url          https://ftp.ensembl.org/pub/release-114/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz \
  --output-dir       ./yeast_hcr_rna_db_output \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz \
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-114/fasta/danio_rerio/cdna/Danio_rerio.GRCz11.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-114/fasta/danio_rerio/ncrna/Danio_rerio.GRCz11.ncrna.fa.gz \
  --gtf-url          https://ftp.ensembl.org/pub/release-114/gtf/danio_rerio/Danio_rerio.GRCz11.114.gtf.gz \
  --output-dir       ~/Desktop/Drer114 \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.UCB_Xtro_10.0.dna.toplevel.fa.gz \
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-114/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.UCB_Xtro_10.0.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-114/fasta/xenopus_tropicalis/ncrna/Xenopus_tropicalis.UCB_Xtro_10.0.ncrna.fa.gz \
  --gtf-url          https://ftp.ensembl.org/pub/release-114/gtf/xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.114.gtf.gz \
  --output-dir       ~/Desktop/Xtro114 \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz \
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz \
  --gtf-url          https://ftp.ensembl.org/pub/release-114/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.114.gtf.gz \
  --output-dir       ~/Desktop/Cele114 \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz \
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-114/fasta/gallus_gallus/cdna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-114/fasta/gallus_gallus/ncrna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.ncrna.fa.gz \
  --gtf-url          https://ftp.ensembl.org/pub/release-114/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.114.gtf.gz \
  --output-dir       ~/Desktop/Ggal114 \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url  \
  --cdna-fasta-url    \
  --ncrna-fasta-url   \
  --gtf-url           \
  --output-dir       ~/Desktop/Ecoli61 \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url  \
  --cdna-fasta-url    \
  --ncrna-fasta-url   \
  --gtf-url           \
  --output-dir       ~/Desktop/ \
  --flanking-nt      40


  