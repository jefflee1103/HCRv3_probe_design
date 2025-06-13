# Guide to Building a Custom BLAST Database

This document provides a detailed walkthrough for using the `build_custom_blast_db.py` script. This script automates the creation of a custom BLAST nucleotide database from genomic data, which is essential for designing HCR probes for species not included by default in this repository.

The script is designed to build a comprehensive RNA-centric database for screening probes against potential off-target binding.

## 1. Core Functionality

The script performs the following key functions:

1.  **Automated Data Retrieval:** Downloads the necessary genomic data directly from URLs. This includes:
    * Genome assembly (DNA FASTA)
    * All cDNA sequences (cDNA FASTA)
    * All non-coding RNA sequences (ncRNA FASTA)
    * Gene annotation file (GTF)
2.  **Intron Extraction:** Extracts intronic sequences from the genome assembly using the GTF file and adds a specified length of flanking exonic sequence to each side of the intron.
3.  **Sequence Concatenation:** Combines the cDNA, ncRNA, and flanked intronic sequences into a single, comprehensive FASTA file (`custom_rna_db.fa`).
4.  **Mapping File Generation:** Creates a `tx2gene.csv` file, which maps transcript identifiers to their corresponding gene identifiers, gene names, and biotypes. This file is critical for interpreting BLAST results in the main probe design workflow.


## 2. Prerequisites

### Software

* **Conda:** It is highly recommended to use Conda for managing the Python environment. You can install Miniconda from [here](https://docs.conda.io/en/latest/miniconda.html).

### Python Environment Setup

A Conda environment file, `environment.yml`, is provided in the `scripts/custom_blast_databases/` directory.

1.  Navigate to the directory containing the environment file:
    ```bash  
    cd scripts/custom_blast_databases/
    ```

2.  Create and activate the Conda environment:
    ```bash  
    conda env create -f environment.yml
    conda activate hcr_custom_blast_db
    ```

## 3. Data Acquisition

The script handles all downloads automatically. Your only task is to find the correct FTP URLs for the required files for your organism of interest. A reliable source for this data is the [Ensembl FTP site](https://ftp.ensembl.org/pub/).

**Crucially,** ensure that all four files (genome, cDNA, ncRNA, and GTF) are from the **exact same release version to maintain consistency between the sequences and their annotations.**

For example, for *Saccharomyces cerevisiae* from Ensembl release 114, you would need the following URLs:

* **Genome:** `https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz`
* **cDNA:** `https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz`
* **ncRNA:** `https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz`
* **GTF:** `https://ftp.ensembl.org/pub/release-114/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz`

## 4. Usage

With the `hcr_custom_blast_db` Conda environment activated, run the script from the command line, providing the URLs you collected.

### Command Structure

```bash  
python build_custom_blast_db.py \
  --genome-fasta-url <URL> \
  --cdna-fasta-url <URL> \
  --ncrna-fasta-url <URL> \
  --gtf-url <URL> \
  --output-dir <PATH> \
  [--flanking-nt <INT>]
```

### Arguments

* `--genome-fasta-url` (Required): URL to the gzipped genome DNA FASTA file.
* `--cdna-fasta-url` (Required): URL to the gzipped cDNA FASTA file.
* `--ncrna-fasta-url` (Required): URL to the gzipped ncRNA FASTA file.
* `--gtf-url` (Required): URL to the gzipped GTF annotation file.
* `--output-dir` (Required): Path to a directory where all output files will be stored. The script will create this directory if it does not exist.
* `--flanking-nt` (Optional): The number of nucleotides to include as flanking regions on each side of an intron. Defaults to `40`.

### Example

Using the yeast URLs from above:

```bash
python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
  --cdna-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  --ncrna-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz \
  --gtf-url https://ftp.ensembl.org/pub/release-114/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz \
  --output-dir ./yeast_hcr_rna_db_output \
  --flanking-nt 40
```

## 5. Script Outputs

Upon successful completion, the specified output directory will be populated with several files:

* **Downloaded & Unzipped Data:**. 
    * `genome.dna.fa`: Uncompressed genomic DNA.
    * `cdna.all.fa`: Uncompressed cDNA sequences.
    * `ncrna.fa`: Uncompressed non-coding RNA sequences.
    * `annotation.gtf`: Uncompressed annotation file.
* **Processed Files:**. 
    * `introns_flank40.fa` (or similar): FASTA file of extracted intronic sequences with 40 nt flanking regions.
    * `tx2gene.csv`: The transcript-to-gene mapping file.
* **Final FASTA for BLAST:**. 
    * `custom_rna_db.fa`: The final, concatenated FASTA file containing all cDNA, ncRNA, and flanked intron sequences.

## 6. Creating and Integrating the Database

1.  **Create the BLAST Database:**
    After the Python script finishes, it will print a command similar to this:

    ```bash
    makeblastdb -in ./yeast_hcr_rna_db_output/custom_rna_db.fa -dbtype nucl -parse_seqids -out ./yeast_hcr_rna_db_output/custom_rna_blast_db -title "Custom_RNA_DB_from_script"
    ```

    Copy and paste this command into your terminal and execute it. This will create the actual BLAST database files (e.g., `custom_rna_blast_db.nhr`, `.nin`, `.nsq`).

1.  **Move Database Files:**
    Move the entire output directory (e.g., `yeast_hcr_rna_db_output/`) into the `data/BLAST/` directory of the HCR probe design repository.

2.  **Update Paths in the HCR Designer:**
    When you are ready to design probes for your custom species:
    * **In the Shiny App:** Select "Custom" from the species dropdown menu. This will reveal fields where you must provide the full path to your new BLAST database name (e.g., `data/BLAST/yeast_hcr_rna_db_output/custom_rna_blast_db`) and your new `tx2gene.csv` file (`data/BLAST/yeast_hcr_rna_db_output/tx2gene.csv`).
    * **In the Script-based Workflow:** Modify the `blast_file` and `tx2gene` variables in the `.qmd` script to point to your new database and mapping file, respectively.
