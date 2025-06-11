# Hybridisation Chain Reaction (HCR) probe designer with Shiny GUI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15537771.svg)](https://zenodo.org/doi/10.5281/zenodo.10491693)

Designs third-generation HCR probes (HCRv3, split-probes) and attaches split initiator sequences.  

* B1-5, Choi 2014 ACS Nano
* B7, B9-10, B13-15, B17, Y Wang 2020 BioRxiv 

Please see user manual PDF file for a detailed guide. 

--- 

## System Requirements
* **BLAST+:** The NCBI BLAST+ suite must be installed and available in the system's PATH. You can download it from the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
* **R:** The scripts are written in R (version 4.1.2 or later is recommended).
* **RStudio (Recommended):** An IDE for R, which simplifies package management and script execution.

---

## Installation

1.  **Clone the repository or download ZIP:**

    ```bash
    git clone https://github.com/jefflee1103/HCRv3_probe_design.git
    
    cd HCRv3_probe_design
    ```

2.  **Install R Packages:**
    Launch R or RStudio and run the following commands to install the required packages from CRAN, Bioconductor and Github.

    ```R
    # Install packages from CRAN
    install.packages(c("tidyverse", "shiny", "shinyFiles", "patchwork", "furrr", "DT", "Rcpp", "scales", "valr"))

    # Install packages from Bioconductor
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("Biostrings")
    
    # Install BLAST helper
    # install.packages("devtools")
    devtools::install_github("HajkD/metablastr", build_vignettes = TRUE, dependencies = TRUE)
    ```

3.  **Check BLAST installation**:
    Check that `R` can find the installed BLAST. If it can't, add the BLAST installation path to R system path using `Sys.setenv(PATH = paste(Sys.getenv("PATH"), path/to/blast/bin, sep = .Platform$path.sep))`. `R` must be able to locate your BLAST binary file. 
    
    ```R
    Sys.which("blastn") # should say sth like "/usr/local/ncbi/blast/bin/blastn" 
    ```
4. **Prepare BLAST database and Tx2gene file**:
    BLAST is used to filter out probes potentially cross-hybridising with non-intended RNA targets. Prepare species-specific transcriptome BLAST databases (unzipped `.fa` or `.fasta`) and Tx2gene file, which is used to map `transcrip_id`s (used in BLAST) to `gene_id`s. *Drosophila melanogaster* BLAST database based on ENSEMBL v99 is provided in this repo. See [Zenodo publication](https://zenodo.org/doi/10.5281/zenodo.10491693) for Human and Mouse databases.  

    To design probes for other species, you will need to build a custom BLAST database. A Python script and instructions for this are provided in the `scripts/custom_blast_databases/` directory.

---

## Usage

This tool offers two workflows for designing HCR probes.

### Workflow 1: Shiny App (Recommended for ease-of-use)

The interactive Shiny web application provides a user-friendly graphical interface for designing probes without writing any code.

1.  **Launch the app:**
    Open RStudio, by clicking the `HCR_probe_designer.Rproj` file, and run the following command:
    ```R
    shiny::runApp("shiny_webapp")
    ```
    
2.  **Using the app:**
    * The application will open in your default web browser.
    * Follow the instructions on the "HCRv3" tab.
    * Upload your target transcript fasta file - only takes one file/gene at a time.
    * Adjust the design parameters (e.g., probe length, GC content) as needed.
    * Follow the instructions in the side panel.
    * Once the design process is complete, you can review the results and download the output files, including the probe sequences, a detailed report, and a plot showing probe locations.

### Workflow 2: Script-based pipeline

This workflow is intended for users who wish to integrate probe design into custom scripts or automated pipelines. The core logic is parameterised in a Quarto (`HCR_probe_design_script-workflow.qmd`) document. You can run the script line-by-line in RStudio or use the `quarto` command-line interface to render the entire document, which will execute the code and generate a report.

---

## Output files

Upon successful completion, following files can be downloaded from the Shiny interface. If using the script workflow, it will generate a new directory within the `output/` folder named with the run date and target name (e.g., `output/YYYY-MM-DD_mygene-exon_HCRB1/`). This tool generates following output files.

    probes.txt          : A table of the final selected probes, ready for ordering.
    probes.pdf          : A plot visualising the location of the selected probes on the target transcript.
    details.csv         : A detailed table of thermodynamics and other calculated values of the final probe set.
    params.txt          : A record of the parameters used for the design run.
    rawblastoutput.csv  : The raw output from the BLAST search for further inspection.
