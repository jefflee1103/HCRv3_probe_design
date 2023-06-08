# HCRv3 split probe design

Designs HCRv3 probes (split probes) and attaches split initiator probes. This workflow is an adapted version of ProbeDesignHD from Raj lab and Sauka-Spengler lab.

## Installation

### Download the entire repo

### BLAST installations

Install BLAST from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/. Use *.dmg file for Mac. Then check that `R` can find the installed BLAST. If it can't, add the BLAST installation path to R system path using `Sys.setenv(PATH = paste(Sys.getenv("PATH"), path/to/blast/bin, sep = .Platform$path.sep))`

  Sys.which("blastn") # should say sth like "/usr/local/ncbi/blast/bin/blastn" 

For **Linux only**, install `libpq-dev` on your linux machine. Not required for MacOS.  

    sudo apt-get install libpq-dev

Install `metablastr` and dependencies (https://github.com/drostlab/metablastr) in `R`. Make sure your Mac has XCode installed in order to install `devtools`.

    # if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")

    BiocManager::install(c("Biostrings", "GenomicFeatures", "GenomicRanges", "Rsamtools", "IRanges", "rtracklayer", "biomaRt"))

    # install.packages("devtools")
    devtools::install_github("HajkD/metablastr", build_vignettes = TRUE, dependencies = TRUE)

Decompress the BLAST database fasta file (from MK) and copy to a desired directory. Dmel BLAST database is included in this repo. Message me for Human and Mouse BLAST databases. 

    # cd to the repo directory
    gzip -cd data/BLAST/Dmel_ens99_from-MK.fa.gz > /path/to/fasta

### Install other R packages

    install.packages(c("tidyverse", "furrr", "patchwork", "valr"))

## Usage

Follow the `HCRv3_probe_design_workbook.qmd` workflow. Input target sequence can be supplied as a FASTA file, or manually paste the sequence to `target_raw` variable.  

Following hairpins are available in Davis lab and from Haram. 

|   B1  |   B2  |   B3  |   B4  |
|:-----:|:-----:|:-----:|:-----:|
| AF647 | AF488 | AF546 | AF594 |

## Output

If ran successfully, there will be five outputs in a named folder:

    probes.txt          : final probe sequences to order
    probes.pdf          : placing of probes along the target RNA sequence
    details.csv         : thermodynamics and other calculated values of the final probes
    params.txt          : parameters used for the design
    rawblastoutput.csv  : Raw blast output of candidate probes 
