# HCRv3 split probe design

Designs HCRv3 probes (split probes) and attaches split initiator probes. This workflow is an adapted version of ProbeDesignHD from Raj lab and Sauka-Spengler lab.

## Installation

### Download the entire repo



### BLAST installation

It is recommended to use conda environment for BLAST installation. Make sure to use x86 version of Miniconda3, not the arm64 version: https://docs.conda.io/en/latest/miniconda.html

Create conda environment for BLAST.

    conda create -n BLAST python=3.7
    conda activate BLAST
    conda install -c bioconda blast
    blastn -h  # check installation

Decompress BLAST database fasta (from MK, Dmel only) and copy to a desired directory and build a database.

    # cd to the repo directory
    gzip -cd data/BLAST/Dmel_ens99_from-MK.fa.gz > /path/to/fasta

    conda activate BLAST
    makeblastdb -in /path/to/fasta -parse_seqids -title Dmel_ens99 -dbtype nucl

### Install rBLAST (R wrapper for BLAST)

Install `rBLAST` (https://github.com/mhahsler/rBLAST) and `Biostrings` in `R`. Make sure your Mac has XCode installed in order to install `devtools`.

    install.packages("devtools")
    devtools::install_bioc("Biostrings")
    devtools::install_bioc("rtracklayer")
    devtools::install_github("mhahsler/rBLAST")

### Add PATH to BLAST in R System

Conda-installed binaries are not found in `R` by default. To append the BLAST PATH in R system, find the installed `BLAST` binary. On macos, it should be under `opt/miniconda3/envs/BLAST/bin`.

    Sys.getenv("PATH") # the current PATH searched in R
    Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/opt/miniconda3/envs/BLAST/bin", sep = .Platform$path.sep))
    Sys.getenv("PATH")

Check `R` can find the `BLAST` Installation

    Sys.which("blastn")
    # should say something like "/opt/miniconda3/envs/BLAST/bin/blastn"

NOTE: the PATH edit is temporary and only applies to the current R session. so the above chunk should be run every time with a new session. The `HCRv3_probe_design_workbook.qmd` has environment chunk that handles this.  

## Usage

Follow the `HCRv3_probe_design_workbook.qmd` workflow. Input target sequence can be supplied as a FASTA file. 

Following hairpins are available in Davis lab and from Haram. 

|   B1  |   B2  |   B3  |   B4  |
|:-----:|:-----:|:-----:|:-----:|
| AF647 | AF488 | AF546 | Af594 |

## Output

If ran successfully, there will be five outputs in a named folder:

    probes.txt          : final probe sequences to order
    probes.pdf          : placing of probes along the target RNA sequence
    details.csv         : thermodynamics and other calculated values of the final probes
    params.txt          : parameters used for the design
    rawblastoutput.csv  : Raw blast output of candidate probes 
