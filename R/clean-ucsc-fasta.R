# ----- Clean UCSC FASTA files in a directory
## Created: 2022.09.23

# Useful for flagging up exon-intron juctions and replacing
# them with character 'N's so that probe designers avoid placing probes over 
# these junctions. 

# Download FASTA files from UCSC `Genomic Sequence`
# Check: `5'UTR Exons`, `CDS Exons`, `3'UTR Exons`
#        `One FASTA record per region (exon, intron, etc.)`
#        `CDS in upper case, UTR in lower case.`

# Place all .fa files in a directory and call the function as follows:
# clean_ucsc_fasta_directory("/path/to/directory")

# ----- Environment

library(tidyverse)

# ----- Function

clean_ucsc_fasta_directory <- function(directory) {
    fastas <- list.files(directory, pattern = "*.fa", full.names = TRUE) %>% set_names()
    fastas_output_full_path <- fastas %>%
        map(~ str_replace_all(.x, ".fa", "_cleaned.fa"))
    cleaned_fastas <- fastas %>%
        map(function(x) {
            ucsc_fasta <- read_lines(x)
            header <- paste0(">", str_replace(basename(x), ".fa", ""))
            fasta_indices <- str_which(ucsc_fasta, ">") %>% .[2:length(.)]
            for (i in fasta_indices) {
                ucsc_fasta[i] <- "N"
            }
            ucsc_fasta[1] <- header
            return(ucsc_fasta)
        })
    cleaned_fastas %>%
        iwalk(~ write_lines(.x, fastas_output_full_path[[.y]]))
}

# ----- Usage 

clean_ucsc_fasta_directory("~/Desktop")


