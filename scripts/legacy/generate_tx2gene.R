# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create tx2gene file to convert tx_id to gene_id for BLAST output
#  when designing HCRv3 probes 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Environment -------------------------------------------------------------

library(dplyr)
library(AnnotationHub) # Install from Bioconductor
ah <- AnnotationHub()

# Fly genes ---------------------------------------------------------------

ens <- query(ah, c("Drosophila melanogaster", "EnsDb"))
tibble(x = ens$ah_id, y = ens$title) # Find release 99
ens <- ens[["AH78755"]] 
genes <- genes(ens, return.type = "data.frame") %>%
  dplyr::select(gene_id, gene_name, gene_biotype)
transcripts <- transcripts(ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
tx2gene <- left_join(transcripts, genes)
# readr::write_csv(tx2gene, "./src/Dmel_tx2gene_ENSEMBL_v99.csv")


# Human genes -------------------------------------------------------------

ens <- query(ah, c("Homo sapiens", "EnsDb"))
tibble(x = ens$ah_id, y = ens$title) %>% head(20) # Find release 99
ens <- ens[["AH78783"]] 
genes <- genes(ens, return.type = "data.frame") %>%
  dplyr::select(gene_id, gene_name, gene_biotype)
transcripts <- transcripts(ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
tx2gene <- left_join(transcripts, genes)
# readr::write_csv(tx2gene, "./src/Hsap_tx2gene_ENSEMBL_v99.csv")


# Mouse genes -------------------------------------------------------------

ens <- query(ah, c("Mus musculus", "EnsDb"))
tibble(x = ens$ah_id, y = ens$title) %>% head(20) # Find release 99
ens <- ens[["AH78811"]] 
genes <- genes(ens, return.type = "data.frame") %>%
  dplyr::select(gene_id, gene_name, gene_biotype)
transcripts <- transcripts(ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
tx2gene <- left_join(transcripts, genes)
# readr::write_csv(tx2gene, "./src/Mmus_tx2gene_ENSEMBL_v99.csv")





