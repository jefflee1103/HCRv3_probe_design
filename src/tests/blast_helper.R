
ws_wrap <- function(paths){
  
  paths <- stringr::str_split(paths,.Platform$file.sep)
  
  ws.paths <- unlist(lapply(paths, function(x) {
    
    stringr::str_c(unlist(sapply(x, function(word) {
      if (stringr::str_detect(word," ")) 
        return(stringr::str_replace(word,word,paste0("'",word,"'"))) else 
          return(word)})),collapse = .Platform$file.sep)
    
  }))
  
  return(ws.paths) 
}

is_outformat <- function(out.format) {
  if (!is.element(out.format, list_outformats()))
    stop("'",out.format, "' is not an available BLAST output format.", call. = FALSE)
  
  return(TRUE)
}

outformat2num <- function(out.format) {
  if (is_outformat(out.format)) {
    return(switch(
      out.format,
      pair = 0,
      qa.ident = 1,
      qa.nonident = 2,
      fq.ident = 3,
      fq.nonident = 4,
      xml = 5,
      tab = 6,
      tab.comment = 7,
      ASN.1.text = 8,
      ASN.1.binary = 9,
      csv = 10,
      ASN.1 = 11,
      json.seq.aln = 12,
      json.blast.multi = 13,
      xml2.blast.multi = 14,
      json.blast.single = 15,
      xml2.blast.single = 16,
      SAM = 17,
      report = 18
    ))
    
  }
}

is_blast_installed <- function(path = NULL) {
  # test if a valid BLAST version is installed
  tryCatch({
    if (is.null(path)) {
      sys_out <-
        system("blastp -version", intern = TRUE)
    } else {
      sys_out <-
        system(paste0(
          'export PATH=$PATH:',
          path, "'; blastp -version '"), intern = TRUE)
    }
    
    
  }, error = function(e)
    stop(
      "It seems like you don't have BLAST installed locally on your machine or the PATH variable to the BLAST program is not set correctly.",
      call. = FALSE
    ))
  
  if (any(stringr::str_detect(sys_out, "blast")))
    return(TRUE)
  
}

blast_nt_to_nt <- function (query, subject, strand = "both", output.path = NULL, 
                            is.subject.db = FALSE, task = "blastn", db.import = FALSE, 
                            postgres.user = NULL, evalue = 0.001, out.format = "csv", 
                            cores = 1, max.target.seqs = 10000L, db.soft.mask = FALSE, 
                            db.hard.mask = FALSE, blast.path = NULL, word_size = word_size) 
{
  if (!is_blast_installed(path = blast.path))
    stop("Please install a valid version of blastn. See Installation Vignette for details.",
         call. = FALSE)
  if (db.import) {
    if (!is.element(out.format, c("xml", "tab", "csv"))) 
      stop("Only output formats: 'xml', 'tab', or 'csv' can be imported.", 
           call. = FALSE)
  }
  if (!is.element(strand, c("both", "plus", "minus"))) 
    stop("Please specify a strand option that is supported by BLAST: strand = 'both', strand = 'plus', strand = 'minus'.")
  file_contains_dna(query, "query")
  file_contains_dna(subject, "subject")
  multi.cores <- parallel::detectCores()
  if (cores > multi.cores) 
    stop("You chose more cores than are available on your machine.", 
         call. = FALSE)
  if (!file.exists(query)) 
    stop("Unfortunately, no query file has been found at ", 
         query, call. = FALSE)
  if (!file.exists(subject)) 
    stop("Unfortunately, no subject file has been found at ", 
         subject, call. = FALSE)
  if (!is.element(task, c("blastn", "blastn-short", "dc-megablast", 
                          "megablast", "rmblastn"))) 
    stop("Please choose a nucleotide-nucleotide comparison task that is supported by BLAST: task = 'blastn', task = 'blastn-short', task = 'dc-megablast', task = 'megablast', task = 'rmblastn'.", 
         call. = FALSE)
  message("Starting 'blastn -task ", task, "' with  query: ", 
          query, " and subject: ", subject, " using ", cores, 
          " core(s) ...")
  blast_call <- paste0("blastn -query ", ws_wrap(query), " -db ", 
                       ws_wrap(subject))
  output_blast <- file.path(ifelse(is.null(output.path), ws_wrap(getwd()), 
                                   ws_wrap(output.path)), paste0(unlist(stringr::str_split(basename(query), 
                                                                                           "[.]"))[1], "_", unlist(stringr::str_split(basename(subject), 
                                                                                                                                      "[.]"))[1], "_", task, "_eval_", evalue, ".blast_tbl"))
  output_read_blast <- file.path(ifelse(is.null(output.path), 
                                        getwd(), output.path), paste0(unlist(stringr::str_split(basename(query), 
                                                                                                "[.]"))[1], "_", unlist(stringr::str_split(basename(subject), 
                                                                                                                                           "[.]"))[1], "_", task, "_eval_", evalue, ".blast_tbl"))
  if (!is.subject.db) {
    if (is.null(blast.path)) {
      system(paste0("makeblastdb -in ", subject, " -input_type fasta -dbtype nucl -hash_index"))
    }
    else {
      system(paste0("export PATH=", blast.path, "; makeblastdb -in ", 
                    subject, " -input_type fasta -dbtype nucl -hash_index"))
    }
  }
  system(paste0(ifelse(is.null(blast.path), blast_call, paste0("export PATH=$PATH:", 
                                                               blast_call)), " -evalue ", evalue, " -strand ", strand, 
                " -max_target_seqs ", as.integer(max.target.seqs), " -out ", 
                output_blast, "-word_size", word_size, " -num_threads ", cores, ifelse(db.soft.mask, 
                                                                                       " -db_soft_mask", ""), ifelse(db.hard.mask, " -db_hard_mask", 
                                                                                                                     ""), paste0(" -task ", task), paste0(" -outfmt \"", 
                                                                                                                                                          outformat2num(out.format = out.format), " qseqid sseqid pident nident length mismatch gapopen gaps positive ppos qstart qend qlen qcovs qcovhsp sstart send slen evalue bitscore score\"")))
  if (!is.null(db.import)) {
    if (db.import) {
      blast_tbl <- read_blast(file = output_read_blast, 
                              out.format = "postgres", postgres.user = postgres.user)
      message("\n")
      message("BLAST search finished! A Postgres  database connection to the BLAST output file has been generated. The BLAST output file can be found at: ", 
              output_blast)
      return(blast_tbl)
    }
    else {
      blast_tbl <- read_blast(file = output_read_blast, 
                              out.format = out.format, postgres.user = NULL)
      message("\n")
      message("BLAST search finished! The BLAST output file was imported into the running R session. The BLAST output file has been stored at: ", 
              output_blast)
      blast_files_output_qry <- list.files(dirname(query))
      blast_files_output_sbj <- list.files(dirname(subject))
      blast_files_to_remove_sbj <- file.path(dirname(subject), 
                                             blast_files_output_sbj[which(stringr::str_detect(blast_files_output_sbj, 
                                                                                              "[.]n.."))])
      message("Removing subject BLAST databases ", paste0(basename(blast_files_to_remove_sbj), 
                                                          collapse = ", "), " from '", dirname(subject), 
              "'.")
      unlink(blast_files_to_remove_sbj, force = TRUE)
      message("All done!")
      return(blast_tbl)
    }
  }
  else {
    message("\n")
    message("BLAST search finished! BLAST output file has been stored at: ", 
            output_blast)
  }
}


run_blastn_short_custom <- function(df, db, tx2gene, cores = n_threads, word_size = word_size){
  # Read the provided tx2gene mapping file
  tx2gene_df <- read_csv(tx2gene)
  
  # Prepare the query sequences in a temporary FASTA file
  candidates <- df
  seqs <- DNAStringSet(candidates$target_sequence)
  names(seqs) <- candidates$unique_id
  
  temp_fasta_path <- file.path(tempdir(), "query.fasta")
  Biostrings::writeXStringSet(seqs, temp_fasta_path, format = "fasta")
  
  # Execute the blastn-short command using the metablastr package
  blast_output <- blast_nt_to_nt(
    query = temp_fasta_path,
    subject = db,
    evalue = 10,
    task = "blastn-short",
    strand = "plus",
    cores = cores,
    word_size = word_size,
    output.path = tempdir()
  )
  
  # Process and annotate the BLAST output
  blast_output <- blast_output %>%
    filter(bit_score > 32) %>%
    tidyr::separate(col = subject_id, into = c("subject_id", "is_intron"), sep = "_") %>%
    mutate(subject_id = str_replace(subject_id, "\\..*", "")) %>%
    left_join(tx2gene_df, by = c("subject_id" = "tx_id")) %>%
    dplyr::select(query_id, subject_id, gene_id, gene_name, everything()) %>%
    mutate(gene_id = if_else(str_detect(subject_id, "FBgn|ENS"), subject_id, gene_id))
  
  return(blast_output)
}

run_blastn_custom <- function(df, db, tx2gene, cores = n_threads, word_size = word_size){
  # Read the provided tx2gene mapping file
  tx2gene_df <- read_csv(tx2gene)
  
  # Prepare the query sequences in a temporary FASTA file
  candidates <- df
  seqs <- DNAStringSet(candidates$target_sequence)
  names(seqs) <- candidates$unique_id
  
  temp_fasta_path <- file.path(tempdir(), "query.fasta")
  Biostrings::writeXStringSet(seqs, temp_fasta_path, format = "fasta")
  
  # Execute the blastn-short command using the metablastr package
  blast_output <- blast_nt_to_nt(
    query = temp_fasta_path,
    subject = db,
    evalue = 10,
    task = "blastn",
    strand = "plus",
    cores = cores,
    word_size = word_size,
    output.path = tempdir()
  )
  
  # Process and annotate the BLAST output
  blast_output <- blast_output %>%
    filter(bit_score > 32) %>%
    tidyr::separate(col = subject_id, into = c("subject_id", "is_intron"), sep = "_") %>%
    mutate(subject_id = str_replace(subject_id, "\\..*", "")) %>%
    left_join(tx2gene_df, by = c("subject_id" = "tx_id")) %>%
    dplyr::select(query_id, subject_id, gene_id, gene_name, everything()) %>%
    mutate(gene_id = if_else(str_detect(subject_id, "FBgn|ENS"), subject_id, gene_id))
  
  return(blast_output)
}