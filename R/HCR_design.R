library(tidyverse) 
library(furrr) 
library(patchwork) 
library(valr) 
library(Biostrings)
library(metablastr)
oligo_length <- 52

clean_pasted_seq <- function(sequence){
  sequence %>%
    str_trim() %>%
    str_replace_all(" ", "") %>%
    str_replace_all("\n", "") %>%
    str_replace_all("\r", "") %>%
    str_to_upper() %>%
    str_replace_all("U", "T") %>%
    return()
}


generate_candidate_probes <- function(target_seq, oligo_length){
  # Generate staggered sequences of probe candidates with fixed length
  df <- map_dfr(1:nchar(target_seq), ~ {
    start <- .x
    end <- .x + (oligo_length - 1)
    sequence <- str_sub(target, start, end)
    tibble(start = start,
           end = end,
           target_sequence = sequence,
           length = nchar(sequence))
  }) %>%
    filter(!str_detect(target_sequence, "N")) %>%
    filter(length == oligo_length) %>%
    mutate(unique_id = paste0("id_", start)) %>%
    mutate(centre_pos = (start + end) / 2) %>%
    dplyr::select(unique_id, centre_pos, everything())
  
  if(nrow(df) == 0){
    return(print("ERROR: Candidate probes could not be generated. Check your FASTA file!"))
  } else {
    return(df)
  }
}


inspect_probe_prep <- function(df){
  ggplot(df,
         aes(x = " ", 
             ymin = start, ymax = end, 
             colour = start)) + 
    geom_linerange(position = position_dodge2(width = 0.75)) + 
    geom_hline(yintercept = as_tibble(str_locate_all(target, pattern = "N")[[1]]) %>% pull(start),
               colour = "gray60") + 
    labs(title = paste0(target_name, 
                        ": all sliding probes - ",
                        target_sequence_total_length, " bp in length"),
         subtitle = "Junction-spanning probes removed (if any present)",
         x = "", y = "") + 
    coord_flip() + 
    scale_colour_viridis_c() +
    theme_minimal() +
    theme(legend.position = "none")
}


str_revcomp <- function(seq){
  # For reverse complementation of character string
  chartr("ATGC", "TACG", seq) %>% stringi::stri_reverse()
}


calculate_thermodynamics <- function(seq, temperature, Na, oligo_conc){
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Calculate thermodymic properties of target **RNA** sequence 
  # Output is a named list of Tm, dG, dH, dS values 
  
  ## RNA-DNA energetics from SantaLucia 98 (Takes RNA sequence dimension as input!!)
  thermodynamics_table <- tibble(
    dimension = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"),
    delH = c(-7.8, -5.9, -9.1, -8.3, -9.0, -9.3, -16.3, -7.0, -5.5, -8.0, -12.8, -7.8, -7.8, -8.6, -10.4, -11.5),
    delS = c(-21.9,-12.3,-23.5,-23.9,-26.1,-23.2,-47.1,-19.7,-13.5,-17.1,-31.9,-21.6,-23.2,-22.9,-28.4,-36.4)
  )
  
  seq_dimension <- map_dfr(1:(nchar(seq)-1), ~ 
                             tibble(dimension = str_sub(seq, .x, .x + 1))) %>%
    left_join(thermodynamics_table, by = "dimension") 
  
  delH <- sum(seq_dimension$delH) ## kcal/mol
  delS <- sum(seq_dimension$delS) ## cal/(mol*Kelvin)
  
  initH <- 1.9
  initS <- -3.9 
  
  dH <- delH + initH
  dS <- delS + initS
  dG <- (dH * 1000 - (temperature + 273.15) * dS) / 1000; ## kcal/mol
  Tm <- (dH * 1000 / (dS + (1.9872 * log(oligo_conc / 4))) ) - 273.15 + 16.6 * log10(Na)
  
  output <- list(Tm, dG, dH, dS) %>% set_names(c("Tm", "dG", "dH", "dS"))
  return(output)
}

calculate_thermodynamics_v <- Vectorize(calculate_thermodynamics)


get_thermodynamic_parameters <- function(candidate_probes, temperature, Na, oligo_conc) {
  candidate_probe_thermodynamics <- future_map(
    candidate_probes$target_sequence,
    ~ calculate_thermodynamics(
      seq = .x,
      temperature = temperature,
      Na = Na,
      oligo_conc = oligo_conc
    )
  )

  candidate_probe_thermodynamics_1st_half <- future_map(
    map_chr(candidate_probes$target_sequence, ~ str_sub(.x, start = 1, end = nchar(.x) / 2 - 1)),
    ~ calculate_thermodynamics(
      seq = .x,
      temperature = temperature,
      Na = Na,
      oligo_conc = oligo_conc
    )
  )

  candidate_probe_thermodynamics_2nd_half <- future_map(
    map_chr(candidate_probes$target_sequence, ~ str_sub(.x, start = nchar(.x) / 2 + 2, end = nchar(.x))),
    ~ calculate_thermodynamics(
      seq = .x,
      temperature = temperature,
      Na = Na,
      oligo_conc = oligo_conc
    )
  )

  output <- list(candidate_probe_thermodynamics, candidate_probe_thermodynamics_1st_half, candidate_probe_thermodynamics_2nd_half) %>%
    set_names(c("Full", "First", "Second"))

  return(output)
}


passed_a_comp <- function(theProbeSeq){

  A_composition <- str_count(theProbeSeq, pattern = "A") / nchar(theProbeSeq)
  theVerdict <- if_else(A_composition < 0.28, TRUE, FALSE)
  return(theVerdict)
}

passed_a_stack <- function(theProbeSeq){

  theVerdict <- if_else(
    str_detect(theProbeSeq, "AAAA"), FALSE, TRUE
  )
  return(theVerdict)
}

passed_c_comp <- function(theProbeSeq){

  C_composition <- str_count(theProbeSeq, pattern = "C") / nchar(theProbeSeq)
  theVerdict <- if_else(
    C_composition > 0.22 & C_composition < 0.28, TRUE, FALSE)
  return(theVerdict)
}

passed_c_stack <- function(theProbeSeq){

  theVerdict <- if_else(
    str_detect(theProbeSeq, "CCCC"), FALSE, TRUE
  )
  return(theVerdict)
}

passed_c_spec_stack <- function(theProbeSeq){
  max_c_comp <- map_dbl(1:7, ~ {
    substring <- str_sub(theProbeSeq, start = .x, end = .x + 5)
    C_composition <- str_count(substring, pattern = "C") / nchar(substring)
  }) %>% max()
  theVerdict <- if_else(max_c_comp <= 0.5, TRUE, FALSE)
  return(theVerdict)
} 


annotate_probes <- function(candidate_probes, thermodynamics){
  # Takes candidate_probes df and adds thermodynamic and nucleotide composition parameters
  
  candidate_probes %>%
    mutate(rev_comp = str_revcomp(target_sequence)) %>%
    mutate(
      GC_content = str_count(target_sequence, pattern = "G|C") / length,
      A_content  = str_count(target_sequence, pattern = "A") / length,
      C_content  = str_count(target_sequence, pattern = "C") / length
    ) %>%
    bind_cols(dG = map_dbl(thermodynamics$Full, pluck("dG"))) %>%
    bind_cols(Tm = map_dbl(thermodynamics$Full, pluck("Tm"))) %>%
    bind_cols(dG_1st_half = map_dbl(thermodynamics$First, pluck("dG"))) %>%
    bind_cols(Tm_1st_half = map_dbl(thermodynamics$First, pluck("Tm"))) %>%
    bind_cols(dG_2nd_half = map_dbl(thermodynamics$Second, pluck("dG"))) %>%
    bind_cols(Tm_2nd_half = map_dbl(thermodynamics$Second, pluck("Tm"))) %>%
    mutate(
      passed_a_comp       = passed_a_comp(rev_comp),
      passed_c_comp       = passed_c_comp(rev_comp),
      passed_a_stack      = passed_a_stack(rev_comp),
      passed_c_stack      = passed_c_stack(rev_comp),
      passed_c_spec_stack = map_lgl(rev_comp, passed_c_spec_stack)
    ) %>%
    return()
}


generate_exploratory_plots <- function(param_vector){
  map(param_vector,
      ~ candidate_probes_annotated %>%
        ggplot(aes_string(x = "start", y = .x, colour = .x)) +
        geom_point() +
        labs(title = paste0(.x, " values along the target RNA sequence")) + 
        scale_colour_viridis_c() +
        theme_minimal() +
        theme(legend.position = "none")
  ) %>%
    wrap_plots()
}


filter_candidate_probes <- function(annotated_df){
  filt_df <- annotated_df %>%
    dplyr::filter(dG > dG_range[[1]] & dG < dG_range[[2]]) %>%
    dplyr::filter(Tm > Tm_range[[1]] & Tm < Tm_range[[2]]) %>%
    dplyr::filter(GC_content > GC_range[[1]] & GC_content < GC_range[[2]])
  if (pass_a_comp == TRUE) {
    filt_df <- dplyr::filter(filt_df, passed_a_comp == TRUE)
  } else {
    filt_df <- filt_df
  }
  if (pass_c_comp == TRUE) {
    filt_df <- dplyr::filter(filt_df, passed_c_comp == TRUE)
  } else {
    filt_df <- filt_df
  }
  if (pass_a_stack == TRUE) {
    filt_df <- dplyr::filter(filt_df, passed_a_stack == TRUE)
  } else {
    filt_df <- filt_df
  }
  if (pass_c_stack == TRUE) {
    filt_df <- dplyr::filter(filt_df, passed_c_stack == TRUE)
  } else {
    filt_df <- filt_df
  }
  if (pass_c_spec_stack == TRUE) {
    filt_df <- dplyr::filter(filt_df, passed_c_spec_stack == TRUE)
  } else {
    filt_df <- filt_df
  }
  filt_df <- filt_df %>%
    mutate(dG_deviation = abs(dG - target_dG)) %>%
    mutate(dG_deviation_halves = abs(target_dG_halves - dG_1st_half) + abs(target_dG_halves - dG_2nd_half))
  
  if(nrow(filt_df) == 0){
    print("ERROR: No candidate passed the filtering. Adjust the filtering parameters!")
  } else {
    print(paste(nrow(filt_df), "potentially overlapping probes passed filtering!"))
    return(filt_df) 
  }
}


run_blastn_short <- function(df, db, tx2gene, cores = n_threads){
  # col_names <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
  tx2gene <- read_csv(tx2gene_csv)
  # predict(db, seqs, 
  #         BLAST_args = paste0("-task blastn-short ",
  #                             "-dust no ",
  #                             "-soft_masking false ",
  #                             "-num_threads 4 ")) %>%
  #   setNames(col_names) %>%
  #   filter(bitscore > 32) %>%
  #   filter(evalue < 50) %>%
  #   separate(col = sseqid, into = c("sseqid", "is_intron"), sep = "_") %>%
  #   mutate(sseqid = str_replace(sseqid, "\\..*", "")) %>%
  #   left_join(tx2gene, by = c("sseqid" = "tx_id")) %>%
  #   dplyr::select(qseqid, sseqid, gene_id, gene_name, everything()) %>%
  #   mutate(gene_id = if_else(str_detect(sseqid, "FBgn|ENSG"), sseqid, gene_id))

  candidates <- df
  seqs <- DNAStringSet(candidates$target_sequence)
  names(seqs) <- candidates$unique_id

  temp_fasta_path <- file.path(tempdir(), "query.fasta")
  Biostrings::writeXStringSet(seqs, temp_fasta_path, format = "fasta")

  blast_output <- blast_nucleotide_to_nucleotide(
    query = temp_fasta_path,
    subject = db,
    evalue = 2,
    task = "blastn-short",
    cores = cores,
    output.path = tempdir()
  )

  blast_output <- blast_output %>%
    filter(bit_score > 32) %>%
    filter(evalue < 50) %>%
    tidyr::separate(col = subject_id, into = c("subject_id", "is_intron"), sep = "_") %>%
    mutate(subject_id = str_replace(subject_id, "\\..*", "")) %>%
    left_join(tx2gene, by = c("subject_id" = "tx_id")) %>%
    dplyr::select(query_id, subject_id, gene_id, gene_name, everything()) %>%
    mutate(gene_id = if_else(str_detect(subject_id, "FBgn|ENSG|ENSMUSG"), subject_id, gene_id))

  return(blast_output)
}

# run_blastn_short_Hsap <- function(db, seqs, tx2gene_csv){
#   col_names <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
#   tx2gene <- read_csv(tx2gene_csv)
#   predict(db, seqs, 
#           BLAST_args = paste0("-task blastn-short ",
#                               "-dust no ",
#                               "-soft_masking false ",
#                               "-num_threads 4 ")) %>%
#     setNames(col_names) %>%
#     filter(bitscore > 32) %>%
#     filter(evalue < 50) %>%
#     separate(col = sseqid, into = c("sseqid", "is_intron"), sep = "_") %>%
#     mutate(sseqid = str_replace(sseqid, "\\..*", "")) %>%
#     left_join(tx2gene, by = c("sseqid" = "tx_id")) %>%
#     dplyr::select(qseqid, sseqid, gene_id, gene_name, everything()) %>%
#     mutate(gene_id = if_else(str_detect(sseqid, "ENSG"), sseqid, gene_id))
# }


summarise_blast_output <- function(blast_output, allowOverlapBreakRegion){
  if(allowOverlapBreakRegion == TRUE){
    blast_summary <- blast_output %>%
      filter(!(alig_length <= 17 & q_start < 25 & q_end > 28)) %>%
      dplyr::select(query_id, gene_id) %>%
      distinct() %>%
      group_by(query_id) %>%
      summarise(n_matches = n(),
                matches = paste(gene_id, collapse = "; ")) %>%
      ungroup()
    
  } else if(allowOverlapBreakRegion == FALSE) {
    blast_summary <- blast_output %>%
      dplyr::select(query_id, gene_id) %>%
      distinct() %>%
      group_by(query_id) %>%
      summarise(n_matches = n(),
                matches = paste(gene_id, collapse = "; ")) %>%
      ungroup()
  }
  
}

screen_with_blast_summary <- function(candidate_probes_filtered, max_blast_matches, blast_summary) {
  ## ---- Screen candidate probes
  screened_ids <- blast_summary %>%
    filter(n_matches <= max_blast_matches) %>%
    pull(query_id)
  candidate_probes_blast_screened <- candidate_probes_filtered %>%
    filter(unique_id %in% screened_ids)

  if (nrow(candidate_probes_blast_screened) == 0) {
    print("ERROR: No probes passed BLAST screen. Relax the filtering parameters!")
  } else {
    ## ---- Merge overlapping intervals for the inspection plot
    candidate_probes_blast_screened <- mutate(candidate_probes_blast_screened, chrom = "valr")
    merge_df <- bed_merge(candidate_probes_blast_screened) %>%
      mutate(region_length = end - start + 1) %>%
      mutate(type = case_when(
        region_length >= 2 * oligo_length ~ 1,
        region_length < 2 * oligo_length ~ 0,
        TRUE ~ 0
      ))

    output <- list(candidate_probes_blast_screened, merge_df)
    return(output)
  }
}

plot_inspection <- function(df, colour_param){
  plot <- ggplot(df) +
    geom_linerange(
      aes(x = " ",
          ymin = start, ymax = end,
          colour = !!sym(colour_param)),
      position = position_dodge2(width = 0.75)
    ) +
    coord_flip() + 
    labs(title = "Probe positions along the target RNA sequence",
         subtitle = paste0(nrow(df), " potentially overlapping probes"),
         x = "", y = "") +
    scale_colour_viridis_c(option = "inferno") + 
    scale_y_continuous(limits = c(0, max(candidate_probes$end))) +
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(candidate_probes$end), fill = "gray80") + 
    theme_gray() +
    theme(legend.position = "bottom")
  return(plot)
}

plot_final_probes <- function(df, colour_param){
  plot <- ggplot(df) +
    geom_linerange(
      aes(x = " ",
          ymin = start, ymax = end,
          colour = !!sym(colour_param)),
      position = position_dodge2(width = 0.75), size = 2
    ) +
    coord_flip() + 
    labs(title = "Probe positions along the target RNA sequence",
         subtitle = paste0(nrow(df), " pairs of HCR probes"),
         x = "", y = "") +
    scale_colour_viridis_c(option = "inferno") + 
    scale_y_continuous(limits = c(0, max(candidate_probes$end))) +
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(candidate_probes$end), fill = "gray80") + 
    theme_gray() +
    theme(legend.position = "bottom")
  return(plot)
}


get_probes_shortregions <- function(candidate_probes_screened, merge_df){
  merge_df_short <- merge_df %>%
    filter(region_length < 2 * oligo_length)
  probes_short_regions <- candidate_probes_screened %>%
    bed_intersect(merge_df_short) %>%
    group_by(start.y, end.y) %>%
    slice_min(dG_deviation.x, n = 1) %>%
    slice_min(dG_deviation_halves.x, n = 1) %>%
    slice_min(start.x, n = 1) %>% ungroup() %>%
    dplyr::select(-contains(".y"), -contains(".overlap")) %>%
    dplyr::rename_with(~ str_replace_all(.x, ".x", "")) %>%
    return()
}


test_distr_combn <- function(probe_midpoints, n_fit, probe_spacing){
  # Distribute non-overlapping "n_fit" number of probes in long regions
  # Find minimal nearest distance between probe combinations (should > 52)
  # Find sum of dG deviation (52nt) - choose minimum
  # Find sum of dG deviation of probe halves (25nt) - choose minimum next in case of ties
  
  # combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
  # future_map_dfr(combinations,
  #         function(combination){
  #           tibble(
  #             midpoint_index = probe_midpoints[combination] %>% paste(collapse = "|"),
  #             nearest_dist = probe_midpoints[combination] %>% dist() %>% min(),
  #             sum_probe_dG_deviations = probe_dG_deviations[combination] %>% sum(),
  #             sum_probe_dG_deviations_halves = probe_dG_deviations_halves[combination] %>% sum()
  #           )
  #         }) %>%
  #   filter(nearest_dist >= 52 + probe_spacing) %>%
  #   slice_min(sum_probe_dG_deviations) %>%
  #   slice_min(sum_probe_dG_deviations_halves) %>%
  #   slice_head(n = 1)
  
  combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
  spread_out_lists <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
  if(length(spread_out_lists) == 0) {
    tibble(
      midpoint_index = character(),
      sum_probe_dG_deviations = numeric(),
      sum_probe_dG_deviations_halves = numeric()
    ) %>% return()
  } else {
    future_map_dfr(spread_out_lists,
                   function(spread_out_list){
                     tibble(
                       midpoint_index = probe_midpoints[spread_out_list] %>% paste(collapse = "|"),
                       sum_probe_dG_deviations = probe_dG_deviations[spread_out_list] %>% sum(),
                       sum_probe_dG_deviations_halves = probe_dG_deviations_halves[spread_out_list] %>% sum()
                     )
                   }) %>%
      slice_min(sum_probe_dG_deviations) %>%
      slice_min(sum_probe_dG_deviations_halves) %>%
      slice_head(n = 1)
  }
}




distribute_probes_longregions <- function(candidate_probes_screened, merge_df, probe_spacing){
  merge_df_long <- merge_df %>%
    filter(region_length >= 2 * oligo_length) %>%
    mutate(longregion_maxfit = floor(region_length / oligo_length))
  probes_long_regions <- candidate_probes_screened %>%
    bed_intersect(merge_df_long) %>%
    mutate(longregion_id = paste0(start.y, ";", end.y)) %>%
    dplyr::select(contains(".x"), longregion_id,  "longregion_maxfit" = longregion_maxfit.y, -contains(".overlap")) %>%
    dplyr::rename_with(~ str_replace_all(.x, "\\.x", ""))
  
  future_map_dfr(
    unique(probes_long_regions$longregion_id) %>% set_names(),
    function(longregion){
      ## subset by "longregion_id"
      subset <- probes_long_regions %>% filter(longregion_id == longregion)
      subset_maxfit <- pull(subset, longregion_maxfit) %>% unique()
      probe_midpoints <- subset %>% pull(centre_pos)
      probe_dG_deviations <- subset %>% pull(dG_deviation)
      probe_dG_deviations_halves <- subset %>% pull(dG_deviation_halves)
      
      ## iterate probe distributions over 1:maxfit number of probes
      future_map_dfr(
        set_names(1:subset_maxfit),
        function(n_fit){
          test_distr_combn(probe_midpoints = probe_midpoints, n_fit = n_fit, probe_spacing = probe_spacing)
        }, .id = "n_fit"
      ) %>%
        mutate(longregion_maxfit = subset_maxfit) %>%
        dplyr::select(longregion_maxfit, everything())
    }, .id = "longregion_id"
  )
}

distribute_overlapping_probes <- function(candidate_probes_screened, merge_df, probe_spacing){
  
  if(sum(merge_df$type) == 0) {
    ## Short regions: can only fit one within overlapping regions
    probes_short_regions <- get_probes_shortregions(candidate_probes_screened, merge_df)
    
    candidate_probes_final <- probes_short_regions
    return(candidate_probes_final)
    
  } else if(sum(merge_df$type) > 0 ) {
    ## Short regions: can only fit one within overlapping regions
    probes_short_regions <- get_probes_shortregions(candidate_probes_screened, merge_df)
    
    ## Long regions: can fit > 1 within overlapping regions -> distribute and select minimum deviation from target Gibbs energy
    probes_long_distributed <- distribute_probes_longregions(
      candidate_probes_screened = candidate_probes_screened,
      merge_df = merge_df, 
      probe_spacing = probe_spacing)
    
    ### Example for fitting any many probes as possible
    probes_long_tokeep <- probes_long_distributed %>%
      group_by(longregion_id) %>%
      slice_max(n_fit) %>%
      ungroup() %>%
      pull(midpoint_index) %>%
      paste(collapse = "|") %>%
      str_split(pattern = "\\|") %>% pluck(1) %>% as.numeric()
    
    probes_long_regions <- candidate_probes_screened %>%
      filter(centre_pos %in% probes_long_tokeep)
    
    ## Combine probes from short and long regions
    candidate_probes_final <- bind_rows(probes_short_regions, probes_long_regions)
    return(candidate_probes_final)
    
  } else {
    print("No probes found. Try again!")
  }
  
  
}

attach_hcr_initiatior <- function(candidate_probes_final, b_identifier){
  # Mutate 1st/2nd split sequences
  df <- candidate_probes_final %>%
    mutate(rev_comp_1st_half_id = seq(from = 1, to = 2 * nrow(candidate_probes_final) - 1, by = 2) %>% str_pad(2, pad = "0")) %>%
    mutate(rev_comp_1st_half = str_sub(rev_comp, start = 1, end = (oligo_length / 2) - 1) %>% tolower()) %>%
    mutate(rev_comp_2nd_half_id = seq(from = 2, to = 2 * nrow(candidate_probes_final), by = 2) %>% str_pad(2, pad = "0")) %>%
    mutate(rev_comp_2nd_half = str_sub(rev_comp, start = (oligo_length / 2) + 2, end = oligo_length) %>% tolower())
  # Add split initiators depending on the "b_identifier" argument
  if(b_identifier == "B1") {
    initiator_I1_a <- "gAggAgggCAgCAAACgg"
    initiator_I1_b <- "gAAgAgTCTTCCTTTACg"
    B1_spacera     <- "AA"
    B1_spacerb     <- "TA"
    
    attachToEndOf1stHalf <- paste0(B1_spacerb, initiator_I1_b)
    attachToStartOf2ndHalf <- paste0(initiator_I1_a, B1_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B2") {
    initiator_I2_a <- "CCTCgTAAATCCTCATCA"
    initiator_I2_b <- "ATCATCCAgTAAACCgCC"
    B2_spacera     <- "AA"
    B2_spacerb     <- "AA"
    
    attachToEndOf1stHalf <- paste0(B2_spacerb, initiator_I2_b)
    attachToStartOf2ndHalf <- paste0(initiator_I2_a, B2_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B3") {
    initiator_I3_a <- "gTCCCTgCCTCTATATCT"
    initiator_I3_b <- "CCACTCAACTTTAACCCg"
    B3_spacera     <- "TT"
    B3_spacerb     <- "TT"
    
    attachToEndOf1stHalf <- paste0(B3_spacerb, initiator_I3_b)
    attachToStartOf2ndHalf <- paste0(initiator_I3_a, B3_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B4") {
    initiator_I4_a <- "CCTCAACCTACCTCCAAC"
    initiator_I4_b <- "TCTCACCATATTCgCTTC"
    B4_spacera     <- "AA"
    B4_spacerb     <- "AT"
    
    attachToEndOf1stHalf <- paste0(B4_spacerb, initiator_I4_b)
    attachToStartOf2ndHalf <- paste0(initiator_I4_a, B4_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B5") {
    initiator_I5_a <- "CTCACTCCCAATCTCTAT"
    initiator_I5_b <- "CTACCCTACAAATCCAAT"
    B5_spacera     <- "AA"
    B5_spacerb     <- "AA"
    
    attachToEndOf1stHalf <- paste0(B5_spacerb, initiator_I5_b)
    attachToStartOf2ndHalf <- paste0(initiator_I5_a, B5_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else {
    print("Please select one of B1, B2, B3, B4, or B5 hairpin identifiers")
    
  }
}


attach_hcr_initiatior_standalone <- function(data, first_half_col, second_half_col, b_identifier, target_name){
  # Mutate 1st/2nd split sequences
  df <- data %>%
    dplyr::rename("rev_comp_1st_half" = !!sym(first_half_col)) %>%
    dplyr::rename("rev_comp_2nd_half" = !!sym(second_half_col)) %>%
    mutate(rev_comp_1st_half = tolower(rev_comp_1st_half)) %>%
    mutate(rev_comp_2nd_half = tolower(rev_comp_2nd_half)) %>%
    mutate(rev_comp_1st_half_id = seq(from = 1, to = 2 * nrow(data) - 1, by = 2) %>% str_pad(2, pad = "0")) %>%
    mutate(rev_comp_2nd_half_id = seq(from = 2, to = 2 * nrow(data), by = 2) %>% str_pad(2, pad = "0")) %>%
    dplyr::select(rev_comp_1st_half_id, rev_comp_1st_half, rev_comp_2nd_half_id, rev_comp_2nd_half)
  # Add split initiators depending on the "b_identifier" argument
  if(b_identifier == "B1") {
    initiator_I1_a <- "gAggAgggCAgCAAACgg"
    initiator_I1_b <- "gAAgAgTCTTCCTTTACg"
    B1_spacera     <- "AA"
    B1_spacerb     <- "TA"
    
    attachToEndOf1stHalf <- paste0(B1_spacerb, initiator_I1_b)
    attachToStartOf2ndHalf <- paste0(initiator_I1_a, B1_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B2") {
    initiator_I2_a <- "CCTCgTAAATCCTCATCA"
    initiator_I2_b <- "ATCATCCAgTAAACCgCC"
    B2_spacera     <- "AA"
    B2_spacerb     <- "AA"
    
    attachToEndOf1stHalf <- paste0(B2_spacerb, initiator_I2_b)
    attachToStartOf2ndHalf <- paste0(initiator_I2_a, B2_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B3") {
    initiator_I3_a <- "gTCCCTgCCTCTATATCT"
    initiator_I3_b <- "CCACTCAACTTTAACCCg"
    B3_spacera     <- "TT"
    B3_spacerb     <- "TT"
    
    attachToEndOf1stHalf <- paste0(B3_spacerb, initiator_I3_b)
    attachToStartOf2ndHalf <- paste0(initiator_I3_a, B3_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B4") {
    initiator_I4_a <- "CCTCAACCTACCTCCAAC"
    initiator_I4_b <- "TCTCACCATATTCgCTTC"
    B4_spacera     <- "AA"
    B4_spacerb     <- "AT"
    
    attachToEndOf1stHalf <- paste0(B4_spacerb, initiator_I4_b)
    attachToStartOf2ndHalf <- paste0(initiator_I4_a, B4_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else if(b_identifier == "B5") {
    initiator_I5_a <- "CTCACTCCCAATCTCTAT"
    initiator_I5_b <- "CTACCCTACAAATCCAAT"
    B5_spacera     <- "AA"
    B5_spacerb     <- "AA"
    
    attachToEndOf1stHalf <- paste0(B5_spacerb, initiator_I5_b)
    attachToStartOf2ndHalf <- paste0(initiator_I5_a, B5_spacera)
    
    df %>%
      mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
      mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
      mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
      mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
      return()
    
  } else {
    print("Please select one of B1, B2, B3, B4, or B5 hairpin identifiers")
    
  }
}



get_final_probe_sequences <- function(probe_details){
  first <- probe_details %>%
    dplyr::select(contains("rev_comp_1st_half")) %>%
    dplyr::select(contains("probe_id"), contains("B")) %>%
    setNames(c("probe_id", "sequence"))
  second <- probe_details %>%
    dplyr::select(contains("rev_comp_2nd_half")) %>%
    dplyr::select(contains("probe_id"), contains("B")) %>%
    setNames(c("probe_id", "sequence"))
  bind_rows(first, second) %>% return()
}

save_params <- function(output){
  tibble(
    creation_date                = Sys.time(),
    n_threads                    = n_threads, 
    target_name                  = target_name,
    HCRv3_B_version              = b_identifier,
    oligo_length                 = oligo_length,
    HCR_probe_hybe_region_length = paste0(oligo_length / 2 - 1, " each"),
    probe_spacing                = probe_spacing,
    total_probe_pairs_generated  = nrow(candidate_probes_final),
    target_dG                    = target_dG,
    target_dG_halves             = target_dG_halves,
    dG_range                     = paste(dG_range, collapse = ":"),
    Tm_range                     = paste(Tm_range, collapse = ":"),
    GC_range                     = paste(GC_range, collapse = ":"),
    pass_a_comp                  = as.character(pass_a_comp),
    pass_c_comp                  = as.character(pass_c_comp),
    pass_a_stack                 = as.character(pass_a_stack),
    pass_c_stack                 = as.character(pass_c_stack),
    pass_c_spec_stack            = as.character(pass_c_spec_stack),
    max_blast_matches            = max_blast_matches,
    allowOverlapBreakRegion      = as.character(allowOverlapBreakRegion),
    BLAST_db_used                = blast_file,
    target_sequence_total_length = target_sequence_total_length,
    input_fasta                  = fasta_file,
    target_sequence_used         = target
  ) %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(everything(), names_to = "parameter", values_to = "value") %>%
    write_tsv(output, col_names = FALSE)
}


export_outputs <- function(output_dir, probe_details, probes, blast_output) {
  ## ---- Probe details
  write_csv(probe_details, paste0(output_dir, target_name, "_", "HCR", b_identifier, "_details.csv"))

  ## ---- Probe sequences
  write_csv(probes, paste0(output_dir, target_name, "_", "HCR", b_identifier, "_probes.csv"))

  ## ---- Raw Blast output
  write_csv(blast_output, paste0(output_dir, target_name, "_", "HCR", b_identifier, "_rawblastoutput.csv"))

  ## ---- Probe design parameters
  save_params(paste0(output_dir, target_name, "_", "HCR", b_identifier, "_params.txt"))
}





















