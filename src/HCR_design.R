library(tidyverse) 
library(patchwork) 
library(valr) 
library(Biostrings)
library(metablastr)
library(Rcpp)
oligo_length <- 52

sourceCpp("./src/thermo_calc.cpp")
sourceCpp("./src/find_optimal_probe_set_dp_cpp.cpp")

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
# 
# inspect_probe_prep <- function(df, sequence, name){
#   ggplot(df,
#          aes(x = " ", 
#              ymin = start, ymax = end, 
#              colour = start)) + 
#     geom_linerange(position = position_dodge2(width = 0.75)) + 
#     geom_hline(yintercept = as_tibble(str_locate_all(sequence, pattern = "N")[[1]]) %>% pull(start),
#                colour = "gray60") + 
#     labs(title = paste0(name, 
#                         ": all sliding probes - ",
#                         nchar(sequence), " bp in length"),
#          subtitle = "Junction-spanning probes removed (if any present)",
#          x = "", y = "") + 
#     coord_flip() + 
#     scale_colour_viridis_c() +
#     theme_minimal() +
#     theme(legend.position = "none")
# }

inspect_probe_prep <- function(df, sequence, name, line_width = 1){
  color_palette <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", 
                     "#EDC948", "#B07AA1")
  df <- df %>%
    dplyr::arrange(start) %>%
    dplyr::mutate(
      probe_group = as.factor(dplyr::row_number() %% 7),
      y_position = (dplyr::row_number() %% 55)
    )
  
  ggplot(df, aes(colour = probe_group)) + 
    geom_segment(aes(x = start, xend = end, y = y_position, yend = y_position),
                 linewidth = line_width) +
    
    geom_vline(xintercept = as_tibble(str_locate_all(sequence, pattern = "N")[[1]]) %>% pull(start),
               colour = "gray20",
               linetype = "dashed", linewidth = 0.85) + 
    
    labs(title = paste0(name, ": All Sliding Window Candidate Probe Positions"),
         subtitle = paste("Sequence Length:", nchar(sequence), "bp. Junction-spanning probes removed. Randomly coloured."),
         y = NULL,
         x = paste0(name, " RNA sequence position (nt)")) + 
    
    scale_y_continuous(expand = expansion(mult = 0.15)) +
    
    scale_colour_manual(values = color_palette, guide = "none") +
    
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none", 
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
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
  
  # Using a standard lapply loop.
  all_results <- lapply(candidate_probes$target_sequence, function(target_seq) {
    
    # Calculate thermodynamics for the full probe and its two halves
    full_thermo <- calculate_thermodynamics_cpp(target_seq, temperature, Na, oligo_conc)
    
    first_half_seq <- substr(target_seq, 1, nchar(target_seq) / 2 - 1)
    first_half_thermo <- calculate_thermodynamics_cpp(first_half_seq, temperature, Na, oligo_conc)
    
    second_half_seq <- substr(target_seq, nchar(target_seq) / 2 + 2, nchar(target_seq))
    second_half_thermo <- calculate_thermodynamics_cpp(second_half_seq, temperature, Na, oligo_conc)
    
    # Return a list containing all results for the probe
    list(
      Full = full_thermo,
      First = first_half_thermo,
      Second = second_half_thermo
    )
  })
  
  # Restructure the list of results into the format expected by the next function
  output <- list(
    Full = purrr::map(all_results, "Full"),
    First = purrr::map(all_results, "First"),
    Second = purrr::map(all_results, "Second")
  )
  
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
    bind_cols(dG = map_dbl(thermodynamics$Full, function(x) x[["dG"]])) %>%
    bind_cols(Tm = map_dbl(thermodynamics$Full, function(x) x[["Tm"]])) %>%
    bind_cols(dG_1st_half = map_dbl(thermodynamics$First, function(x) x[["dG"]])) %>%
    bind_cols(Tm_1st_half = map_dbl(thermodynamics$First, function(x) x[["Tm"]])) %>%
    bind_cols(dG_2nd_half = map_dbl(thermodynamics$Second, function(x) x[["dG"]])) %>%
    bind_cols(Tm_2nd_half = map_dbl(thermodynamics$Second, function(x) x[["Tm"]])) %>%
    mutate(
      passed_a_comp       = passed_a_comp(rev_comp),
      passed_c_comp       = passed_c_comp(rev_comp),
      passed_a_stack      = passed_a_stack(rev_comp),
      passed_c_stack      = passed_c_stack(rev_comp),
      passed_c_spec_stack = map_lgl(rev_comp, passed_c_spec_stack)
    ) %>%
    return()
}


#' Generate Exploratory Thermodynamic Plots
#'
#' @param annotated_df A data frame of annotated candidate probes.
#' @param param_vector A character vector of parameters to plot (e.g., "Tm", "dG").
#' @param dG_range_in A numeric vector of length 2 for the dG filter range. Optional.
#' @param Tm_range_in A numeric vector of length 2 for the Tm filter range. Optional.
#' @param GC_range_in A numeric vector of length 2 for the GC content filter range. Optional.
#'
#' @return A patchwork object of ggplot plots.
generate_exploratory_plots <- function(annotated_df, param_vector, dG_range_in = NULL, Tm_range_in = NULL, GC_range_in = NULL){
  
  # Create a list to hold the plots
  plot_list <- map(param_vector,
                   ~ {
                     p <- ggplot(annotated_df, aes_string(x = "start", y = .x, colour = .x)) +
                       geom_point(alpha = 0.7) +
                       labs(title = paste0(.x, " along the target sequence")) + 
                       scale_colour_viridis_c() +
                       theme_minimal() +
                       theme(legend.position = "none")
                     
                     # Conditionally add dotted lines for the selected range using annotate()
                     if (.x == "dG" && !is.null(dG_range_in)) {
                       p <- p + 
                         annotate("segment", x = -Inf, xend = Inf, y = dG_range_in[1], yend = dG_range_in[1], colour = "indianred2", linetype = "dotted", linewidth = 0.75) +
                         annotate("segment", x = -Inf, xend = Inf, y = dG_range_in[2], yend = dG_range_in[2], colour = "indianred2", linetype = "dotted", linewidth = 0.75)
                     } else if (.x == "Tm" && !is.null(Tm_range_in)) {
                       p <- p + 
                         annotate("segment", x = -Inf, xend = Inf, y = Tm_range_in[1], yend = Tm_range_in[1], colour = "indianred2", linetype = "dotted", linewidth = 0.75) +
                         annotate("segment", x = -Inf, xend = Inf, y = Tm_range_in[2], yend = Tm_range_in[2], colour = "indianred2", linetype = "dotted", linewidth = 0.75)
                     } else if (.x == "GC_content" && !is.null(GC_range_in)) {
                       p <- p + 
                         annotate("segment", x = -Inf, xend = Inf, y = GC_range_in[1], yend = GC_range_in[1], colour = "indianred2", linetype = "dotted", linewidth = 0.75) +
                         annotate("segment", x = -Inf, xend = Inf, y = GC_range_in[2], yend = GC_range_in[2], colour = "indianred2", linetype = "dotted", linewidth = 0.75)
                     }
                     
                     return(p)
                   }
  )
  
  # Arrange the plots using patchwork
  wrap_plots(plot_list)
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
  # Read the provided tx2gene mapping file
  tx2gene_df <- read_csv(tx2gene)
  
  # Prepare the query sequences in a temporary FASTA file
  candidates <- df
  seqs <- DNAStringSet(candidates$target_sequence)
  names(seqs) <- candidates$unique_id
  
  temp_fasta_path <- file.path(tempdir(), "query.fasta")
  Biostrings::writeXStringSet(seqs, temp_fasta_path, format = "fasta")
  
  # Execute the blastn-short command using the metablastr package
  blast_output <- blast_nucleotide_to_nucleotide(
    query = temp_fasta_path,
    subject = db,
    evalue = 1,
    task = "blastn",
    strand = "plus",
    cores = cores,
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

run_blastn <- function(df, db, tx2gene, cores = n_threads, task){
  # Read the provided tx2gene mapping file
  tx2gene_df <- read_csv(tx2gene)
  
  # Prepare the query sequences in a temporary FASTA file
  candidates <- df
  seqs <- DNAStringSet(candidates$target_sequence)
  names(seqs) <- candidates$unique_id
  
  temp_fasta_path <- file.path(tempdir(), "query.fasta")
  Biostrings::writeXStringSet(seqs, temp_fasta_path, format = "fasta")
  
  # Execute the blastn-short command using the metablastr package
  blast_output <- blast_nucleotide_to_nucleotide(
    query = temp_fasta_path,
    subject = db,
    evalue = 1,
    task = task,
    strand = "plus",
    cores = cores,
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





summarise_blast_output <- function(blast_output, allowOverlapBreakRegion, evalue_cutoff, consider_pseudogene){
  
  # First, filter by pseudogene consideration, which is common to both paths
  filtered_blast <- blast_output %>%
    filter(
      if(consider_pseudogene == FALSE){
        !str_detect(gene_biotype, "pseudogene")
      } else {
        # This condition effectively does nothing, keeping all rows
        TRUE 
      }
    )
  
  # Next, apply the conditional filtering for the break region
  if(allowOverlapBreakRegion == TRUE){
    filtered_blast <- filtered_blast %>%
      # Exclude alignments that are short AND span the 25/25 break region
      filter(!(alig_length <= 25 & q_start < 25 & q_end > 28))
  }
  
  # Finally, perform the summarisation on the fully filtered data
  blast_summary <- filtered_blast %>%
    filter(evalue < evalue_cutoff) %>%
    dplyr::select(query_id, gene_id) %>%
    distinct() %>%
    group_by(query_id) %>%
    summarise(n_matches = n(),
              matches = paste(gene_id, collapse = "; ")) %>%
    ungroup()
  
  return(blast_summary)
}

screen_with_blast_summary <- function(candidate_probes_filtered, max_blast_matches, blast_summary) {
  ## ---- Screen candidate probes
  screen_out_ids <- blast_summary %>%
    filter(n_matches > max_blast_matches) %>%
    pull(query_id)
  candidate_probes_blast_screened <- candidate_probes_filtered %>%
    filter(!(unique_id %in% screen_out_ids))

  if (nrow(candidate_probes_blast_screened) == 0) {
    print("ERROR: No probes passed BLAST screen. Relax the filtering parameters!")
  } else {
    ## ---- Merge overlapping intervals for the inspection plot
    candidate_probes_blast_screened <- mutate(candidate_probes_blast_screened, chrom = "valr")
    merge_df <- bed_merge(candidate_probes_blast_screened) %>%
      mutate(region_length = end - start + 1)
      # mutate(type = case_when(
      #   region_length >= 2 * oligo_length ~ 1,
      #   region_length < 2 * oligo_length ~ 0,
      #   TRUE ~ 0
      # ))

    output <- list(candidate_probes_blast_screened, merge_df)
    return(output)
  }
}

plot_inspection <- function(df, all_candidates, colour_param, target_name){
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
         x = "", y = paste0(target_name, " RNA sequence position (nt)")) +
    scale_colour_viridis_c(option = "inferno") + 
    scale_y_continuous(limits = c(0, max(all_candidates$end))) +
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(all_candidates$end), fill = "gray80") + 
    theme_gray() +
    theme(legend.position = "bottom")
  return(plot)
}

plot_final_probes <- function(df, all_candidates, colour_param, target_name){
  plot <- ggplot(df) +
    # geom_linerange(
    #   aes(x = " ",
    #       ymin = start, ymax = end,
    #       colour = !!sym(colour_param)),
    #   # position = position_dodge2(width = 0.2), 
    #   size = 2
    # ) +
    geom_crossbar(
      aes(x = " ",
          y = (start + end) / 2, # Calculate midpoint for the 'y' aesthetic on the fly.
          ymin = start, ymax = end,
          fill = !!sym(colour_param)),
      fatten = 0, # Setting fatten=0 removes the central horizontal line, making it look like a linerange.
      width = 0.03, # Controls the width of the Tufte-style end caps.
      linewidth = 0.3,
      position = position_jitter(width = 0.1, height = 0) # Jitters the 'x' aesthetic.
    ) +
    coord_flip() + 
    labs(title = "Probe positions along the target RNA sequence",
         subtitle = paste0(nrow(df), " pairs of HCR probes"),
         x = "", y = paste0(target_name, " RNA sequence")) +
    scale_fill_viridis_c(option = "inferno") + 
    scale_y_continuous(limits = c(0, max(all_candidates$end))) +
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(all_candidates$end), fill = "gray80") + 
    theme_gray() +
    theme(legend.position = "bottom")
  return(plot)
}


plot_final_probes <- function(df, all_candidates, colour_param, target_name){
  plot <- ggplot(df) +
    # Draw a background rectangle for context.
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(all_candidates$end), fill = "gray80", colour = "black", linewidth = 0.3) +
    # Use geom_crossbar, which is compatible with position_jitter.
    geom_crossbar(
      aes(x = " ",
          y = (start + end) / 2, # Calculate midpoint for the 'y' aesthetic on the fly.
          ymin = start, ymax = end,
          fill = !!sym(colour_param)),
      fatten = 0, # Setting fatten=0 removes the central horizontal line, making it look like a linerange.
      width = 0.025, # Controls the width of the Tufte-style end caps.
      linewidth = 0.3,
      position = position_jitter(width = 0.2, height = 0) # Jitters the 'x' aesthetic.
    ) +
    coord_flip() + 
    labs(title = "Probe positions along the target RNA sequence",
         subtitle = paste0(nrow(df), " pairs of HCR probes"),
         x = "", y = paste0(target_name, " RNA sequence position (nt)")) +
    scale_fill_viridis_c(option = "inferno") + 
    scale_y_continuous(limits = c(0, max(all_candidates$end))) +
    theme_gray() +
    theme(
      legend.position = "bottom",
      # Hide the original x-axis text/ticks as they are just a single category.
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  return(plot)
}

# get_probes_shortregions <- function(candidate_probes_screened, merge_df){
#   merge_df_short <- merge_df %>%
#     filter(region_length < 2 * oligo_length)
#   probes_short_regions <- candidate_probes_screened %>%
#     bed_intersect(merge_df_short) %>%
#     group_by(start.y, end.y) %>%
#     slice_min(dG_deviation.x, n = 1) %>%
#     slice_min(dG_deviation_halves.x, n = 1) %>%
#     slice_min(start.x, n = 1) %>% ungroup() %>%
#     dplyr::select(-contains(".y"), -contains(".overlap")) %>%
#     dplyr::rename_with(~ str_replace_all(.x, ".x", "")) %>%
#     return()
# }
# 
# 
# test_distr_combn <- function(probe_midpoints, n_fit, probe_spacing){
#   # Distribute non-overlapping "n_fit" number of probes in long regions
#   # Find minimal nearest distance between probe combinations (should > 52)
#   # Find sum of dG deviation (52nt) - choose minimum
#   # Find sum of dG deviation of probe halves (25nt) - choose minimum next in case of ties
#   
#   # combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
#   # future_map_dfr(combinations,
#   #         function(combination){
#   #           tibble(
#   #             midpoint_index = probe_midpoints[combination] %>% paste(collapse = "|"),
#   #             nearest_dist = probe_midpoints[combination] %>% dist() %>% min(),
#   #             sum_probe_dG_deviations = probe_dG_deviations[combination] %>% sum(),
#   #             sum_probe_dG_deviations_halves = probe_dG_deviations_halves[combination] %>% sum()
#   #           )
#   #         }) %>%
#   #   filter(nearest_dist >= 52 + probe_spacing) %>%
#   #   slice_min(sum_probe_dG_deviations) %>%
#   #   slice_min(sum_probe_dG_deviations_halves) %>%
#   #   slice_head(n = 1)
#   
#   combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
#   combinations_nonoverlapping <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
#   if(length(combinations_nonoverlapping) == 0) {
#     tibble(
#       midpoint_index = character(),
#       sum_probe_dG_deviations = numeric(),
#       sum_probe_dG_deviations_halves = numeric()
#     ) %>% return()
#   } else {
#     future_map_dfr(combinations_nonoverlapping,
#                    function(combination){
#                      tibble(
#                        midpoint_index = probe_midpoints[combination] %>% paste(collapse = "|"),
#                        sum_probe_dG_deviations = probe_dG_deviations[combination] %>% sum(),
#                        sum_probe_dG_deviations_halves = probe_dG_deviations_halves[combination] %>% sum()
#                      )
#                    }) %>%
#       slice_min(sum_probe_dG_deviations) %>%
#       slice_min(sum_probe_dG_deviations_halves) %>%
#       slice_head(n = 1)
#   }
# }
# 
# 
# 
# 
# distribute_probes_longregions <- function(candidate_probes_screened, merge_df, probe_spacing){
#   merge_df_long <- merge_df %>%
#     filter(region_length >= 2 * oligo_length) %>%
#     mutate(longregion_maxfit = floor(region_length / oligo_length))
#   probes_long_regions <- candidate_probes_screened %>%
#     bed_intersect(merge_df_long) %>%
#     mutate(longregion_id = paste0(start.y, ";", end.y)) %>%
#     dplyr::select(contains(".x"), longregion_id,  "longregion_maxfit" = longregion_maxfit.y, -contains(".overlap")) %>%
#     dplyr::rename_with(~ str_replace_all(.x, "\\.x", ""))
#   
#   future_map_dfr(
#     unique(probes_long_regions$longregion_id) %>% set_names(),
#     function(longregion){
#       ## subset by "longregion_id"
#       subset <- probes_long_regions %>% filter(longregion_id == longregion)
#       subset_maxfit <- pull(subset, longregion_maxfit) %>% unique()
#       probe_midpoints <- subset %>% pull(centre_pos)
#       probe_dG_deviations <- subset %>% pull(dG_deviation)
#       probe_dG_deviations_halves <- subset %>% pull(dG_deviation_halves)
#       
#       ## iterate probe distributions over 1:maxfit number of probes
#       future_map_dfr(
#         set_names(1:subset_maxfit),
#         function(n_fit){
#           test_distr_combn(probe_midpoints = probe_midpoints, n_fit = n_fit, probe_spacing = probe_spacing)
#         }, .id = "n_fit"
#       ) %>%
#         mutate(longregion_maxfit = subset_maxfit) %>%
#         dplyr::select(longregion_maxfit, everything())
#     }, .id = "longregion_id"
#   )
# }
# 
# # distribute_overlapping_probes <- function(candidate_probes_screened, merge_df, probe_spacing){
# #   
# #   if(sum(merge_df$type) == 0) {
# #     ## Short regions: can only fit one within overlapping regions
# #     probes_short_regions <- get_probes_shortregions(candidate_probes_screened, merge_df)
# #     
# #     candidate_probes_final <- probes_short_regions
# #     return(candidate_probes_final)
# #     
# #   } else if(sum(merge_df$type) > 0 ) {
# #     ## Short regions: can only fit one within overlapping regions
# #     probes_short_regions <- get_probes_shortregions(candidate_probes_screened, merge_df)
# #     
# #     ## Long regions: can fit > 1 within overlapping regions -> distribute and select minimum deviation from target Gibbs energy
# #     probes_long_distributed <- distribute_probes_longregions(
# #       candidate_probes_screened = candidate_probes_screened,
# #       merge_df = merge_df, 
# #       probe_spacing = probe_spacing)
# #     
# #     ### Example for fitting any many probes as possible
# #     probes_long_tokeep <- probes_long_distributed %>%
# #       group_by(longregion_id) %>%
# #       slice_max(n_fit) %>%
# #       ungroup() %>%
# #       pull(midpoint_index) %>%
# #       paste(collapse = "|") %>%
# #       str_split(pattern = "\\|") %>% pluck(1) %>% as.numeric()
# #     
# #     probes_long_regions <- candidate_probes_screened %>%
# #       filter(centre_pos %in% probes_long_tokeep)
# #     
# #     ## Combine probes from short and long regions
# #     candidate_probes_final <- bind_rows(probes_short_regions, probes_long_regions)
# #     return(candidate_probes_final)
# #     
# #   } else {
# #     print("No probes found. Try again!")
# #   }
# # }
# 
# ## Overlapping probe distributor 
# distribute_overlapping_probes <- function(candidate_probes_screened, merge_df, probe_spacing){
#   
#   merged_regions <- merge_df %>%
#     mutate(
#       potential_max_fit = floor((region_length + probe_spacing) / (oligo_length + probe_spacing))
#     )
#   
#   ### Distribute probes
#   if(nrow(merged_regions) == 0){
#     ### Error message when no probes are found
#     print("No probes found. Check blast outputs for any issues.")
#     
#   } else {
#     ### Loop around potential max fit of probes per merged regions
#     potential_max_fits <- unique(merged_regions$potential_max_fit) %>% sort()
#     distributed_probes_list <- list()
#     
#     for (i in potential_max_fits){
#       ### single fit
#       if(i == 1){
#         merged_subregions <- filter(merged_regions, potential_max_fit == i)
#         probes_subregions <- candidate_probes_screened %>%
#           bed_intersect(merged_subregions) %>%
#           group_by(start.y, end.y) %>%
#           slice_min(dG_deviation.x, n = 1) %>%
#           slice_min(dG_deviation_halves.x, n = 1) %>%
#           slice_min(start.x, n = 1) %>% ungroup() %>%
#           dplyr::select(-contains(".y"), -contains(".overlap")) %>%
#           dplyr::rename_with(~ str_replace_all(.x, ".x", ""))
#         distributed_probes_list[[i]] <- probes_subregions
#       } else if (i > 1 & i < 5) {
#         merged_subregions <- filter(merged_regions, potential_max_fit == i)
#         probes_subregions_overlapping <- candidate_probes_screened %>%
#           bed_intersect(merged_subregions) %>%
#           mutate(longregion_id = paste0(start.y, ";", end.y)) %>%
#           dplyr::select(contains(".x"), longregion_id, "potential_max_fit" = potential_max_fit.y, -contains(".overlap")) %>%
#           dplyr::rename_with(~ str_replace_all(.x, "\\.x", ""))
#         
#         probes_subregions_nonoverlapping <- unique(probes_subregions_overlapping$longregion_id) %>% 
#           purrr::set_names() %>%
#           future_map_dfr(
#             function(multifit_region){
#               ## subset by "longregion_id"
#               subset <- filter(probes_subregions_overlapping, longregion_id == multifit_region)
#               subset_maxfit <- i
#               probe_midpoints <- subset %>% pull(centre_pos)
#               probe_dG_deviations <- subset %>% pull(dG_deviation)
#               probe_dG_deviations_halves <- subset %>% pull(dG_deviation_halves)
#               
#               ## iterate probe distributions from potential maxfit number to single fit for each longregion_id
#               n_fit <- subset_maxfit
#               
#               combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
#               combinations_nonoverlapping <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
#               
#               while(length(combinations_nonoverlapping) == 0 & n_fit > 0){
#                 n_fit <- n_fit - 1
#                 combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
#                 combinations_nonoverlapping <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
#               }
#               
#               if(length(combinations_nonoverlapping) == 0 & n_fit == 0){
#                 tibble(
#                   midpoint_index = character(),
#                   sum_probe_dG_deviations = numeric(),
#                   sum_probe_dG_deviations_halves = numeric()
#                 ) %>% return()
#               } else {
#                 future_map_dfr(combinations_nonoverlapping,
#                                function(combination){
#                                  tibble(
#                                    midpoint_index = probe_midpoints[combination] %>% paste(collapse = "|"),
#                                    sum_probe_dG_deviations = probe_dG_deviations[combination] %>% sum(),
#                                    sum_probe_dG_deviations_halves = probe_dG_deviations_halves[combination] %>% sum()
#                                  )
#                                }) %>%
#                   slice_min(sum_probe_dG_deviations) %>%
#                   slice_min(sum_probe_dG_deviations_halves) %>%
#                   slice_head(n = 1) %>% return()
#               }
#             }, .id = "longregion_id"
#           )
#         
#         probes_subregions_nonoverlapping_tokeep <- probes_subregions_nonoverlapping %>%
#           group_by(longregion_id) %>%
#           # slice_max(n_fit) %>%
#           ungroup() %>%
#           pull(midpoint_index) %>%
#           paste(collapse = "|") %>%
#           str_split(pattern = "\\|") %>% pluck(1) %>% as.numeric()
#         
#         probes_subregions <- candidate_probes_screened %>%
#           filter(centre_pos %in% probes_subregions_nonoverlapping_tokeep)
#         
#         distributed_probes_list[[i]] <- probes_subregions
#       }
#       else if (i >= 5){
#         ### greedy algorithm if i>=5. Memory constrains :p 
#         merged_subregions <- filter(merged_regions, potential_max_fit == i)
#         probes_subregions_overlapping <- candidate_probes_screened %>%
#           bed_intersect(merged_subregions) %>%
#           mutate(longregion_id = paste0(start.y, ";", end.y)) %>%
#           dplyr::select(contains(".x"), longregion_id, "potential_max_fit" = potential_max_fit.y, -contains(".overlap")) %>%
#           dplyr::rename_with(~ str_replace_all(.x, "\\.x", ""))
#         
#         probes_subregions_overlapping_greedy <- probes_subregions_overlapping %>%
#           mutate(end = end + probe_spacing) %>%
#           arrange(end)
#         
#         greedy_df <- probes_subregions_overlapping_greedy
#         nonoverlapping_probes <- c()
#         while (nrow(greedy_df) > 0){
#           probe_slice <- slice_head(greedy_df, n = 1) 
#           probe_id <- pull(probe_slice, unique_id)
#           probe_end <- pull(probe_slice, end)
#           nonoverlapping_probes <- c(nonoverlapping_probes, probe_id)
#           greedy_df <- filter(greedy_df, start > probe_end)
#           
#           probes_subregions <- candidate_probes_screened %>%
#             filter(unique_id %in% nonoverlapping_probes)
#           
#           distributed_probes_list[[i]] <- probes_subregions
#         }
#       }
#     }
#   }
#   output <- purrr::reduce(distributed_probes_list, bind_rows) 
# }

#' Find the Optimal Set of Non-Overlapping Probes Prioritizing Probe Count
#'
#' This function uses a modified dynamic programming algorithm to find the
#' non-overlapping probe set that first maximizes the number of probes, and
#' second, maximizes the score (i.e., minimizes dG deviation).
#'
#' @param probes_df A dataframe of candidate probes for a single contiguous region.
#' @param probe_spacing The minimum number of nucleotides between probes.
#' @return A tibble containing the optimal set of non-overlapping probes.

find_optimal_probe_set_dp <- function(probes_df, probe_spacing) {
  # First, filter out any probes that might have bad data (NA coordinates)
  probes_df <- probes_df %>%
    filter(!is.na(start) & !is.na(end))
  
  # --- FIX ---
  # Add a guard clause to handle the simple cases of 0 or 1 probe.
  # If there is only one probe, it is by definition the optimal set.
  # This prevents errors by avoiding the main algorithm.
  if (nrow(probes_df) < 2) {
    return(probes_df)
  }
  
  # For clarity in maximization, we use a positive score.
  probes_df <- probes_df %>%
    mutate(score = 5000 - dG_deviation)
  
  # Sort probes by their end position.
  probes_df <- probes_df %>%
    arrange(end)
  
  n <- nrow(probes_df)
  p <- integer(n)
  ends <- probes_df$end
  
  # Predecessor calculation loop
  for (i in 2:n) {
    required_end <- probes_df$start[i] - probe_spacing
    if (required_end < ends[1]) {
      p[i] <- 0
      next
    }
    low <- 1
    high <- i - 1
    j <- 0
    while (low <= high) {
      mid <- floor((low + high) / 2)
      if (ends[mid] <= required_end) {
        j <- mid
        low <- mid + 1
      } else {
        high <- mid - 1
      }
    }
    p[i] <- j
  }
  
  # Dynamic programming main loop
  dp <- vector("list", n + 1)
  dp[[1]] <- c(num_probes = 0, score = 0)
  
  for (i in 1:n) {
    pred_solution <- dp[[p[i] + 1]]
    solution_with_i <- c(
      num_probes = 1 + pred_solution[["num_probes"]],
      score = probes_df$score[i] + pred_solution[["score"]]
    )
    
    solution_without_i <- dp[[i]]
    
    if (solution_with_i[["num_probes"]] > solution_without_i[["num_probes"]]) {
      dp[[i + 1]] <- solution_with_i
    } else if (solution_with_i[["num_probes"]] == solution_without_i[["num_probes"]] &&
               solution_with_i[["score"]] > solution_without_i[["score"]]) {
      dp[[i + 1]] <- solution_with_i
    } else {
      dp[[i + 1]] <- solution_without_i
    }
  }
  
  # Backtracking to find the selected probes
  selected_indices <- c()
  i <- n
  while (i > 0) {
    pred_solution <- dp[[p[i] + 1]]
    solution_with_i <- c(
      num_probes = 1 + pred_solution[["num_probes"]],
      score = probes_df$score[i] + pred_solution[["score"]]
    )
    
    if (isTRUE(all.equal(dp[[i + 1]], solution_with_i))) {
      selected_indices <- c(i, selected_indices)
      i <- p[i]
    } else {
      i <- i - 1
    }
  }
  
  return(probes_df[selected_indices, ])
}

#' #' Distribute Probes Across Overlapping Regions using Dynamic Programming
#' #'
#' #' This function takes all candidate probes, identifies contiguous regions of
#' #' overlapping probes, and then uses a dynamic programming approach to find the
#' #' optimal non-overlapping set for each region.
#' #'
#' #' @param candidate_probes_screened A dataframe of all candidate probes that passed previous filters.
#' #' @param merge_df A dataframe of merged genomic intervals from `valr::bed_merge`.
#' #' @param probe_spacing The minimum number of nucleotides between probes.
#' #' @return A tibble containing the final, globally optimal set of non-overlapping probes.
#' 
#' distribute_overlapping_probes <- function(candidate_probes_screened, merge_df, probe_spacing){
#'   
#'   # Use bed_intersect to assign each probe to its corresponding merged region
#'   probes_in_regions <- candidate_probes_screened %>%
#'     valr::bed_intersect(merge_df) %>%
#'     # Create a unique ID for each merged region to group by
#'     mutate(region_id = paste(start.y, end.y, sep="-")) %>%
#'     # Clean up column names after the intersect operation
#'     dplyr::select(
#'       -contains(".y"),
#'       -contains(".overlap")
#'     ) %>%
#'     dplyr::rename_with(~ str_replace_all(.x, ".x", ""))
#'   
#'   # Split the dataframe into a list of dataframes, one for each region.
#'   # This prepares the data for parallel processing.
#'   probes_by_region <- probes_in_regions %>%
#'     group_by(region_id) %>%
#'     group_split()
#'   
#'   # Use furrr to apply the DP function to each region in parallel.
#'   # This significantly speeds up the process if there are many distinct regions.
#'   optimal_probes_list <- furrr::future_map(
#'     probes_by_region,
#'     ~ find_optimal_probe_set_dp(.x, probe_spacing)
#'   )
#'   
#'   # Combine the lists of optimal probes from each region into a single dataframe.
#'   final_probes <- bind_rows(optimal_probes_list)
#'   
#'   return(final_probes)
#' }

#' Distribute Probes Across Overlapping Regions using Dynamic Programming
#'
#' This function takes all candidate probes, identifies contiguous regions of
#' overlapping probes, and then uses a dynamic programming approach to find the
#' optimal non-overlapping set for each region. It now calls the fast C++
#' implementation.
#'
#' @param candidate_probes_screened A dataframe of all candidate probes that passed previous filters.
#' @param merge_df A dataframe of merged genomic intervals from `valr::bed_merge`.
#' @param probe_spacing The minimum number of nucleotides between probes.
#' @return A tibble containing the final, globally optimal set of non-overlapping probes.
distribute_overlapping_probes_cpp <- function(candidate_probes_screened, merge_df, probe_spacing){
  
  # Use bed_intersect to assign each probe to its corresponding merged region
  probes_in_regions <- candidate_probes_screened %>%
    valr::bed_intersect(merge_df) %>%
    # Create a unique ID for each merged region to group by
    mutate(region_id = paste(start.y, end.y, sep="-")) %>%
    # Clean up column names after the intersect operation
    dplyr::select(
      -contains(".y"),
      -contains(".overlap")
    ) %>%
    dplyr::rename_with(~ str_replace_all(.x, ".x", ""))
  
  # Split the dataframe into a list of dataframes, one for each region.
  # This prepares the data for parallel processing.
  probes_by_region <- probes_in_regions %>%
    group_by(region_id) %>%
    group_split()
  
  # Not using furrr future functions (incompatible).
  optimal_probes_list <- map(
    probes_by_region,
    ~ find_optimal_probe_set_dp_cpp(.x, probe_spacing) 
  )
  
  # Combine the lists of optimal probes from each region into a single dataframe.
  final_probes <- bind_rows(optimal_probes_list) %>% as_tibble()
  
  return(final_probes)
}

cull_excess_pairs <- function(df, max_probe_pairs = max_probe_pairs){
  if(nrow(df) <= max_probe_pairs){
    print("Designed probe pairs are less than the max probe pairs parameter. No change applied.")
    return(df)
  } else {
    df %>%
      slice_min(dG_deviation, n = max_probe_pairs) %>%
      slice_head(n = 33)
  }
} 


# attach_hcr_initiatior <- function(candidate_probes_final, b_identifier){
#   # Mutate 1st/2nd split sequences
#   df <- candidate_probes_final %>%
#     mutate(rev_comp_1st_half_id = seq(from = 1, to = 2 * nrow(candidate_probes_final) - 1, by = 2) %>% str_pad(2, pad = "0")) %>%
#     mutate(rev_comp_1st_half = str_sub(rev_comp, start = 1, end = (oligo_length / 2) - 1) %>% tolower()) %>%
#     mutate(rev_comp_2nd_half_id = seq(from = 2, to = 2 * nrow(candidate_probes_final), by = 2) %>% str_pad(2, pad = "0")) %>%
#     mutate(rev_comp_2nd_half = str_sub(rev_comp, start = (oligo_length / 2) + 2, end = oligo_length) %>% tolower())
#   # Add split initiators depending on the "b_identifier" argument
#   if(b_identifier == "B1") {
#     initiator_I1_a <- "gAggAgggCAgCAAACgg"
#     initiator_I1_b <- "gAAgAgTCTTCCTTTACg"
#     B1_spacera     <- "AA"
#     B1_spacerb     <- "TA"
#     
#     attachToEndOf1stHalf <- paste0(B1_spacerb, initiator_I1_b)
#     attachToStartOf2ndHalf <- paste0(initiator_I1_a, B1_spacera)
#     
#     df %>%
#       mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
#       mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
#       mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
#       mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
#       return()
#     
#   } else if(b_identifier == "B2") {
#     initiator_I2_a <- "CCTCgTAAATCCTCATCA"
#     initiator_I2_b <- "ATCATCCAgTAAACCgCC"
#     B2_spacera     <- "AA"
#     B2_spacerb     <- "AA"
#     
#     attachToEndOf1stHalf <- paste0(B2_spacerb, initiator_I2_b)
#     attachToStartOf2ndHalf <- paste0(initiator_I2_a, B2_spacera)
#     
#     df %>%
#       mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
#       mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
#       mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
#       mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
#       return()
#     
#   } else if(b_identifier == "B3") {
#     initiator_I3_a <- "gTCCCTgCCTCTATATCT"
#     initiator_I3_b <- "CCACTCAACTTTAACCCg"
#     B3_spacera     <- "TT"
#     B3_spacerb     <- "TT"
#     
#     attachToEndOf1stHalf <- paste0(B3_spacerb, initiator_I3_b)
#     attachToStartOf2ndHalf <- paste0(initiator_I3_a, B3_spacera)
#     
#     df %>%
#       mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
#       mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
#       mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
#       mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
#       return()
#     
#   } else if(b_identifier == "B4") {
#     initiator_I4_a <- "CCTCAACCTACCTCCAAC"
#     initiator_I4_b <- "TCTCACCATATTCgCTTC"
#     B4_spacera     <- "AA"
#     B4_spacerb     <- "AT"
#     
#     attachToEndOf1stHalf <- paste0(B4_spacerb, initiator_I4_b)
#     attachToStartOf2ndHalf <- paste0(initiator_I4_a, B4_spacera)
#     
#     df %>%
#       mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
#       mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
#       mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
#       mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
#       return()
#     
#   } else if(b_identifier == "B5") {
#     initiator_I5_a <- "CTCACTCCCAATCTCTAT"
#     initiator_I5_b <- "CTACCCTACAAATCCAAT"
#     B5_spacera     <- "AA"
#     B5_spacerb     <- "AA"
#     
#     attachToEndOf1stHalf <- paste0(B5_spacerb, initiator_I5_b)
#     attachToStartOf2ndHalf <- paste0(initiator_I5_a, B5_spacera)
#     
#     df %>%
#       mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
#       mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
#       mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
#       mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
#       return()
#     
#   } else {
#     print("Please select one of B1, B2, B3, B4, or B5 hairpin identifiers")
#     
#   }
# }


attach_hcr_initiatior <- function(candidate_probes_final, b_identifier, target_name, oligo_length = 52){
  # Mutate 1st/2nd split sequences
  df <- candidate_probes_final %>%
    mutate(rev_comp_1st_half_id = seq(from = 1, to = 2 * nrow(candidate_probes_final) - 1, by = 2) %>% str_pad(2, pad = "0")) %>%
    mutate(rev_comp_1st_half = str_sub(rev_comp, start = 1, end = (oligo_length / 2) - 1) %>% tolower()) %>%
    mutate(rev_comp_2nd_half_id = seq(from = 2, to = 2 * nrow(candidate_probes_final), by = 2) %>% str_pad(2, pad = "0")) %>%
    mutate(rev_comp_2nd_half = str_sub(rev_comp, start = (oligo_length / 2) + 2, end = oligo_length) %>% tolower())
  
  # Add split initiators depending on the "b_identifier" argument
  if(b_identifier == "B1") {
    initiator_a <- "gAggAgggCAgCAAACgg"
    initiator_b <- "gAAgAgTCTTCCTTTACg"
    spacera     <- "AA"
    spacerb     <- "TA"
    
  } else if(b_identifier == "B2") {
    initiator_a <- "CCTCgTAAATCCTCATCA"
    initiator_b <- "ATCATCCAgTAAACCgCC"
    spacera     <- "AA"
    spacerb     <- "AA"
    
  } else if(b_identifier == "B3") {
    initiator_a <- "gTCCCTgCCTCTATATCT"
    initiator_b <- "CCACTCAACTTTAACCCg"
    spacera     <- "TT"
    spacerb     <- "TT"
    
  } else if(b_identifier == "B4") {
    initiator_a <- "CCTCAACCTACCTCCAAC"
    initiator_b <- "TCTCACCATATTCgCTTC"
    spacera     <- "AA"
    spacerb     <- "AT"
    
  } else if(b_identifier == "B5") {
    initiator_a <- "CTCACTCCCAATCTCTAT"
    initiator_b <- "CTACCCTACAAATCCAAT"
    spacera     <- "AA"
    spacerb     <- "AA"
    
  } else if(b_identifier == "B7") {
    initiator_a <- "CTTCAACCTCCACCTACC"
    initiator_b <- "TCCAATCCCTACCCTCAC"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B9") {
    initiator_a <- "CACGTATCTACTCCACTC"
    initiator_b <- "TCAGCACACTCCCAACCC"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B10") {
    initiator_a <- "CCTCAAGATACTCCTCTA"
    initiator_b <- "CCTACTCGACTACCCTAG"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B11") {
    initiator_a <- "CGCTTAGATATCACTCCT"
    initiator_b <- "ACGTCGACCACACTCATC"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B13") {
    initiator_a <- "AGGTAACGCCTTCCTGCT"
    initiator_b <- "TTATGCTCAACATACAAC"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B14") {
    initiator_a <- "AATGTCAATAGCGAGCGA"
    initiator_b <- "CCCTATATTTCTGCACAG"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B15") {
    initiator_a <- "CAGATTAACACACCACAA"
    initiator_b <- "GGTATCTCGAACACTCTC"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else if(b_identifier == "B17") {
    initiator_a <- "CGATTGTTTGTTGTGGAC"
    initiator_b <- "GCATGCTAATCGGATGAG"
    spacera     <- "WW"
    spacerb     <- "WW"
    
  } else {
    print("Please select one of B1-5, B7, B9-10, B13-15, B17 hairpin identifiers")
    
  }
  
  attachToEndOf1stHalf <- paste0(spacerb, initiator_b)
  attachToStartOf2ndHalf <- paste0(initiator_a, spacera)
  
  df %>%
    mutate(rev_comp_1st_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_1st_half_id)) %>%
    mutate(rev_comp_1st_half_B1 = paste0(rev_comp_1st_half, attachToEndOf1stHalf)) %>%
    mutate(rev_comp_2nd_half_probe_id = paste0(target_name, "_HCR", b_identifier, "_", rev_comp_2nd_half_id)) %>%
    mutate(rev_comp_2nd_half_B1 = paste0(attachToStartOf2ndHalf, rev_comp_2nd_half)) %>%
    return()
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
    evalue_cutoff                = evalue_cutoff,
    max_blast_matches            = max_blast_matches,
    allowOverlapBreakRegion      = as.character(allowOverlapBreakRegion),
    consider_pseudogene.         = as.character(consider_pseudogene),
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





















