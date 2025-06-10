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
    # CORRECTED LINE: Use the 'target_seq' argument passed to the function
    sequence <- str_sub(target_seq, start, end) 
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
    # This return statement will cause issues in Shiny. It's better to stop with an error.
    stop("Candidate probes could not be generated. Check your FASTA file for sequence content and length.")
  } else {
    return(df)
  }
}


inspect_probe_prep <- function(df, sequence, name){
  ggplot(df,
         aes(x = " ", 
             ymin = start, ymax = end, 
             colour = start)) + 
    geom_linerange(position = position_dodge2(width = 0.75)) + 
    geom_hline(yintercept = as_tibble(str_locate_all(sequence, pattern = "N")[[1]]) %>% pull(start),
               colour = "gray60") + 
    labs(title = paste0(name, 
                        ": all sliding probes - ",
                        nchar(sequence), " bp in length"),
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

# # Corrected generate_exploratory_plots function
# generate_exploratory_plots <- function(annotated_df, param_vector){
#   map(param_vector,
#       ~ annotated_df %>%
#         ggplot(aes_string(x = "start", y = .x, colour = .x)) +
#         geom_point() +
#         labs(title = paste0(.x, " along the target sequence")) + 
#         scale_colour_viridis_c() +
#         theme_minimal() +
#         theme(legend.position = "none")
#   ) %>%
#     wrap_plots()
# }

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

# 1. Corrected `run_blastn_short` function
# ------------------------------------------
# Issue: The original function contained a reference to an undefined variable 
# ('tx2gene_csv') which would cause a runtime error.
# Fix: The 'read_csv' call has been corrected to use the 'tx2gene' argument,
# which correctly passes the file path provided by the user in the UI.

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
    evalue = 10,
    task = "blastn-short",
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


# 2. Corrected `summarise_blast_output` function
# -----------------------------------------------
# Issue: The original function contained a syntactically incorrect and duplicated line
# ('blast_summary <- blast_output %>%') in the 'else if' block, which would
# prevent the code from executing.
# Fix: The duplicated line has been removed. The logic has also been slightly
# streamlined to improve readability without altering the scientific outcome.

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

plot_inspection <- function(df, all_candidates, colour_param){
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
    scale_y_continuous(limits = c(0, max(all_candidates$end))) +
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(all_candidates$end), fill = "gray80") + 
    theme_gray() +
    theme(legend.position = "bottom")
  return(plot)
}

# Corrected plot_final_probes function
plot_final_probes <- function(df, all_candidates, colour_param){
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
    # Corrected to use the 'all_candidates' argument for correct y-axis scaling
    scale_y_continuous(limits = c(0, max(all_candidates$end))) +
    annotate(geom = "rect", xmin = 0.5, xmax = 0.54, ymin = 0, ymax = max(all_candidates$end), fill = "gray80") + 
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
  
  combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
  combinations_nonoverlapping <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
  if(length(combinations_nonoverlapping) == 0) {
    tibble(
      midpoint_index = character(),
      sum_probe_dG_deviations = numeric(),
      sum_probe_dG_deviations_halves = numeric()
    ) %>% return()
  } else {
    future_map_dfr(combinations_nonoverlapping,
                   function(combination){
                     tibble(
                       midpoint_index = probe_midpoints[combination] %>% paste(collapse = "|"),
                       sum_probe_dG_deviations = probe_dG_deviations[combination] %>% sum(),
                       sum_probe_dG_deviations_halves = probe_dG_deviations_halves[combination] %>% sum()
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



## Overlapping probe distributor 
distribute_overlapping_probes <- function(candidate_probes_screened, merge_df, probe_spacing){
  
  merged_regions <- merge_df %>%
    mutate(
      potential_max_fit = floor((region_length + probe_spacing) / (oligo_length + probe_spacing))
    )
  
  ### Distribute probes
  if(nrow(merged_regions) == 0){
    ### Error message when no probes are found
    print("No probes found. Check blast outputs for any issues.")
    
  } else {
    ### Loop around potential max fit of probes per merged regions
    potential_max_fits <- unique(merged_regions$potential_max_fit) %>% sort()
    distributed_probes_list <- list()
    
    for (i in potential_max_fits){
      ### single fit
      if(i == 1){
        merged_subregions <- filter(merged_regions, potential_max_fit == i)
        probes_subregions <- candidate_probes_screened %>%
          bed_intersect(merged_subregions) %>%
          group_by(start.y, end.y) %>%
          slice_min(dG_deviation.x, n = 1) %>%
          slice_min(dG_deviation_halves.x, n = 1) %>%
          slice_min(start.x, n = 1) %>% ungroup() %>%
          dplyr::select(-contains(".y"), -contains(".overlap")) %>%
          dplyr::rename_with(~ str_replace_all(.x, ".x", ""))
        distributed_probes_list[[i]] <- probes_subregions
      } else if (i > 1 & i < 5) {
        merged_subregions <- filter(merged_regions, potential_max_fit == i)
        probes_subregions_overlapping <- candidate_probes_screened %>%
          bed_intersect(merged_subregions) %>%
          mutate(longregion_id = paste0(start.y, ";", end.y)) %>%
          dplyr::select(contains(".x"), longregion_id, "potential_max_fit" = potential_max_fit.y, -contains(".overlap")) %>%
          dplyr::rename_with(~ str_replace_all(.x, "\\.x", ""))
        
        probes_subregions_nonoverlapping <- unique(probes_subregions_overlapping$longregion_id) %>% 
          purrr::set_names() %>%
          future_map_dfr(
            function(multifit_region){
              ## subset by "longregion_id"
              subset <- filter(probes_subregions_overlapping, longregion_id == multifit_region)
              subset_maxfit <- i
              probe_midpoints <- subset %>% pull(centre_pos)
              probe_dG_deviations <- subset %>% pull(dG_deviation)
              probe_dG_deviations_halves <- subset %>% pull(dG_deviation_halves)
              
              ## iterate probe distributions from potential maxfit number to single fit for each longregion_id
              n_fit <- subset_maxfit
              
              combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
              combinations_nonoverlapping <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
              
              while(length(combinations_nonoverlapping) == 0 & n_fit > 0){
                n_fit <- n_fit - 1
                combinations <- combn(length(probe_midpoints), n_fit, simplify = FALSE)
                combinations_nonoverlapping <- keep(combinations, ~ probe_midpoints[.x] %>% dist() %>% min() >= 52 + probe_spacing)
              }
              
              if(length(combinations_nonoverlapping) == 0 & n_fit == 0){
                tibble(
                  midpoint_index = character(),
                  sum_probe_dG_deviations = numeric(),
                  sum_probe_dG_deviations_halves = numeric()
                ) %>% return()
              } else {
                future_map_dfr(combinations_nonoverlapping,
                               function(combination){
                                 tibble(
                                   midpoint_index = probe_midpoints[combination] %>% paste(collapse = "|"),
                                   sum_probe_dG_deviations = probe_dG_deviations[combination] %>% sum(),
                                   sum_probe_dG_deviations_halves = probe_dG_deviations_halves[combination] %>% sum()
                                 )
                               }) %>%
                  slice_min(sum_probe_dG_deviations) %>%
                  slice_min(sum_probe_dG_deviations_halves) %>%
                  slice_head(n = 1) %>% return()
              }
            }, .id = "longregion_id"
          )
        
        probes_subregions_nonoverlapping_tokeep <- probes_subregions_nonoverlapping %>%
          group_by(longregion_id) %>%
          # slice_max(n_fit) %>%
          ungroup() %>%
          pull(midpoint_index) %>%
          paste(collapse = "|") %>%
          str_split(pattern = "\\|") %>% pluck(1) %>% as.numeric()
        
        probes_subregions <- candidate_probes_screened %>%
          filter(centre_pos %in% probes_subregions_nonoverlapping_tokeep)
        
        distributed_probes_list[[i]] <- probes_subregions
      }
      else if (i >= 5){
        ### greedy algorithm if i>=5. Memory constrains :p 
        merged_subregions <- filter(merged_regions, potential_max_fit == i)
        probes_subregions_overlapping <- candidate_probes_screened %>%
          bed_intersect(merged_subregions) %>%
          mutate(longregion_id = paste0(start.y, ";", end.y)) %>%
          dplyr::select(contains(".x"), longregion_id, "potential_max_fit" = potential_max_fit.y, -contains(".overlap")) %>%
          dplyr::rename_with(~ str_replace_all(.x, "\\.x", ""))
        
        probes_subregions_overlapping_greedy <- probes_subregions_overlapping %>%
          mutate(end = end + probe_spacing) %>%
          arrange(end)
        
        greedy_df <- probes_subregions_overlapping_greedy
        nonoverlapping_probes <- c()
        while (nrow(greedy_df) > 0){
          probe_slice <- slice_head(greedy_df, n = 1) 
          probe_id <- pull(probe_slice, unique_id)
          probe_end <- pull(probe_slice, end)
          nonoverlapping_probes <- c(nonoverlapping_probes, probe_id)
          greedy_df <- filter(greedy_df, start > probe_end)
          
          probes_subregions <- candidate_probes_screened %>%
            filter(unique_id %in% nonoverlapping_probes)
          
          distributed_probes_list[[i]] <- probes_subregions
        }
      }
    }
  }
  output <- purrr::reduce(distributed_probes_list, bind_rows) 
}

cull_excess_pairs <- function(df, max_probe_pairs = max_probe_pairs){
  if(nrow(df) <= max_probe_pairs){
    print("Designed probe pairs are less than the max probe pairs parameter. No change applied.")
    return(df)
  } else {
    df %>%
      slice_min(dG_deviation, n = max_probe_pairs)
  }
} 




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













