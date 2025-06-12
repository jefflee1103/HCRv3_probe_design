# server.R
# Defines the server-side logic for the HCR Probe Designer.
# This script contains all the reactive logic that processes user inputs
# and generates the outputs displayed in the UI.

server <- function(input, output, session) {
  
  # --- 1. Reactive Values for Staged Pipeline ---
  # These reactive values store the output of each major step in the pipeline,
  # allowing subsequent steps to access the data without re-computation.
  # They are initialized to NULL and are updated as each step completes.
  
  fasta_info <- reactiveVal(NULL)
  target_seq <- reactiveVal(NULL)
  candidate_probes <- reactiveVal(NULL)
  inspection_plot <- reactiveVal(NULL)
  annotated_candidate_probes <- reactiveVal(NULL)
  exploratory_plots <- reactiveVal(NULL)
  thermo_filtered_probes <- reactiveVal(NULL)
  filtered_probes_plot <- reactiveVal(NULL)
  blast_results <- reactiveVal(NULL)
  blast_summary_data <- reactiveVal(NULL)
  blast_summary_plot <- reactiveVal(NULL)
  blast_screened_probes <- reactiveVal(NULL)
  blast_screened_plot <- reactiveVal(NULL)
  final_probe_set <- reactiveVal(NULL)
  probe_map_plot <- reactiveVal(NULL)
  
  
  # --- 2. Process Sequence on File Upload ---
  # --- This observer triggers automatically when a user uploads a FASTA file.
  observeEvent(input$target_fasta, {
    req(input$target_fasta, input$target_name) # Ensure file and name are provided
    
    tryCatch({
      showNotification("FASTA file detected. Processing sequence...", id = "seq_proc", type = "message")
      
      # Store FASTA file info and read the sequence
      fasta_info(input$target_fasta)
      datapath <- input$target_fasta$datapath
      
      # Read and validate the FASTA file content. Collapse multi-FASTA entries.
      target_fasta_set <- Biostrings::readBStringSet(datapath)
      target_raw <- paste(target_fasta_set, collapse = "N")
      target <- clean_pasted_seq(target_raw)
      target_seq(target) 
      
      # Generate all possible candidate probes from the target sequence
      all_candidates <- generate_candidate_probes(target_seq = target, oligo_length = 52)
      candidate_probes(all_candidates)
      
      # Generate and store an initial plot showing all candidate probes
      insp_plot <- inspect_probe_prep(all_candidates, target, input$target_name)
      inspection_plot(insp_plot)
      
      removeNotification("seq_proc")
      showNotification(paste("Sequence processed. Found", nrow(all_candidates), "candidate probes."), type = "default", duration = 8)
      
      # Automatically switch to the first results tab
      updateTabsetPanel(session, "results_tabs", selected = "Candidate Probes")
      
      # Reset all downstream reactive values to clear results from any previous runs
      annotated_candidate_probes(NULL); exploratory_plots(NULL); thermo_filtered_probes(NULL); 
      filtered_probes_plot(NULL); blast_results(NULL); blast_summary_data(NULL); 
      blast_summary_plot(NULL); blast_screened_probes(NULL); blast_screened_plot(NULL); 
      final_probe_set(NULL); probe_map_plot(NULL)
      
    }, error = function(e) {
      removeNotification("seq_proc")
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  
  # --- 3. Calculate Thermodynamics ---
  # --- This observer triggers when the "Calculate Thermodynamics" button is clicked.
  observeEvent(input$run_thermo_calc, {
    req(candidate_probes(), input$n_threads) # Requires candidate probes to exist
    
    tryCatch({
      showNotification("Calculating thermodynamic parameters (this may take a moment)...", id = "thermo_calc", type = "message")
      
      # Set up parallel processing and run thermodynamic calculations
      thermodynamics <- get_thermodynamic_parameters(candidate_probes(), 37, 0.3, 5e-5)
      
      # Annotate the candidate probes with the calculated thermo properties
      annotated_probes <- annotate_probes(candidate_probes(), thermodynamics)
      annotated_candidate_probes(annotated_probes)

      # Generate and store exploratory plots for the thermo data
      explore_plots <- generate_exploratory_plots(annotated_probes, c("Tm", "dG", "GC_content"))
      exploratory_plots(explore_plots)
      
      removeNotification("thermo_calc")
      showNotification("Thermodynamics calculation complete.", type = "default")
      updateTabsetPanel(session, "results_tabs", selected = "Thermo Filtering")
      
    }, error = function(e) {
      removeNotification("thermo_calc")
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  
  # --- 4. Thermodynamics Filtering ---
  
  # --- Slider and Numeric Input Synchronization
  observeEvent(c(input$Tm_min, input$Tm_max), {
    req(input$Tm_min, input$Tm_max, input$Tm_min < input$Tm_max)
    updateSliderInput(session, "Tm_range", value = c(input$Tm_min, input$Tm_max))
  })
  
  observeEvent(input$Tm_range, {
    updateNumericInput(session, "Tm_min", value = input$Tm_range[1])
    updateNumericInput(session, "Tm_max", value = input$Tm_range[2])
  })
  
  observeEvent(c(input$dG_min, input$dG_max), {
    req(input$dG_min, input$dG_max, input$dG_min < input$dG_max)
    updateSliderInput(session, "dG_range", value = c(input$dG_min, input$dG_max))
  })
  
  observeEvent(input$dG_range, {
    updateNumericInput(session, "dG_min", value = input$dG_range[1])
    updateNumericInput(session, "dG_max", value = input$dG_range[2])
  })
  
  observeEvent(c(input$GC_min, input$GC_max), {
    req(input$GC_min, input$GC_max, input$GC_min < input$GC_max)
    updateSliderInput(session, "GC_range", value = c(input$GC_min, input$GC_max))
  })
  
  observeEvent(input$GC_range, {
    updateNumericInput(session, "GC_min", value = input$GC_range[1])
    updateNumericInput(session, "GC_max", value = input$GC_range[2])
  })
  
  
  # --- This observer triggers when the "Run Thermo Filter" button is clicked.
  observeEvent(input$run_thermo_filter, {
    req(annotated_candidate_probes()) # Requires annotated probes
    showNotification("Filtering probes based on parameters...", type = "message")
    
    # Regenerate exploratory plots to show the user-selected filter ranges as lines
    explore_plots_with_filters <- generate_exploratory_plots(
      annotated_df = annotated_candidate_probes(),
      param_vector = c("Tm", "dG", "GC_content"),
      dG_range_in = input$dG_range,
      Tm_range_in = input$Tm_range,
      GC_range_in = input$GC_range
    )
    exploratory_plots(explore_plots_with_filters) # Update the reactive value for the plot
    
    # Filter the probes based on the user-defined thermodynamic and composition criteria
    filt_df <- annotated_candidate_probes() %>%
      dplyr::filter(
        dG > input$dG_range[1] & dG < input$dG_range[2], 
        Tm > input$Tm_range[1] & Tm < input$Tm_range[2], 
        GC_content > input$GC_range[1] & GC_content < input$GC_range[2]
      )
    if (input$pass_a_comp) { filt_df <- dplyr::filter(filt_df, passed_a_comp == TRUE) }
    if (input$pass_c_comp) { filt_df <- dplyr::filter(filt_df, passed_c_comp == TRUE) }
    if (input$pass_a_stack) { filt_df <- dplyr::filter(filt_df, passed_a_stack == TRUE) }
    if (input$pass_c_stack) { filt_df <- dplyr::filter(filt_df, passed_c_stack == TRUE) }
    if (input$pass_c_spec_stack) { filt_df <- dplyr::filter(filt_df, passed_c_spec_stack == TRUE) }
    
    # Calculate deviation from target Gibbs free energy values
    filt_df <- filt_df %>% 
      mutate(
        dG_deviation = abs(dG - input$target_dG), 
        dG_deviation_halves = abs(input$target_dG_halves - dG_1st_half) + abs(input$target_dG_halves - dG_2nd_half)
      )
    
    # Handle case where no probes pass the filter
    if(nrow(filt_df) == 0){ 
      showNotification("Filtering complete. 0 probes passed.", type = "warning", duration = 8)
      thermo_filtered_probes(NULL)
      filtered_probes_plot(NULL)  
    } else {
      showNotification(paste("Filtering complete.", nrow(filt_df), "probes passed."), type = "default")
      thermo_filtered_probes(filt_df)
      
      # Generate a plot showing the positions of the filtered probes
      filtered_plot <- plot_inspection(df = filt_df, all_candidates = candidate_probes(), colour_param = "Tm", target_name = input$target_name)
      filtered_probes_plot(filtered_plot)
    }
    
    # Reset downstream reactives
    blast_results(NULL); blast_summary_data(NULL); blast_summary_plot(NULL); 
    blast_screened_probes(NULL); blast_screened_plot(NULL); final_probe_set(NULL); 
    probe_map_plot(NULL)
  })
  
  
  # --- 5. Run BLAST ---
  
  # # --- shinyFiles Logic for BLAST 
  # 
  # # Define root volumes for file access. For local use, 'Home' is sufficient.
  # # C: is added for Windows compatibility.
  # volumes <- c(Home = fs::path_home(), "C:" = "C:/")
  # 
  # # Server-side logic for the BLAST database file selection
  # shinyFileChoose(input, "blast_db_path", roots = volumes, session = session, filetypes = c('fa', 'fasta', ''))
  # 
  # # Server-side logic for the tx2gene map file selection
  # shinyFileChoose(input, "tx2gene_path", roots = volumes, session = session, filetypes = c('csv'))
  # 
  # # Reactive expressions to store the selected paths
  # # These reactives will hold the parsed file path strings.
  # blast_db_filepath <- reactive({
  #   req(input$blast_db_path)
  #   parseFilePaths(volumes, input$blast_db_path)$datapath
  # })
  # 
  # tx2gene_filepath <- reactive({
  #   req(input$tx2gene_path)
  #   parseFilePaths(volumes, input$tx2gene_path)$datapath
  # })
  # 
  # # Observers to display the selected paths in the UI
  # # Display selected BLAST DB path
  # output$blast_db_display <- renderText({
  #   req(blast_db_filepath())
  #   blast_db_filepath()
  # })
  # 
  # # Display selected tx2gene map path
  # output$tx2gene_display <- renderText({
  #   req(tx2gene_filepath())
  #   tx2gene_filepath()
  # })
  
  # --- This observer triggers when the "Run BLAST" button is clicked.
  observeEvent(input$run_blast, {
    req(thermo_filtered_probes(), input$blast_db_path, input$tx2gene_path) # Requires filtered probes and paths
    showNotification("Running BLASTn against filtered probes (this may take a long time)...", id = "blast_run", type = "message", duration = NULL)
    
    tryCatch({
      # Run BLAST search using the filtered probes as a query
      blast_output <- run_blastn_short(df = thermo_filtered_probes(), db = input$blast_db_path, tx2gene = input$tx2gene_path, cores = input$n_threads)
      blast_results(blast_output)
      
      # Summarize the raw BLAST output to count off-target hits per probe
      summary_data <- summarise_blast_output(blast_output = blast_output, allowOverlapBreakRegion = input$allowOverlapBreakRegion, evalue_cutoff = input$blast_evalue, consider_pseudogene = input$consider_pseudogene)
      blast_summary_data(summary_data)
      
      # Create a histogram of the number of BLAST matches
      summary_hist <- ggplot(summary_data, aes(x = n_matches)) + 
        geom_histogram(binwidth = 1, fill = "royalblue3", alpha = 0.7, colour = "black") + 
        labs(title = "Distribution of BLAST Matches per Probe", x = "Number of Matches", y = "Probe Count") + 
        scale_x_continuous(breaks = scales::breaks_width(1)) +
        theme_minimal()
      blast_summary_plot(summary_hist)
      
      removeNotification("blast_run")
      showNotification("BLAST run and summarisation complete.", type = "default")
      updateTabsetPanel(session, "results_tabs", selected = "BLAST Screening")
      
    }, error = function(e) {
      removeNotification("blast_run")
      showNotification(paste("BLAST Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # observeEvent(input$run_blast, {
  #   req(thermo_filtered_probes(), blast_db_filepath(), tx2gene_filepath()) # Requires filtered probes and paths
  #   showNotification("Running BLASTn against filtered probes (this may take a long time)...", id = "blast_run", type = "message", duration = NULL)
  #   
  #   tryCatch({
  #     # The 'run_blastn_short' function now uses the reactive filepaths
  #     blast_output <- run_blastn_short(
  #       df = thermo_filtered_probes(), 
  #       db = blast_db_filepath(), 
  #       tx2gene = tx2gene_filepath(), 
  #       cores = input$n_threads
  #     )
  #     blast_results(blast_output)
  #     
  #     # Summarize the raw BLAST output to count off-target hits per probe
  #     summary_data <- summarise_blast_output(blast_output = blast_output, allowOverlapBreakRegion = input$allowOverlapBreakRegion, evalue_cutoff = input$blast_evalue, consider_pseudogene = input$consider_pseudogene)
  #     blast_summary_data(summary_data)
  #     
  #     # Create a histogram of the number of BLAST matches
  #     summary_hist <- ggplot(summary_data, aes(x = n_matches)) + 
  #       geom_histogram(binwidth = 1, fill = "royalblue3", alpha = 0.7, colour = "black") + 
  #       labs(title = "Distribution of BLAST Matches per Probe", x = "Number of Matches", y = "Probe Count") + 
  #       scale_x_continuous(breaks = scales::breaks_width(1)) + 
  #       theme_minimal()
  #     blast_summary_plot(summary_hist)
  #     
  #     removeNotification("blast_run")
  #     showNotification("BLAST run and summarisation complete.", type = "default")
  #     updateTabsetPanel(session, "results_tabs", selected = "BLAST Screening")
  #     
  #   }, error = function(e) {
  #     removeNotification("blast_run")
  #     showNotification(paste("BLAST Error:", e$message), type = "error", duration = 10)
  #   })
  # })
  
  
  # --- 6. Screen with BLAST Results ---
  # --- This observer triggers when the "Run BLAST Screen" button is clicked.
  observeEvent(input$run_blast_screen, {
    req(blast_summary_data(), thermo_filtered_probes()) # Requires BLAST summary
    showNotification("Screening probes based on BLAST results...", type = "message")
    
    # Filter out probes that have more than the maximum allowed number of BLAST hits
    screened_probes_list <- screen_with_blast_summary(
      candidate_probes_filtered = thermo_filtered_probes(), 
      max_blast_matches = input$max_blast_matches, 
      blast_summary = blast_summary_data()
    )
    blast_screened_probes(screened_probes_list)
    
    screened_probes_df <- screened_probes_list[[1]]
    
    # Handle case where no probes pass the screen
    if(nrow(screened_probes_df) == 0){
      showNotification("Screening complete. 0 probes passed. Please adjust screening parameters.", type = "warning", duration = 8)
      blast_screened_plot(NULL)
    } else {
      showNotification(paste("Screening complete.", nrow(screened_probes_df), "probes passed."), type = "default")
      
      # Generate a plot of the BLAST-screened probes
      screened_plot <- plot_inspection(df = screened_probes_df, all_candidates = candidate_probes(), colour_param = "Tm", target_name = input$target_name)
      blast_screened_plot(screened_plot)
    }
  })
  
  
  # --- 7. Configure Final Probe Set ---
  # This observer triggers when the "Configure Final Probe Set" button is clicked.
  observeEvent(input$run_probe_config, {
    req(blast_screened_probes(), input$initiator_set) # Requires BLAST-screened probes
    showNotification("Configuring final probe set...", id = "probe_config", type = "message")
    
    tryCatch({
      plan(multisession, workers = input$n_threads) # Use parallel processing
      
      # Resolve overlaps to ensure probes are spaced correctly
      screened_probes_list <- blast_screened_probes()
      candidate_probes_screened <- screened_probes_list[[1]]
      merge_df <- screened_probes_list[[2]]
      distributed_probes <- distribute_overlapping_probes_cpp(candidate_probes_screened, merge_df, input$probe_spacing)
      
      # Cull the probe set down to the maximum number requested by the user
      culled_probes <- cull_excess_pairs(distributed_probes, input$num_probes)
      plan(sequential)
      
      # Attach the HCR initiator sequences to the final probes
      final_probes_with_initiators <- attach_hcr_initiatior(
        candidate_probes_final = culled_probes, 
        b_identifier = input$initiator_set, 
        target_name = input$target_name
      )
      final_probe_set(final_probes_with_initiators)
      
      # Generate the final probe map plot
      map_plot <- plot_final_probes(df = culled_probes, all_candidates = candidate_probes(), colour_param = "Tm", target_name = input$target_name)
      probe_map_plot(map_plot)
      
      removeNotification("probe_config")
      showNotification(paste("Successfully configured", nrow(culled_probes), "probe pairs."), type = "default")
      updateTabsetPanel(session, "results_tabs", selected = "Final HCR Probes")
      
    }, error = function(e) {
      plan(sequential)
      removeNotification("probe_config")
      showNotification(paste("Configuration Error:", e$message), type = "error", duration = 10)
    })
  })
  
  
  # --- 8. Render Outputs & Downloads ---
  
  # Render all the plots generated in the steps above
  output$inspection_plot_output <- renderPlot({ req(inspection_plot()); inspection_plot() }, res = 96) 
  output$exploratory_plots_output <- renderPlot({ req(exploratory_plots()); exploratory_plots() }, res = 96) 
  output$filtered_probes_plot_output <- renderPlot({ req(filtered_probes_plot()); filtered_probes_plot() }, res = 96)
  output$blast_summary_histogram <- renderPlot({ req(blast_summary_plot()); blast_summary_plot() }, res = 96)
  output$blast_screened_plot <- renderPlot({ req(blast_screened_plot()); blast_screened_plot() }, res = 96)
  output$probe_map_plot <- renderPlot({ req(probe_map_plot()); probe_map_plot() }, res = 96)
  
  # Render all the data tables
  output$blast_summary_table <- DT::renderDataTable({ req(blast_summary_data()); DT::datatable(blast_summary_data(), options = list(pageLength = 5, scrollX = TRUE), rownames = FALSE) })
  final_probes_for_display <- reactive({ req(final_probe_set()); get_final_probe_sequences(final_probe_set()) })
  output$final_probes_table <- DT::renderDataTable({ DT::datatable(final_probes_for_display(), options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE, caption = "Final probe sequences ready for ordering.") })
  output$details_table <- DT::renderDataTable({ req(final_probe_set()); DT::datatable(final_probe_set(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE, caption = "Detailed properties of the final selected probes.") })
  output$raw_blast_table <- DT::renderDataTable({ req(blast_results()); DT::datatable(blast_results(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE, caption = "Raw BLAST output for all thermodynamically filtered probes.") })
  
  # Define download handlers for all downloadable data
  output$download_probes <- downloadHandler(
    filename = function() { paste0(input$target_name, "_HCR-", input$initiator_set, "_probes.csv") }, 
    content = function(file) { readr::write_csv(final_probes_for_display(), file) }
  )
  
  output$download_details <- downloadHandler(
    filename = function() { paste0(input$target_name, "_HCR-", input$initiator_set, "_details.csv") }, 
    content = function(file) { readr::write_csv(final_probe_set(), file) }
  )
  
  output$download_raw_blast <- downloadHandler(
    filename = function() { paste0(input$target_name, "_HCR-", input$initiator_set, "_raw_blast.csv") }, 
    content = function(file) { req(blast_results()); readr::write_csv(blast_results(), file) }
  )
  
  # Reactive that compiles all run parameters into a data frame for display and download
  run_parameters_reactive <- reactive({
    req(final_probe_set(), target_seq(), fasta_info())
    params <- list(
      "Creation Date" = as.character(Sys.time()),
      "Threads Used" = input$n_threads,
      "Probe Set Name" = input$target_name,
      "HCRv3 Initiator Set" = input$initiator_set,
      "Oligo Length" = 52,
      "HCR Hybridization Region Length" = "25 nt each",
      "Minimum Probe Spacing" = paste(input$probe_spacing, "nt"),
      "Total Probe Pairs Generated" = nrow(final_probe_set()),
      "Target Probe-Pair dG" = input$target_dG,
      "Target Probe-Halves dG" = input$target_dG_halves,
      "dG Range Filter" = paste(input$dG_range, collapse = " : "),
      "Tm Range Filter" = paste(input$Tm_range, collapse = " : "),
      "GC Range Filter" = paste(input$GC_range, collapse = " : "),
      "Filter by A-composition" = as.character(input$pass_a_comp),
      "Filter by C-composition" = as.character(input$pass_c_comp),
      "Filter by A-stack" = as.character(input$pass_a_stack),
      "Filter by C-stack" = as.character(input$pass_c_stack),
      "Filter by C-spec stack" = as.character(input$pass_c_spec_stack),
      "BLAST E-value Cutoff" = input$blast_evalue,
      "Max BLAST Matches" = input$max_blast_matches,
      "Allow Overlap in Break Region" = as.character(input$allowOverlapBreakRegion),
      "Consider Pseudogenes" = as.character(input$consider_pseudogene),
      "BLAST Database Used" = input$blast_db_path,
      "Target Sequence Length" = nchar(target_seq()),
      "Input FASTA File" = fasta_info()$name,
      "Target Sequence Used (first 60 nt)" = paste0(substr(target_seq(), 1, 60), "...")
    )
    tibble(parameter = names(params), value = as.character(params))
  })
  
  # Render the parameters as formatted text
  output$params_output <- renderPrint({ 
    params_df <- run_parameters_reactive()
    cat("--- HCR Probe Design Run Parameters ---\n\n")
    for (i in 1:nrow(params_df)) { 
      cat(sprintf("%-35s: %s\n", params_df$parameter[i], params_df$value[i])) 
    }
  })
  
  # Handler for downloading the parameters as a text file
  output$download_params <- downloadHandler(
    filename = function() { paste0(input$target_name, "_HCR-", input$initiator_set, "_parameters.txt") }, 
    content = function(file) { readr::write_tsv(run_parameters_reactive(), file) }
  )
}
