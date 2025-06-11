# ui.R
# Defines the user interface for the FISH Probe Designer, now with a multi-tab navigation bar.

# Load necessary libraries for the UI
library(shiny)
library(shinythemes)

# --- UI Definition with Navbar ---
ui <- navbarPage(
  title = "FISH Probe Designer",
  theme = shinytheme("lumen"),
  
  # --- Tab 1: Home (Landing Page) ---
  # This panel serves as the landing page and displays content from an external markdown file.
  tabPanel("Home",
           fluidRow(
             column(width = 8, offset = 2,
                    # For this to work, place home.md in the app's root directory.
                    includeMarkdown("markdowns/home.md") 
             )
           )
  ),
  
  # --- Tab 2: HCRv3 Probe Designer ---
  # This panel contains the entire HCRv3 application.
  # The instructions are now the first tab in the main results panel.
  tabPanel("HCRv3",
           sidebarLayout(
             
             # --- Sidebar Panel for HCRv3 Inputs ---
             sidebarPanel(
               width = 3,
               
               # --- Section 1: Target Sequence & Processing ---
               h4("1. Target Sequence & Processing"),
               wellPanel(
                 textInput("target_name", "Probe-set name:", value = "Give probe-set a name"),
                 fileInput("target_fasta", "Upload FASTA File:", accept = c(".fasta", ".fa")),
                 hr(),
                 numericInput("n_threads", "Number of Cores:", value = 2, min = 1, max = parallel::detectCores()),
                 actionButton("run_thermo_calc", "Calculate Thermodynamics", icon = icon("cogs"), width = "100%")
               ),
               
               # --- Section 2: Thermodynamics Filtering ---
               h4("2. Thermodynamics Filtering"),
               wellPanel(
                 sliderInput("Tm_range", "Melting Temperature (°C):", min = 20, max = 100, value = c(60, 90), step = 1),
                 sliderInput("dG_range", "Probe-Target dG (kcal/mol):", min = -120, max = -30, value = c(-100, -60), step = 1),
                 sliderInput("GC_range", "GC Content (Proportion):", min = 0.1, max = 0.9, value = c(0.4, 0.6), step = 0.01),
                 numericInput("target_dG", "Target dG (Pair):", value = -80),
                 numericInput("target_dG_halves", "Target dG (Halves):", value = -38),
                 checkboxInput("pass_a_comp", "Filter by A-composition", value = TRUE),
                 checkboxInput("pass_c_comp", "Filter by C-composition", value = FALSE),
                 checkboxInput("pass_a_stack", "Filter by A-stack", value = TRUE),
                 checkboxInput("pass_c_stack", "Filter by C-stack", value = TRUE),
                 checkboxInput("pass_c_spec_stack", "Filter by C-spec stack", value = FALSE),
                 actionButton("run_thermo_filter", "Run Thermo Filter", icon = icon("fire"), width = "100%")
               ),
               
               # --- Section 3: BLAST Running ---
               h4("3. BLAST Running"),
               wellPanel(
                 textInput("blast_db_path", "Path to BLAST Database:", value = "~/Documents/Data/BLAST/Dmel/Dmel_BLASTdb_ens99.fa"),
                 textInput("tx2gene_path", "Path to Tx2gene map (.csv):", value = "../data/Dmel_tx2gene_ENSEMBL_v99.csv"),
                 numericInput("blast_evalue", "BLAST E-value Cutoff:", value = 0.1, min = 0),
                 checkboxInput("allowOverlapBreakRegion", "Allow overlap spanning the 25/25 break region", value = TRUE),
                 actionButton("run_blast", "Run BLAST", icon = icon("searchengin"), width = "100%")
               ),
               
               # --- Section 4: BLAST Screening ---
               h4("4. BLAST Screening"),
               wellPanel(
                 numericInput("max_blast_matches", "Max number of BLAST hits", value = 1, min = 0),
                 checkboxInput("consider_pseudogene", "Consider pseudogene", value = TRUE),
                 actionButton("run_blast_screen", "Run BLAST Screen", icon = icon("filter"), width = "100%")
               ),
               
               # --- Section 5: Probe Set Configuration ---
               h4("5. Probe Set Configuration"),
               wellPanel(
                 selectInput("initiator_set", "Initiator Set:", choices = c("B1", "B2", "B3", "B4", "B5", "B7", "B9", "B10", "B13", "B14", "B15", "B17")),
                 numericInput("probe_spacing", "Minimum Spacing (nt):", value = 3, min = 2, step = 1),
                 numericInput("num_probes", "Maximum Probes to Select:", value = 30, min = 1, max = 200, step = 1),
                 actionButton("run_probe_config", "Configure Final Probe Set", icon = icon("sitemap"), width = "100%")
               )
             ),
             
             # --- Main Panel for HCRv3 Outputs ---
             mainPanel(
               width = 9,
               tabsetPanel(
                 id = "results_tabs",
                 
                 # --- Instructions Tab ---
                 tabPanel("Instructions",
                          # This will look for the file in a 'markdowns' sub-folder.
                          # You must create this folder in your app's root directory.
                          includeMarkdown("markdowns/hcrv3_instruction.md")
                 ),
                 
                 # --- Subsequent Results Tabs ---
                 tabPanel("Candidate Probes", 
                          h4("Initial Candidate Probe Map"),
                          h5("Splice junctions indicated by vertical lines"),
                          plotOutput("inspection_plot_output", height = "800px")
                 ),
                 tabPanel("Thermo Filtering",
                          h4("Exploratory Thermodynamic Plots"),
                          plotOutput("exploratory_plots_output", height = "400px"),
                          hr(),
                          h4("Thermodynamically Filtered Probes"),
                          plotOutput("filtered_probes_plot_output", height = "400px")
                 ),
                 tabPanel("BLAST Screening",
                          h4("BLAST Match Distribution"),
                          plotOutput("blast_summary_histogram", height = "400px"),
                          hr(),
                          h4("BLAST Summary Table"),
                          DT::dataTableOutput("blast_summary_table"),
                          hr(),
                          h4("BLAST Screened Probes"),
                          plotOutput("blast_screened_plot", height = "400px")
                 ),
                 tabPanel("Final HCR Probes", 
                          h4("Probe Map"),
                          plotOutput("probe_map_plot", height = "600px"),
                          hr(),
                          h4("Final Probe Sequences to Order"),
                          DT::dataTableOutput("final_probes_table"),
                          downloadButton("download_probes", "Download Probes (.csv)")
                 ),
                 tabPanel("Run Parameters & Properties",
                          h4("Parameters Used for This Design Run"),
                          verbatimTextOutput("params_output"),
                          downloadButton("download_params", "Download Parameters (.txt)"),
                          hr(),
                          h4("Detailed Probe Properties and Thermodynamics"),
                          DT::dataTableOutput("details_table"),
                          downloadButton("download_details", "Download Details (.csv)"),
                          hr(),
                          h4("Raw BLAST Output"),
                          DT::dataTableOutput("raw_blast_table"),
                          downloadButton("download_raw_blast", "Download Raw BLAST Output (.csv)")
                 )
               )
             )
           )
  ),
  
  # --- Tab 3: FLARIM ---
  tabPanel("FLARIM",
           fluidRow(
             column(width = 8, offset = 2,
                    h4("Coming SOON!") 
             )
           )
  ),
  
  
  
  # --- Tab 4: Contact (Placeholder) ---
  tabPanel("Contact",
           fluidRow(
             column(width = 8, offset = 2,
                    h4("jeff.lee{at}glasgow.ac.uk") 
             )
           )
  )
)
