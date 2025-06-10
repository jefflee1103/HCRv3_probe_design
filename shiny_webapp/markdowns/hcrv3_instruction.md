### **Instructions for the HCRv3 Probe Designer**

[![DOI](https://zenodo.org/badge/511863830.svg)](https://zenodo.org/doi/10.5281/zenodo.10491693)

This guide provides a step-by-step walkthrough for using the HCRv3 Probe Designer tool.

<img src="images/HCRv3_schematics.png" alt="HCR (version 3) in situ schematics." width="50%"/>

<img src="images/HCRv3_probe-design_schematics.png" alt="HCR probe design workflow." width="50%"/>

#### Step 1: Target Sequence & Processing

1.  **Probe-set name**: Provide a unique and descriptive name for your target. This name will be used in the filenames of the downloaded results.
2.  **Upload FASTA File**: Click the "Browse..." button and select a valid FASTA file (.fasta or .fa) containing the transcript sequence for your target gene. The tool will automatically process the sequence upon upload.
3.  **Number of Cores**: Specify the number of CPU cores to use for parallelizable tasks like thermodynamic calculations and BLAST searches. Using more cores will significantly speed up these processes.
4.  **Calculate Thermodynamics**: After the sequence is processed, click this button to calculate the thermodynamic properties for all potential probe candidates.

#### Step 2: Thermodynamics Filtering

This panel allows you to filter the candidate probes based on their thermodynamics properties. Adjust the sliders and checkboxes to meet your experimental requirements. Probes outside these ranges will be discarded.

-   **Melting Temperature (Tm)**: The temperature at which 50% of the probe-target duplex is dissociated.
-   **Probe-Target dG**: The Gibbs free energy of the probe binding to its target. More negative values indicate more stable binding.
-   **GC Content**: The proportion of Guanine and Cytosine bases in the probe sequence.
-   **Composition and Stack Filters**: These checkboxes allow you to filter out probes with undesirable sequence features, such as stretches of a single nucleotide (e.g., "AAAA"), which can cause non-specific binding.

Click **"Run Thermo Filter"** to apply your selections.

#### Step 3: BLAST Running

This step screens the remaining probes for potential off-target binding.

1.  **Path to BLAST Database**: Provide the full file path to your local BLAST database (e.g., a genome or transcriptome).
2.  **Path to Tx2gene map**: Provide the path to a CSV file that maps transcript IDs to gene IDs. This is used to properly attribute BLAST hits.
3.  **BLAST E-value Cutoff**: Set the expectation value cutoff for BLAST hits. Lower values are more stringent.

Click **"Run BLAST"** to begin the search. This may take a considerable amount of time depending on the size of your database and the number of probes.

#### Step 4: BLAST Screening

After the BLAST search is complete, use this panel to filter probes based on the number of significant off-target hits they produced.

-   **Max number of BLAST hits**: Probes with more than this number of hits will be discarded. A value of 1 is stringent and typically means the probe only hits the intended target transcript.

Click **"Run BLAST Screen"** to apply the filter.

#### Step 5: Probe Set Configuration

The final step configures the probe set from the pool of specific, thermodynamically stable candidates.

1.  **Initiator Set**: Select the HCRv3 initiator pair (e.g., B1, B2) you will be using in your experiment. The tool will append the correct initiator sequences to the probes.
2.  **Minimum Spacing**: The minimum number of nucleotides between adjacent probes. This ensures that probes do not sterically hinder each other.
3.  **Maximum Probes to Select**: The maximum number of non-overlapping probe pairs to include in the final set.

Click **"Configure Final Probe Set"** to generate the final results, which can be viewed and downloaded in the "Final HCR Probes" tab.
