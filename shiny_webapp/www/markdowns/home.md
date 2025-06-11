# **Welcome to the FISH Probe Designer**

<img src="images/home_funny-image.png" alt="Nonsensical funny image." width="400"> 

<br>  
<br>  

This application provides a suite of tools for designing fluorescent *in situ* hybridization (FISH) probes for various molecular imaging techniques.

## **Available Tools**

### **HCRv3 Probe Designer**

The **Hybridization Chain Reaction v3 (HCRv3)** tool implements a comprehensive pipeline for designing split-initiator probesets. The workflow is as follows:

1. **Sequence Upload**: Upload a FASTA file containing your target transcript sequence.  
2. **Candidate Generation**: The tool generates all possible 52 nt probe candidates across the transcript.  
3. **Thermodynamic Calculation & Filtering**: Key thermodynamic parameters (Tm, dG, GC-content) are calculated for each candidate probe. Probes that do not meet the specified criteria are filtered out.  
4. **BLAST Screening**: Remaining candidates are screened against a specified BLAST database to identify and remove probes with potential off-target binding.  
5. **Probe Set Configuration**: The final set of non-overlapping probes is selected, and the appropriate HCR initiator sequences are appended.

### **FLARIM Probe Designer**

*(This section is a placeholder for future development).*

Please select a tool from the navigation bar to begin.