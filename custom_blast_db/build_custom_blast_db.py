#!/usr/bin/env python3

"""
build_custom_blast_db.py

Description:
This script automates the creation of a custom BLAST nucleotide database tailored for
filtering Hybridisation Chain Reaction (HCR) probes against potential cross-reactivity.
It focuses on RNA sequences, including cDNAs, non-coding RNAs (ncRNAs), and
intronic sequences with specified flanking regions. The script downloads the necessary
FASTA and GTF files from Ensembl FTP, extracts intronic sequences, concatenates
all relevant sequences into a single FASTA file, and generates a transcript-to-gene
mapping file (tx2gene.csv). This final FASTA file can then be
used with NCBI's makeblastdb tool to generate the BLAST database.

Usage:
python build_custom_blast_db.py \
  --genome-fasta-url <URL_TO_GENOME_FASTA_GZ> \
  --cdna-fasta-url <URL_TO_CDNA_FASTA_GZ> \
  --ncrna-fasta-url <URL_TO_NCRNA_FASTA_GZ> \
  --gtf-url <URL_TO_GTF_GZ> \
  --output-dir <PATH_TO_OUTPUT_DIRECTORY> \
  [--flanking-nt <INTEGER>]

Examples:
python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz \
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz \
  --gtf-url          https://ftp.ensembl.org/pub/release-114/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.114.gtf.gz \
  --output-dir       ./yeast_hcr_rna_db_output \
  --flanking-nt      40

python build_custom_blast_db.py \
  --genome-fasta-url https://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa.gz\
  --cdna-fasta-url   https://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz \
  --ncrna-fasta-url  https://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz\
  --gtf-url          https://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz \
  --output-dir       ./mmus_hcr_rna_db_output \
  --flanking-nt      40

Expected Outputs (in the specified --output-dir):
1.  genome.dna.fa: Uncompressed genomic DNA FASTA file.
2.  cdna.all.fa: Uncompressed cDNA FASTA file.
3.  ncrna.fa: Uncompressed ncRNA FASTA file.
4.  annotation.gtf: Uncompressed GTF annotation file.
5.  introns_flank<N>.fa: FASTA file containing extracted intronic sequences with N
    nucleotides of flanking sequence on each side (e.g., introns_flank40.fa).
6.  custom_rna_db.fa: The final concatenated FASTA file containing sequences from
    cdna.all.fa, ncrna.fa, and introns_flank<N>.fa. This file is ready to be
    used with 'makeblastdb'.
7.  tx2gene.csv: A CSV file mapping transcript IDs (as tx_id) to gene IDs,
    gene names, and gene biotypes.
8.  (Temporary) temp_annotation.gtf.db: A temporary SQLite database created by
    gffutils during intron extraction. This file is removed upon successful completion.

The script will also print the 'makeblastdb' command to the console, which can be
used to create the BLAST database from 'custom_rna_db.fa'.
"""

import argparse
import os
import sys
import requests
import gzip
import shutil
import gffutils
from pyfaidx import Fasta
import numpy as np
import csv # Added for tx2gene generation

def download_and_unzip(url, final_unzipped_path, download_dir):
    """
    Downloads a file from a URL, unzips if it's a .gz file.
    The final_unzipped_path is the desired path for the uncompressed file.
    download_dir is where the initial download (potentially gzipped) will be stored.
    """
    # Determine a temporary name for the downloaded file (could be gzipped or not)
    base_name = os.path.basename(final_unzipped_path)
    if url.endswith(".gz") and not base_name.endswith(".gz"):
        temp_downloaded_filename = os.path.join(download_dir, base_name + ".gz")
    else:
        temp_downloaded_filename = os.path.join(download_dir, base_name)
    
    # Ensure download_dir exists
    os.makedirs(download_dir, exist_ok=True)

    print(f"  Downloading {url} to {temp_downloaded_filename}...")
    try:
        response = requests.get(url, stream=True, timeout=300) # 5 min timeout
        response.raise_for_status()  # Raise an exception for HTTP errors
        with open(temp_downloaded_filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"  Successfully downloaded {temp_downloaded_filename}")

        if temp_downloaded_filename.endswith(".gz"):
            print(f"  Unzipping {temp_downloaded_filename} to {final_unzipped_path}...")
            with gzip.open(temp_downloaded_filename, 'rb') as f_in:
                with open(final_unzipped_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f"  Successfully unzipped to {final_unzipped_path}")
            os.remove(temp_downloaded_filename)  # Remove the .gz file after unzipping
        elif temp_downloaded_filename != final_unzipped_path:
            # If not gzipped, move/rename the downloaded file to the final_unzipped_path
            shutil.move(temp_downloaded_filename, final_unzipped_path)
            print(f"  File moved to {final_unzipped_path}")
        else:
            # If not gzipped and paths are the same, do nothing extra
            print(f"  File is not gzipped. Final path: {final_unzipped_path}")
        
        return final_unzipped_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {url}: {e}", file=sys.stderr)
        if os.path.exists(temp_downloaded_filename):
            os.remove(temp_downloaded_filename)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred during download/unzip of {url}: {e}", file=sys.stderr)
        if os.path.exists(temp_downloaded_filename):
            os.remove(temp_downloaded_filename)
        sys.exit(1)

def extract_intronic_sequences(gtf_filepath, genome_fasta_filepath, output_introns_fasta_path, flanking_nt, temp_db_path):
    """
    Extracts intronic sequences with flanking regions from a genome FASTA based on GTF annotations.
    Also creates and populates the gffutils database used by other functions.
    Args:
        gtf_filepath (str): Path to the GTF annotation file.
        genome_fasta_filepath (str): Path to the genome FASTA file.
        output_introns_fasta_path (str): Path to write the output FASTA file of intronic sequences.
        flanking_nt (int): Number of nucleotides to include on each side of the intron.
        temp_db_path (str): Path for the temporary gffutils SQLite database.
    """
    print(f"    Loading genome from {genome_fasta_filepath}...")
    try:
        genome = Fasta(genome_fasta_filepath, sequence_always_upper=True) 
    except Exception as e:
        print(f"Error loading genome FASTA with pyfaidx: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"    Creating/loading gffutils database from {gtf_filepath} (temp DB at {temp_db_path})...")
    try:
        db = gffutils.create_db(gtf_filepath, temp_db_path, force=True,
                                merge_strategy='create_unique',
                                disable_infer_genes=True, 
                                disable_infer_transcripts=True) 
        
        print("    Creating intron features from GTF...")
        introns_list = list(db.create_introns())
        
        if not introns_list:
            print("    No introns were derived by db.create_introns(). Check GTF content and format.")
        else:
            print(f"    Found {len(introns_list)} potential intron features. Updating database with these intron features...")
            db.update(introns_list, merge_strategy='create_unique',
                      disable_infer_genes=True, disable_infer_transcripts=True, make_backup=False)

    except Exception as e:
        print(f"Error creating or processing gffutils database: {e}", file=sys.stderr)
        if os.path.exists(temp_db_path): os.remove(temp_db_path)
        if os.path.exists(temp_db_path + "-journal"): os.remove(temp_db_path + "-journal")
        sys.exit(1)

    print(f"    Extracting intronic sequences with {flanking_nt}nt flanking regions to {output_introns_fasta_path}...")
    count_extracted_sequences = 0
    with open(output_introns_fasta_path, 'w') as out_f:
        for gene in db.features_of_type('gene'):
            gene_id_attr = gene.attributes.get('gene_id', [gene.id])[0]
            
            current_gene_introns = list(db.children(gene, featuretype='intron', order_by='start'))

            if not current_gene_introns:
                continue

            intron_coords_for_gene = []
            for intron_feature in current_gene_introns:
                flanked_start_0based = intron_feature.start - 1 - flanking_nt
                flanked_end_0based = intron_feature.end + flanking_nt 

                try:
                    chrom_len = len(genome[intron_feature.chrom])
                except KeyError:
                    print(f"Warning: Chromosome '{intron_feature.chrom}' for intron in gene '{gene_id_attr}' not found in FASTA. Skipping this intron.", file=sys.stderr)
                    continue
                
                actual_start = max(0, flanked_start_0based)
                actual_end = min(chrom_len, flanked_end_0based)

                if actual_start < actual_end: 
                    intron_coords_for_gene.append([actual_start, actual_end, intron_feature.chrom, intron_feature.strand])
            
            if not intron_coords_for_gene:
                continue
            
            unique_intervals_details = []
            seen_coords_details = set()

            for start_coord, end_coord, chrom_val, strand_val in intron_coords_for_gene:
                coord_tuple = (chrom_val, start_coord, end_coord, strand_val)
                if coord_tuple not in seen_coords_details:
                    unique_intervals_details.append(coord_tuple)
                    seen_coords_details.add(coord_tuple)
            
            intron_counter_for_gene = 1
            for chrom_val, interval_0based_start, interval_0based_end, strand_val in unique_intervals_details:
                try:
                    seq_record = genome[chrom_val][interval_0based_start:interval_0based_end]
                except KeyError: 
                    print(f"Warning: Chromosome '{chrom_val}' not found when extracting sequence for gene '{gene_id_attr}'. Skipping.", file=sys.stderr)
                    continue
                
                sequence_str = seq_record.seq
                if strand_val == '-':
                    sequence_str = seq_record.complement.reverse.seq 
                
                if sequence_str:
                    safe_gene_id = str(gene_id_attr).replace(" ", "_").replace(";", "_").replace(":", "_").replace("|","_")
                    header = f">{safe_gene_id}_I{intron_counter_for_gene} {chrom_val}:{interval_0based_start+1}-{interval_0based_end}({strand_val})"
                    out_f.write(f"{header}\n{sequence_str}\n")
                    intron_counter_for_gene += 1
                    count_extracted_sequences += 1
    
    print(f"    Extracted {count_extracted_sequences} unique intronic sequences (with flanking regions).")
    if count_extracted_sequences == 0:
        print("    Warning: No intronic sequences were extracted. The output intron FASTA will be empty or incomplete.", file=sys.stderr)

def generate_tx2gene_file(gff_db_path, output_csv_path):
    """
    Generates a CSV file mapping transcript IDs to gene IDs, gene names, and gene biotypes.
    Args:
        gff_db_path (str): Path to the gffutils SQLite database file.
        output_csv_path (str): Path to write the output CSV file.
    """
    print(f"  Generating transcript-to-gene map: {output_csv_path}...")
    try:
        # Load the existing database; do not recreate it here.
        db = gffutils.FeatureDB(gff_db_path, keep_order=True)
    except Exception as e:
        print(f"Error loading GFF database {gff_db_path}: {e}", file=sys.stderr)
        sys.exit(1)

    records = set() # Use a set to store tuples of (tx_id, gene_id, gene_name, gene_biotype) to ensure uniqueness

    # Iterate through all transcript features in the database
    for transcript_feature in db.features_of_type('transcript', order_by='id'):
        # Extract transcript_id, fallback to feature.id if 'transcript_id' attribute is missing
        tx_id_list = transcript_feature.attributes.get('transcript_id')
        tx_id = tx_id_list[0] if tx_id_list else transcript_feature.id

        # Extract gene_id from transcript attributes
        gene_id_from_attrs_list = transcript_feature.attributes.get('gene_id')
        gene_id = gene_id_from_attrs_list[0] if gene_id_from_attrs_list else None

        gene_name = ""
        gene_biotype = ""

        if gene_id:
            try:
                # Fetch the corresponding gene feature using the gene_id.
                # This assumes gene_id attribute of transcript refers to the ID of a gene feature.
                gene_feature = db[gene_id] 
                
                gene_name_list = gene_feature.attributes.get('gene_name')
                gene_name = gene_name_list[0] if gene_name_list else ""
                
                # Try 'gene_biotype' first, then 'gene_type' as a fallback (common in Ensembl GTFs)
                gene_biotype_list = gene_feature.attributes.get('gene_biotype')
                if not gene_biotype_list:
                    gene_biotype_list = gene_feature.attributes.get('gene_type') # Fallback
                gene_biotype = gene_biotype_list[0] if gene_biotype_list else ""

            except gffutils.exceptions.FeatureNotFoundError:
                # This occurs if a transcript's gene_id attribute points to an ID not present as a gene feature's primary ID.
                print(f"    Warning: Gene ID '{gene_id}' from transcript '{tx_id}' not found as a gene feature. Corresponding gene_name and gene_biotype will be empty.", file=sys.stderr)
            except Exception as e:
                # Catch any other errors during attribute extraction for a specific gene
                print(f"    Warning: Error processing gene attributes for gene ID '{gene_id}' (transcript '{tx_id}'): {e}", file=sys.stderr)
        
        # Only add record if both essential IDs are present
        if tx_id and gene_id:
            records.add((tx_id, gene_id, gene_name, gene_biotype))
        else:
            if not tx_id:
                 print(f"    Warning: Transcript feature {transcript_feature.id} missing 'transcript_id' attribute or could not be determined. Skipping for tx2gene.", file=sys.stderr)
            if not gene_id: # This implies gene_id_from_attrs_list was None or empty
                 print(f"    Warning: Transcript feature {tx_id or transcript_feature.id} missing 'gene_id' attribute. Skipping for tx2gene.", file=sys.stderr)


    if not records:
        print(f"    Warning: No transcript-to-gene records were generated. CSV file '{output_csv_path}' will be empty or contain only headers.", file=sys.stderr)
    
    try:
        with open(output_csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['tx_id', 'gene_id', 'gene_name', 'gene_biotype'])
            # Sort records for consistent output order, though set itself is unordered.
            for record in sorted(list(records)):
                writer.writerow(record)
        print(f"  Successfully wrote {len(records)} unique records to {output_csv_path}")
    except IOError as e:
        print(f"Error writing tx2gene CSV file {output_csv_path}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Build a custom BLAST database from Ensembl cDNA, ncRNA, and extracted intronic sequences. Also generates a tx2gene mapping file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--genome-fasta-url", required=True, help="URL to the genome DNA FASTA file (e.g., .dna.toplevel.fa.gz)")
    parser.add_argument("--cdna-fasta-url", required=True, help="URL to the cDNA FASTA file (e.g., .cdna.all.fa.gz)")
    parser.add_argument("--ncrna-fasta-url", required=True, help="URL to the ncRNA FASTA file (e.g., .ncrna.fa.gz)")
    parser.add_argument("--gtf-url", required=True, help="URL to the GTF annotation file (e.g., .gtf.gz)")
    parser.add_argument("--output-dir", required=True, help="Directory to store downloaded files, intermediate files, and the final concatenated FASTA.")
    parser.add_argument("--flanking-nt", type=int, default=40, help="Number of nucleotides to include as flanking regions for introns.")

    args = parser.parse_args()

    # --- 1. Setup Output Directory ---
    print(f"Creating output directory: {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define local file paths
    genome_fasta_local_unzipped_path = os.path.join(args.output_dir, "genome.dna.fa")
    cdna_fasta_local_unzipped_path = os.path.join(args.output_dir, "cdna.all.fa")
    ncrna_fasta_local_unzipped_path = os.path.join(args.output_dir, "ncrna.fa")
    gtf_local_unzipped_path = os.path.join(args.output_dir, "annotation.gtf")
    
    # Intermediate and final files
    introns_fasta_path = os.path.join(args.output_dir, f"introns_flank{args.flanking_nt}.fa")
    concatenated_fasta_path = os.path.join(args.output_dir, "custom_rna_db.fa")
    # temp_gtf_db_path is where gffutils stores its SQLite DB for the GTF.
    temp_gtf_db_path = os.path.join(args.output_dir, "temp_annotation.gtf.db") 
    tx2gene_csv_path = os.path.join(args.output_dir, "tx2gene.csv") # Path for the new CSV file

    # --- 2. Download and Unzip Files ---
    print("\n--- Downloading and Unzipping Files ---")
    print("Downloading genome FASTA...")
    download_and_unzip(args.genome_fasta_url, genome_fasta_local_unzipped_path, args.output_dir)
    
    print("\nDownloading cDNA FASTA...")
    download_and_unzip(args.cdna_fasta_url, cdna_fasta_local_unzipped_path, args.output_dir)

    print("\nDownloading ncRNA FASTA...")
    download_and_unzip(args.ncrna_fasta_url, ncrna_fasta_local_unzipped_path, args.output_dir)

    print("\nDownloading GTF Annotation File...")
    download_and_unzip(args.gtf_url, gtf_local_unzipped_path, args.output_dir)

    # --- 3. Extract Intronic Sequences (also creates the gffutils DB at temp_gtf_db_path) ---
    print("\n--- Extracting Intronic Sequences ---")
    extract_intronic_sequences(
        gtf_filepath=gtf_local_unzipped_path,
        genome_fasta_filepath=genome_fasta_local_unzipped_path,
        output_introns_fasta_path=introns_fasta_path,
        flanking_nt=args.flanking_nt,
        temp_db_path=temp_gtf_db_path # This path will contain the gffutils DB
    )

    # --- 4. Generate Tx2Gene Mapping File ---
    # This step requires the gffutils database (temp_gtf_db_path) created by extract_intronic_sequences.
    print("\n--- Generating Tx2Gene Map ---")
    if os.path.exists(temp_gtf_db_path):
         generate_tx2gene_file(temp_gtf_db_path, tx2gene_csv_path)
    else:
        # This case should ideally not be reached if extract_intronic_sequences succeeded and wasn't skipped.
        print(f"  Warning: Temporary GTF database ({temp_gtf_db_path}) not found. Skipping Tx2Gene map generation.", file=sys.stderr)


    # --- 5. Concatenate FASTA Files ---
    print("\n--- Concatenating FASTA Files ---")
    files_to_concatenate = [cdna_fasta_local_unzipped_path, ncrna_fasta_local_unzipped_path, introns_fasta_path]
    
    with open(concatenated_fasta_path, 'wb') as outfile: 
        for fasta_file_path in files_to_concatenate:
            if os.path.exists(fasta_file_path) and os.path.getsize(fasta_file_path) > 0:
                print(f"  Adding {os.path.basename(fasta_file_path)} to {os.path.basename(concatenated_fasta_path)}...")
                with open(fasta_file_path, 'rb') as infile:
                    shutil.copyfileobj(infile, outfile)
            else:
                print(f"  Warning: {os.path.basename(fasta_file_path)} is missing or empty. Skipping its concatenation.", file=sys.stderr)

    print(f"\nSuccessfully created concatenated FASTA for BLAST database: {concatenated_fasta_path}")

    # --- 6. Clean up temporary files ---
    print("\n--- Cleaning Up Temporary Files ---")
    if os.path.exists(temp_gtf_db_path):
        os.remove(temp_gtf_db_path)
        print(f"  Removed temporary GTF database: {temp_gtf_db_path}")
    if os.path.exists(temp_gtf_db_path + "-journal"): 
        os.remove(temp_gtf_db_path + "-journal")
        print(f"  Removed temporary GTF database journal file.")
    
    # Optional: Clean up .fai index files created by pyfaidx
    # if os.path.exists(genome_fasta_local_unzipped_path + ".fai"): os.remove(genome_fasta_local_unzipped_path + ".fai")

    print("\n--- Process Complete ---")
    print(f"The combined FASTA file for BLAST database creation is: {concatenated_fasta_path}")
    print(f"The transcript-to-gene mapping file is: {tx2gene_csv_path}") # Added message for the new file
    print("You can now create the BLAST database using a command like:")
    # Example makeblastdb command
    db_title = "Custom_RNA_DB_from_script" # Make title shell-friendly
    blast_db_name = os.path.join(args.output_dir, "custom_rna_blast_db") # Place output DB in output_dir
    print(f"makeblastdb -in {concatenated_fasta_path} -dbtype nucl -parse_seqids -out {blast_db_name} -title \"{db_title}\"")


if __name__ == '__main__':
    main()
