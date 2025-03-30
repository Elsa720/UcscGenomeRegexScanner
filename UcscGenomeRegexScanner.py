#!/usr/bin/env python3

import os
import re
import gzip
import argparse
import csv
from Bio import SeqIO

VERSION='0.2.1'

'''
0.2.1 Update: The update output bed contains sequence length information.
'''

'''
version = 0.2.0 Compared to version = 0.1.0, the antisense chain filtering option for regular expressions is added. 

Main changes:
1. Add `--reverse` argument and reverse complement function to add antisense chain filtering information. Support IUPAC fuzzy base code and preserve case. 
2. Add the `--quiet` argument, to not print any output information, which is disabled by default. 
3. Add the ignore-case condition of the input regular expression. 
4. Increase output BED file information, add sense,antisense chain and sequence information; 
5. The start and end of the BED file are absolute location information, considering the two input cases of genome fasta file and fragment sequence. 
    For example: 
    >chr1 (default is genome fa, starting from 0) 
    >chr1:1-100 (Check whether the header includes: -, if there is, it is considered to be fragmented sequence fa, and the start-1 is converted to 0-based before processing)

Multi-threading has not been implemented yet.
'''

FASTA_DIR = "fasta"
BED_DIR = "bed_files"

def download_genome(species, quiet=False):
    """Download the UCSC species genome from the UCSC server."""
    os.makedirs(FASTA_DIR, exist_ok=True) 
    local_file = os.path.join(FASTA_DIR, f"{species}.fa.gz")

    if os.path.exists(local_file):
        if not quiet:
            print(f"{species} genome already downloaded.")
        return local_file

    if species.startswith("GCF_") and re.match(r"GCF_\d+\.\d+",species):
        # GCF format：GCF_xxxxxxxxx.x
        gcf_parts = species.split("_")[1]
        gcf_path = f"{gcf_parts[:3]}/{gcf_parts[3:6]}/{gcf_parts[6:9]}/{species}"
        fasta_url = f"ftp://hgdownload.soe.ucsc.edu/hubs/GCF/{gcf_path}/{species}.fa.gz"

    else:
        fasta_urls = [
            f"ftp://hgdownload.cse.ucsc.edu/goldenPath/{species}/bigZips/{species}.fa.gz",
            f"ftp://hgdownload.soe.ucsc.edu/goldenPath/{species}/bigZips/{species}.fa.gz"
        ]

        for fasta_url in fasta_urls:
            if not quiet:
                print(f"Trying {fasta_url}...")
            if os.system(f"wget --timestamping '{fasta_url}' -O {local_file}") == 0:
                return local_file

        if not quiet:
            print(f"The fasta file for this species {species} was not found in either the cse or soe paths.")
        return None

    # GCF species download
    if not quiet:
        print(f"Downloading {species} genome from {fasta_url}...")
    if os.system(f"wget --timestamping '{fasta_url}' -O {local_file}") == 0:
        return local_file
    else:
        if not quiet:
            print(f"The fasta file for this bacterial species {species} was not found.")
        return None

def load_fasta_from_directory(fastas, quiet=False):
    """Look for the .fa/.fa.gz file in the FASTA DIR directory"""
    for ext in [".fa", ".fa.gz"]:
        local_file = os.path.join(FASTA_DIR, f"{fastas}{ext}")
        if os.path.exists(local_file):
            if not quiet:
                print(f"Found local FASTA file: {local_file}")
            return local_file
    if not quiet:
        print(f"Warning: No FASTA file found for {fastas} in {FASTA_DIR}")
    return None

def process_fasta(species_or_file, patterns, mode, reverse=False, quiet=False):
    """
    Processing species according to mode: 
    -mode ="species" : download UCSC reference genome 
    -mode ="fasta" : find local FASTA file
    """
    if mode == "species":
        fasta_file = download_genome(species_or_file, quiet=quiet)
    elif mode == "fasta":
        fasta_file = load_fasta_from_directory(species_or_file, quiet=quiet)
    else:
        raise ValueError("Invalid mode. Expected 'species' or 'fasta'.")

    if fasta_file and os.path.exists(fasta_file):
        return process_fasta_file(fasta_file, patterns, reverse=reverse, quiet=quiet)
    
    if not quiet:
        print(f"Skipping {species_or_file} due to missing FASTA file.")
    return []

def reverse_complement(seq):
    compdict = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'R': 'Y', 'Y': 'R', 'S': 'W', 'W': 'S',
        'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
        'H': 'D', 'V': 'B', 'N': 'N',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
        'r': 'y', 'y': 'r', 's': 'w', 'w': 's',
        'k': 'm', 'm': 'k', 'b': 'v', 'd': 'h',
        'h': 'd', 'v': 'b', 'n': 'n'
    }
    # use list comprehensions to generate complementary sequences and reverse them
    return ''.join(compdict[n] for n in seq[::-1])

def process_fasta_file(fasta_file, patterns, reverse=False, quiet=False):
    if not os.path.exists(fasta_file):
        if not quiet:
            print(f"Skipping {fasta_file} due to missing file.")
        return []
    
    if not quiet:
        print(f"{fasta_file} is on the processing...")

    total_fasta_length = 0
    results = []
    compiled_patterns = {pat: re.compile(pat, re.IGNORECASE) for pat in patterns} # ignore case
    match_stats = {pat: {'count': 0, 'bases': 0, 'regions': []} for pat in patterns}

    # displays the regular expression being processed at the beginning
    if not quiet:
        for pat in patterns:
            print(f"Searching for pattern: {pat} on the {fasta_file}...")

    # parse the fasta file
    with gzip.open(fasta_file, "rt") if fasta_file.endswith(".gz") else open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            header = record.id
            seq = str(record.seq)
            total_fasta_length += len(seq)
            
            # check that the header contains fragmented information
            if ":" in header and "-" in header:
                chrom, range_info = header.split(":")
                fragment_start = int(range_info.split("-")[0]) - 1  # converts 1-based from chrx:x-x format to 0-based
            else:
                chrom = header
                fragment_start = 0


            for pat, regex in compiled_patterns.items():
                matches = list(regex.finditer(seq))
                match_stats[pat]['count'] += len(matches)
                for match in matches:
                    start, end = match.start(), match.end() #python 0-based
                    matched_seq = seq[start:end]
                    match_stats[pat]['bases'] += (end - start)
                    if fragment_start > 0:
                        bed_start = fragment_start + start  
                        bed_end = bed_start + len(matched_seq)  
                    else:
                        bed_start = start
                        bed_end = end
                    match_stats[pat]['regions'].append(f"{chrom}\t{bed_start}\t{bed_end}\t{len(matched_seq)}\t+\t{matched_seq}\n")
                    
            # If --reverse is specified, the antisense chain is handled
            if reverse:
                rev_seq = reverse_complement(seq)
                for pat, regex in compiled_patterns.items():
                    matches = list(regex.finditer(rev_seq))
                    match_stats[pat]['count'] += len(matches)
                    for match in matches:
                        start, end = match.start(), match.end()
                        matched_seq = rev_seq[start:end]
                        match_stats[pat]['bases'] += (end - start)
                        if fragment_start > 0:
                            rev_start = fragment_start + (len(seq) - end)  
                            rev_end = fragment_start + (len(seq) - start)  
                        else:
                            rev_start = len(seq) - end
                            rev_end = len(seq) - start
                        match_stats[pat]['regions'].append(f"{chrom}\t{rev_start}\t{rev_end}\t{len(matched_seq)}\t-\t{matched_seq}\n")


    # BED files saved
    os.makedirs(BED_DIR, exist_ok=True)
    for pat in patterns:
        bed_filename = os.path.join(BED_DIR, f"{os.path.basename(fasta_file).replace('.fa', '').replace('.gz', '')}_{re.sub(r'[^a-zA-Z0-9]', '', pat)}.bed")
        with open(bed_filename, "w") as bed_file:
            bed_file.writelines(match_stats[pat]['regions'])
        occupancy_rate = 100 * match_stats[pat]['bases'] / total_fasta_length if total_fasta_length > 0 else 0
        if not quiet:
            print("="*50)
            print(f"FASTA File: {fasta_file}")
            print(f"    Pattern:  {pat}")
            print(f"    Matches Counts: {match_stats[pat]['count']}, Matched Bases: {match_stats[pat]['bases']}, Total FASTA File Bases: {total_fasta_length}, Bases Occupancy Rate(%): {occupancy_rate}")
            print(f"BED files saved:{bed_filename}")
            print("="*50)
        results.append((os.path.basename(fasta_file), pat, match_stats[pat]['count'], match_stats[pat]['bases'], total_fasta_length, occupancy_rate))

    return results

def load_list_from_file(filename):
    if not os.path.exists(filename):
        print(f"Error: File '{filename}' not found!")
        return []

    with open(filename, "r") as file:
        fasta_list = [line.strip() for line in file.readlines() if line.strip()]
    return fasta_list

def load_patterns_from_file(filename):
    if not os.path.exists(filename):
        print(f"Error: Pattern file '{filename}' not found!")
        return []

    with open(filename, "r") as file:
        patterns = [line.strip() for line in file.readlines() if line.strip()]
    return patterns

def main():
    parser = argparse.ArgumentParser(description="Scan UCSC genomes or FASTA files for motif patterns and generate BED files.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--species","-s", nargs="+", help="List of species (e.g., hg38 mm10)")
    group.add_argument("--species-list","-sl", type=str, help="Path to a text file with species names (one per line)")
    group.add_argument("--fasta", "-f", type=str, help="Path to a specific FASTA file")
    group.add_argument("--fasta-list", "-fl", type=str, help="Path to a text file containing FASTA file names (one per line) in ./fasta/")


    group_patterns = parser.add_mutually_exclusive_group(required=True)
    group_patterns.add_argument("--patterns","-p", nargs="+", help="List of regular expressions for searching  motifs")
    group_patterns.add_argument("--patterns-list","-pl", type=str, help="Path to a text file with regular expressions (one per line)")

    # parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads (default: 4)")
    parser.add_argument("--quiet", "-q", action="store_true", help="Do not print any progress report during the task running.")
    parser.add_argument("--reverse","-r", action="store_true", help="Enable reverse complement strand search")
    parser.add_argument("--output-csv","-o", type=str, help="Path to output CSV file(optional)")
    args = parser.parse_args()

    patterns = load_patterns_from_file(args.patterns_list) if args.patterns_list else args.patterns
    if not patterns:
        print("Error: No patterns found to process!")
        return

    results = []
    if args.species or args.species_list:
        species_list = load_list_from_file(args.species_list) if args.species_list else args.species
        if not species_list:
            print("Error: No species found to process!")
            return

        for species in species_list:
            results.extend(process_fasta(species, patterns, mode="species", reverse=args.reverse, quiet=args.quiet))

    elif args.fasta or args.fasta_list:
        fasta_list = load_list_from_file(args.fasta_list) if args.fasta_list else [args.fasta]
        if not fasta_list:
            print("Error: No FASTA files found to process!")
            return
        for fasta_name in fasta_list:
            results.extend(process_fasta(fasta_name, patterns, mode="fasta", reverse=args.reverse, quiet=args.quiet))

    if not args.quiet:
        print("\n=== Final Summary ===")
        print("Fasta file\tPattern\tTotal Matches Counts\tMatched Bases\tTotal FASTA File Bases\tBases Occupancy Rate (%)")
        for fasta_file, pat, count, bases, fasta_size, rate in results:
            print(f"{fasta_file}\t{pat}\t{count}\t{bases}\t{fasta_size}\t{rate:.6f}")

    if args.output_csv:
        with open(args.output_csv, "w", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["Fasta file", "Pattern", "Total Matches Counts", "Matched Bases", "Total FASTA File Bases", "Bases Occupancy Rate (%)"])
            csv_writer.writerows(results)
        print(f"Results saved to {args.output_csv}")

if __name__ == "__main__":
    main()