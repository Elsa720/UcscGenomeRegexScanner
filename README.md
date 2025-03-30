# UcscGenomeRegexScanner

## Description

`UscsGenomeRegexScanner.py` is a script used to scan multiple genome or user-defined FASTA files in the specified regular expression pattern. It supports both sense and anti-sense strand matching and generates a standard BED file to record the matching results.

### Features

- Support simultaneous retrieval of multiple genomes from UCSC or user-defined FASTA files.
- Supports simultaneous matching of multiple regular expression modes, in which regular expressions are case-insensitive.
- Support antisense chain search, and IUPAC fuzzy base code is considered.
- Support UCSC genome files downloaded automatically from UCSC website.
- Generate a standard BED file with matching sequence information.
- Output the matching statistics, including matching counts, matching bases, and base occupancy rate. The results can be output to a csv file.

## Usage

### Command line options

```bash
-h, --help            show this help message and exit
  --species SPECIES [SPECIES ...], -s SPECIES [SPECIES ...]
                        Name list of species (e.g., hg38 mm10)
  --species-list SPECIES_LIST, -sl SPECIES_LIST
                        Path to a text file with species names (one per line)
  --fasta FASTA, -f FASTA
                        Name list to a specific FASTA file in ./fasta/ path
  --fasta-list FASTA_LIST, -fl FASTA_LIST
                        Path to a text file containing FASTA file names (one per line)
  --patterns PATTERNS [PATTERNS ...], -p PATTERNS [PATTERNS ...]
                        List of regular expressions for searching motifs, the case is ignored
  --patterns-list PATTERNS_LIST, -pl PATTERNS_LIST
                        Path to a text file with regular expressions (one per line)
  --quiet, -q           Do not print any progress report during the task running(optional, default FALSE).
  --reverse, -r         Enable reverse complement strand search(optional, default FALSE).
  --output-csv OUTPUT_CSV, -o OUTPUT_CSV
                        Path to output CSV file(optional)

Optional arguments notes:
1. If the argument of --species/-s or --species-list/-sl is chosen,
    - the species name involved must be a species name that can be retrieved from the genome fasta file on the ucsc website to successfully access genome files.
    - if run successfully, downloaded species genome .fa.gz files will be stored in the./fasta/ directory.
2. In order to expand the usability of this script, the argument of --fasta/-f and --fasta-list/-fl to allow users to customize the upload fasta file for sequence retrieval.
    - For custom fasta files, you have to create a new ./fasta/ directory and store the interested fasta files in it first to let the script recognize successfully.
    - If the argument of --fasta/-f and --fasta-list/-fl is chosen, you only need to follow the argument with the name of the file in .fa or .fa.gz format instead of the entire file path.
    - To run the script successfully, one of the four arguments --species/-s, --species-list/-sl, --fasta/-f and --fasta-list/-fl must be selected.
3. If the argument of --reverse or -r is chosen,
    - the coordinate information of the antisense chain in bed file maps the coordinate information of the justice chain in the direction of 5'-3', and the sequence information of the antisense chain in bed file is output in the direction of 5'-3'.
4. If the argument of --output-csv or -o is chosen,
    - the final summary included the matching counts, matching bases, whole genome/fasta file bases and base occupancy rate(%) will be exported to a csv file.

```

### Output file

Bed files are automatically saved in the `./bed_files/` directory. Each regular expression pattern generates a BED file in the format of `[genome or FASTA file name]_[simplified regular expression pattern].bed`,the columns contains:

```bash
1. Chromosome
2. Start of the match sequence
3. End of the match sequence
4. Length of the match sequence
5. Strand
6. Matched sequence as it appears on the forward strand
```

In order to correctly identify the position coordinates of the regular expression, the input fasta file must be in `>chrX:X-X` format; if `>chrX`, it is recognized as genomic information, and the default count starts from 0.
The coordinate information of the antisense chain maps the coordinate information of the justice chain in the direction of 5'-3', and the sequence information of the antisense chain is output in the direction of 5'-3'.

The csv file will be exported automatically after selecting the `--output-csv` or `-o` arguments + `output file path`. The file will summarize all the matching results of each fasta file, which contains:

```bash
1. Fasta file: genome file or fasta file user-defined.
2. Pattern: regular expressions user-defined.
3. Total Matches Counts: the number of each regular expressions retrieved in each fasta file
4. Matched Bases: the number of sequence bases that match the regular expression
5. Total FASTA File Bases: the number of FASTA file sequence bases
6. Bases Occupancy Rate (%): 100 * Matched Bases/Total FASTA File Bases, the occupancy rate of the regular expression
```

### Basic Example

1. Scan multiple UCSC species genomes

``` bash
python ./UcscGenomeRegexScanner.py --species hg38 mm39 --patterns "C{3,}" "G{4,}" --reverse
```

At present, most of the genomes of eukaryotes and bacteria can be successfully downloaded.
2. Scan local FASTA files

``` bash
python ./GenomeRegexScanner.py --fasta example --patterns "C{3,}" "G{4,}" --reverse
```

Make sure the `example.fa` or `example.fa.gz` file can be found in `./fasta/` directory before running this command.
3. Load species and regular expressions from file (if a lot of searching is required)

``` bash
python  ./GenomeRegexScanner.py  --species-list  ./species.txt --patterns-list ./patterns.txt  --reverse
```

`./species.txt` or `./patterns.txt` is a text file with species names or regular expressions, one per line.
4. Save the results to a CSV file

``` bash
python  ./GenomeRegexScanner.py --species hg38 --patterns "C{3,}" --output-csv ./results.csv
```

The result contained the bases occupancy infomation will be save to the `./results.csv`.

Other Matters Needing Attention:

- **Network connection:** If scanning the UCSC genome, make sure the network connection is accessable.
- **File path:** Ensure that the fasta and bed files directories have write permission.
- **Regular expression:** Supports standard regular expression syntax, ignoring case.

## Download, Installation and Requirements

Make sure you have the following dependencies installed in your environment:

- Python 3.6 or later
- Biopython (for parsing the fa file's sequence information)
- wget (for downloading genome files)
  
You can install the biopython package with the following instructions:

``` bash
pip install biopython
or
conda install -c conda-forge biopython
```

## Acknowledgment and Contribution

This script refers to the previous script [fastaRegexFinder.py](https://github.com/dariober/bioinformatics-cafe/tree/master/fastaRegexFinder) on the regular expression of the reverse chain. Thanks for the contribution to this script!

If you have any other questions or suggestions, please feel free to open an issue. We appreciate everyone's contribution!
