Folder description: This folder contains the test files and the results shown if you run this script successfully.

file：
`fastalist_test.txt`: contain the name lists of .fa files, the .fa files are pre-stored in `./fasta/` folders.

`regexlist_test.txt`: contain the list of the regular expressions, this file is an example of regular expression retrieval with two different i-motif structures.

`specieslist_test.txt`:contain the genomic name list of species in UCSC.

`./fasta/`: store the downloaded genomic .fa.gz files from UCSC and pre-saved .fa files.

`./bed_files/` and the `./*_test_result.csv`: the directory saved the bed files and summary results that get after running the following command:

```bash
python ./UcscGenomeRegexScanner.py -sl specieslist_test.txt -p "([cC]{3,}\w{1,12}){3,}[cC]{3,}" -r -o specieslist_test_result.csv
```

```bash
python ./UcscGenomeRegexScanner.py -fl fastalist_test.txt -p "([cC]{3,}\w{1,12}){3,}[cC]{3,}" -r -o fastalist_test_result.csv
```

