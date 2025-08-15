# *A. junii* genomic analysis #
This repository outlines the analysis pipeline for the paper *Population Genomics and Host Adaptation of Acinetobacter junii: Integrated Genomic and Transcriptomic Analyses Reveals Novel Virulence Factors and Bloodstream Infection Mechanisms （in review)*

## Genome assembly and QC assessment ##
High-quality reads were assembled using SPAdes (https://github.com/ablab/spades)  
```
python spades.py --pe1-1 file1 --pe1-2 file2 -o assmebly --careful -k 21,33,55,77,99,127
```
Species confirmation by GTDB-Tk (https://github.com/Ecogenomics/GTDBTk)  
```
