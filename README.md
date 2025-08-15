# *A. junii* genomic analysis #
This repository outlines the analysis pipeline for the paper *Population Genomics and Host Adaptation of Acinetobacter junii: Integrated Genomic and Transcriptomic Analyses Reveals Novel Virulence Factors and Bloodstream Infection Mechanisms （in review)*


## Remove adapter sequence ##
java -jar /home/jason/Bioinfo/Training/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE raw_read_1.fq.gz raw_read_2.fq.gz -baseout clean_read.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## Genome assembly and QC assessment ##
High-quality reads were assembled using SPAdes (https://github.com/ablab/spades)  
```
python spades.py --pe1-1 file1 --pe1-2 file2 -o assmebly --careful -k 21,33,55,77,99,127
```


