# *A. junii* genomic analysis #
This repository outlines the analysis pipeline for the paper *Population Genomics and Host Adaptation of Acinetobacter junii: Integrated Genomic and Transcriptomic Analyses Reveals Novel Virulence Factors and Bloodstream Infection Mechanisms （in review)*


## Remove adapter sequence ##
```
java -jar /home/jason/Bioinfo/Training/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE raw_read_1.fq.gz raw_read_2.fq.gz -baseout clean_read.fq.gz ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

## Genome assembly and QC assessment ##
High-quality reads were assembled using SPAdes (https://github.com/ablab/spades)  
```
python spades.py --pe1-1 file1 --pe1-2 file2 -o assmebly --careful -k 21,33,55,77,99,127
```

## Genome annotation ##   
Virulence factors (VFs) identification using core setA VF data ```VFDB_setA_pro.fas``` by abriate (https://github.com/tseemann/abricate)  

Antibiotic resistance genes (ARGs) identification using RGI (https://github.com/arpcard/rgi?tab=readme-ov-file)
```
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
rgi load --card_json /path/to/card.json --local
rgi main --input_sequence nucleotide_input.fasta --output_file output_file --local --clean
```

CPS region prediction using Kaptive (https://github.com/klebgenomics/Kaptive)
```
kaptive assembly database assemblies/*.fasta -o kaptive_results.tsv -f
```

Pangenome analysis using prokka (https://github.com/tseemann/prokka) and roary (https://sanger-pathogens.github.io/Roary/)  
```
prokka /path/to/"$sample".fasta --quiet --outdir /path/to/prokka_output/"$sample" --force --prefix $sample
roary –f output_dir *.gff
```

Secretion systems annotation by MacSyFinder

```
macsyfinder --db-type ordered_replicon --sequence-db seqence.faa --models TXSScan -o txss/outdir -w 8
```

## Genome population ##
```
library(fastbaps)
library(ape)

fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
sparse.data <- import_fasta_sparse_nt(fasta.file.name)
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
baps.hc <- fast_baps(sparse.data)
clusters <- best_baps_partition(sparse.data, as.phylo(baps.hc))
```

## GWAS analysis ##
The GWAS was conducted simultaneously on the pangenome and k-mer association.
Pangenome was retrieved from roary and the gene_presence_absence.csv and traits file were input into Scoary (https://github.com/AdmiralenOla/Scoary).
```
scoary -g gene_presence_absence.csv -t traits.csv -c BH --no_pairwise
```
Pyseer was used for K-mer-based analysis and the tutorial can be accessed in (https://pyseer.readthedocs.io/en/master/tutorial.html#k-mer-association-with-mixed-effects-model)

## Phylogenetic analysis ##
All public A. junii sequence was deposited in genome_fna directory.
core genome alignment was generated using snippy (https://github.com/tseemann/snippy), and recombination sites were removed with Gubbins (https://github.com/nickjcroucher/gubbins). A maximum-likelihood phylogenetic tree was then constructed using IQ-TREE (http://www.iqtree.org/) based on clean core genome SNP alignments.
```
snippy --outdir mut1 --ref ref.gbk --ctgs fna/seq1.fasta
run_gubbins.py -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
iqtree -s clean.core.aln --boot-trees --wbtl -m GTR+I+G -B 1000 -nt 18
```
## Plasmid analysis ##
All public plasmid sequence was deposited in plasmid_fna directory
Plasmid sequence distances are calculated using Mash (https://github.com/marbl/Mash).
```
mash triangle -E -s 5000 -k 13 /path/to/fasta.fa > /path/to/edgelist.tsv
```

The community detection was generated based on similarity using the Louvain algorithm (https://github.com/taynaud/python-louvain)
```
usage: louvain_community.py [-h] -i INPUT -o OUTPUT [--resolution RESOLUTION]

Calculate Louvain communities from Mash distance results.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input Mash distance results file (tab-delimited).
  -o OUTPUT, --output OUTPUT
                        Output file for Louvain community results.
  --resolution RESOLUTION
                        Resolution parameter for Louvain algorithm (default:
                        1.0)
```
Acinetobacter plasmid typing was based on three replicons families: Rep_3 (PF01051), replicase (PF03090), Rep_1 (PF01446), and RepC (PF06504).
```
hmmscan --tblout /replicon/ourput.tbl -E 1e-5 Database/00.PfamDB/pfam.hmm seq.faa
```

## Transcriptome analysis ##

Different enriched sequence (DES) analysis
```
bwa indx ref.fna
bwa mem ref.fna R1.fq.gz R2.fq.gz > file.sam

hisat2-build reference_fna output_basename
hisat2 hisat_index -1 RNA_1.fq.gz -2 RNA_2.fq.gz -S output.sam

samtools view -S file.sam -b -o file.bam
samtools sort file.bam -o file.sort.bam
featureCounts -a ref.gtf -p -t gene -o Deq_all file.sorted.bam
```
RNA sequencing figure plot
```
library(ggplot2)

data <- read.table("data.txt", sep = "\t")


colnames(data) <- c("location", "log2fold", "a1", "a2", "a3")


ggplot(data, aes(x = location, y = log2fold)) +
  geom_point(aes(color = log2fold > 0), size = 1, shape = 21) +
  scale_color_manual(values = c("cadetblue2", "brown2")) +
  theme_minimal() +
  geom_point(data = subset(data, log2fold > 2 | log2fold < -2),
             aes(x = location, y = log2fold, color = log2fold > 0),
             size = 1, shape = 21, fill = NA, stroke = 0.5, color = "black") +
  geom_text(data = subset(data, log2fold > 2 | log2fold < -2),
            aes(label = a1),
            vjust = -0.5, hjust = 0.5, size = 1, color = "black") +
  labs(x = "Location", y = "Log2 Fold Change", color = "Direction") +
  coord_cartesian(ylim = c(-5, 7)) +
  scale_y_continuous(breaks = seq(-5, 7, by = 2)) +
  theme_bw() +
  theme(aspect.ratio = 0.4)
```

