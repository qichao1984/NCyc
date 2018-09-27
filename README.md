# NCyc：a curated integrative database for fast and accurate metagenomic profiling of nitrogen cycling genes
Nitrogen (N) cycle is one of the most important biogeochemical cycles in the Earth ecosystem and has gained extensive foci in ecology and environmental studies. During the past decades, extensive efforts have been applied to characterize N cycle pathways in various ecosystems using different approaches. Until recently, shotgun metagenome sequencing have been applied to explore N cycle related gene families and link them with environmental processes. However, there are problems in applying publically available orthology databases to profile N cycle gene families in shotgun metagenomes, such as inefficient database searching, unspecific orthology groups, and low coverage of N cycle genes and/or gene (sub)families in large-scale orthology databases. A curated high quality reference database is therefore in urgent need for metagenomic analysis of N cycle gene families. To solve these issues, we built a manually curated integrative database (NCyc) for fast and accurate profiling of N cycle gene (sub)families from shotgun metagenome sequencing data. 

Here, protein sequences for NCyc gene families were recruited from multiple public databases such as UniProt, KEGG, COG, eggNOG, and the SEED. The NCyc database contains a total of 68 gene (sub)families and covers  eight N cycle processes with 84,759 and 219,146 representative sequences at 95% and 100% identity cutoffs, respectively. We also identified 1,958 homologous orthology groups and included corresponding sequences in the database to avoid false positive assignments due to “small database” issues. 

Three files are generated for NCyc:

<b>1. NCyc_100.faa.gz</b>: fasta format representative sequences obtained by clustering curated sequences at 100% sequence identity. This file can be used for "BLAST" searching NCyc genes in shotgun metagenomes.

<b>2. NCyc_95.faa.gz</b>: fasta format representative sequences obtained by clustering curated sequences at 95% sequence identity. This file can be used for "BLAST" searching NCyc genes in shotgun metagenomes.

<b>3. id2gene.map</b>: a mapping file that maps sequence IDs to gene names, only sequences belonging to NCyc gene families are included. Sequences for NCyc homologs are not included. This file is used to generate NCyc profiles from BLAST-like results against the NCyc database. 



<b>Example for using NCycProfiler.PL:</b>

perl NCycProfiler.pl -d \<workdir\> -m \<diamond|usearch|blast\> -f \<filetype\> -s \<seqtype\> -si \<sample size info file\> -rs \<random sampling size\> -o \<outfile\>
  
Detailed explanations: 

-d : specify the directory where your fasta/fastq (or gzipped) files are located. 

-m : specify the database searching program you plan to use, currently diamond, usearch and blast are supported. 

-f : specify the extensions of your sequence files, e.g. fastq, fastq.gz, fasta,fasta.gz, fq, fq.gz, fa, fa.gz

-s : sequence type, nucl or prot

-si: a tab delimited file containing the sample/file name and the number of sequences they have

-rs: specify the number of sequences for random subsampling, if not specified, the lowest number in -si will be used.

-o : the output file for N cycle gene profiles.   
  
