This directory contains data to test figeno. There are 5 _config.json files which can be used to generate example figures, using data contained in this directory. Please note that the input files were subsetted to only contain regions of interest (in order to save space), and that relative paths were used in the config files, so in order to generate the figures you either have to start figeno from this directory, or change the paths to absolute paths. In order to generate the figures, you can either
- run `figeno make GDM1_config.json` (or any other of the config files) while in this directory.
- run `figeno gui` from this directory, click "Load config", select the config file that you want, and click "Generate figure".

## GDM1 (nanopore)
GDM1_subset.bam contains nanopore reads for the leukemic GDM-1 cell line, subsetted to regions around the t(6;7) breakpoint (6:135,505,353 to 7:156,812,311; hg19). The full data can be downloaded from the SRA (accession SRR28257102). GDM1_config.json will show allele-specific DNA methylation at the MNX1 locus, and GDM1_splitread_config.json will show both sides of the breakpoints, with lines connecting different alignments coming from the same splitread.

## LNCaP (HiC)
HiC data for the LNCaP cell line, subsetted around the t(7;14). This data is from the [NeoLoopFinder manuscript (Wang et al. 2021, Nature Methods)](https://www.nature.com/articles/s41592-021-01164-w), and the full data can be downloaded from GEO accession [GSE161493](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161493). LNCaP_config.json will display HiC data around the breakpoint, as well as a DNase-seq bigwig track.

## THP1 (WGS)
Processed whole genome sequencing data from the THP-1 cell line (WGS data from the cancer cell line encyclopedia (SRA accession SRR8670675). THP1_SV.vcf contains structural variants in vcf format, THP1_ratio.txt and THP1_CNVs contain CNA information produced by Control-FREEC. THP1_circos_config.json will generate a circos plot with all chromosomes, while THP1_symmetrical_config.json will produce a symmetrical plot with a few chromosomes.
