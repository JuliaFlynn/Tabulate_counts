# BolonLab

This is the pipeline used in the paper "Comprehensive fitness maps of Hsp90 show widespead envrionmental dependence."  It takes raw fastq files generated on Illumina deep sequencing platforms and counts Hsp90 alleles.

01_doit_full.sh
01bc_freq_parse.pl
01_Hsp90_full.input

The 01_doit_full.sh script calls the 01bc_freq_parse.pl script and the input file and outputs all the Hsp90 reads:

#read in single end fastq file
#read in file containing information on how to parse samples including sequence for constant region for sequences that match constant region and have PHRED>cutoff
#output bc and index sequences to parsed files (one for each constant region
#NOTE: match updated to tolerate up to 2 snp - output as separate files *.snp0 *.snp1 and *.snp2
#output unparsed sequences as leftover.fastq

02bcfreq_parse2.pl calls 02_Hsp90_full.input and parses the Hsp90 reads into identifiers

#read in single 01seq file containing bc,index for samples matching constant region (e.g. parsed based on parent plasmid/construct)
#read in input file containing information on how to parse indexes and output a separate file for each index with bc only

bash 02b_sort_uniq.sh
sort2uniq.pl
02_Hsp90_full.input

The 02b_sort_uniq.sh script calls the sort2uniq.pl scrit and the input file and compiles all reads to unique barcodes.

#read in sorted sequence file (sequences only)
#count number of each unique sequence
#output each unique sequence and its number of occurences

03_tabulate_counts.pl calls the 03_hsp90_full.input file and the barcode-ORF assembly file created from the paired-end scripts and outputs the tabulation of reads of each barcode

#read in uniq files of barcodes for different timepoints
#read in barcode-ORF assembly files
#write counts file organized by position, amino acid, codon

