# FishTEA: Zebrafish transposable element analyzer pipeline

Analysis of RNA-seq and DNA-seq data for comparison of overlapping differentially expressed genes and transposable elements.
Also contains a script for chromatin accessibility to display open/closed regions of DanRer11 genome.
All scripts use Python v 3.11.


# Scripts and files
All scripts and files can be found here >> https://github.com/alicegodden/FishTEA/tree/main/scripts 

Or to clone:

$ git clone https://github.com/alicegodden/FishTEA/

$ cd FishTEA # navigate to FishTEA directory

$ python FishTEA.py # to start using FishTEA 




# Preparatory work (DNA_FishTEA.py)
Analysis of TE's in the genome to be able to plot a phenogram. Can also take the chrom and pos columns of the retroseq vcf file input to this, and put this through Biomart online and use the co-ordinates tool, to generate a list of genes that overlap with the TE insertion sites. 
Biomart > Filters > Region > Co-ordinates

![Screenshot 2024-05-22 at 16 16 54](https://github.com/alicegodden/FishTEA/assets/136358959/ea7a7399-5f64-4e55-a8f6-fc4e6ff9b84a)


# FishTEA.py
Analysis of lists of genes and TE's to analyse for danRer11. This script has several steps-
1. Adding TE location co-ordinates
2. Matching overlapping genes and TEs
3. Phenogram plotting of overlapping genes and TEs
4. & 5. Plotting bar charts counting TE family/class counts

 Note-  This pipeline needs file "GRCz11_Ensembl_rmsk_TE.gtf" to run. You can make your own TE gtf file, or use TETranscripts pre-made repository from the Hammell lab:
 https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/ 
 You may need to change chromosome numbering from chr1 > 1 if using Ensembl reference genomes. 
 The file of TEs (sigTEs_positions.csv) to match needs to look like this and be .csv:

TE
ERV3_DR-I:ERV1:LTR
L1-2_DR:L1-Tx1:LINE
Gypsy-21-LTR_DR:Gypsy:LTR

The file of genes (sigGenes.csv) needs to look like this and be .csv:
Gene_ID,Gene_start,Gene_end,Gene_name,Gene_chromosome
ENSDARG00000000563,42732992,42873700,ttn.1,9
ENSDARG00000001154,43941313,44015210,rimbp2,8

Missing information can be generated through use of Biomart. 

Data for centromere (chrcen.txt) and chromosome length (chrom_end.txt) are annotated at the end of the pipe and attached here as files.

