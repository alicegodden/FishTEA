# FishTEA pipeline

Analysis of RNA-seq and DNA-seq data for comparison of overlapping differentially expressed genes and transposable elements.
Also contains a script for chromatin accessibility to display open/closed regions of DanRer11 genome.
All scripts use Python v 3.11.


# Files

DNA_FishTEA.py - Analysis of TE's in the genome, plot a phenogram. Can also take the chrom and pos columns of the retroseq vcf file input to this, and put this through Biomart online and use the co-ordinates tool, to generate a list of genes that overlap with the TE insertion sites. 
Biomart > Region > Co-ordinates

![Screenshot 2024-05-22 at 16 16 54](https://github.com/alicegodden/FishTEA/assets/136358959/ea7a7399-5f64-4e55-a8f6-fc4e6ff9b84a)


FishTEA.py - Analysis of lists of genes and TE's to analyse. This script has several steps-
1. Adding TE location co-ordinates
2. Matching overlapping genes and TEs
3. Phenogram plotting of overlapping genes and TEs
4. & 5. Plotting bar charts counting TE family/class counts


chromatin_phenogram.py - Analysis of open regions in danRer11 genome

socstress_phenogram_TE_GENES.py - Phenogram plotting of TEs that overlap with genes by class/family.
