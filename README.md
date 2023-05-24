# ProteinSequenceAssembly

This open-sourced sequence assembly pipeline enable site-based sequence reconstruction to generate consensus sequences directly from MaxQuant´s output tables. It provides faster, more reproducible, and transparent data analysis processes. The pipeline was developed by Patrick L. Rüther.

# Input
1. MaxQuant´s txt folder (summary.txt, msms.txt, evidence.txt)
2. Aligned version of the the FASTA file used for the MaxQuant´s search (aligned at the gene level)
3. Exclusion list as a .csv file (facultative)

# Principle

1.	The aligned FASTA file is loaded and expanded from a protein level table to a site level table where each row corresponds to an amino acid position in a protein sequence. The amino acids are thus referenced using the “aligned” and “gapless” positions. Therefore it is easily possible to go from the position referenced in MaxQuant (“gapless”) to a gene position “aligned”. 
2.	Experimental annotations are extracted from the summary table and added to the evidence and msms.txt files. It is recommended to set up the Max Quant search with the experiment and fraction information. 
3.	Unique peptide sequences are extracted from the msms.txt table and mapped to all possible protein matches in the database. The table is expanded and each peptide sequence - protein pair is referenced with the start and end position of the peptide at the protein (gapeless) and gene (aligned) levels. The razor principle is then applied to peptide sequences matching multiple genes, where the referred peptides are assigned to the gene with the most peptide identifications. 
Output: PepSeqs.csv → gather unique peptide sequences with unique gene mapping. 
4.	The delta score is recalculated from the msms.txt table for the PSMs where the second best identification differs from the first one by a I/L, G/E or N/D. In those cases the delta score is reprocessed as the difference between the score of the first best identification with the score of the third best identification. 
5.	Facultative step: Identification of the optimal delta score cut off. The algorithm will loop with different delta score cutoffs and apply the subsequent FDR at 1%. A plot representing the number of PSMs as a function of the delta score will be displayed. The optimal delta score cutoff has to be chosen as the maximum of the curve to maximise the identifications. 
6.	The chosen delta score cutoff is applied to the msms.txt table, and the identifications are filtered to reach a FDR of 1% per sample. The delta score has been optimised for the Paranthropus samples at 15, that value can be modified in the script. 
7.	The evidence table is filtered based on the re-calculated delta score and FDR. The reverse sequences, contaminants and sequences found in the exclusion list are filtered out. Identifications (peptide sequences, modified peptide sequences and PSMs) are calculated per experiment and per sample. 
Output: IDbyExp.csv → table with the identifications per experiment. IDbySample.csv → table with the identifications per sample.
8.	Entries without genes are filtered out. Entries from the evidence file are collapsed per charge state to obtain unique modified peptides per row. The table is then expanded to the site level. Modified and unmodified amino acids are then collapsed to obtain unique unmodified amino acids. Information from overlapping peptides is collapsed on the site level. Absolute amino acid coverage per sample and per experiment can be calculated. 
Output: AACoverageRawFile.csv → table with the total amino acid coverage per experiment. AACoverageSample.csv → table with the total amino acid coverage per sample. sites.csv → table gathering all the covered positions per sample with the median and maximum score, the average mass error, the maximum log10 intensity of all peptides covering that site in that sample, the PSM count and the list of peptide sequences covering that site. 
9.	For a given sample, more than one amino acid can be identified covering the same position on the gene sequence. Those positions have to be detected to generate a consensus by only keeping the amino acid with the highest coverage. The relative coverage (PSM count) and the relative intensity of the conflicting sites are generated through the ratio between the site count/intensity of the given amino acid and the total count/intensity of the site for one sample. The consensus is based on the most common variant. 
Output: VarSites.csv → Display the list of conflicting sites with their relative coverage. 
10.	An empty canva is generated where the gaps are filled with dashes. The sites are then pasted together into a consensus sequence by sample.
11.	Reference sequences from the database are added to the consensus. The sequences can be exported as FASTA files. 
Output: Hominid-Consensus_Ref_MQ.fasta → Consensus sequences for all the samples aligned with the reference sequences.  Paranthropus-Consensus_Ref_MQ.fasta → Consensus sequences of given samples aligned with the reference sequences. Paranthropus-Consensus_MQ.fasta → Consensus sequences of given samples only. 

# How to get started ? 
1. Create a folder on the computer including the full MaxQuant´s txt folder, the aligned version of the database used for the search, the exclusion list (facultative) and the R-project folder (1). The R-Project folder contains the R-Project and 3 distinct R scripts: Functions.R, Prepare.R and Reshape.R (2).

 ![image](https://github.com/ClaireKoenig/ProteinSequenceAssembly/assets/134442809/de8a2517-b511-4cf1-882f-e8b9e56bde1c)


2. The exclusion list is a csv file with 4 columns: Sequence, Sample, Gene, Comment. The sequence and the sample columns are necessary for the script to run. The gene and comment column were designed for traceability in the data analysis. 

![image](https://github.com/ClaireKoenig/ProteinSequenceAssembly/assets/134442809/182e2547-7a89-48b0-8f58-d157a1ef88be)
 
3. The database used as an input has to be aligned on the gene level. 

![image](https://github.com/ClaireKoenig/ProteinSequenceAssembly/assets/134442809/ce7758aa-c40a-44ec-9b5d-4eefe8c21fef)

# Set the R environment

1. Open the 3 scripts: Function.R, Prepare.R and Reshape.R. 
Function.R is used for loading functions that will be used in the other files. Prepare.R is loading the relevant MaxQuant tables for further analysis as well as the database, which is then indexed with gene level and protein level positions. Reshape.R is processing the tables and generating the output files. 
2. In the Prepare.R script, change the name of the database (line 30). 
3. In the Reshape.R scprit, run line 8 to line 14. 
All the relevant packages are loaded. The functions used in the script are imported and the MaxQuant´s output tables are loaded in the R environment. The aligned FASTA file is imported and expanded from a protein level table to a site level table where each row corresponds to an amino acid position in a protein sequence. The amino acids are thus referenced using the “aligned” and “gapless” positions. Therefore, it is easily possible to go from the position referenced in MaxQuant (“gapless”) to a gene position “aligned”. 

 # Perform the sequence assembly
 
 1. Experimental annotation from summary table: 
The script is creating an ExpAnnot table. It is recommended to set up the Max Quant search with the experiment and fraction information. The structure of the ExpAnnot table should be the following: 

![image](https://github.com/ClaireKoenig/ProteinSequenceAssembly/assets/134442809/68be242e-f9a2-4f62-8bb3-6280b221a2f6)

2. Annotate tables:
Run lines 39 and 40 to add the annotations to the msms and evidence tables. 

3. Annotate peptide sequences:
Run lines 44 to 102 to generate the PepSeqs table. 

4. Delta Score recalculation and cutoff optimization:
Run lines 107 to 124 to recalculate the delta score. 

! FACULTATIVE STEP: 
Run lines 128 to 161 to assess the optimal delta score cutoff for the given dataset. Line 141 should be modified based on the sample the delta score cut-off has to be optimized for. The optimal delta score cut-off has to be chosen as the maximum of the curve to maximise the identifications. 
Change line 167 to indicate the optimal delta score to apply. Run lines 166 to 176 to apply the new delta score cut-off and the PSM FDR at 1%. 
Line 178 is project specific and should not be run.  

5. Transfer the filtered sequences from the MSMS table to the evidence table and collapse to site level to combine information from overlapping peptides on 1 site:
Run from lines 183 to 192 to filter the evidence table based on the previous filtering and remove reverse hits and contaminants. 
FACUTLATIVE STEP: 
Filtering based on the exclusion list is facultative. If you have an exclusion list, run lines 195 to 198. 

If you want to skip that step, just go straight to line 200 to export the filtered evidence table. 
Run from lines 203 to 256 to calculate the number of identifications per samples and per experiment.
Lines 158 to 260 are specific to the example project and should not be run on another data set.
Run lines 263 to 275 to collapse the peptide identifications to unique amino acid positions. 

6. Generate the consensus sequence
Run lines 280 to 320 to generate the consensus sequences.

7. Generate FASTA format outputs 
Run lines 326 to 343 to generate and export the consensus sequences as FASTA files. 
The outputs can be renamed line 335 for option 1, 339 for option 2 and 343 for option 3. Option 2 and 3 are selectively only exporting the Paranthropus sequences from the example dataset. The text string to extract can be changed lines 339 and 343 for option 2 and 3 respectively. 

# Test data

The test data used here represents the first successful mass spectrometric sequencing of enamel proteomes from four ca. 2 million years old dental specimens morphologically identified as P. robustus from the site of Swartkrans. The Paranthropus raw files were searched with the corresponding lab blanks, modern homo sapiens, gorilla and pongo reference files, as well as raw files corresponding to homo antecessor and gigantopithecus samples. 

The raw files were searched on MaxQuant (V. 1.6.0.17). Unspecific digestion was selected. Cysteine trioxidation was set as a fixed modification, and Arg->Ornithine, Gln->pyro-Glu, Glu->pyro-Glu, Deamidation (NQ), Phospho (ST) and Oxidation (MPW) were set as variable modifications. The minimum peptide length was set at 6 amino acids for unspecific search and the maximum was set at 20 amino acids. The peptide mass was limited to 3500 Da. The minimum Andromeda score for both modified and unmodified peptides was fixed at 25. Protein, PSM and site FDR were set at 10% and were manually adjusted during the data analysis. The minimum delta score for modified and unmodified peptides was set at 0 and was manually adjusted afterwards. 

# Setting up the search on MaxQuant

1. Use a database without gaps (non-aligned).
2. Set the protein, PSM and site FDR to 10%.
3. Set the minimum delta score for modified and unmodified peptides to 0.
4. Set the experiment names for each raw file.
5. Set the fraction number for each raw file.
