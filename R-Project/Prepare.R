##########################
#  Date: 190404
#  User: PR
#  Purpose: Load files
##########################

library(data.table)
library(stringr)
library(bit64)
library(seqinr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
source("R-Project/Functions.R")

##LOAD TABLES

# 1. Load Evidence
evidence <- fread("txt\\evidence.txt")

# 2. Load Summary
Summary <- fread("txt\\summary.txt")
Summary <- Summary[Experiment != ""]

# 3. Load MSMS 
MSMS <- fread("txt\\msms.txt")


##LOAD FASTA

FASTA <- fasta_to_datatable(paths = "PR220121_GreatApes_Enamelome_aligned.fasta")
#Generate gapless sequences
FASTA[, Gapless := str_remove_all(Sequence, "-")]
# Add identifier (number 1 to N)
FASTA[,id := 1:.N]
# generate site level fasta site by cutting the sequence after each a.a. and making it a new row in a new table
FASTAsites <- FASTA[, .(AA = unlist(str_split(Sequence, "")),
                        AlignedPos = 1:nchar(Sequence)),
                    by = .(Protein, Annot, Gene, Species, id)]
#Add gapless position (add +1 for all "_")
FASTAsites[, GaplessPos := cumsum(AA != "-"), by = .(Annot)]
#Add identifier for each amino acid: gene_aligned position
FASTAsites[, GeneSite := paste(Gene, AlignedPos, sep = "_")]

