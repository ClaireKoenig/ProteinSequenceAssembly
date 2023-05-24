##########################
#  Date: 190404
#  User: PR
#  Purpose: Reshape
##########################

library(data.table)
library(stringr)
library(stringi)
library(seqinr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("..")
source("R-Project\\Prepare.R")

## 1. Experimental annotation from summary table 

ExpAnnot <- Summary[nzchar(`Enzyme mode`), .(Experiment, `Raw file`)]
setorder(ExpAnnot, `Raw file`)
ExpAnnot[, Sample := str_extract(Experiment, "^[^_]+_[^_]+")]
ExpAnnot[grepl("MSblank", Sample), Sample := paste("MSblank", 1:.N, sep = "_")]
ExpAnnot[grepl("HomoSapiens", Sample), Sample := "HomoSapiens_1"]
ExpAnnot[grepl("Pogo", Sample), Sample := "Pongo_1"]
ExpAnnot[grepl("(1848)|(14132)", Experiment), Sample := "Paranthropus_SK14132"]
ExpAnnot[grepl("(1847)|(SK830)", Experiment), Sample := "Paranthropus_SK850"]
ExpAnnot[grepl("(1846)|([^K]830)", Experiment), Sample := "Paranthropus_SK830"]
ExpAnnot[grepl("Paranthropus_SKX2", Sample), Sample := "Paranthropus_SK835"]
ExpAnnot[grepl("Paranthropus_SKX", Sample), Sample := "Paranthropus_SK835"]
ExpAnnot[grepl("Paranthropus_SKX2_single", Experiment), Experiment := "Paranthropus_SK835_single"]
ExpAnnot[grepl("Paranthropus_SKX_single", Experiment), Experiment := "Paranthropus_SK835_single"]
ExpAnnot[grepl("Paranthropus_SKX_frc", Experiment), Experiment := "Paranthropus_SK835_frc"]
ExpAnnot[, Sample := as.factor(Sample)]
ExpAnnot[, Fractionated := as.factor(grepl("frc", `Raw file`, ignore.case = T))]
levels(ExpAnnot$Fractionated) <- c("single", "frc")
ExpAnnot[Fractionated == "frc", Fraction := as.factor(as.integer(str_extract(`Raw file`, "\\d$")))]
ExpAnnot[, Experiment := as.factor(Experiment)]

## 2. Annotate tables

MSMS <- merge(MSMS, ExpAnnot, by = "Raw file")
evidence <- merge(evidence[, !c("Experiment", "Fraction")], ExpAnnot, by = c("Raw file"))

## 3. Annotate peptide sequences

if(length(list.files(".", pattern = "PepSeq")) == 1) {
  PepSeqs <- fread(list.files(".", pattern = "PepSeq"))
} else {
#Get all unique peptide sequence from MSMS.txt -> list of unique peptide sequences
PepSeqs <- unique(MSMS[, .(Sequence)])
# Find all proteins in the fasta (gapless sequences) that contain this peptide -> add second column with all the proteins containing that sequence
PepSeqs[, Proteins := lapply(Sequence, function(x) FASTA[stri_detect(Gapless,  regex = x), ]$Protein)]
# Take only peptides found in the FASTA and Make 1 row for each PepSequence-Protein pair (= expand table to long format)
PepSeqs <- PepSeqs[, .(Protein = unlist(Proteins)), by = .(Sequence)]
# Locate the peptide sequence in the respective FASTA gapless sequence -> add start and end position of the peptide from the gapless sequence
PepSeqs[, GaplessPos := mapply(function(x, y) str_locate_all(FASTA[match(x, Protein)]$Gapless, y), x = Protein, y = Sequence)]
# Only keep the start position of each peptide match
PepSeqs[, GaplessPos := lapply(GaplessPos, function(x) x[,1])]
# Repeat peptide with multiple matches in the same protein sequence (probably irrelevant for enamel)
PepSeqs <- PepSeqs[, .(GaplessPos = unlist(GaplessPos)), by = .(Sequence, Protein)]
# Add an identifier before conversion to site level
PepSeqs[, id := 1:.N]
# Transform to site level for mapping aligned positions -> expand table: 1 line per amino acid for each peptide-prot sequence
PepSeqs <- PepSeqs[, .(GaplessPos = GaplessPos:(GaplessPos+nchar(Sequence)-1)), by = .(Sequence, Protein, id)]
# Add aligned positions and gene
PepSeqs <- merge(PepSeqs, FASTAsites[AA != "-", .(Gene, Protein, AlignedPos, GaplessPos)], by = c("Protein", "GaplessPos"), all.x = T)
# Sanity check: Any peptides without aligned position?
PepSeqs[is.na(AlignedPos)]
# Go back to peptide level
PepSeqs <- PepSeqs[, .(GaplessPos = paste(GaplessPos, collapse = ","),
                       AlignedPos = paste(AlignedPos, collapse = ",")),
                   by = .(Gene, Protein, Sequence, id)]
# Filter down to unique peptide sequence and Gene
PepSeqs <- PepSeqs[, .(Proteins = paste(Protein, collapse = ";"),
                       GaplessPositions = paste(GaplessPos, collapse = ";")),
                   .(Sequence, Gene, AlignedPos)]
# Check for peptides that match multiple proteins per gene
PepSeqs[, Multimatch := .N, by = .(Sequence, Gene)]
# Get rid of unique weird alignment artifact
PepSeqs <- PepSeqs[order(Sequence, Gene, str_count(Proteins, ","))]
PepSeqs <- PepSeqs[, .(Proteins = paste(Proteins, collapse = ";"),
                       AlignedPos = AlignedPos[1],
                       GaplessPositions = GaplessPositions[1]), by = .(Gene, Sequence)]


# Calculate how often a peptide sequence is repeated in the table
PepSeqs[, Multiplicity := .N, by = .(Sequence)]
PepSeqs[Multiplicity > 2, .(Genes = paste(Gene, collapse = ";")), by = .(Sequence, Multiplicity)]
# Apply razor principle (winner takes all) to peptides mapping to multiple genes
PepSeqs[, AnyMultiPep := any(Multiplicity == 2), by = .(Gene)]
# Subset only peptides matching multiple genes and count how many peptides per gene in the entire dataset
MultiPeps <- merge(PepSeqs[Multiplicity == 2], PepSeqs[AnyMultiPep == T, .(PepCount = .N), by = .(Gene)],
                   by = "Gene")
# Find and keep "razor" gene with most peptides in the entire dataset
MultiPeps[, MaxPepCount := max(PepCount), by = .(Sequence)]
MultiPeps <- MultiPeps[PepCount == MaxPepCount]
MultiPeps[, .N, by = .(Sequence)][N > 1]
MultiPeps[, .N, by = .(Gene)]
# Combine the peptide sequences with unique gene mapping
PepSeqs <- rbind(PepSeqs[Multiplicity == 1], MultiPeps[, !c("MaxPepCount", "PepCount")])
#Export site table as csv
fwrite(PepSeqs, "PepSeqs.csv")
}


## 4. Delta Score recalculation and cutoff optimization

# Recalculate delta scores
MSMS[, DeltaBackup := copy(`Delta score`)]
# Split top 3 sequences
MSMS[, Sequences := str_split(`All sequences`, ";")]
# Split top 3 modified sequences
MSMS[, ModSequences := str_split(`All modified sequences`, ";")]
# Split top 3 scores
MSMS[, AllScores := str_split(`All scores`, ";")]
MSMS[, AllScores := lapply(AllScores, as.numeric)]

# Identify rank 1 and 2 hits that are identical except 1 leucine/isoleucine change
MSMS[, IL_hit := unlist(lapply(Sequences, function(x) str_replace_all(x[1], "I", "L") == str_replace_all(x[2], "I", "L")))]
MSMS[IL_hit == T, `Delta score` := unlist(lapply(AllScores, function(x) x[2] - x[3]))]

# Identify rank 1 and 2 hits that are identical except 1 Q/E or N/D
MSMS[, QE_hit := unlist(lapply(ModSequences, 
                               function(x) str_replace_all(x[1], c("Q\\(de\\)" = "E", "N\\(de\\)" = "D")) == str_replace_all(x[2], c("Q\\(de\\)" = "E", "N\\(de\\)" = "D"))))]
MSMS[QE_hit == T, `Delta score` := unlist(lapply(AllScores, function(x) x[2] - x[3]))]


# ! FACULTATIVE! Delta Score optimization plots - Sample-wise FDR at different delta score cutoffs - peptides, evidence hits, PSMs
library(ggplot2)
MSMS[, Subset := Sample]
setorder(MSMS, Subset, -Score)
DeltaOpt <- data.table(MinDelta = c(0:30, 10 * (4:10)))
DeltaOpt[, IDs := lapply(MinDelta, function(x) {
  tmp <- MSMS[`Delta score` >= x]
  setorder(tmp, Subset, -Score)
  tmp[, qValue := cumsum(Reverse == "+")/1:.N, by = .(Subset)]
  evidence[,.(EvidenceHits = .N,
              PSMs = sum(`MS/MS count`),
              Peptides = length(unique(Sequence)),
              Sequences = list(unique(Sequence))), 
           keyby = .(PassedMSMS_FDR = id %in% tmp[qValue <= FDR]$`Evidence ID`,
                     IsParanthropus = str_extract(Sample, "Paranthropus_SK\\d+"))][PassedMSMS_FDR == T]
})]
DeltaOpt[, Samples := lapply(IDs, function(x) x$IsParanthropus)]
DeltaOpt[, Peptides := lapply(IDs, function(x) x$Peptides)]
DeltaOpt[, EvidenceHits := lapply(IDs, function(x) x$EvidenceHits)]
DeltaOpt[, PSMs := lapply(IDs, function(x) x$PSMs)]
DeltaOptPlotting <- DeltaOpt[, .(Sample = unlist(Samples),
             Peptides = unlist(Peptides),
             EvidenceHits = unlist(EvidenceHits),
             PSMs = unlist(PSMs)),
         by = .(MinDelta)]
ggplot(melt(DeltaOptPlotting[!is.na(Sample), ], measure.vars = c("Peptides", "EvidenceHits", "PSMs"),
            variable.name = "Identification", value.name = "IDs-at-1%-FDR"),
       aes(x = MinDelta, y = `IDs-at-1%-FDR`, color = Identification)) +
  geom_line() +
  facet_grid(Identification~Sample, scales = "free") +
  theme_bw()

ggsave("Paranthropus_DeltaScore-Opt.pdf",
       width = 10, height = 8)


# Apply PSM FDR control with optimized delta score cutoff

FDR <- 0.01
MSMS_F <- copy(MSMS[`Delta score` >= 15, ])
MSMS_F[, Subset := Sample]
# Order the MSMS table by 1. sample 2. decreasing score
setorder(MSMS_F, Subset, -Score)
# Calculate q-Value
MSMS_F <- MSMS_F[, qValue := cumsum(Reverse == "+")/1:.N, by = .(Subset)]
# Select the last PSM that's under the q-value threshold
MSMS_F[, PassedFDR := 1:.N <= which.max(which(qValue <= FDR)), by = .(Subset)]
MSMS_F[, .N, by = .(Subset, PassedFDR)]
MSMS_F <- MSMS_F[PassedFDR == T]
# Addition of the manually recovered collagen17A1 sequence
MSMS_F <- rbind(MSMS_F, MSMS[46412,], fill=TRUE)


## 5. Transfer the filtered sequences from the MSMS table to the evidence table and collapse to site level to combine information from overlapping peptides on 1 site

#Filter the evidence table based on re-calculated delta scores and FDR
evidenceF <- evidence[id %in% MSMS_F$`Evidence ID`]
#FDR from the 1% applied on the MSMS table. Still need to remove the 1% reverse hits.
evidenceF <- evidenceF[Reverse != "+"]
evidenceF <- evidenceF[!all.contain(Proteins, pattern = "CON"),]

# Annotate evidence with pep locations and filter based on exclusion list
evidenceF <- merge(evidenceF[, !c("Proteins")], 
                      unique(PepSeqs[, .(Sequence, Gene, Proteins, AlignedPos, GaplessPositions)]), 
                      "Sequence", all.x = T, allow.cartesian = T)

# Facultative - if exclusion list is available
ExclusionList <- fread("CK_ExclusionListSample.csv")
ExclusionList[,Excl := paste0(Sequence, "_", Sample)]
evidenceF [, Excl := paste0(Sequence, "_", Sample)]
evidenceF <- evidenceF[!Excl %in% ExclusionList$Excl]

fwrite(evidenceF, "EvidenceFiltered.csv")

# IDs - generate table with PSMs, mod peptides and peptides count per experiment (ie. sample + frac or SS)
IDbyExperiment <- evidenceF[, .(PSMs = sum(`MS/MS count`),
                                      ModPeps = length(unique(`Modified sequence`)),
                                      Peptides = length(unique(Sequence))),
                                  keyby = .(Sample, Fractionated, Experiment)]
# IDs to long - generate long table with PSMs, mod peptides and peptides count per experiment (ie. sample + frac or SS)
IDbyExperimentLong <- melt(IDbyExperiment, 
                           measure.vars = c("PSMs", "ModPeps", "Peptides"), 
                           variable.name = "Identifications", 
                           value.name = "IDs")

# IDs by sample - generate table with PSMs, mod peptides and peptides count per sample
IDbySample <- evidenceF[, .(PSMs = sum(`MS/MS count`),
                                ModPeps = length(unique(`Modified sequence`)),
                                Peptides = length(unique(Sequence))),
                            keyby = .(Sample)]
# IDs to long - generate long table with PSMs, mod peptides and peptides count per sample
IDbySampleLong <- melt(IDbySample, 
                           measure.vars = c("PSMs", "ModPeps", "Peptides"), 
                           variable.name = "Identifications", 
                           value.name = "IDs")
#Export tables
fwrite(IDbyExperiment, "IDbyExp.csv")
fwrite(IDbyExperimentLong, "IDbyExpLong.csv")
fwrite(IDbySample, "IDbySample.csv")
fwrite(IDbySampleLong, "IDbySampleLong.csv")

# Collapse unique modified peptides = Collapse charge state
evidenceF[is.na(Intensity), Intensity := 0]
modpepsUnique <- evidenceF[!is.na(Gene), 
                           .(PSMcount = sum(`MS/MS count`),
                               Charges = paste(unique(Charge), collapse = "/"),
                               MedScore = median(Score, na.rm = T),
                               MaxScore = max(Score, na.rm = T),
                               AvgMass = mean(`m/z`),
                               AvgMassErrorPPM = mean(`Mass error [ppm]`),
                               MaxLogInt = max(log10(Intensity), na.rm = T)),
                           by = .(Sequence, Proteins, AlignedPos, `Modified sequence`, Gene, 
                                  Sample, Fractionated, Experiment)]
modpepsUnique[MaxLogInt < 0, MaxLogInt := 0]

# Expand to site level - add identifier
modpepsUnique[, id := 1:.N]
# For (i) all peptides mapped to the alignment, split the peptide sequence into AAs and split the aligned positions
sites <- modpepsUnique[!is.na(AlignedPos), 
                       .(AA = unlist(str_split(Sequence, "")),
                         Site = as.integer(unlist(str_split(AlignedPos, ",")))),
                       by = .(Proteins, Gene, MedScore, MaxScore, AvgMassErrorPPM, MaxLogInt, 
                              PSMcount, id, Sequence, `Modified sequence`, Sample, Fractionated, Experiment)]
# Annotate individual sites
sites[, GeneSite := paste(Gene, Site, sep = "_")]
# Calculate sequence coverage per raw file
sites[, .(AbsCoverage = length(unique(GeneSite))), by = .(Experiment)]
AACoverageRawFile<-sites[, .(AbsCoverage = length(unique(GeneSite))), by = .(Experiment)]
fwrite(AACoverageRawFile, "AACoverageRawFile.csv")
# Rescue PTM-based false identifications
sites[GeneSite == "ENAM_174" & AA == "D", AA := "W"]
sites[GeneSite == "AMELX_33" & AA == "D", AA := "N"]
sites[GeneSite == "AMELX_110" & !grepl("Pongo", Sample), AA := "I"]

# Collapse sites to have unique unmodified amino acids
sites <- sites[, .(MedScore = median(MedScore, na.rm = T),
                   MaxScore = max(MaxScore, na.rm = T),
                   AvgMassErrorPPM = mean(AvgMassErrorPPM, na.rm = T),
                   MaxLogInt = max(MaxLogInt, na.rm = T),
                   PSMcount = sum(PSMcount, na.rm = T),
                   Sequences = paste(unique(Sequence), collapse = ";")),
               by = .(Gene, Site, GeneSite, AA, Sample)]

# Absolute sequence coverage per sample
sites[, .(AbsCoverage = length(unique(GeneSite))), by = .(Sample)]
AACoverageSample <- sites[, .(AbsCoverage = length(unique(GeneSite))), by = .(Sample)]

fwrite(AACoverageSample, "AACoverageSample.csv")

## 6. Generate the consensus sequence

# Calculate relative variant coverage
sites[, SiteCoverage := sum(PSMcount), by = .(GeneSite, Sample)]
sites[, SiteInt := sum(MaxLogInt), by = .(GeneSite, Sample)]

fwrite(sites, "Site table.csv")

# Consensus based on most common variant
sites[, RelCoverage := PSMcount/SiteCoverage, by = .(GeneSite, Sample)]
sites[, RelInt := MaxLogInt/SiteInt, by = .(GeneSite, Sample)]
sites[, MaxCoverage := max(RelCoverage), by = .(GeneSite, Sample)]
sites[, Multiplicity := length(unique(AA[RelCoverage > 0.1 & RelInt > 0.1 & PSMcount > 1])), by = .(GeneSite, Sample)]
VarSites <- sites[Multiplicity > 1, .(Sample, GeneSite, RelCoverage, RelInt, 
                          PSMcount, SiteCoverage, AA)][order(Sample, GeneSite, RelCoverage, RelInt, 
                                                                      PSMcount, SiteCoverage, AA)]
# Export table
fwrite(VarSites, "VarSites.csv")

# Generate consensus table with the most commun variant
consensus <- sites[RelCoverage == MaxCoverage]
consensus[, .N, by = .(GeneSite, Sample)][N > 1]
consensus[, MaxSiteInt := max(MaxLogInt), by = .(Sample, GeneSite)]
consensus <- consensus[MaxLogInt == MaxSiteInt]
consensus[, .N, by = .(GeneSite, Sample)][N > 1]
#All covered position per sample
AllPos <- copy(FASTA)
AllPos <- AllPos[Gene %in% sites$Gene, .(Sequence = Sequence[1]), by = .(Gene)]
AllPos <- AllPos[, .(Site = 1:nchar(Sequence)), by = .(Gene)]
AllPos <- AllPos[, .(Sample = unique(sites$Sample)), by = .(Gene, Site)]

# Combine the identified sites/AAs with all possible alignment positions
consensus <- merge(consensus, AllPos, by = c("Sample", "Gene", "Site"), all = T)
# Fill up gaps with dashes
consensus[is.na(AA), AA := "-"]
consensus[, .N, by = .(Site, Gene, Sample)][N > 1]
# Sort for sequence assembly
setorder(consensus, Sample, Gene, Site)
# Paste together the sites into a consensus sequence
consensus <- consensus[, .(Consensus = paste(AA, collapse = "")), by = .(Sample, Gene)]
consensus[, AbsCoverage := str_count(Consensus, "[^-]")]
consensus[, RelCoverage := AbsCoverage/nchar(Consensus)]
consensus[, Annot := paste(Sample, Gene, AbsCoverage, sep = "_")]
consensus <- merge(consensus, unique(ExpAnnot[, .(Sample)]), by = "Sample")


## 7. Generate output

# Add reference sequences - 000 --> add reference sequences with sample name "000"
consensus <- rbind(consensus[, .(Consensus, Annot, Gene, Sample = as.character(Sample))], 
                   FASTA[Gene %in% consensus$Gene, .(Consensus = Sequence, Annot, Gene, Sample = "000")])
sort(c("000", as.character(unique(sites$Sample))))
setorder(consensus, Gene, Sample)

# Export as FASTA File
#Option 1: Consensus sequences for all the samples aligned with the reference sequences
write.fasta(as.list(consensus[]$Consensus),
            consensus[]$Annot,
            paste0("PR", format(Sys.Date(), "%y%m%d"), "_", "Hominid-Consensus_Ref_MQ.fasta"), as.string = T)
#Option 2: Consensus sequences of given samples aligned with the reference sequences
write.fasta(as.list(consensus[grepl("(Paranthropus)|000", Sample)]$Consensus),
            consensus[grepl("(Paranthropus)|000", Sample)]$Annot,
            paste0("PR", format(Sys.Date(), "%y%m%d"), "_", "Paranthropus-Consensus_Ref_MQ.fasta"), as.string = T)
#Option 3: Consensus sequences of given samples only
write.fasta(as.list(consensus[grepl("Paranthropus", Sample)]$Consensus),
            consensus[grepl("(Paranthropus)", Sample)]$Annot,
            paste0("PR", format(Sys.Date(), "%y%m%d"), "_", "Paranthropus-Consensus_MQ.fasta"), as.string = T)




