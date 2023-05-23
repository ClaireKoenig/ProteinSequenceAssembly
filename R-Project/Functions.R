##################################################
# Project:
# Author(s): Patrick Ruether
# Date: 03.03.2020
# Version:
# Script purpose: Function collection
##
##
##################################################

# Packages

require(data.table)
require(seqinr)
require(bit64)

# Load FASTA(s)
# Transfor FASTA file into data table, well annotated

fasta_to_datatable <- function(path = ".", paths = as.character()) {
  if(length(paths > 0)){
    FASTApaths <- paths
  } else {
    if(length(list.files(path, pattern = ".fasta")) == 0) stop("No FASTA file in specified path")
    FASTApaths <- list.files(path, pattern = ".fasta", full.names = T)
  }
  for(i in FASTApaths){ #read all the lines in the FASTA?
    iFASTA <- read.fasta(i, seqtype = "AA", as.string = T, strip.desc = T) #strip.desc:	if TRUE the '>' at the beginning of the description lines is removed in the annotations of the sequences
    iFASTAmatrix <- matrix(unlist(lapply(iFASTA, attributes)), byrow = T,
                           ncol = length(names(attributes(iFASTA[[1]])))) #generate matrix with prot sequences
    colnames(iFASTAmatrix) <- names(attributes(iFASTA[[1]]))
    iFASTAmatrix <- data.table(iFASTAmatrix, Sequence = unlist(iFASTA, use.names = F)) #convert matrix to data table
    iFASTA <- data.table(iFASTAmatrix, str_split_fixed(iFASTAmatrix$name, "\\|", 3))
    names(iFASTA)[grep("V\\d", names(iFASTA))] <- c("Database", "Protein", "Gene_Species")
    iFASTA$File <- i
    if(exists("output")){
      output <- rbind(output, iFASTA)
    } else {
      output <- iFASTA
    }
  }
  output[, Gene := str_extract(Annot, "(?<=GN\\=)[^\\s]+")]
  output[, Gene := toupper(Gene)]
  output[, Species := str_extract(Annot, "(?<=OS\\=)[^\\s]+\\s[^\\s]+")]
  output[, SpeciesID := str_extract(Annot, "(?<=OX\\=)[^\\s]+")]
  return(output)
}

# MQ proteinGroups Intensities to Int64

MQintensity64 <- function(DT = x) {
  if(any(!sapply(DT[, grepl("(Intensity)|(iBAQ)|(LFQ)", 
                            names(DT)), with = F], is.integer64))) {
    keycol <- names(DT)
    
    DT <- DT[, lapply(.SD, as.integer64), 
             .SDcols = names(DT)[grepl("(Intensity)|(iBAQ)|(LFQ)", names(DT)) & 
                                   (!sapply(DT, is.integer64))],
             by = c(names(DT)[!(grepl("(Intensity)|(iBAQ)|(LFQ)", names(DT)) & 
                                  (!sapply(DT, is.integer64)))])]
    if(sum(grepl("%", names(DT)) & 
           (!sapply(DT, is.double))) > 0) {
      DT <- DT[, lapply(.SD, as.double), 
             .SDcols = names(DT)[grepl("%", names(DT)) & 
                                   (!sapply(DT, is.double))],
             by = c(names(DT)[!(grepl("%", names(DT)) & 
                                  (!sapply(DT, is.double)))])] }
    setcolorder(DT, keycol)
    remove(keycol)
    return(DT)
  } else {
    return(DT)
  }
}

# Most common element in a Vector

MostCommon <- function(x) {
  # Check if type character
  if(!is.character(x)) {
    stop("Input has to be type Character")
  }
  # Check if empty character vector
  if(length(x[x != ""]) == 0){
    x <- NA
  } else {
    x <- sort(table(x), decreasing = T)
    x <- x[nchar(names(x)) > 0]
    x <- x[x == max(x, na.rm = T)]
    # Check if more than 1 most common gene name
    if(length(x) > 1){
      # Take the shortest gene name
      x <- x[nchar(names(x)) == min(nchar(names(x)))]
      # If multiple most common and shortest gene names, take first one
      x <- names(x[1])
    } else {
      x <- names(head(x, 1))
    }
  }
  return(x)
}

# All entries in a pasted list are e.g. contaminants

all.contain <- function(x, pattern = "CON__", sep = ";") {
  require(stringr)
  if(length(grep(sep, x, value = T)) == 0){
    warning("Separator not found in input values")
  }
  x <- str_split(x, sep)
  x <- lapply(x, grepl, pattern = pattern)
  x <- lapply(x, all)
  x <- unlist(x)
  return(x)
}

# Paste all list columns in a data table

# Paste list-type columns

pasteListCols <- function(dt, sep = ";"){
  if(!is.data.table(dt)){
    stop("Input is not a data.table")
  } else {
    dt[, names(dt)[dt[, sapply(.SD, typeof)] == "list"] := lapply(.SD, function(x) {
      unlist(lapply(x, paste, collapse = ";"))}), 
      .SDcols = names(dt)[dt[, sapply(.SD, typeof)] == "list"]]
  }
}


# Melt Maxquant proteingroups table to long format

MeltMQprot <- function(dt = x, appendix = "_long") {
  if(!is.data.table(dt)){
    stop("Input has to be a data.table. Generate by reading proteinGroups.txt with fread function")
  }else{
    iBAQ <- any(grepl("^iBAQ", names(dt)))
    LFQ <- any(grepl("^LFQ", names(dt)))
    MSMScount <- any(grepl("^MS\\/MS count\\s.+", names(dt)))
    print(paste0("Melted Peptides, Razor + unique peptides, Unique peptides, Sequence coverage, Intensity",
          ifelse(MSMScount == T, ", MSMScount", ""), ifelse(LFQ == T, ", LFQ", ""), ifelse(iBAQ == T, ", iBAQ", "")))
    # List all column groups for melting
    
    meltlist <- list(grep("^Peptides .+", names(dt)), 
                     grep("^Razor \\+ unique peptides .+", names(dt)), 
                     grep("^Unique peptides .+", names(dt)), 
                     grep("^Sequence coverage [^\\[]", names(dt)), 
                     grep("^Intensity .+", names(dt)))
    if(MSMScount == T){meltlist <- c(meltlist, list(grep("^MS\\/MS count\\s.+", names(dt))))}
    if(LFQ == T){meltlist <- c(meltlist, list(grep("^LFQ .+", names(dt))))}
    if(iBAQ == T){meltlist <- c(meltlist, list(grep("^iBAQ .+", names(dt))))}  
    
    meltcols <- c("Peptides",
                  "Razor + unique peptides",
                  "Unique peptides",
                  "Sequence coverage",
                  "Intensity")
    if(MSMScount == T){meltcols <- c(meltcols, "MS/MS count")}
    if(LFQ == T){meltcols <- c(meltcols, "LFQ")}
    if(iBAQ == T){meltcols <- c(meltcols, "iBAQ")} 
    
    # Stop if meltlist contains column groups of different sizes
    
    if(length(unique(unlist(lapply(meltlist, length)))) != 1){
      print(lapply(meltlist, function(x) names(dt)[x]))
      stop("Something wrong with your column names")
    }
    
    # Melt that data table
    
    longDT <- melt(dt, measure.vars = meltlist, variable.name = "Experiment", 
                          value.name = paste0(meltcols, appendix),
                          value.factor = T)
    
    levels(longDT$Experiment) <- str_extract(names(dt), "(?<=Peptides ).+")[grep("Peptides ", names(dt))]

    # Order columns
    
    setcolorder(longDT, c("Experiment", "Protein IDs", "Majority protein IDs", "Protein names",
                          "Gene names", "Fasta headers", "Number of proteins", 
                          paste0(meltcols, appendix))[c("Experiment", "Protein IDs", "Majority protein IDs", "Protein names",
                          "Gene names", "Fasta headers", "Number of proteins", paste0(meltcols, appendix)) %in% names(longDT)])
    
  }
}

