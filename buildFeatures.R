library('protr')
library('Peptides')

#' A set of helper functions for the buildFeature() function
#' Helper functions include: divideSeq(), buildSAAC, and seqStats
#' @param string peptide sequence
#' @param prop1 First split (left-side) of the peptide sequence
#' @param prop2 Second split (middle) of the peptide sequence

## simple function to break the peptide sequence into 3 parts
divideSeq <- function(string, prop1, prop2) {
  length <- nchar(string)
  
  if (length %% 3 != 0) {
    part1_length <- floor(length * prop1)
    part2_length <- floor(length * prop2)
    part1 <- substr(string, 1, part1_length)
    part2 <- substr(string, part1_length + 1, part1_length + part2_length)
    part3 <- substr(string, part1_length + part2_length + 1, length)
    return(c(part1[1], part2[1], part3[1]))
  } else {
    part_length <- length / 3
    part1 <- substr(string, 1, part_length)
    part2 <- substr(string, part_length + 1, 2 * part_length)
    part3 <- substr(string, 2 * part_length + 1, length)
    return(c(part1[1], part2[1], part3[1]))
  }
}

## build SAAC feature dataframe
buildSAAC <- function(seq, prop1, prop2) {
  SAAC <- data.frame()
  
  for (i in seq) {
    ## break into 3 parts based on the proportion specificed
    parts <- divideSeq(i, prop1, prop2)
    
    ## go through each part and calculate AAC for each part
    composition <- lapply(parts, function(part) {
      comp <- extractAAC(part)
      temp <- as.data.frame(t(comp))
      temp <- temp[, order(names(temp))]
      return(temp)
    })
    ## combine everything
    SAAC <- rbind(SAAC, do.call(cbind, composition))
  }
  
  rownames(SAAC) <- make.names(seq, unique=TRUE)
  return(SAAC)
}

## get peptide sequence summary info
seqStats <- function(seq,AMP=T){
  ## get the simple statistics of sequences
  count <- length(seq)
  if(AMP==T){
    cat('#########AMP-Sequence Summary############\n')
  }else{
    cat('#########NonAMP-Sequence Summary############\n')
  }
  cat(count,"peptide sequences are found after preprocessing!\n")
  cat('Minimum sequence length is',min(nchar(seq)),'\n')
  cat('Maximum sequence length is',max(nchar(seq)),'\n')
  cat('Average sequence length is', mean(nchar(seq)),'\n')
  cat('##########################################\n')
}

preProcess <- function(seq){
  seq <- seq[!grepl("^>", seq)]
  seq <- seq[!grepl("[BJOUXZ]", seq)] ## remove non-standard rare residues
  seq <- Filter(function(x) nchar(x) >= 10, seq) ## filter out the ones less than 10 AAs
  seq <- Filter(function(x) nchar(x) <= 500, seq) ## filter out the ones greater than 500 AAs
  return(seq)
}

#' Build the features for peptide sequences
#' Features include: AAC + AAC.group + AAC.class + PAAC + NAAC + SAAC
#'
#' This function computes the necessary features for later model prediction
#'
#' @param ampfile Path to a fasta format AMP file, or a text file in fasta format
#' @param nonampfile Path to a fasta format Non-AMP file, or a text file in fasta format
#' @param out Name of output file
#' @param split1_prop First split (left-side) fraction of the peptide sequence 
#' @param split2_prop Second split (middle) fraction of the peptide sequence
#' @return The feature dataframe of peptide sequences
#'
#' @examples
#' features <- buildFeatures('testAMP_data.fasta','testNonAMP_data.fasta',
#'                            'final_features',0.2,0.6)
#' 
#'
#' @author Junxi Feng <junxifeng@hsph.harvard.edu>
#' @import protr
#' @import Peptides
#' @export
buildFeatures <- function(ampFile,nonampFile,out,split1_prop=0.2,split2_prop=0.6){
  ## INPUT: fasta files
  amp_seq <- readLines(ampFile)
  nonamp_seq <- readLines(nonampFile)
  
  ## Preprocess
  amp_seq <-preProcess(amp_seq)
  nonamp_seq <- preProcess(nonamp_seq)
  
  cat('Preprocessing done!\n')
  
  seq <- unlist(c(amp_seq,nonamp_seq))
  ## Get stats
  seqStats(amp_seq,AMP=T)
  seqStats(nonamp_seq,AMP=F)
  
  ## Calculate Amino Acid Composition (AAC)
  AAC <- data.frame(do.call(rbind, lapply(seq, function(i) t(data.frame(extractAAC(i))))))
  rownames(AAC) <- make.names(seq, unique=TRUE)
  AAC <- AAC[, order(names(AAC))]
  
  ## AAC Group
  AAC$GP <- AAC$G + AAC$P
  AAC$TSYQN <- AAC$`T` + AAC$S + AAC$Y + AAC$Q + AAC$N
  AAC$DE <- AAC$D + AAC$E 
  AAC$HKR <- AAC$H + AAC$K + AAC$R
  cat('AAC done!\n')
  
  ## PHYC (full name)
  AAC$Iso_electric <- sapply(seq, function(i) pI(i))
  AAC$Hydrophobicity <- sapply(seq, function(i) hydrophobicity(i))
  AAC$Net_charge <- sapply(seq, function(i) charge(i))
  cat('PHYC done!\n')
  
  ## Peptide
  #Peptide <- data.frame(do.call(rbind, lapply(seq, function(i) t(data.frame(aaComp(i)[[1]]))[1,])))
  #rownames(Peptide) <- make.names(seq, unique=TRUE)
  #colnames(Peptide) <- rownames(aaComp(seq[2])[[1]])
  #cat('Peptide done!\n')
  
  ## PAAC
  PAAC_list <- lapply(seq, function(i) {
    PAAC1 <- extractPAAC(i, lambda = 9)
    PAAC2 <- extractAPAAC(i, lambda = 9)
    cbind(t(data.frame(PAAC1)), t(data.frame(PAAC2)))
  })
  
  PAAC <- do.call(rbind, PAAC_list)
  rownames(PAAC) <- seq
  cat('PAAC done!\n')
  
  ## NAAC
  NAAC <- AAC[1:20] / sqrt(rowSums(AAC[1:20]^2))
  cat('NAAC done!\n')
  
  ## SAAC
  SAAC <- buildSAAC(seq,split1_prop,split2_prop)
  cat('SAAC done!\n')
  
  ## Assemble
  res <- cbind(AAC,PAAC,NAAC,SAAC)
                           
  ## Labels
  labels <- c(rep(1, length(amp_seq)), rep(0, length(nonamp_seq)))
  res <- as.data.frame(scale(as.data.frame(res)))
  res$labels <- labels
  
  cat('Finished!\n')

  output_path <- 'feature_matrix'
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  file_path <- file.path(output_path, out)

  write.csv(res,file_path)
  cat('Feature matrix stored in',output_path,'\n')
}



args <- commandArgs(trailingOnly = TRUE)

param1 <- as.character(args[1])
param2 <- as.character(args[2])
param3 <- as.character(args[3])
param4 <- as.numeric(args[4])
param5 <- as.numeric(args[5])

buildFeatures(ampFile = param1,nonampFile = param2,out = param3,
              split1_prop = param4,split2_prop = param5)
