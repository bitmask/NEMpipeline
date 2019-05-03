#!/usr/bin/env Rscript

# pipeline
# go from aligned bam files to log fold change counts

source("config.R")

#' Compute lfc from bam files
#'
#' @param project
#' @param aligners
#' @param sampleTable
#' @param base_input_dir
#' @param lfc_dir
#' @return Number representing the distance between x and y
#' @export

step_010_lfc <- function(project, aligners, sampleTable, base_input_dir, lfc_dir) {
    for (aligner in aligners) {
        # set up the input/output locations
        align_dir <- file.path(base_input_dir, project, "aligned", aligner)

        # loading bam.files
        filenames <- as.character(sampleTable$Bam.File)
        filenames <- paste(align_dir, filenames, sep="/")
        print(filenames)
        file.exists(filenames)
        bamfiles <- Rsamtools::BamFileList(filenames)
        seqinfo(bamfiles[1])
              
        # set to single core
        #register(SerialParam())
              
        # counting reads
        fc <- Rsubread::featureCounts(filenames,annot.inbuilt="hg38",isPairedEnd=TRUE)

        saveRDS(fc, file.path(lfc_dir, paste(aligner, project, "Rds", sep=".")))
    }
}
