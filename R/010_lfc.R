#!/usr/bin/env Rscript

# pipeline
# go from aligned bam files to log fold change counts


#' Compute log fold change from bam files
#'
#' Output of this function is also stored in data/hisat2.er.rda
#'
#' @param project
#' @param aligners
#' @param sampleTable
#' @param base_input_dir
#' @param lfc_dir
#' @return nothing. writes output to file
#' @export

step_010_lfc <- function(project, aligners, sampleTable, base_input_dir, lfc_dir) {
    for (aligner in aligners) {
        # set up the input/output locations
        align_dir <- file.path(base_input_dir, project, "aligned", aligner)

        # loading bam.files
        filenames <- as.character(sampleTable$Bam.File)
        filenames <- paste(align_dir, filenames, sep="/")
        print(filenames)
        if (file.exists(filenames)) {
            bamfiles <- Rsamtools::BamFileList(filenames)
            seqinfo(bamfiles[1])
                  
            # counting reads
            fc <- Rsubread::featureCounts(filenames,annot.inbuilt="hg38",isPairedEnd=TRUE)

            return fc
        }
        else {
            print("using cached file. To run these steps, Download BAM files from ... and put in ...")
            data(fc, package="NEMpipeline")
        }
        saveRDS(fc, file.path(lfc_dir, paste(aligner, project, "Rds", sep=".")))
    }
}
