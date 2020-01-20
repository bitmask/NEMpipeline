#!/usr/bin/env Rscript

#library(GenomicRanges)
#library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)
library(ChIPpeakAnno)
library(org.Hs.eg.db)


library(dplyr)

get_all_chipseq_data <- function(project, chip_dir) {
    annoData <- ChIPpeakAnno::toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene)

    project_chip_dir <- file.path(chip_dir, project)

    # get all projects 
    sgene_peaks <- NULL
    for (sgene in selected.genes) {
        if (!dir.exists(file.path(project_chip_dir, sgene))) {
            print(paste0("skipping ", sgene))
            next
        }
        print(paste0("starting ", sgene))
        mapping <- data.frame()
        projects <- dir(file.path(project_chip_dir, sgene))
        reps <- list()
        i <- 1
        # push all data in all projects onto reps, and we'll take the union of the peaks that are found in all
        for (project in projects) {
            samples <- dir(file.path(project_chip_dir, sgene, project))
            for (sam in samples) {
                if (sam == "url" | sam == "orig") {
                    next
                }
                samfile <- file.path(project_chip_dir, sgene, project, sam)
                if (!file_test("-f", samfile)) {
                    # skip directories at this level, they may contain orig, unprocessed data
                    next
                }
                print(sam)
                chip <- NULL
                format <- tail(strsplit(sam, "[.]")[[1]], n=1)
                tryCatch({
                    #this_file <- system.file(project_chip_dir, sgene, project, sam, package="ChIPpeakAnno")
                    #chip <- toGRanges(this_file, format=format)
                    chip <- rtracklayer::import(samfile, format=format)
                }, warning = function(w) {
                  message(w)
                }, error = function(e) {
                  message(e)
                })
                # TODO: handle other weird formats

                reps[[i]] <- chip
                i <- i+1


            }
        }
        
        peaks <- Reduce(union, reps)
        
        # look through all the data that was read for this sgene, and turn this into a list of affected genes for this sgene
        anno <- NULL

        #anno <- annotatePeakInBatch(peaks, AnnotationData=annoData)
        #anno <- addGeneIDs(anno, orgAnn="org.Hs.eg.db", feature_id_type="entrez_id", IDs2Add=c("symbol"))

        anno <- annotatePeakInBatch(peaks, AnnotationData=annoData, output="overlapping", FeatureLocForDistance="TSS", bindingRegion=c(-10000, 10000)) 
        anno$symbol <- xget(anno$feature, org.Hs.egSYMBOL)


        mapping <- rbind(cbind(sgene, unique(anno$symbol)), mapping)
        colnames(mapping) <- c("sgene", "peak")
        sgene_peaks <- rbind(mapping, sgene_peaks)
    }
    out_file_name <- file.path(peaks_dir, "peaks")
    write.table(sgene_peaks, out_file_name, quote=FALSE, row.names=FALSE, sep="\t")
    return(sgene_peaks)
}

step_090_chip <- function(project, aligner, diffexp_method, prep_method, egenes_dir, peaks_dir, chip_dir, sgenes) {
    # peaks_dir is standard output dir
    # chip_dir is where the chipseq data is stored (large)

    # get chip data
    sgene_peaks <- NULL
    if (file.exists(file.path(peaks_dir, "peaks"))) {
        sgene_peaks <- read.csv(file.path(peaks_dir, "peaks"), sep="\t")
    }
    else {
        sgene_peaks <- get_all_chipseq_data(project, chip_dir)
    }
    lengthannoData <- 26034

    # get egenes for each file matching
    matching <- dir(egenes_dir, pattern="egenes")
    matching <- Filter(function(x) grepl(paste("\\.", prep_method, sep=""), x), matching)
    matching <- Filter(function(x) grepl(paste("\\.", diffexp_method, "\\.", sep=""), x), matching)
    matching <- Filter(function(x) grepl(paste("\\.", aligner, "\\.", sep=""), x), matching)
    matching <- Filter(function(x) grepl(paste("\\.", project, "\\.", sep=""), x), matching)
    for (input_file_name in matching) {
        attachments <- read.table(file.path(egenes_dir, input_file_name), header=TRUE)

        #shuffle attachments randomly
        attachments$egene <- sample(attachments$egene, length(attachments$egene), replace=FALSE)
        # should do this so that egenes per sgene are unique

        attachments$peak <- 0
        for (sgene in sgenes) {
            these_egenes <- attachments[attachments$sgene == sgene,]$egene
            these_peaks <- sgene_peaks[sgene_peaks$sgene == sgene,]$peak
            these_logic <- attachments$egene %in% these_peaks
            if (sum(these_logic) > 1) {
                attachments[(attachments$sgene == sgene) & these_logic,]$peak <- 1
            }
        }

        se <- attachments %>% tbl_df %>% group_by(sgene, peak) %>% summarize(count=n())
        se$peak <- factor(se$peak)

        sep <- attachments %>% tbl_df %>% group_by(sgene) %>% summarize(validated=sum(peak[peak==1]), total=n()) %>% mutate(percent = validated*100/total)

        # how many genes are peaks for each chip experiment?
        in_the_urn <- sgene_peaks %>% tbl_df %>% group_by(sgene) %>% summarize(count=n()) %>% unique()
        sep$nlog10p <- -1
        chip_enrichment <- list()
        for (sg in sgenes) {
            q <- sep[sep$sgene == sg,]$validated # number of white balls drawn
            m <- in_the_urn[in_the_urn$sgene == sg,]$count # number of white balls in the urn
            k <- sep[sep$sgene == sg,]$total - q    # number of black balls drawn
            n <- lengthannoData - k   # number of black balls in the urn
            if (q > 0) {
                print(paste0("sg ", sg))
                print(paste0("q ", q))
                print(paste0("m ", m))
                nlog10p <- -log10(phyper(q, m, n, k, lower.tail=FALSE))
                if (is.infinite(nlog10p)) {
                    nlog10p <- 100 #XXX if input data changes, check this is larger than max
                }
                chip_enrichment[[sg]] <- nlog10p
                sep[sep$sgene == sg,]$nlog10p <- nlog10p
            }
        }

        ggplot2::ggplot(sep %>% filter(nlog10p > -1), ggplot2::aes(x=reorder(sgene, percent), y=percent, size=nlog10p)) + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylim(c(0,100)) + ggplot2::scale_size(labels=function(x) {x[length(x)] <- expression("">=100); return(x)}) + ggplot2::xlab("Co-factor") + ggplot2::ylab("Percent") + ggplot2::labs(size="-log(P value)") 
        out_file_name <- file.path(peaks_dir, paste("validation", input_file_name, "pdf", sep="."))
        ggplot2::ggsave(out_file_name)
    }
}


