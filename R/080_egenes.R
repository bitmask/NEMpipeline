
step_080_egenes <- function(project, aligner, diffexp_method, prep_method, nem_method, nems_dir, egenes_dir, sgenes) {
    # read in the egene attachments from the nems
    matching <- dir(nems_dir, pattern=nem_method)
    matching <- Filter(function(x) grepl(paste("\\.", prep_method, sep=""), x), matching)
    matching <- Filter(function(x) grepl(paste("\\.", diffexp_method, "\\.", sep=""), x), matching)
    matching <- Filter(function(x) grepl(paste("\\.", aligner, "\\.", sep=""), x), matching)
    matching <- Filter(function(x) grepl(paste("\\.", project, "\\.", sep=""), x), matching)
    for (input_file_name in matching) {
        tryCatch({
            in_file <- file.path(nems_dir, input_file_name)
            if (file_test("-f", in_file)) {
                print(input_file_name)
                nem_model <- readRDS(in_file)
                tocsv <- data.frame()
                for (sgene in sgenes) {
                    egene_list <- unique(unlist(nem_model$mappos[[sgene]]))
                    tocsv <- rbind(tocsv, cbind(sgene, egene_list))
                }
                tocsv <- as.data.frame(tocsv)
                colnames(tocsv) <- c("sgene", "egene")
                write.table(tocsv, file.path(egenes_dir, paste0("egenes.", input_file_name)), quote=FALSE, row.names=FALSE, sep="\t")
            }
        })
    }
}

