

#' Correlation between s-genes
#'
#' @param project
#' @return nothing
#' @export
step_054_proportionality <- function(project, perturbed_genes, prepared_dir, nems_dir) {

    #from prepared data, find egenes == sgenes

    matching <- dir(prepared_dir, pattern=project)
    matching <- Filter(function(x) grepl(paste("\\.", "Rds", sep=""), x), matching)
    for (input_file_name in matching) {
        expr_data <- readRDS(file.path(prepared_dir, input_file_name))
        x <- expr_data[perturbed_genes, perturbed_genes]

        # proportionality measure phi_s
        x[x==1] <- 100
        phs <- propr::propr(x+1, metric = "phi", ivar = 0)
        propr::plot(phs)
        propr::dendrogram(phs)
        net <- propr::getAdjacency(phs)

        # put object in nem format
        l <- list()
        l$graph <- net
        output_file_name <- paste("prop", input_file_name, sep=".")
        saveRDS(l, file.path(nems_dir, output_file_name))
    }

}

