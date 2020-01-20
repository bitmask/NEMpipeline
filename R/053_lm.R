#' Elastic net linear model
#'
#' @param project
#' @return nothing
#' @export

library(glmnet)

step_053_lm <- function( project, perturbed_genes, replicates, prepared_dir, nems_dir) {
    #from prepared data, find egenes == sgenes
    matching <- dir(prepared_dir, pattern=project)
    matching <- Filter(function(x) grepl(paste("\\.", "Rds", sep=""), x), matching)
    for (input_file_name in matching) {
        expr_data <- readRDS(file.path(prepared_dir, input_file_name))


        X <- t(expr_data)
        y <- rep(perturbed_genes, replicates)

        fit <- glmnet(X, y, family="multinomial", alpha=1)
        #plot(fit, xvar="dev", label=T)
        #plot(fit, xvar="lambda", label=T)

        p <- predict(fit, X, s=0.001, type="response")
        pmat <- matrix(p, nrow=replicates*length(perturbed_genes), ncol=length(perturbed_genes))
        colnames(pmat) <- perturbed_genes
        rownames(pmat) <- rownames(X)
        heatmap(pmat, scale="none", Rowv=NA, Colv=NA)

        cvfit <- cv.glmnet(X, y, family="multinomial", alpha=0.5, type.measure="class")


        # correction from perturb seq paper to account for cells that have a guideRNA but that have not been perturbed


        # put object in nem format
        l <- list()
        l$graph <- pmat
        output_file_name <- paste("lm", input_file_name, sep=".")
        saveRDS(l, file.path(nems_dir, output_file_name))
    }

}

