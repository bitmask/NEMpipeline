
# functions
trim_bootstrap <- function(graphs, bootstrap_threshold=bootstrap_threshold) {
    lapply( graphs, 
        function(gr) {
            values <- lapply(gr, function(x) { ifelse(x > bootstrap_threshold, 1, 0) })
            trimmed <- matrix( unlist(values), nrow=nrow(gr), ncol=ncol(gr), byrow=FALSE)
            dimnames(trimmed) <- dimnames(gr)
            return(trimmed)
        }
    )
}


# TODO tryCatch reading all files



step_070_plot <- function(distance_method, prep_method, project, nems_dir, plots_dir, draw_nets_max_nodes, draw_nets_max_count) {
    draw_networks <- TRUE
    # compare all the data (at first)
    compare <- list()
    names_list <- list()
    idx <- 1
    matching <- dir(nems_dir, pattern=prep_method)
    # filter for project name
    matching <- Filter(function(x) grepl(paste("\\.", project, "\\.", sep=""), x), matching)
    for (input_file_name in matching) {
        in_file <- file.path(nems_dir, input_file_name)
        if (file_test("-f", in_file)) {
            this_nem_method <- strsplit(input_file_name, "[.]")[[1]][1] # first element
            names_list[[idx]] <- input_file_name
            message("idx ", idx, " model ", input_file_name)
            tryCatch({
                nem_model <- readRDS(in_file)
                compare[[idx]] <- nem_model
                idx <- idx+1
            }, warning = function(w) {
              message(w)
            }, error = function(e) {
              message(e)
            })
        }
    }

    graphs <- lapply(compare, function(x) {as(x$graph, "matrix")})
    names(graphs) <- names_list

    # write out a text version of the adjacency matrix
    #for (idx in 1:length(graphs)) {
    #    name <- names(graphs)[[idx]]
    #    name <- gsub("Rds", "data", name)
    #    write.table(graphs[[idx]], file.path(plots_dir, name), quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    #}

    # compare bootstrap values in steps of 10%
    for (bootstrap_threshold in rev(seq(0.0,0.9,0.1))) {
        # bootstrap values less than the bootstrap_threshold -> 0, above -> 1
        trimmed_graphs <- trim_bootstrap(graphs, bootstrap_threshold=bootstrap_threshold)
        # should we display the networks?
        if (nrow(graphs[[1]]) > draw_nets_max_nodes || length(graphs) > draw_nets_max_count) {
            draw_networks <- FALSE
        }
        # TODO name this better
        if (requireNamespace("labnetmet")) {
            if (distance_method == "transitive") {
                output_pdf <- file.path(plots_dir, paste(paste("bootstrapgraphs_", prep_method, "_", bootstrap_threshold*100, sep=""), project, "pdf", sep="."))
                labnetmet::plot_dist(trimmed_graphs, labnetmet::trans_dist, output_pdf, draw_networks=draw_networks)
                
            } else if (distance_method == "intersection") {
                output_pdf <- file.path(plots_dir, paste(paste("intersect_bootstrapgraphs_", prep_method, "_", bootstrap_threshold*100, sep=""), project, "pdf", sep="."))
                labnetmet::plot_dist_list(trimmed_graphs, labnetmet::intersect_dist_list, output_pdf)
                # generate null models from the data
                params <- labnetmet::find_parameters_for_nullmodel(trimmed_graphs)
                n <- params[[1]]
                k <- params[[2]]
                p <- params[[3]]
                p <- signif(p, digits=2)

                output_pdf <- file.path(plots_dir, paste(paste("nullmodel_n-", n, "k-", k, "p-", p, "_", prep_method, "_", bootstrap_threshold*100, sep=""), project, "pdf", sep="."))
                labnetmet::calculate_random_dist_list(labnetmet::intersect_dist_list, output_pdf, n, k, p)
            } else {
                print(paste0("unknown distance method ", distance_method))
            }
        } else {
            print("To print the graphs")
            print("devtools::install_github('bitmask/labnetmet')")
        }
    }
}
