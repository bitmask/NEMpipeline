
# execute the analysis pipeline for the simulated data
# each step will write intermediate output files to disk

run_simulated_pipeline <- function() {

    ################################################################################
    #
    # Definitions
    #
    ################################################################################

    # file locations
    base_input_dir <- "~/NEMpipelineBAM"

    # output locations
    base_output_dir <- "~/NEMpipelineoutput"
    source_dir <- file.path(base_output_dir, "000_source")
    lfc_dir <- file.path(base_output_dir, "010_lfc")
    diffexp_dir <- file.path(base_output_dir, "030_diffexp")
    #e2_regulon_dir <- file.path(base_output_dir, "039_e2_regulon")
    prepared_dir <- file.path(base_output_dir, "040_prepared")
    heatmap_dir <- file.path(base_output_dir, "045_heatmap")
    nems_dir <- file.path(base_output_dir, "050_nems")
    consensus_dir <- file.path(base_output_dir, "060_consensus")
    plots_dir <- file.path(base_output_dir, "070_plots")
    benchmark_dir <- file.path(base_output_dir, "075_benchmark")
    egenes_dir <- file.path(base_output_dir, "080_egenes")
    peaks_dir <- file.path(base_output_dir, "090_chip")

    # and ensure they exist
    for (output_dir in c(base_output_dir, lfc_dir, diffexp_dir, prepared_dir, heatmap_dir, nems_dir, consensus_dir, plots_dir, benchmark_dir, egenes_dir, peaks_dir)) {
        if( ! dir.exists(output_dir)) {
            dir.create(output_dir)
        }
    }

    benchmark_file <- "timing_simulated.log"

    # prep_methods: binary, lfc, pvalue, boolean, fgnem, decompose, progressive, bootstrap, random, geneset, cluster_bin
    # binary - discretize expression data to 0 or 1 (there is an effect or not)
    # lfc - use log fold change of differential expression
    # pvalue - use pvalue of differential expression
    # boolean - model inhibition and activation with bnem
    # fgnem - use factor graph nem
    # decompose - subset sgenes into a set of smaller nems, then evaluate with search to see if the models are compositional
    # progressive - take N most significant egenes, then 2N most significant, then 3N etc
    # bootstrap - randomly sample from the egenes
    # random - randomly generate data, to compare the bootstraps against
    # geneset - use gsea to group the egenes into genesets to run the nems on
    # cluster_bin - cluster binary expression data to generate ensemble egenes to run nems on
    prep_method <- "binary"

    # Which NEM methods should be applied?
    # The way the data is prepared determines which NEM methods are compatible with which prepared data
    #
    # nem_methods: search, triples, pairwise, BN.exhaustive, BN.greedy, ModuleNetwork.orig, ModuleNetwork, mc.eminem, depn, lem
        # search -- Markowetz 2005 "Non-transcriptional pathway features reconstructed from secondary effects of RNA interference"
        # triples -- Markowetz 2007 "Nested effects models for high-dimensional phenotyping screens"
        # pairwise -- Markowetz 2007 "Nested effects models for high-dimensional phenotyping screens"
        # BN.exhaustive -- Zeller 2008 "A Bayesian Network View on Nested Effects Models"
        # BN.greedy -- Zeller 2008 "A Bayesian Network View on Nested Effects Models"
        # ModuleNetwork.orig -- Frohlich 2007 "Large scale statistical inference of signaling pathways from RNAi and microarray data"
        # ModuleNetwork -- Frohlich 2008, "Estimating large-scale signaling networks through nested effect models with intervention effects from microarray data"
        # nem.greedyMAP -  Tresch 2008 "Structure Learning in Nested Effects Models"
        # nem.greedy -- H. Fröhlich, A. Tresch, T. Beißbarth, Nested Effects Models for Learning Signaling Networks from Perturbation Data, Biometrical Journal, 2(51):304 - 323, 2009
        # mc.eminem -- Niederberger 2012 "MC EMiNEM Maps the Interaction Landscape of the Mediator"
        # depn -- Frohlich 2009 "Deterministic Effects Propagation Networks for reconstructing protein signaling networks from multiple interventions"
        # lem -- Szczurek 2016 "Linear effects models of signaling pathways from combinatorial perturbation data"
    nem_method_compat <- list("binary" = c("greedy", "search", "nem.greedy", "triples", "pairwise", "ModuleNetwork", "ModuleNetwork.orig"), 
                         "cluster_bin" = c("greedy", "search", "nem.greedy", "triples", "pairwise", "ModuleNetwork", "ModuleNetwork.orig"), 
                         "lfc" = c("mnem", "bnem", "lem", "mnem", "mc.eminem"),
                         "pvalue" = c("mc.eminem", "bnem", "lem", "mnem", "mc.eminem"),
                         "boolean" = c("bnem"), 
                         "geneset" = c("mc.eminem"),
                         "bootstrap" = c("nem.greedy"),
                         "random" = c("nem.greedy"),
                         "decompose" = c("search"),
                         "progressive" = c("nem.greedy"),
                         "sc" = c("lem")
                         )

    nem_method <- "triples"

    report_attached_egenes <- FALSE

    ################################################################################
    #
    # Pipeline
    #
    ################################################################################

    print("Running simulated pipeline")
    print(paste("project: ", project))
    
    perturbed_genes <- letters[1:7]

    for (alpha in seq(0.05, 0.5, 0.1)) {
        for (beta in seq(0.05, 0.5, 0.1)) {
            step_042_simulate_data(project, alpha, beta, perturbed_genes, prepared_dir, lfc=TRUE)
        }
    }



    # nems
    step_050_nems(project, aligner, diffexp_method, prep_method, nem_method, nem_method_compat, prepared_dir, nems_dir, egenes_dir, benchmark_file, report_attached_egenes)

    # run HIS matlab code here
    r <- readline(prompt="Run HIS matlab code in matlab now and press enter ")

    # read in matlab output and save as an R object so we can generate network comparison plots
    step_051_his(project, "his", perturbed_genes, prepared_dir, nems_dir)


}

