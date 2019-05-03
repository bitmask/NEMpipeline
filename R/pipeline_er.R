#!/usr/bin/env Rscript

# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

################################################################################
#
# Definitions
#
################################################################################

# file locations
base_input_dir <- "../input"

project <- "er"

projects_definition <- list("er" = list("regulon" = "Jagannathan&Robinson-Rechavi2011.csv",
                                        "experiment" = "experiment_definitions.csv"))

# read in the sample.table excel sheet specifying the experiment details
#csv.file <- file.path(base_input_dir, project, projects_definition[[project]]$experiment)
#if (file.exists(csv.file)) {
    #experiment_definitions <- read.csv(csv.file)
    #sampleTable <- experiment_definitions
#}


data(experiment_definitions, package="NEMpipeline")

# aligners: hisat2, bowtie
aligners <- c("hisat2") 

# diffexp_methods: DESeq, edgeR
diffexp_methods <- c("DESeq") 

################################################################################
#
# Pipeline
#
################################################################################


step_010_lfc(project, aligners, experiment_definitions, base_input_dir, lfc_dir)


