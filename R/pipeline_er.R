
# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

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

# project name and definitions
project <- "er"
projects_definition <- list("er" = list("regulon" = "Jagannathan&Robinson-Rechavi2011.csv",
                                        "experiment" = "experiment_definitions.csv"))
data(ER_experiment_definitions, package="NEMpipeline")

# read in the sample.table excel sheet specifying the experiment details
#csv.file <- file.path(base_input_dir, project, projects_definition[[project]]$experiment)
#if (file.exists(csv.file)) {
    #experiment_definitions <- read.csv(csv.file)
    #sampleTable <- experiment_definitions
#}



# aligners: hisat2, bowtie
aligner <- "hisat2" 

# diffexp_methods: DESeq, edgeR
diffexp_method <- "DESeq" 

################################################################################
#
# Pipeline
#
################################################################################


step_010_lfc(project, aligner, experiment_definitions, base_input_dir, lfc_dir)
step_030_diffexp(project, aligner, diffexp_method, lfc_dir, diffexp_dir, experiment_definitions)


