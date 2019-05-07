
# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

run_ER_pipeline <- function() {

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

# successful screens
samples <- c("EP300 12BR3 A", "EP300 12BR3 B", "EP300 12BR3 C", "ESRRA 2C A", "ESRRA 2C B", "ESRRA 2C C", "NCOA3 1DR3 A", "NCOA3 1DR3 B", "NCOA3 1DR3 C", "NCOA3 2BR A", "NCOA3 2BR B", "NCOA3 2BR C", "NR2F2 17AR3 A", "NR2F2 17AR3 B", "NR2F2 17AR3 C", "NRIP1 5C A", "NRIP1 5C B", "NRIP1 6A C", "RARA 7B A", "RARA 8C A", "RARA 8C C", "SUMO1 M1 A", "SUMO1 M1 B",     "SUMO1 M1 C",    "SUMO1 M35 A",    "SUMO1 M35 B",    "SUMO1 M35 C", "SUMO3 24AR3 A", "SUMO3 24AR3 B", "SUMO3 24AR3 C", "SUMO3 M16 A", "SUMO3 M16 C", "TRIM33 13CR3 A", "TRIM33 13CR3 B", "TRIM33 13CR3 C",     "TRIM33 M2 B", "TRIM33 M2 C", "ZMIZ1 19AR3 A", "ZMIZ1 19AR3 B", "ZMIZ1 19AR3 C", "ZMIZ1 22m35b A", "ZMIZ1 22m35b B", "ZMIZ1 22m35b C")



# 2 experiments must be above this lfc
expr.cutoff <- 0.5

################################################################################
#
# Pipeline
#
################################################################################

print("Running ER pipeline")
print(paste("project: ", project))
print(paste("aligner: ", aligner))
    step_010_lfc(project, aligner, ER_experiment_definitions, base_input_dir, lfc_dir)
    step_030_diffexp(project, aligner, diffexp_method, lfc_dir, diffexp_dir, ER_experiment_definitions, expr.cutoff, samples)
}

