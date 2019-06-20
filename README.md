==NEM Pipeline

To generate simulated data with false positive and false negative rates as specified, run methods, and plot results.  Results will be written into the directory specified.

```
> set.seed(42)
> alpha <- 0.15
> beta <- 0.05
> output_dir <- "~/projects/NEMpipelineoutput")
> run_simulated_pipeline(alpha, beta, output_dir)
```
