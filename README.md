
# Instructions to regenerate figures in _Wassenaar, Gorelick, et al. 2024_.

## Step 1. Install dependencies.

**Note:** Most of the dependencies for this project can be installed from the Conda environment included in this repository, `environment.yml`. All additional dependencies can be manually installed into the same conda environment following the steps below. This has the benefit that all newly downloaded software will not conflict with previously installed software. Additionally, all software installed during these instructions can be removed without permanent changes to your computer.

### 1.1 Clone this GitHub repository to your computer

On a computer with `git` installed, open a terminal and clone this repository to your current directory. 
```
git clone git@github.com:agorelick/peritoneal_metastasis.git
```

### 1.2 Install a version of the Conda package manager

If you already have Conda installed, proceed to the next step (1.2). Otherwise, install a version of Conda. I prefer **miniforge**, which is a minimal version of Mamba, a C++ re-implementation of Conda for improved speed. Instructions for download here: https://github.com/conda-forge/miniforge#mambaforge.

**Note**: This has so far only been tested on macOS (Big Sur and Sonoma), but it should work with little/no modification on any UNIX-based operating system.

### 1.3 Install the `peritoneal_metastasis` Conda environment from the .yml file.

Change directory to the cloned peritoneal_metastasis repo, then install the conda environment. 

```
cd peritoneal_metastasis
conda env create -n "peritoneal_metastasis" -f environment.yml
```
**Tips:** I ran into trouble creating this conda environment on a mac which already had R installed globally. Most of the conda packages installed, but then in R I was unable to install new packages such as the Quartet package (which is not on Conda) due to dyld files not being found. What worked for me was the following steps:
- completely remove the `peritoneal_metastasis` conda environment
- remove the files/directories for r-base-4.3.3: `rm -r ~/miniforge3/pkgs/r-base-4.3.3-*`
- then rerun `conda env create -n "peritoneal_metastasis" -f environment.yml`

Once the installation is complete, activate the peritoneal_metastasis environment.

```
conda activate peritoneal_metastasis
```

### 1.4 Install additional dependencies

With the peritoneal_metastasis environment activated, run the following to install some additional dependencies into R. 
```
Rscript R/get_prerequisites.R
```

## Step 2. (Optional) Process raw polyguanine genotype data using the _polyG_ R package

This GitHub repo contains pregenerated angular distance data (distance matrices and bootstrapped distance matrix R objects) which will be used to recreate figures from the manuscript. These were generated from raw polyguanine fingerprint data (i.e. "marker files") using the _polyG_ R package. It is not necessary to run this step if you only wish to recreate the manuscript's figures. However, if you are interested in regenerating the angular distance data, or you would just like to try out the _polyG_ package, that can be accomplished in this step.

### 2.1 Try out the polyG R package 

The _polyG_ R package will be automatically installed during **Section 1.4**. With the `peritoneal_metastasis` Conda environment attached, open R from the top directory of this GitHub repository. Then load the _polyG_ library and run the poly-G processing pipeline on the small example dataset (patients C12, C31 and C36 from Naxerova, _Science_ 2017) using the code chunk below. This only takes a few minutes to run on a 2020 MacBook Pro. **Please note:** The _polyG_ package is available as a stand-alone repository and can be installed independently of the instructions shown here. Please see https://github.com/agorelick/polyG for more detailed information. 

```r
# load the polyG package
library(polyG)

# run the pipeline on the example data
polyG(input_dir='original_data/polyG/example_input', results_dir='processed_data/polyG/example_output', seed=42)
```

The above command will populate output data into the directory `processed_data/polyG/example_output/results/sample_exclusion_0.3_rep_cut_0.11/`. This directory will contain numerous intermediate files for each patient. The key output files are the angular distance data, which will be generated in the following subdirectories:

- `results_angular_distance_representativeReplicates/angular_dist_trees_w_root_usedmarkers` - This directory contains PDF files of angular distance phylogenies for each patient (multiple versions will be generated for each patient).
- `results_angular_distance_representativeReplicates/angular_dist_matrix_w_root_usedmarkers` - This directory contains angular distance matrices as .txt files for each patient. It also includes .rds R objects containing 1,000 bootstrap replicates of the angular distance matrices.
- `results_angular_distance_representativeReplicates/angular_dist_heatmap_usedmarkers` - This directory contains PDFs of normalized polyguanine genotype heatmaps, as in **SI Figs 7-32, panel d**. 


### 2.2 Generate angular distance data for all patients

Raw poly-G marker files for all patients in this study are included in this GitHub repo under `original_data/polyG`, split up by patient cohort. Our pipeline to process these raw poly-G data and regenerate the angular distance data for each patient can be executed with the following command. This may take 30 minutes to an hour to complete. 

```
Rscript R/run_polyG_pipeline.R
```

Note: the above code will generate the same key output files as described in **Section 2.1** for each patient. These will be populated in the `processed_data/polyG` directory, organized by patient cohort (e.g. "peritoneal") and polyG processing workflow parameters. For example, the output data for patient C146 will be populated into `processed_data/polyG/peritoneal/results/sample_exclusion_0.1_rep_cut_0.11/results_angular_distance_representativeReplicates/`. (This will also generate the normalized polyguanine genotype heatmaps shown in **SI Figs 7-32, panel d**.)



## Step 3. Make all figures from the paper

### 3.1 Generate all figures except SCNA data (SI Figs 1-6) and angular distance heatmaps (SI Figs 7-32 panel d)
All data-based figures (with the exception of the SCNA data in SI Figs 1-6, and the angular distance heatmaps in SI Figs 7-32, panel d) can be regenerated by running the `R/make_figures.R` script with the following command. The generated figure panels will be populated in the directory `figures_and_tables/`.

**NB** In a few cases, it would be impractical to fully regenerate _all_ the input data (either because it would take a very long time to run, or because it would require installing software with dependencies that directly conflict with the dependencies for rest of the figures). In these cases, I have included pre-generated input data in this repo which will be automatically used. These cases are documented in the `R/make_figures.R` file with instructions to regenerate them if desired.  

```
Rscript R/make_figures.R
```

### 3.2 (Optional) Process lpWGS copy number data and generate SI Figs 1-6

The `R/make_figures.R` script uses preprocessed copy number data which are already included in this github repo. However, the complete set of copy number data data for each of the six patients with lpWGS can be processed with the following command (note, this may 5-10 minutes to complete). This will also generate the components of SI Figs 1-6.
```
Rscript R/process_scna_data.R
```
The components of SI Figs 1-6 will be generated into the following locations (similarly for the other 5 lpWGS patients):
* Panel A: `figures_and_tables/copynumber/heatmaps/C146_cnv_segment_heatmap.pdf`
* Panel B: `figures_and_tables/copynumber/distance_matrix_comparisons/C146_cnv_bins_euclidean_matrix_comparison.pdf`
* Panel C: `figures_and_tables/copynumber/tree_comparisons/C146_cnv_bins_euclidean_nj_tree_comparison.pdf`



## Step 4. Remove all newly-added software

All newly software installed can be removed by removiing the peritoneal_metastasis conda environment. First, detach the conda environment if you have not already:
```
conda deactivate
```

Next, remove the environment:
```
conda remove --name peritoneal_metastasis --all
```

If you also want to remove the conda application, see the instructions for the version of Conda you installed. If you used _miniforge_ with the link I suggested above, see the **Uninstallation** directions on the same website: https://github.com/conda-forge/miniforge#mambaforge.

## Contact

If you have any trouble running the code here, or would like more information, please feel free to contact me at alexander_gorelick@hms.harvard.edu. Thanks for looking!


