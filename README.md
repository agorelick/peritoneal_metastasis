
# Instructions to regenerate figures in _Wassenaar, Gorelick, et al. 2024_.

## Step 1. Install dependencies.

**Note:** Most of the dependencies for this project can be installed from the Conda environment included in this repository, `environment.yml`. All additional dependencies can be manually installed into the same conda environment following the steps below. This has the benefit that all newly downloaded software will not conflict with previously installed software. Additionally, all software installed during these instructions can be removed without permanent changes to your computer.

### 1.1 Clone this GitHub repository to your computer

On a computer with `git` installed, open a terminal and clone this repository to your current directory. 
```
git clone git@github.com:agorelick/peritoneal_metastasis.git
```

### 1.1 Install a version of the Conda package manager

If you already have Conda installed, proceed to the next step (1.2). Otherwise, install a version of Conda. I prefer **miniforge**, which is a minimal version of Mamba, a C++ re-implementation of Conda for improved speed. Instructions for download here: https://github.com/conda-forge/miniforge#mambaforge.

**Note**: This has so far only been tested on macOS (Big Sur and Sonoma), but it should work with little/no modification on any UNIX-based operating system.

### 1.2 Install the `peritoneal_metastasis` Conda environment from the .yml file.

Change directory to the cloned peritoneal_metastasis repo, then install the conda environment. 

```
cd peritoneal_metastasis
conda env create -n "peritoneal_metastasis" -f environment.yml
```
**Tips:** I ran into trouble creating this conda environment on a mac which already had R installed globally. Most of the conda packages installed, but then in R I was unable to install new packages such as the Quartet package (which is not on Conda) due to dyld files not being found. What worked for me was the following steps:
- temporarily give ~/.R/Makevars and ~/.Rprofile new names
- completely remove the `peritoneal_metastasis` conda environment
- remove the files/directories for r-base-4.3.3: `rm -r ~/miniforge3/pkgs/r-base-4.3.3-*`
- then rerun `conda env create -n "peritoneal_metastasis" -f environment.yml`

Once the installation is complete, activate the peritoneal_metastasis environment.

```
conda activate peritoneal_metastasis
```

### 1.3 Install additional dependencies

With the peritoneal_metastasis environment activated, run the install some additional dependencies into R. 
```
Rscript R/get_prerequisites.R
```

## Step 2. Make all figures from the paper

### 2.1 Generate all figures except SCNA data (SI Figs 1-6)
All data-based figures (with the exception of the SCNA data in SI Figs 1-6) can be regenerated by running the `R/make_figures.R` script with the following command. The generated figure panels will be populated in the directory `figures_and_tables/`.

**NB** In a few cases, it would be impractical to fully regenerate _all_ the input data (either because it would take a very long time to run, or because it would require installing software with dependencies that directly conflict with the dependencies for rest of the figures). In these cases, I have included pre-generated input data in this repo which will be automatically used. These cases are documented in the `R/make_figures.R` file with instructions for to regenerate them if desired.  

```
Rscript R/make_figures.R
```

### 2.2 (Optional) Process lpWGS copy number data and generate SI Figs 1-6

The `R/make_figures.R` script uses preprocessed copy number data which are already included in this github repo. However, the complete set of copy number data data for each of the six patients with lpWGS can be processed with the following command (note, this may 5-10 minutes to complete). This will also generate the components of SI Figs 1-6.
```
Rscript R/process_scna_data.R
```

### 2.2 (Optional) Process raw polyG data to generate angular distance matrices and heatmaps (SI Figs 7-32, panel d)

The `R/make_figures.R` script uses pregenerated angular distance data (distance matrices and bootstrapped distance matrix R objects) which are already included in this github repo. These were generated from raw polyguanine fingerprint data (i.e. "marker files") using the _polyG_ R package (https://github.com/agorelick/polyG). Raw poly-G marker files are included in this `peritoneal_metastasis` GitHub repo under `original_data/polyG`, split up by patient cohort. Our complete pipeline to process the raw poly-G data and re-generate the angular distance matrices can be executed with the following commands. Note, the above step to install pre-requisites will already have installed the _polyG_ R package.

```
Rscript R/process_scna_data.R
```

## Step 3. Remove all newly-added software

All newly software installed can be removed by removiing the peritoneal_metastasis conda environment. First, detach the conda environment if you have not already:
```
conda deactivate peritoneal_metastasis
```

Next, remove the environment:
```
conda remove --name peritoneal_metastasis --all
```

If you also want to remove the conda application, see the instructions for the version of Conda you installed. If you used _miniforge_ with the link I suggested above, see the **Uninstallation** directions on the same website: https://github.com/conda-forge/miniforge#mambaforge.

## Contact

If you have any trouble running the code here, or would like more information, please feel free to contact me at alexander_gorelick@hms.harvard.edu. Thanks for looking!


