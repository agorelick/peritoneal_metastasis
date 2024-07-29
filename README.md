
# Instructions to regenerate analyses in the Main Figures and Extended Data Figures of _Wassenaar, Gorelick, et al. 2024_.

## 1. Install dependencies.

**Note:** Most of the dependencies for this project can be installed from the Conda environment included in this repository, `environment.yml`. All additional dependencies can be manually installed into the same conda environment following the steps below. This has the benefit that all newly downloaded software will not conflict with previously installed software. Additionally, all software installed during these instructions can be removed without permanent changes to your computer.

### 1.1 Clone this GitHub repository to your computer

On a computer with `git` installed, open a terminal and clone this repository to your current directory. 
```
git clone git@github.com:agorelick/peritoneal_metastasis.git
```

### 1.1 Install a version of the Conda package manager

If you already have Conda installed, proceed to the next step (1.2). Otherwise, install a version of Conda. I prefer **miniforge**, which is a minimal version of Mamba, a C++ re-implementation of Conda for improved speed. Instructions for download here: https://github.com/conda-forge/miniforge#mambaforge.

**Note**: This has so far only been tested on macOS (v11.5.2), but it should work with little/no modification on any UNIX-based operating system.

### 1.2 Install the `peritoneal_metastasis` Conda environment from the .yml file.

Change directory to the cloned peritoneal_metastasis repo, then install the conda environment. 

```
cd peritoneal_metastasis
conda env create -n "peritoneal_metastasis" -f environment.yml
```

Once the installation is complete, activate the peritoneal_metastasis environment.

```
conda activate peritoneal_metastasis
```

### 1.3 Install additional dependencies

With the peritoneal_metastasis environment activated, run the install some additional dependencies into R. 
```
Rscript R/get_prerequisites.R
```

## 2. Make all figures from the paper

```
Rscript R/make_figures.R
```


### (Optional) Process lpWGS copy number data

The `make_figures.R` script uses preprocessed copy number data which are already included in this github repo. However, the complete set of copy number data data for each of the six patients with lpWGS can be processed with the following command (note, this may 5-10 minutes to complete).
```
Rscript R/process_scna_data.R
```

## Remove all newly-added software

All newly software installed can be removed by removiing the peritoneal_metastasis conda environment. First, detach the conda environment if you have not already:
```
conda deactivate peritoneal_metastasis
```

Next, remove the environment:
```
conda remove --name peritoneal_metastasis --all
```

If you also want to remove the conda application, see the instructions for the version of Conda you installed. If you used _miniforge_ with the link I suggested above, see the **Uninstallation** directions on the same website: https://github.com/conda-forge/miniforge#mambaforge.



