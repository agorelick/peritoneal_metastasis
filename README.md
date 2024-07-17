
# Instructions to regenerate analyses in the Main Figures and Extended Data Figures of _Wassenaar, Gorelick, et al. 2024_.

## 1. Install dependencies.

**Note:** Most of the dependencies for this project can be installed from the Conda environment included in this repository, `environment.yml`. All additional dependencies can be manually installed into the same conda environment following the steps below. This has the benefit that all newly downloaded software will not conflict with previously installed software. Additionally, all software installed during these instructions can be removed without permanent changes to your computer.

### 1.1 Clone this GitHub repository to your computer

On a computer with `git` installed, open a terminal and clone this repository to your current directory. 
```
git clone git@github.com:agorelick/peritoneal_metastasis.git
```

### 1.1 Install a version of the Conda package manager

If you already have Conda installed, proceed to the next step (1.2). Otherwise, install a version of Conda. I prefer miniforge, which is a minimal version of Mamba, a C++ re-implementation of Conda for improved speed. Instructions for download here: https://github.com/conda-forge/miniforge#mambaforge.

**Note**: This has so far only been tested on macOS (v11.5.2), but it should work with little/no modification on any UNIX-based operating system.

### 1.2 Install the `peritoneal_metastasis` Conda environment from the .yml file.

Change directory to the cloned peritoneal_metastasis repo, then install the conda environment. 

```
cd peritoneal_metastasis
conda env create -f environment.yml
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



### XXXX Process lpWGS copy number data

With the peritoneal_metastasis environment activated, run the install some additional dependencies into R. 
```
Rscript R/process_scna_C146.R
Rscript R/process_scna_C146.R
Rscript R/process_scna_C146.R
Rscript R/process_scna_C146.R
Rscript R/process_scna_C146.R
Rscript R/process_scna_C146.R
```

