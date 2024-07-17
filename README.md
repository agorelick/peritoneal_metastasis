
# Instructions to regenerate analyses in the Main Figures and Extended Data Figures of _Wassenaar, Gorelick, et al. 2024_.

## 1. Install dependencies.


### 1.1 Clone this GitHub repository to your computer

On a computer with `git` installed, open a terminal and clone this repository to your current directory, then `cd` into this directory. 
```
git clone git@github.com:agorelick/peritoneal_metastasis.git
cd peritoneal_metastasis
```

### 1.1 Install a version of the Conda package manager

Most of the dependencies for this project can be installed from the included Conda environment, `environment.yml`. If you already have Conda installed, proceed to the next step (1.2). Otherwise, install a version of Conda. I prefer miniforge, which is a minimal version of Mamba, a C++ re-implementation of Conda for improved speed. Instructions for download here: https://github.com/conda-forge/miniforge#mambaforge.

**Note**: This has so far only been tested on macOS (v11.5.2), but it should work with little/no modification on any UNIX-based operating system.

### 1.2 Install the `peritoneal_metastasis` Conda environment from the .yml file.

Install the conda environment. 
```
conda env create -f environment.yml
```
