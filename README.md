# Project name: Directory Template

## Description

A short sentence describing the project. For example, this repository serves as a template for the directory structure of future projects.

### Goal

A paragraph discussing the goal and results of the project. For example, in order to increase the computational reproducibility of my projects, I have created this template repository for future projects. I believe that creating a directory structure that is adaptable and scales up well is critical to creating reproducible code, and hope this will be helpful. I find this structure useful because it has space for both simple R scripts (load, analyze, and visualize data) and it has a simple mechanism to source those R scripts for sharing with collaborators in a bookdown format. Specifically, the R scripts live within the `src` directory and can be separated into directories for discrete analyses. Each script will be numbered in the order in which they are intended to be sourced. For sharing the entire analyses workflow, the `docs` directory contains a simple file (`index.Rmd`) that can be `knit`ed to a bookdown format by sequentially sourcing each R script within separate code chunks. Analyses that take a non-trivial amount of time do not need to be evaluated during the `knit` process. Simply include the chunk option `eval = FALSE`. The code from the script will be printed into the document, but not executed. This way, only scripts that contain visuals need to be executed. Finally once the document is `knit`, GitHub Pages will build the a bookdown website from the files within the `doc` directory which is available at `username.github.io/repo_name` (_e.g._ [https://tsoleary.github.io/template_Rproject/](https://tsoleary.github.io/template_Rproject/])). All changes made to project files will be logged with sufficient detail in sentence form so they will be easily read and understood within the commit log to track the progress overtime.

### Project structure

This project is organized into the following directory structure:

- `data/` - all data used in code
  - `raw/` - all raw data
  - `processed/` - all processed data
- `src/` - source code used for analyses -- files saved with numerical prefix
  - `starwars/` - source code related to starwars data
  - `iris/` - source code related to iris data
- `output/` - all output files 
  - `figs/` - all figures generated
- `docs/` - a directory to host a `bookdown::gitbook` sourcing of all analyses
- `scratch/` - any temporary or extraneous files and information
  - `writing` - early writing scraps (_e.g._ methods sections) to be compiled later
  
    
