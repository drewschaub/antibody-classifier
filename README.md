# Antibody Classifier

## Overview
This repository contains tools for the RMSD alignment and analysis of antibody structures. The pipeline includes dataset curation, structure alignment, RMSD computation, and result visualization.

## Installation
To set up the environment, follow these steps:

1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/antibody-classifier.git
    cd antibody-classifier
    ```

2. Create a virtual environment and install dependencies:
    ```bash
	conda create -n antibody-classifier
	conda activate antibody-classifier
	conda install pandas pymol-open-source scipy -c conda-forge
    ```

## Usage

### Curating Datasets
Store your pdb files in `pdb.gz` or `pdb` format in the `data` folder. You can separate projects. To match with the manuscript, this has folders for `hiv1`, `influenza`, and `sarscov2`. To gzip on MacOS, simply copy your `pdb` files to your folder and run:
```bash
gzip *.pdb
```

### Defining the epitopes
Define each epitope in the epitope.csv file which will be stored in your `data_dir`. The first column will be the name of the `pdb.gz` or `pdb` file. The second column wll be the epitope defined as a dictionary object:

```
{'B': [512,513,514,515,516,517,518,519]}
```