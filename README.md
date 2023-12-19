# Stage-dependent_cell_differentiation_in_multicellularity_with_specialized_cell_types
This repository includes the codes, data, and other related materials for the work "Dynamic cell differentiation in multicellularity with specialized cell types".

Overview

This repository contains the simulation of "Stage-dependent cell differentiation in multicellularity with specialized cell types" by Yuanxiao et al. The aim is to investigate the population growth rate of multicellular organisms under different dynamic cell differentiation strategies.

Organization

Code.

This folder includes the Python scripts for calculating the growth rate of organisms under different cell differentiation strategies.

Data.

This folder stores the data (population growth rate and according to differentiation probabilities across cell divisions) of each optimal dynamic cell differentiation strategy. Since the data size, here we show the compressed data rather than the original data.

Figure.

This folder contains the Python scripts that read the relevant data stored in the folder data and generate corresponding figures.

Usage

dynamic_heatmap.py contains the source code for calculating population growth rates of dynamic cell differentiation strategies.

make.py is the values of parameters. The file generates a run.sh file, which is executed on the cluster.

data_* stores the optimal dynamic cell differentiation strategy and its corresponding population growth rate under different parameter spaces.

fig*.py is the Python script that generates the figures in the main text by reading the corresponding data in the data folder.

Requirement

All the codes were developed in Python 3.6.
