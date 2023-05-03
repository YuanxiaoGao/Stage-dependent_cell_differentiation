# Dynamic_cell_differentiation_in_multicellularity_with_specialized_cell_types
This repository includes the codes, data and other related materials for the work "Dynamic cell differentiation in multicellularity with specialized cell types".

Overview

This repository contains the simulation of "Dynamic cell differentiation in multicellularity with specialized cell types" by Yuanxiao et al. The aim is to investigate the population growth rate of multicellular organisms under different dynamic cell differentiation strategies.

Organization

Code.

This folder includes the python scripts for calculating the growth rate of organisms under different cell differentiation strategies.

Data.

This folder stores the data (population growth rate and according differentiation probabilities acrocess cell divisions) of each optimal dynamic cell differentiation strategies. Since the data size, here we show the compressed data rather than the original data.

Figure.

This folder contains the python scripts that read the relevant data strored in the folder data and generate corresponding figures.

Usage

dynamic_heatmap.py contains the source code for calculating population growth rates of dynamic cell differentiation strategies.

make_*.py are the values of parameters. The file generates a run.sh file, which is executed on cluster.

Requirement

All the codes were developed in Python 3.6.
