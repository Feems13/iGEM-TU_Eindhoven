In this project, a kinetic model of a generalized extracellular molecule sensor (GEMS) platform that activates the Jak-STAT3 pathway is developed and analyzed. The goal is to explore receptor activation and signal transduction in the GEMS system. This gives the opportunity to explore optimal design parameters which can be implemented experimentally. 


## Table of Contents
* [General Info](#general-information)
* [Technologies Used](#technologies-used)
* [Features](#features)
* [Usage](#usage)
* [Acknowledgements](#acknowledgements)


## General Information
iGEM-TU_Eindhoven contains a model that describes a synthetic biological system of GEMS system and the JAK-STAT3 pathway. The model consist of 22 ordinary differential equations (ODEs), that are solved with the ode15s solver.
The purpose of this project is to investigate the behavior of the system (the receptor activation, signal transduction and the output concentration) and to explore the effect of changing some constants.

## Technologies Used
- MATLAB - version 2018a

## Features
The following four scripts simulate the ODE model:
- ODEs_IL10_deg.m 			- Contains the ODEs and different cases.
- Defaults_deg.m			- The default settings and initial concentrations.
- OdeSolver_IL10_deg.m			- Solves the ODEs and plot the results.
- paramIL10_deg.m 			- Contains the constants.

The following five scripts run the MPSA model:
- MPSA_process_results_dummies.m	- Processes the results of the dummy values.
- MPSA_process_results_final.m		- Processes the results of the parameters.
- MPSA_run_results.m			- Simulate MPSA.
- Simulate_EP.m				- Simulate MPSA with specific features.
- find_characteristics.m		- MPSA is done for this feature. 
- lhsdesign_modified.m			- Design latin hypercube sample.
- defaults_MPSA.m			- Defaults for the MPSA.

The following scripts can be used to investigate the system:
- TEST_ODE_deg.m 			- To execute some simple tests to check whether the model statisfy predefined boundary conditions.
- Vary_Kd.m				- The effect of variatins in the KD value on the required times for STAT3npd to stabalize.
- ratio_RJL.m				- To plot the effect of variations of the ligand-receptor ratio on the signal transduction. 
- features.m				- Features determined for the ODE.
- features_RJ.m				- Features for receptor-ligand ratio.
- features_kd.m				- Features for difference in Kd.

## Usage
To simulate the model and to plot the results, run the file OdeSolver_IL10_deg.m.
To run the MPSA run MPSA_run_results.m

## Acknowledgements
- This project was inspired by the work of prof. Tom de Greef.
- Many thanks to Tom de Greef, Bryan, Alex, and Anna.
