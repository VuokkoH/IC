# Script for iterating the IC tool for multiple input points

Purpose: InterfaceCatchment analysis for STYLE-project walkability analysis.
Author: VH 4 July 2022 / updated 11 Jan 2023

## Details of the latest run:

- Final analysis was coducted using QGIS via CSC Puhti.
- Run settings: 300m walkable distance, no dead-end removal.
- Final analysis was conducted for 1578 starting points in the city centre, and 1014 starting points in the Kaarela suburb

## Contents
--------

IC_inputs 		- Input files used in the analysis (input points, blocks layer) in  EPSG:3067.
IC_results 		- Result files
IC_script_drafts	- Older script versions /testing
icsimple.py		- This file contains the IC analysis in function run() -modified from original so that references to GUI were removed by removing all references to 'self'. Function returns the IC value.
icsetup.py	 	- This file imports the icsimple tool run-function and handles the input filepaths and final result file writing. There are two versions; one that uses pandas for output file (used locally), and one that uses pure python (used in CSC)