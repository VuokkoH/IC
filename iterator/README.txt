VH 4 July 2022

InterfaceCatchment analysis for STYLE-project walkability analysis.

The tool was run with default settings: 400m walkable distance, no dead-end removal.

Folders:
--------

IC 			- original IC tool forked from https://github.com/Awapic/IC
IC_inputs 		- Input files used in the analysis, modified from Maija's files in here:  T:\STYLE\Paper 4\Walkability index\Kortteleiden huokoisuus ja saavutettava julkisivu\IC 
	 		- Re-projected all input data to EPSG:3067 prior to analysis, and split input point files into multiple files in order to run in subsets.

IC_results 		- Result files
IC_script_drafts	- Older script drafts /testing
iciterator.py		- This file contains the IC analysis in function run() -modified from original so that references to GUI were removed by removing all references to 'self'. Function returns the IC value.
import_iterate.py 	- This file imports iciterator and handles the input filepaths and final result file writing. There are two versions; one that uses pandas for output file (used locally), and one that uses pure python (used in CSC)