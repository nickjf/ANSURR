
# ANSURR | Accuracy of NMR Structures Using RCI and Rigidity 
[![DOI](https://zenodo.org/badge/234519929.svg)](https://zenodo.org/badge/latestdoi/234519929)

ANSURR uses backbone chemical shifts to validate the accuracy of NMR protein structures as described here https://www.nature.com/articles/s41467-020-20177-1. This repository contains the code required to install and run ANSURR on a Linux (and possibly Mac) machine. Please let me know if you have any issues. 

## Installation

To install run the shell script `install.sh` and follow the instructions (see below for an overview of these instructions):

`./install.sh`

## Running ANSURR

ANSURR requires two input files, a NMR protein structure in PDB format and a shifts file in NMR Star v3 format. To re-reference chemical shifts using PANAV before running ANSURR (recommended):

`ansurr -p xxxx.pdb -s xxxx.str -r`

To run without re-referencing chemical shifts:

`ansurr -p xxxx.pdb -s xxxx.str`

Options:

`-b` (under development) read in extra covalent bonds from CONECT recoords in the PDB file. FIRST will identify most covalent bonds but not metal-ligand bonds. This option allows you to add these in by specifying them in CONECT records (for details on how to add these to your PDB file see here: https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html). Note that modelling metal-ligand bonds as covalent bonds in this way may not be reliable, hence this option is under development. Currently, I would not recommend using ANSURR to validate proteins containing bound metals. 

`-h` print the help message 

`-l` include free ligands when computing flexibility. Free ligands are defined as HETATMs that appear after the TER record in the PDB file. Note that metals will generally not be bonded unless specified in CONECT records (see option -b). 

`-n` include non-standard residues when computing flexibility. Non-standard residues are defined as HETATMs that appear before the TER record in the pdb file. Note that RCI will not be calculated for non-standard residues and so they will not be used to compute validation scores. Regardless, including non-standard residues is a good idea to avoid breaks in the protein structure which would otherwise make those regions too floppy.

`-o` combine chains into a single structure when calculating flexibility. This is useful when the structure is an oligomer as oligomerisation will often result in changes in flexibility. Currently this option combines all chains present in the pdb file. A future release will allow you to choose which of the chains should be combined. 

`-r` re-reference chemical shifts using PANAV before running ANSURR


## Output

A directory called `<yourpdbfile>_<yourshiftfile>` is made to save the output generated. This directory will be overwritten if you run ANSURR again with input files with the same names as before. This directory contains two directories called  `ANSURR_output` and `other_output`. `ANSURR_output` contains:  

* `scores.out` - a text file with the validation scores for each model 
* `<yourpdbfile>_<yourshiftfile>.png` - a graphical summary of the validation scores for each model 
* `out/` - text files with protein flexibility predicted by RCI and FIRST for each model
* `figs/` - plots of protein flexibility predicted by RCI and FIRST for each model

`other_output` contains output from various programs run as part of ANSURR:

* `PANAV/` - re-referenced chemical shifts
* `RCI/` - flexibility predicted from chemical shifts using RCI
* `extracted_pdbs/` - PDB files for each model extracted from the NMR structure
* `DSSP/` - secondary structure for each model according to the program DSSP
* `FIRST/` - flexibility predicted for each model using FIRST


## Overview of installation instructions

1. You are asked whether you want to install ANSURR for all users:<br><br> 
`would you like to install ANSURR for all users? [recommended - requires root privileges] (y/n) `<br><br>
Answering yes is recommended because then all users on the machine will be able to run ANSURR (obvious), but also you will be able to run ANSURR anywhere by simply typing `ansurr`, without having to place it in your `$PATH`, ammending your `$PATH` or copying the program into the directory you want to run it each time. In order to do this you need root/superuser privileges and so you will likely need to provide your password.

2. The program will check your default version of python and ask if you want ANSURR to use that one:<br><br>
`default version of python on this machine is: python x.x`<br>
&nbsp;` -> would you like ANSURR to use this version? [python3 is recommended] (y/n)`<br><br> 
ANSURR will work with either python2 or python3 but as python2 is no longer supported, it is recommended that you use python3. If you don't want to use the default version, the program will look for other python versions installed on your machine and list them:<br><br>
`other versions of python installed on this machine include: python2 python3 pythonx.x`<br>
&nbsp;` -> please type which to use from the list above (or type the path if not in the list):`<br><br>
To select a python version just type it when prompted. If you don't see the version of python you want but you know where it lives you can type the path to it. 

3. The program will check if the python version you selected has the required modules: `numpy scipy matplotlib`. If any are missing, it will suggest how you might go about installing them using the program `pip` and offer to run this command for you:<br><br>
`checking pythonx.x has required modules (numpy, scipy, matplotlib)`<br>
`could not find scipy matplotlib, you could try using pip to install by running: pythonx.x -m pip install scipy matplotlib`<br>
&nbsp;` -> would you like to try this now? (y/n)`<br><br>
If module installation fails then you will need to find another way to install these modules in order to run ANSURR (see https://docs.python.org/3/installing/index.html).

4. ANSURR can optionally use PANAV to re-reference chemical shifts, which requires Java. The program checks to see if Java is installed and if not will offer to try and install it using `APT` on Linux machines or `Homebrew` on Macs:<br><br>
`ANSURR can optionally use PANAV to re-reference chemical shifts, which requires Java, checking to see if Java is installed`<br>
`Java not found, you could try to install Java by running: sudo apt install default-jre`<br>
&nbsp;` -> would you like to try this now? (y/n)`<br><br>
Java is only required to run PANAV, so ANSURR will run fine without it - but make sure your chemical shifts are referenced correctly or re-reference them using some other software.


## Help

Contact Nick Fowler (njfowler.com) for support, queries or suggestions.

## Known Issues

- PANAV (a program for re-referencing chemical shifts) has some kind of GUI interface built in to the code and will not run if it can't find a display. This means it won't work when running ANSURR on a remote machine without X forwarding.  


## Acknowledgements

Random Coil Index (RCI) | Berjanskii, M.V. &amp; Wishart, D.S. A simple method to predict protein flexibility using secondary chemical shifts. Journal of the American Chemical Society 127, 14970-14971 (2005).

Floppy Inclusions and Rigid Substructure Topography (FIRST) | Jacobs, D.J., Rader, A.J., Kuhn, L.A. &amp; Thorpe, M.F. Protein flexibility predictions using graph theory. Proteins-Structure Function and Genetics 44, 150-165 (2001).

Probabilistic Approach to NMR Assignment and Validation (PANAV) | Bowei Wang, Yunjun Wang and David S. Wishart. "A probabilistic approach for validating protein NMR chemical shift assignments". Journal of Biomolecular NMR. Volume 47, Number 2 / June 2010: 85-99

DSSP | A series of PDB related databases for everyday needs. Wouter G Touw, Coos Baakman, Jon Black, Tim AH te Beek, E Krieger, Robbie P Joosten, Gert Vriend. Nucleic Acids Research 2015 January; 43(Database issue): D364-D368. | Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Kabsch W, Sander C, Biopolymers. 1983 22 2577-2637.

adjustText - automatic label placement for matplotlib | https://github.com/Phlya/adjustText











