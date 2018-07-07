#Created by Alaa Abdel Latif on 25/08/2016. 
#Copyright © 2016 Alaa Abdel Latif. All rights reserved.

This is a currently working version of SELEXIM. This program was initially developed as part of a research project at University of Cambridge under supervision of Dr. Lucy Colwell. 

SELEXIM is a numerical simulation program for Systematic Evolution of Ligands through eXponential enrichment (SELEX) experiments. SELEX experiments are performed for aptamer discovery. Aptamers are small biomolecules 
(oligonucleotides) that can be used as candidates for therapeutics, diagnostics, and biomarkers. 

This program allows the user to simulate the dynamic changes in sequence populations during SELEX experiments through probabilistic modelling. The simulation can be carried out using either options that take into account primary or secondary structures or a combination of both through energy-based affinity measures. 

Installation
============================
Download or clone this repository.

Extract using:
$tar -czvf selex_simulation.tar.gz

PACKAGE DEPENDENCies
============================
Ensure you have the following Python (>2.7) module dependencies:

- sys
- time
- math
- linecache
- itertools
- operator
- collections
- scipy
- numpy
- matplotlib
- gc
- ConfigParser

If any of these modules are missing they can be installed using:
$sudo pip install [modulename]
If you do not have root priviledge and cannot use sudo, then you can install the missing module(s) locally using:
$pip install --user --upgrade [modulename]
Any further issues with installation should be informed to respective system administrators (or if it's a toughy just email us).

The program also requires the ViennaRNA package to be installed. This can be downloaded from:
https://www.tbi.univie.ac.at/RNA/#download
Extract the downloaded ViennaRNA package using:
$tar -czvf ViennaRNA-[version].tar.gz
This can then be installed using:
$cd ViennaRNA-[version]/
$./configure --without-perl --with-python
$make
$sudo make install
Again, if user does not have root privileges, then a local installation can be done using:
$./configure --prefix=[LOCALPATH] --without-perl --with-python
$make
$make install

USAGE
============================
All of the parameters for the simulation can be specified from the 'settings.init' file. Descriptions for each parameter is given inside the file. The default values in the settings file correspond to the conditions used to report the results in the corresponding thesis. 
After specifying the parameters, save the settings file and then run the simulation from the command-line using:

$python sim_.py

The software outputs updates to the command-line through out the simulation. It is recommended to save the updates 
to a log file. This will allow the user to fully understand the steps that were carried out during the simulation. 

Please note that under the default parameters, the simulation run takes almost 4 hours on an Intel(R)Core(TM) Quad CPU Q9400 machine. Using a large scale parameter or a large number of pcr cycles can result in excessive CPU time and memory use. 

Please report any issues to al.a.latif.94@gmail.com

Applications
============================
The SELEXIM software tool can be used to study the evolutionary dynamics of SELEX experiments. These can include 
the effects of round numbers, PCR cycles, nucleotide-bias, polymerase error rates, selection pressure, and sequence length. An important aspect about SELEXIM is that it has a modular structure. This means that it can be used as general framework for sequence evolution during SELEX. A user can change the models that are used to simulate each 
step of the experiment. Here is an example figure of sequence populations during a simulation:

![HE4_LOOP_DIST_FREQS](plots/he4_loop_small_SELEX_Analytics_distFreqs.pdf)

Another application that SELEXIM can be used for is aptamer discovery. Real SELEX experiments can require significant time and cost. At the end of an experiment, scientists usually identify a small number of aptamer candidates that have good affinity and specificity to the target ligand. SELEXIM can be used to generate a large number of alternative aptamer candidates that possess similar primary and secondary structures to the aptamers that were discovered from the real SELEX experiment. Here are some examples of alternative aptamers generated using SELEXIM against a candidate that was discovered during a SELEX experiment against a HE4 target ligand:

[HE4_InSilico_Aptamer_Discovery](plots/InSilico_Apt_Discovery.pdf)
