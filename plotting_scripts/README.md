# CURRENT INSTRUCTIONS #

These scripts take the layer ntuples from the alignment code in this repo and plot either a single chamber or a ring.

`mualvisualMUAL.py` plots a single chamber.
Usage: `python mualvisualMUAL.py MEp_1_3_17_MC_TEST.root`
The script has convenient labeling features if the filename is given in a certain format.
` MEp_1_3_17_MC_TEST.root ` tells the script which chamber to plot, if it is MC or data, and a label string ("TEST")

For data, no "_MC" is required, so a filename should look like `MEp_1_3_1_EXAMPLE.root`.


`ringvisualmultiple.py` plots the ME13 ring.
Usage: `python mualvisualMUAL.py filename1.root filename2.root filename3.root ...`


# OLD INSTRUCTIONS #
Files:

indexbase.php
mualvisualMUAL.py
MEp_1_3_20_MCRING.root //example data

Example output:

http://muon-alignment-rymuelle.web.cern.ch/muon-alignment-rymuelle/nick_python_script/MEp_1_3_20_MCRING_plots/_index.php

Here are some basic instructions written by Nick:

 really only have one script file. All the functionality has been built in and would need to be modified for use with the alignment procedure. Here are some things that would need to be removed (I'll work on this):
 - "Specialsuffix" stuff for organization
 - Automatic index page creation
 - Cut implementation
 - Alignment/analyzer toggle hack
 - Autodetection of chamber from filename
 - Clean up comments/"old" code
I suppose when uploading to github, I'll keep all of these features for standalone use

I have placed 3 files (excluding the folder) in http://namin.web.cern.ch/namin/alignment/plots/final/ryan/ 
As an example, MEp_1_3_17_MCTEST.root is the output of running over monte carlo using single muon gun shooting into the ME13 ring.
indexbase.php is a php file that collects all picture files in a folder and places them in a grid in one html page (so we don't have to click on all the pictures individually to look at them)
mualvisualMUAL.py is the main python script

Instructions for running:
1) Ssh into lxplus.cern.ch (note, not lxplus5, since their python version was too low)
2) Copy these 3 files from above to some area
3) Run: `python mualvisualMUAL.py MEp_1_3_17_MCTEST.root -b`
4) A folder called "MEp_1_3_17_MCTEST_plots" will be produced. Inside you will find an _index.php file. 
5) Clicking on _index.php, for example, here http://namin.web.cern.ch/namin/alignment/plots/final/ryan/MEp_1_3_17_MCTEST_plots/ will neatly show all the produced plots
Tested this procedure just now and it works for me, so let me know if you run into any issues.
