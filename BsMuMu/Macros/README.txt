The macros in this folder can be used to filter the tree and produce a version with preselection cuts applied and with all the cut-based analysis cuts applied.
Usage:
./all.py

At this moment the file must be edited to change wether to apply the analysis cuts or only the preselection.
The variable cut_based is a boolen, if true it means the analysis cuts are applied.

Then, run
./runTMVA.sh

To produce all the tables and plots, run
python TMVAPlots.py
