The cfg and macros in this dir can be used to evaluate trigger efficiency from MC samples (MC truth efficiency).
The HLTrigReport.py cfg can be run to get the trigger report using CMSSW.

The rest are pyRoot macros that require CMSSW environment (cmsenv) to run.
The genVertexRadius.py produces control plots for the generator level quantities, including decay radius (with respect to the center of CMS).
The L2Muons.py computes the efficiency for the trigger used in the displaced leptons analysis and plots of L2 muons.
The L3Muons.py computes the efficiency for a standard dimuon (prompt) trigger and plots of L3 muons.

IMPORTANT: the triggers are selected by their index in the L2Muons and L3Muons macros. These numbers can change if the HLT menu changes.
To get the correct numbers run the HLTrigReport.py cfg to get the list of triggers and their indeces. Replace the numbers with the correct
ones if needed in the L2Muons.py and L3Muons.py.
