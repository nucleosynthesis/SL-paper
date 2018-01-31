# SL-paper
Repo for work on simplified likelihood paper

To produce the "inputs" - i.e how ATLAS/CMS would make the analysis run 

  `root -l -b -q makeSimpleLikelihoodToy.c`

Then convert this into a RooWorkspace (for fitting etc)

  `python letsUseRooFit.py`

From here, we can generate the toys (sampling from the nuisance parameters of the RooWorkspace) and then,
using those toys we can calculate the m1,m3 moments and the covariance matrix (set  makePlots=True to also make some 
comparison plots of the SL form vs toys)
 
  `python genToys`  # makes a ROOT TTree with hat{b} for each bin, use this output for the next step
  `python toys2ModelFile.py toy_trees.root model-90.py` 

The `model-90.py` can be used with the simplified likelihood code here: https://github.com/nucleosynthesis/work-tools/tree/master/stats-tools/SL 
