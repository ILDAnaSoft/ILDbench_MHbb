## scripts directory


Put your Marlin steering files here and explain to which processor they refer (if the file name is not obvious) !

.) the steering file is zhll2j.xml

.) steps to follow to run the processors:

a.) prepare a link to the directory of MC samples, something like (on kekcc) data_ildopt_500 -> /group/ilc/soft/samples/mc-opt-3/ild/dst-merged/500-TDR_ws

b.) prepare a link LCFIPlusConfig -> $ILDConfig/v02-00-01/LCFIPlusConfig/

c.) prepare a link to gear file gear_ILD_l5_o1_v02.xml

d.) prepare links to weight files for isolated lepton tagging processors (each for electron and muon)

e.) run the master script "goSteve all"; jobs will be submitted; and then wait for completion

d.) that's it!

.) the outcome from above steps is a collection of root files under "rootfiles/" for all the signal and background samples;
in each root file, there are a few NTuples which contain various useful variables at generator level and reconstruction level

.) after all jobs finished successfully, run the goMerge script, all the root files from a same process will be merged.

.) go to "macros/Readme.md" to see steps for running analysis and making plots
