This repository contains the code for simulating the intercalation of LiMnO and Na in hard carbon.

The following provides a desctiption for each file in this repo:

**results**
_Description:_ This folder contains all of the input arguments and csv files for the data. All simulations will save in
here with a random id. I suggest changing the name after a simulation has finished to keep a logical order of results 
and some way to distinguish them.

**hard_carbon_animation.py**
_Description:_ This was an attempt to visualise the sodium insertion in hard carbon.

**hard_carbon_MC_clean.py**
_Description:_ This code simulates the sodium insertion in hard carbon and outputs a csv file for the results. This will
then have to be run on 'overlay_results.py' to display the results visually. This version is suitable for the HEC.

**LiMnO_MMC_Simulation.py**
_Description:_ This simulates the insertion of Li in LiMnO but does not use numpy arrays. 

**LMO_MC_efficient.py**
_Description:_ This is the most clear and easiest to understand version of the LiMnO system.

**MC_clean.py**
_Description:_ Michaels code for the LiMnO system including the effects of defects. 

**Na_hard_carbon_MC.py**
_Description:_ Na insertion in hard carbon with the plotting script built in which will display the results after the 
code has finished. 

**overlay_results.py**
_Description:_ This is the plotting script for all of the Na results.
