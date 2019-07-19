# sim
FOXSI simulation codes to share with Säm and others

What it does: Sets up X-ray sources and combines them into one map.  The main perks are that it is very easy to add as many sources 
              as you want; each new source is a line or two of code.  Sources can be Gaussian (specified by parameters) or can be 
              maps.  Everything is coregistered.
              
What it doesn't do (yet): Folding through the FOXSI response. This is coming soon (along with documentation).

foxsi_sim_image__define.pro defines the simulation object that contains all the machinery to set up the sources.

sample_script.pro shows an example of using this object to simulate some of the things we put in the FOXSI proposal.

