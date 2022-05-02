# MC_veto
Repository for vetoed monte carlo 
This is a quick-and-dirty how-to-use functionality review of the MC integrator I've made.

Step 1) Move the files 'MC_MT_veto.cc', 'Makefile', and all .sh files into the src directory of dShower, which I assume you will have set up with openloops, have PDF.cc in the polarised/unpolarised basis, and have generated all the 3flv DPDF grids. The makefile is identical to the dShower makefile but includes instructions on how to make the MC_MT_veto program.

Step 2) Move the 'Integrals' folder also into src - this has the file structure already set up to compute both plots of dPDF distributions as a function of the rapidity product, and plots of the asymmetry values as functions of a minimum rapidity.

Step 3) open MC_MT_veto.cc and change the path from my '/home/lthsmith/Desktop/' section to wherever you have stored dShower. Similarly change the paths to dshower and OpenLoops in all the 'runner_[].sh' files corresponding to the DPDs you wish to use, and also in the 'linker.sh' file.

Step 4) run 'source linker.sh' and then 'make MC_MT_veto' to compile the integrator.

_____ASIDE_____
if you wish to use the integrator on its own, the syntax is as follows:

./MC_MT_veto [dpd selection (int, 0-3)] [iterations (int, >0)] [quark flavour 1] [quark flavour 2] [quark flavour 3] [quark flavour 4] [minimum rapidity (double, >0 & <sqrt(19))]

where the quark flavours ints labelled as follows: up=2,down=1,strange=3,ubar=-2,dbar=-1,sbar=-3
______________

NOTE: if you select dpd=2, this is the pospol prescription and you *MUST* open PDF.cc comment out lines 312 and 313, and uncomment lines 314 and 315 to access the correct DPD grid.

To calculate total DPD distributions, you MUST edit MC_MT_veto's file allocation lines - (206-221 at time of writing) by, taking PDGS as an example, uncommenting line 207 and commenting away line 208. All other DPDs can be changed in this manner. This is due to the code being set up for asymmetry calculations at this time. After this, you can run 'source runner_[DPDF].sh', where [DPDF] can be either PDGS,Mixpol,Pospol, or MSTW. This will run a bash file I've set up to loop over all the quark flavours, stored in combs_WW.txt, with minimum rapidity=0 in batches of fourteen (I wrote these on a 16-thread PC), and combine (with the relevant multiplicity in flavour as there are some symmetries) all the binned data in the 'Integrals/[DPDF]_Plot' folder. These can then be plotted by use of the 'DPD_plotter.py' file in the 'Integrals' folder.

In order to plot the asymmetries, you start by running 'source runner_veto_[DPDF].sh', which computes the same process as above, but also loops over minimum rapidity from 0-3 in steps of 0.5. These will be combined for each minimum rapidity in the same fashion as the prior case, but the generated combined .txt files that start with 'Combined' will themselves be evaluated by the python file 'Asym_plotter_[DPDF].py' which will generate a .txt and a pyplot of the asymmetry over the range of minimum rapidity.

If you run all of the [DPDF] options for the above, you can move each of the Asyms_[DPDF].txt files from their respective folders into the 'Asyms' folder, and run 'Asymplotter.py', which will plot the distribution of the asymmetries as a function of the min. rapidity, and as a ratio against the DGS DPDF.

It is worth noting that the .sh files can be easily edited to change the number of concurrent processes, the max and min rapidities, the number of iterations, etc.

If you have any questions, comments, suggestions or gripes, please feel free to get in touch with me: lawrie.smith@postgrad.manchester.ac.uk
