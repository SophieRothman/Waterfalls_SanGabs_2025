This repository contains codes used to produce the results for "Waterfalls perturb channel 
and hillslope dynamics in the San Gabriel Mountains, California" by Rothman, Scheingross, 
and McCoy. It includes:

width_calc_tribs_interp_nochute.m - which is a function that produces cross sectional profiles 
as well as channel slope, drainage area, and distance from a waterfall across a drainage

wf_finder.m - which is a function called by width_calc_tribs_interp_nochute.m, and identifies 
waterfalls in a drainage

width_modeling.m - which uses Ferguson et al., (2007)'s equation 20 to model flow depth, width, 
and velocity using cross sectional form, channel slope, dischrage and grain size.