## Array Imaging Sidewalls - Matlab

A Matlab library based on the Python 3 `arim` language and the NDT Library. Implements more complex geometry for TFM and sensitivity images for a range of scatterers.

There are three groups of functions provided, found in the following subdirectories:
* `engine` :
	* All functions which are used to run the simulations and imaging steps.
* `home` :
	* The `main` functions (i.e. run to produce TFM or Sensitivity images) for use on a personal Matlab installation. Developed with Matlab version r2020b.
* `bp` :
	* The `main` functions for use on BluePebble. Designed for Matlab version r2019a.
	
There are two additional folders provided:
* `output` :
	* Default save location for files when `main` functions within the home folder are run.
* `support` :
	* Folder for miscellaneous other files which do not contribute to this program (i.e. Python testing, etc.).
	
To do:
- [ ] Fix multi frequency model.
- [ ] Adapt code to allow for multiple types of scatterers to be simulated at once.
- [x] Move away from specifying wall indices in ray weights function, so that general walls will be defined.
	- [ ] This may require probe and scatterer definitions as walls, so find a way to make this work also.
- [x] Align sensitivity with arim to improve runtime.
- [x] Update TFM back and frontwall views to make sure they still work with new geometry structures.
- [x] Go back and fix backwall views and frontwall view functions with respect to new method of calculating ray_weights.
- [x] Create a function fn_scat_info, comparable to fn_path_info.
- [x] Make sure that the function still works for the immersion case.
- [x] Vectorize functions as much as possible, to try and speed it up.
	- [x] Vectorize lookup of FMC data.
	- [x] Vectorize scatterer loop (i.e. for xpt = 1:total_xpts, for zpt ...)
- [x] When adding complicated geometry, make sure that conv_angle in fn_inc_out_angles is shifted with respect to the definition of the normal of the interface, rather than the z-component simply being +ve or -ve.
- [x] In fn_make_views, have it sort out views more elegantly like arim does using a function to sort, rather than just reference a big ol' list of valid names.