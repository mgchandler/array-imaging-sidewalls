## Array Imaging Sidewalls - Matlab

A Matlab library based on the Python 3 `arim` language and the NDT Library. Implements more complex geometry for TFM and sensitivity images for a range of scatterers.

There are three groups of functions provided, found in the following subdirectories:
* engine
All functions which are used to run the simulations and imaging steps.
* home
The `main` functions (i.e. run to produce TFM or Sensitivity images) for use on a personal Matlab installation. Developed with Matlab version r2020b.
* bp
The `main` functions for use on BluePebble. Designed for Matlab version r2019a.
	
There are two additional folders provided:
* output
Default save location for files when `main` functions within the home folder are run.
* support
Folder for miscellaneous other files which do not contribute to this program (i.e. Python testing, etc.).
	
To do:
- [ ] Go back and fix backwall views and frontwall view functions with respect to new method of calculating ray_weights.
- [ ] Create a function fn_scat_info, comparable to fn_path_info.
- [ ] Make sure that the function still works for the immersion case.
- [ ] Adapt code to allow for multiple types of scatterers to be simulated.
- [ ] Adapt functions to allow for matrix multiplication as much as possible, to try and speed it up.