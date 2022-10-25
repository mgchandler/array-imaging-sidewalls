close all
clear
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution";
cd(dirname)

yaml_name = "scat_vs_artefact_dist.yml";
yaml_options = yaml.loadFile(yaml_name);