# SmB6_XYmap_analysis
Code updated and uploaded on Jan 5 2022.
Used in analysis published in Kotta et al., "Metallic chemical potentials in an insulating topological state." \
For questions, or to request a set of demo ARPES data (.itx or .mat), email Erica Kotta at eck325@nyu.edu. 

## Instructions on running quick demo: 
Save all files in this repository into the same folder, then run XYmap_gui. \
Demo dataset "XYSmB6_demo_data.mat" should be listed in drop-down menu near top of gui. \
Click "Import mat" to import the mat file data, then click "Import params" to import the pre-saved set of parameters (FL E, BG E, etc). \
Click "Convert axes" to convert theta to momenta. New momenta axis will be plotted at the top of the main figure. \
Now parameters including "th" will be ignored for parameters including "k" for both plotting and calculating feature maps (except for tilt map, which searches for L/R d-band crossing in raw-theta to find surface tilt angle and uses that value to convert individual image k-axis). \
If desired, adjust parameters, and click "Save" to overwrite parameter set mat file. \
Click "Run" at the bottom of the left-side panel to calculate feature maps. \
Click "Map error" to estimate the error of feature maps. (Necessary step for making binning series.) \
\
Select maps from the drop-down menus in bottom panel to choose axes of binning series, and set corresponding number of bins for each axis. Check "k-shift" to apply k-axis correction, then click "Make BS" to plot a binning series. \
\
A new figure will pop up. Click "Plot It" to plot binned spectra, with H-axis distribution plotted at the top and V-axis distribution plotted to the right. \
Numbers above each binned spectrum indicate how many individual spectra contribute to that bin. H/V bin edges and widths may be adjusted.\
Click "Export BS" button at top-right to export the current series of spectra and corresponding energy/momentum axes to your current workspace for further analysis. 

## Instructions on how to run the code:

Save all files in this repository into the same folder, then run XYmap_gui.

#### If you don't already have an XY scan dataset saved into a mat file: 
You can take a set of ARPES images, saved in the .itx format, saved into a folder that you specify within XYmap_gui.m script as variable data.itx_foldername and corresponding .csv file specified with data.csv_filename. \
The .csv file contains XY scan metadata, listing the manipulator X,Y coordinates of each ARPES image in the folder of itx files. \
Functions itx_csv_to_spec.m and csv_to_XY_list.m convert the dataset into a set of variables. You can now save these into a .mat file using the commented-out code located in the convertItx2Mat function inside the gui script. \
Restart the gui with the new .mat file saved in the same folder. 

#### If you already have the XY scan dataset saved into a mat file: 
Select the mat file from the drop-down menu near the top of the gui and click the "Import mat" button. 
##### If you need to set new parameters: 
There is a set of default parameter values (FL E, BG E, ...) listed on the left hand side with default values for SmB6 spectral feature analysis. 
- All parameters with "E" are energies, in units of kinetic energy. 
- All paremeters with "th" are theta, units of degrees. 
- All paremeters with "k" are momenta, units of inverse angstrom. 

See XYmap_gui.m script for more descriptions of the parameters. \
Click "Save" at bottom left to save the current set of parameter values into a new .mat file (name set same as mat file of spectral data with added "\_params") \
Click "Run" to run the feature map analysis code as defined in run_the_params.m script and produce a set of feature maps. 

##### If you have a saved set of parameters from a previous run:
Click "Import params" to import the parameter set included in that .mat file. 

##### To convert theta-->k or adjust geometry parameters: 
Click "Convert axes" to convert to momenta with the specified experimental geometry parameters. 

#### If you already have set of feature maps output from gui: 
Click "Map error" button to calculate the standard deviation of each feature map, defined in the calculate_map_error.m script. \
Click "Replot maps" to replot a figure of all the feature maps. 

#### If you already have feature maps and map error calculations: 
Select a feature in the drop-down menu at the bottom under "H axis" to set the spectral feature histogram that will be used to bin the dataset along the horizontal axis. \
Type the number of bins you want to divide the horizontal axis histogram into in the window "Num bins" next to it. \
If you want to further bin along a second axis, select another feature map in the drop-down menu under "V axis" and a number >1 in the "Num bins" next to it. \
Click "Make BS" to plot a series of binned spectra partitioned along the horizontal and/or vertical axis as selected. \
Select the "k-shift" option if you want to also add a k-correction to each individual spectrum before including into the binning series figure. 
