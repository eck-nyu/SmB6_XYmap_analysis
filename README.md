# SmB6_XYmap_analysis
Code updated and uploaded on Jan 5 2022.
Used in analysis published in Kotta et al., "Metallic chemical potentials in an insulating topological state." \
Direct questions to Erica Kotta at eck325@nyu.edu.

## Instructions on how to run the code:

Save all files in this repository into the same folder, then run XYmap_gui.m.

#### If you don't already have an XY scan dataset saved into a mat file: 
You can take a set of ARPES images, saved in the .itx format, saved into a folder that you specify within XYmap_gui.m script as variable data.itx_foldername and corresponding .csv file specified with data.csv_filename. \
The .csv file contains XY scan metadata, listing the manipulator X,Y coordinates of each ARPES image in the folder of itx files. \
Functions itx_csv_to_spec.m and csv_to_XY_list.m convert the dataset into a set of variables. You can now save these into a .mat file using the commented-out code located in the convertItx2Mat function inside the gui script. \

#### If you already have the XY scan dataset saved into a mat file: 
Select the mat file from the drop-down menu near the top of the gui, click the "Import mat" button. \
There is a set of default parameter values (FL E, BG E, ...) listed on the left hand side with default values for SmB6 spectral feature analysis. \
All parameters with "E" are energies, in units of kinetic energy. \
All paremeters with "th" are theta, units of degrees. \
All paremeters with "k" are momenta, units of inverse angstrom. \
See XYmap_gui.m script for more descriptions of the parameters. \
Click "Save" at bottom left to save the current set of parameter values into the selected .mat file along with the spectra dataset. \
Click "Run" to run the feature map analysis code as defined in run_the_params.m script and produce a set of feature maps. \
Default set of features included in this demo are: 
- total intensity (tot int)
- sample tilt (tilt)
- J=5/2 f-band energy (f E)
- J=5/2 f-band width (f S)
- J=7/2 f-band energy (f2 E) 
- J=7/2 f-band width (f2 S) 
- J=7/2:J=7/2 f-band amplitude ratio (f21 rat)
- Sm-multiplet feature intensity (Ism)
- Sm-multiplet feature energy (Ism E)

#### If you already have set of feature maps output from gui: 
Click "Map error" button to calculate the standard deviation of each feature map, defined in the calculate_map_error.m script. \
Click "Replot maps" to replot a figure of all the feature maps. \

#### If you already have feature maps and map error calculations: 
Select a feature in the drop-down menu at the bottom under "H axis" to set the spectral feature histogram that will be used to bin the dataset along the horizontal axis. \
Type the number of bins you want to divide the horizontal axis histogram into in the window "Num bins" next to it. \
If you want to further bin along a second axis, select another feature map in the drop-down menu under "V axis" and a number >1 in the "Num bins" next to it. \
Click "Make BS" to plot a series of binned spectra partitioned along the horizontal and/or vertical axis as selected. \
Select the "k-shift" option if you want to also add a k-correction to each individual spectrum before including into the binning series figure. \
