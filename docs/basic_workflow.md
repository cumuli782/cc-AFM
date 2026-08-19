This workflow shows the general calculation of cc-AFM images. It results in a single cc-AFM image based on the pure isovalue. For a more refined calculation (including tip deflections) and the calculation of several images at once, refer to the advanced example.

# Requirements

In order to calculate cc-AFM images the electronic density (and if applicable the partial electronic density) of the system is needed. This density has to be in the xsf format and on a surface-cell grid: 
The first two grid vectors have to be parallel to the xy plane (surface), the third one has to be parallel to the z axis (perpendicular to the surface). For this workflow it will be assumed that the charge density file
has the name CHG.xsf. Also some familiarity with the ppafm program package is assumed. *Note:* The ppafm programs are parallelized with OpenMp, the program extract_isosurface_height_map_parallel.py is parallelized using the multiprocess package. Therefore, they are best run on machines with several cores.

# Step 1: AFM calculation

## 1a) Choosing ppafm model
In order for the cc-AFM calculation to yield reasonable results, the height of the isosurface has to fit the changes in frequency. For example, if the radii in the Lennard-Jones model are chosen to large, the frequencies change at a higher z coordinate than it should be the case. In the best case, this results in an offset between charge density and frequencies, in the worst case the result becomes entirely unusable. Therefore, one has to choose a model where frequencies and density fit to each other. There are two options for this:
1) (recommended): Use the full density model, where the Pauli forces are calculated via a convolution of tip density and sample density. This requires some experimentation with the parameters A and beta of the model (see https://github.com/Probe-Particle/ppafm/wiki/Forces#full-density-based-model).
2) Fit the atomic radii of the Lennart-Jones model. This can be done using the script "fitPauli.py" found in "ppafm/cli/utilities".

## 1b) Preparing input

First, the 3D stack of constant height AFM images has to be calculated using the ppafm program. For this, the parameter file "params.ini" has to be created. For the structure of this file and the choice of parameters, refer to the ppafm wiki https://github.com/Probe-Particle/ppafm/wiki. Most of the parameters can be chosen as in an ordinary AFM calculation, with two exceptions: 
1) The scan area (defined by ScanMin and ScanMax) has to be smaller or equal to the cell size of CHG.xsf -- the scan area has to fit in the charge density cell. This means, the CHG cell has to be big enough, that the area of interest (the molecule and its surroundings) can be covered by the rectangular cell defined by scanMin and scanMax. The z coordinate of scanMin and scanMax can be chosen as large as the CHG cell, but can be chosen smaller in order to save on output data.
2) The cc-AFM images are sensitive to changes in z direction. Therefore a finer grid is recommended in z direction. Change ScanStep accordingly, for example ScanStep 0.1 0.1 0.01

## 1c) Doing the calculations
For the force field calculation the usage of the full density model is assumed, the file names are assumed to be CHG_tip.xsf for the tip density and LOCPOT.xsf for the electrostatic potential. For other models, refer to the ppafm wiki. 
Generate forcefields:
```
ppafm-conv-rho -s CHG.xsf -t CHG_tip.xsf -A 1.0 -B 1.1
ppafm-generate-dftd3 -i CHG.xsf --df_name PBE
ppafm-generate-elff -i LOCPOT.xsf -t dz2
```
Replace the chosen Beta (here 1.1), the chosen DFT functional (here PBE) and the tip model (here dz2) to your needs. The Amplitude for the Pauli forces are changed in the next step.

Relax probe particle:
```
ppafm-relaxed-scan --noLJ -A 18.0
```
Here the Amplitude of the Pauli forces can be changed to your needs.

Calculate the frequency shifts:
```
ppafm-plot-results --save_df
```
This saves the frequency shifts as a single xsf file. Though it is not neccessary for the next steps, one can also add the option "--df" to save the frequency shifts as png images in order to check the calculation.

# Step 2: (effective) isosurface calculation
Now the isosurface has to be calculated. This workflow covers only the pure isosurface, in order to include the probe-particle deflections and the probe-tip interaction refer to the advanced example. The isosurface corresponds to a surface of constant electronic density, which translates (in the Tersoff-Hamann approximation) to a constant tunnel current. This is done in the script "extract_isosurface_height_map_parallel.py"
```
python extract_isosurface_height_map_parallel.py CHG.xsf output_iso.txt 0.00001 -s 0.6
```
output_iso.txt is the name of the output file, the third value is the isovalue. Depending on the program used to generate and extract the density, it can have different units (for example charge/cellvolume, charge/A^3, charge/cell). In the case of charge/cell, the isovalue is typicall in the range of 10^-4 to 10^-6. The option "-s" defines the start of the scan for the isosurface, given as fraction of the cell height (default 0.9). In principle this can be left out, the code can scan both upwards and downwards. However, sometimes there are artifacts in the density in the form of small areas of heightend or lowered density, especially near the atom cores. Therefore a start in the vacuum region is preferable.

# Step 3: cc-AFM image calculation
Finally by using the calculated frequency shifts and the isosurface, the cc-AFM image is calculated. This is done in the script "Afm_at_Isosurface.py".
```
python Afm_at_Isosurface.py Q-0.05K0.20/Amp1.00/df.xsf output_iso.txt cc-AFM_image.png
```
df.xsf contains the frequency shifts calculated in step 1. The name of the folder is dependent on your choice of parameters. output_iso.txt is the isosurface file calculated in step 2. cc-AFM_image.png is the output image. Additionally, the height profile is printed in the file cc-AFM_image_height.png
