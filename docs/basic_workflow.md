This workflow shows the general calculation of cc-AFM images. It results in a single cc-AFM image based on the pure isovalue. For a more refined calculation (including tip deflections) and the calculation of several images at once, refer to the advanced example.

# Requirements

In order to calculate cc-AFM images the electronic density (and if applicable the partial electronic density) of the system is needed. This density has to be in the xsf format and on a surface-cell grid: 
The first two grid vectors have to be parallel to the xy plane (surface), the third one has to be parallel to the z axis (perpendicular to the surface). For this workflow it will be assumed that the charge density file
has the name CHG.xsf. Also some familiarity with the ppafm program package is assumed. 

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

