This example demonstrates the calculation of a series of cc-AFM images of the 2-Iodotriphenylene on Ag(111) from the full electronic density. The calculations include the effects of the probe-particle deflection and the probe-tip interaction. Starting point are the charge density of the sample (CHG_sample), the charge density of the probe-particle (CHG_tip) and the electrostatic potential (LOCPOT) as they are calculated by VASP (all can be found in the example folder). The resulting images should be similar to those shown in ... figure ... (they will not be the same, since the resolution of the files has been reduced to save space)

# Step 0: File preparation
The output files of VASP are not in the xsf format, they have to be converted first. This can be done using the xsfconvert program, which can be found here: https://github.com/jenskunstmann/xsfconvert/tree/master
```
v2xsf CHG_sample -d -o CHG_sample.xsf
v2xsf CHG_tip -d -o CHG_tip.xsf
v2xsf LOCPOT -d -o LOCPOT.xsf
```
With this the files have the proper format and scaling.

# Step 1: AFM calculations

## 1a) params.ini
A full parameter file is present in the example directory. While nothing has to be changed here, some parameters should be disussed. As discussed in the basic workflow, most of the parameters are set as in a regular constant height AFM calculation. In this case, the setup is chosen to emulate a CO molecule adsorbed on a silver tip. The parameters which have to be chosen specially for the cc-AFM calculation are the scan range (ScanMin and ScanMax) and the scan resolution (ScanStep). Since the unit cell of the system is rectangular, the scan range was chosen to be as large as the unit cell in x and y direction. While the z direction can also ecompass the whole cell, the scan range has been reduced, since the values above and below the thresholds will not contribute to the cc-AFM images. This saves on input and output operations. The resolution has been increased in z-direction to accomodate for the sensetivity of the cc-AFM method in this direction.

## 1b) calculate the force fields
The force fields are computed with the full density method. We choose a beta of 1.12 and an A of 18 for the Pauli forces, the dz2 model for the tip for the electrostatic forces and the PBE functional parameters for the van-der-Vaals forces.
```
ppafm-conv-rho -s CHG_sample.xsf -t CHG_tip.xsf -B 1.12
ppafm-generate-elff -i LOCPOT.xsf -t dz2
ppafm-generate-dftd3 -i CHG_sample.xsf --df_name PBE
```
## 1c) relax probe particle
Using the force fields, the probe-particle has to be relaxed. Since the probe-particle deflections will be used later in the generation of the isosurface, they have to be saved using the --pos option.
```
ppafm-relaxed-scan --noLJ -A 18 --pos
```
The amplitude for the Pauli forces also have to be provided, henceforth the option -A 18.

## 1d) calculate frequency shifts
Lastly the frequency shifts have to be calculated and saved in the xsf format:
```
ppafm-plot-results --df --save_df
```
With this command the frequency shifts are outputted in the xsf format and as actual images. The images are not needed for the further calculations, but can be useful to check if something goes wrong in the calculation. Remove --df if no image output is wanted.

# Step 2: Generate the effective iso surfaces
Now the effective charge iso surface has to be calculated. The effective isosurface differs from a regular charge isosurface, since it includes position changes due to the probe-particle deflections, potential decreases or increases in tunnel current due to the probe-tip interaction, and an averaging over the oscillation amplitude. For this, several parameters have to be known from the AFM calculation:
  1) The particle deflections
  2) The equilibrium distance between probe particle and tip
  3) The oscillation amplitude

The probe-tip interaction is modeled as an exponential function, for which the decay constant has to be set. This is dependent on the tip and adsorbed molecule, and generally not known beforehand. It is therefore best to test several decay constants. Furthermore, the isovalue has to be set. It is again useful to consider several isovalues.
In order to avoid repeating the reading of the files (which can take a considerable amount of time), a range of those parameters can be considered in a single call of the program:
```
python extract_isosurface_height_map_parallel.py CHG_sample.xsf Test_calculation 0.1 -ppd Q-0.05K0.15 -R 3.0 -A 1.0 -dr "0 6 5" -ir "0.000002 0.00002 0.000002" -s 0.65
```
The first three define the charge density file, the prefix for the output files, and the isovalue. However, the isovalue is overwritten by options later on.

The option -ppd sets the directory where the particle deflections are stored. Those are called PPpos_x.xsf, PPpos_y.xsf and PPpos_z.xsf and are in the directory created by ppafm (directory name is dependent on parameter choice). The options -R defines the tip-probe equilibrium distance, and -A defines the oscillation amplitude. Both parameters have to be chosen to be the same as in params.ini. The parameters -dr and -ir define the different decay constants and isovalues for which the isosurface is computed. The syntax is similar to the python range function: The first value within the quotes is the starting point of the range, the second is the end point of the range (excluded), and the last value is the step size. E.g. with -dr "0 6 5" the decay constants of 0 and 5 are considered. Lastly the option -s defines the starting point of the scan, relative to the z-axis. The starting point should idealy be within the vacuum region.

With the option -dr set, a subdirectory is created for each decay constant. Each directory contains the isosurfaces for each isovalue. The first line of those files indicate the number of x and y grid points, the second and third define the x and y vector of the cell. The other lines represent the isosurface in the xyz format.
