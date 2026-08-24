This example demonstrates the calculation of a series of cc-AFM images of the 2-Iodotriphenylene on Ag(111). The calculations include the effects of the probe-particle deflection and the probe-tip interaction. Starting point are the charge density of the sample (CHG_sample), the charge density of the probe-particle (CHG_tip) and the electrostatic potential (LOCPOT) as they are calculated by VASP (all can be found in the example folder). The resulting images should be similar to those shown in ... figure ... (they will not be the same, since the resolution of the files has been reduced to save space)

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

