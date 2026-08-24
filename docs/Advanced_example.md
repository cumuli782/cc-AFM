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
A full parameter file is present in the example directory. As discussed in the basic workflow, most of the parameters are set as in a regular constant height AFM calculation. In this case, the setup is chosen to emulate a CO molecule adsorbed on a silver tip. The parameters which have to be chosen specially for the cc-AFM calculation are the scan range and the scan resol
