This example demonstrates the calculation of a series of cc-AFM images of the 2-Iodotriphenylene on Ag(111). The calculations include the effects of the probe-particle deflection and the probe-tip interaction. Starting point are the charge density of the sample (CHG_sample), the charge density of the probe-particle (CHG_tip) and the electrostatic potential (LOCPOT) as they are calculated by VASP (all can be found in the example folder). The resulting images should be similar to those shown in ... figure ... (they will not be the same, since the resolution of the files has been reduced to save space)

# Step 0: File preparation
The output files of VASP are not in the xsf format, they have to be converted first. This can be done using the xsfconvert program, which can be found here: https://github.com/jenskunstmann/xsfconvert/tree/master
```
v2xsf CHG_sample -d -o CHG_sample.xsf
v2xsf CHG_tip -d -o CHG_tip.xsf
v2xsf LOCPOT -d -o LOCPOT.xsf
```

