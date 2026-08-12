# cc-AFM
Code for constant current atomic force microscopy simulations
# Overview
This code allows for the generation of constant current atomic force mircroscopy (cc-AFM) images, meaning AFM images, where the height is controlled by a constant current scanning tunnel microscopy (STM) feedback loop. This allows the imaging of bulky or corrugated molecules on surfaces, in contrast to the regular constant height AFM images. The code takes (partial) charge densities calculated by density functional theory, as well as the AFM frequency shifts calculated by the probe-particle model (ppafm) as input, and generates the constant current AFM and STM images as output.

# Installation

The code consists of two python scripts and thus only requires the installation of some libraries (listed below), which can be installed via pip. A virtual environment for the installation is recommended. The code was tested with python 3.9.21.

Set up environment (optional):
```
python -m venv cc-AFM
source cc-AFM/bin/activate
```
Install libraries:
```
pip install matplotlib
pip install numpy
pip install multiprocess
```

# Basic usage

The generation of cc-AFM images happens in two steps: First, the script extract_isosurface_parallel.py extracts a height profile from the charge density dependent
 on a provided iso value -- the profile essentially follows this iso surface. A higher iso value corresponds to a higher STM tunnel current. The provided density 
has to be in the xcrysden format (.xsf). The basic usage of this script is
```
python extract_isosurface_parallel.py <density.xsf> <output_file> <iso_value>
```

In a second step, the frequency shifts are collected according to the height profile. This occurs in the script Afm_at_Isosurface.py. This script takes the freque
ncy shifts at different heights calculated by the ppafm program (also in xsf format), and interpolates the frequency at the height according to the height profile
. The basic usage of this script is
```
python Afm_at_Isosurface.py <Frequency_shifts.xsf> <height_profile> <output_file.png>
```
The whole workflow and a more avanced example are given in the wiki.

# Citation
This code is free to use and free to modify. However, we appreaciate the citation of the paper, where this code has been introduced:
