# cc-AFM
Code for constant current atomic force microscopy simulations
# Overview
This code allows for the generation of constant current atomic force mircroscopy (cc-AFM) images, meaning AFM images, where the height is controlled by a constant current scanning tunnel microscopy (STM) feedback loop. This allows the imaging of bulky or corrugated molecules on surfaces, in contrast to the regular constant height AFM images. The code takes (partial) charge densities calculated by density functional theory, as well as the AFM frequency shifts calculated by the probe-particle model (ppafm) as input, and generates the constant current AFM and STM images as output.

# Installation

The code consists of two python scripts and thus only requires the installation of some libraries (listed below), which can be installed via pip. A virtual environment for the installation is recommended. The code was tested with python 3.9.21.

Required libraries:
pip install matplotlib
pip install numpy
pip install multiprocess
