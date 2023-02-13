# Kerr Vector Map
Kerr Microscope Vector Mapping

This code is intended for use with Evico made Kerr microscopes.
The aim of this code is to take two images obtained simultaneously along the transverse and longitudinal directions and conver them into one vector image (see: https://link.aps.org/doi/10.1103/PhysRevB.95.014426).

Images can be obtained in this way through the use of LED-selective sensitivity. The longitudinal and transverse light levels are set to be the same and then, in the KerrLab software that accompanies an Evico microscope, it is possible to image longitudinal and transverse sensitives at the same time from the same background. 

This code will then read in both sets of magnetic images in the form outputted from KerrLab and perform the analysis necessary to obtain a 2D vector plot. To do this it converts the images into some vector information and then calculates the angle at which the magnetization is pointing within a pixel cell. 
Included are two images (labelled x and y) with which to test the code before use. 
Also included is some basic code to produce a line scan through the angular data.
