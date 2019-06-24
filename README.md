# CLAMS

CLAMS is a Computational Library for Analysis of Mass Spectra, or in short, a peak picker developed in R.
In this manner, it accepts MZML files and identifies a discrete set of peaks with their associated experimental data. 

Authors: Chris Dejong, Michael Cannon


This package works best with high resolution mass spectrometers such as a microTOF, which is sensitive to at least 2 decimal places.

Utilizes bioconductor package "xcms" for reading the mass spec and drawing EICs

Functioning:
Will Group ranges of peaks that are really a single peak
Will find isotope peaks
Chlorine checker, to see if the compound likely has a chlorine based on the isotopic distribution

Other dependencies:
Rpackages: tools, data.table, arg.parser, mzR, xcms

