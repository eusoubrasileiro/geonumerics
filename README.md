Geonumerics
===========

Numerical tools for geology and geophysics in python and C. 
This repo is my *sandbox*.

Most of the new code (work related created tools) are jupyter notebooks at `ipython_notebooks` folder. Some of those notebooks are prototypes to the **Fatiando a Terra** project where I focus on seismic modelling. Some are directly related to my work, for example, seismic and digital signal processing,  vibroseis field parameter test, vibroseis field data qc and tools for segy manipulation: coordinates fix, amplitude normalize and plotting. Finally a few notebooks are result of free time effort trying to implement reverse time depth migration. 

Here are some of the more recent (2012+) prototypes based on [fatiando style](https://github.com/fatiando/prototypes):
They use some of the branches of `fatiando_seismic` repository:

* [Analytic solution for scalar wave equation](http://nbviewer.ipython.org/github/eusoubrasileiro/geonumerics/blob/master/ipython_notebooks/Fatiando%20-%20F.D.%20vs%20Analytic%20Solution.ipynb) based on Alford 1974 Gephy. and testing of explicit finite diferences wave solution
* [Scalar wave equation in 3D](http://nbviewer.ipython.org/github/eusoubrasileiro/geonumerics/blob/master/ipython_notebooks/Fatiando%20Scalar3.ipynb) simple test for Reynolds 1D absorbing boundaries.
* [Nmo correction](http://nbviewer.ipython.org/github/eusoubrasileiro/geonumerics/blob/master/ipython_notebooks/Geonumerics%20-%20Nmo%20and%20Rms%20velocity.ipynb) base on Oz Yilmaz in a simple sythentic data example.
* [Pos-stack reverse time depth migration](http://nbviewer.ipython.org/github/eusoubrasileiro/geonumerics/blob/master/ipython_notebooks/Fatiando%20RTM%20zero-offset.ipynb) for zero-offset sections a simple example

This repo also includes a little DSP library I developed in C.

**Historical background:**
This is a collection of some tools I developed along the years. I started with C and now I only use python. 
Some tools are related to projects at college (long time ago), work or simple for learning. 
This repo is a mix between my old google-code repo, work related code, fatiando under development tools (not in master branch) and cool stuff I am trying to learn.

