# Janson.Ichter.Pavone.IJRR15
Deterministic Sampling-Based Motion Planning: Optimality, Complexity, and Performance

This repo contains the code used for results in the paper "Deterministic Sampling-Based Motion Planning: Optimality, Complexity, and Performance" submitted to IJRR in 2015 by Janson, Ichter, and Pavone of Stanford's Autonomous Systems Lab.

All the C code was run through the Open Motion Planning Library (OMPL) at http://ompl.kavrakilab.org/. The MATLAB planning is run through the runFMT.m file. The Julia code is run through the iPython notebook.

## Dependences
- For the Julia Code (kinodynamic planning), the iPython notebooks use the code from https://raw.githubusercontent.com/schmrlng/MotionPlanning.jl.
- For the C++ code, the Open Motion Planning Library (OMPL) was used and can be found at http://ompl.kavrakilab.org/.

## Disclaimers
This code is fairly rough and most likely useful as a reference only. It is also subject to changes, bugs, etc. 
