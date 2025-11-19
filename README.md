# PSP Liquids 6-DOF  Guide

This readme details the use and the implementation of the 6-DoF model for PSP Liquids. This code is meant to be versatile and to be able to model a variety of rockets, from small scale solids to large scale liquid systems. The primary use case of the code is to run the modeling for the PSP Liquids launch vehicles, but is more generally capable.

This code is intended to be dynamic and model a variety of different rocket characteristics in a general way. This code is also intended to model sensors and other characteristics neccesary for GNC of the rocket.

This document will not dive deeply into the methodology and mathematics behind the model. For that, please refer to our 'Flight Dynamics Bible' book series:
- [Volume I](https://www.overleaf.com/read/cgjkwxwwzxmc#1818cf)
- [Volume II](https://www.overleaf.com/project/67c9db5f55304a0119638f83)
- [Volume III](https://www.overleaf.com/read/ttdvwvvznjgb#da02d8)

## Getting Started

The main user interaction with the 6-DoF in the MATLAB app. This can be accessed via the `RocketGUI.mlapp` file within MATLAB. The GUI is designed for operation with MATLAB 2025a or 2025b. Compatibility with other versions of MATLAB is not gauranteed.

To run this script, you **must** have `TheSixDoF` and all subfolders on your MATLAB path. The current directly should be the top level folder, `PSPL-6DOF`. 

### Creating a custom simulation

The main simulation is comprised of a few main components. All simulations start by defining a ```Rocket``` object, which contains information on the  
structure of the rocket system. A default rocket is created upon calling ```Rocket```, but custom parameters can be added upon initialization of the object.

Any simulation also needs to define an ```Environment```, which contains data on the launch location, as well as weather, elevation, local gravity, and launch date. Most of this data comes for free with the specification of the launch location! A default location and time are automatically set for this object.


If you have any further questions, feel free to reach out to [Hudson Reynolds](mailto:reyno250@purdue.edu) or refer to any of our books.
