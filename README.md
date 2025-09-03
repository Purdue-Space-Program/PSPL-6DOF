# PSP Liquids 6-DOF  Guide

This readme details the use and the implementation of the 6-DoF model for PSP Liquids. This code is meant to be versatile and to be able to model a variety of rockets, from small scale solids to large scale liquid systems. The primary use case of the code is to run the modeling for the PSP Liquids launch vehicles, but is more generally capable.

This code is intended to be dynamic and model a variety of different rocket characteristics in a general way. This code is also intended to model sensors and other characteristics neccesary for GNC of the rocket.

This document will not dive deeply into the methodology and mathematics behind the model. For that, please refer to our 'Flight Dynamics Bible' book series:
- [Volume I](https://www.overleaf.com/read/cgjkwxwwzxmc#1818cf)
- [Volume II](https://www.overleaf.com/project/67c9db5f55304a0119638f83)
- [Volume III](https://www.overleaf.com/read/ttdvwvvznjgb#da02d8)

## Code Structure

The PSPL 6-DoF takes an Object-Oriented Programming (OOP) approach. As a result, each of the primary components of the simulation are called through a class structure in MATLAB. To access this class structure, the folder `TheSixDoF` must be on the MATLAB file path. 

The OOP approach helps us to maintain a clean code base and store related objects within a single MATLAB object. As our codebase expands, this class structure becomes increasingly important. However, most MATLAB users are unfamilar with its class structure. A desciption of the MATLAB class structure may be found in the final chapter of Volume I.

## Getting Started

An example script ```MainRK4``` is included in the repo as an example of a simulation. Simply running this simulation will output the results of the simulation. This model defaults to a run with the CMS rocket. Various settings in the model may be adjusted to the liking of the user.

Simulation can start off very simple and complexity can be added thanks to the class structure design. Every simulation will require a `Rocket` object, as well as some basic `Sim` parameters to be set.

### Creating a custom simulation

The main simulation is comprised of a few main components. All simulations start by defining a ```Rocket``` object, which contains information on the  
structure of the rocket system. A default rocket is created upon calling ```Rocket```, but custom parameters can be added upon initialization of the object.

Any simulation also needs to define an ```Environment```, which contains data on the launch location, as well as weather, elevation, local gravity, and launch date. Most of this data comes for free with the specification of the launch location! A default location and time are automatically set for this object.


If you have any further questions, feel free to reach out to [Hudson Reynolds](mailto:reyno250@purdue.edu) or refer to any of our books.
