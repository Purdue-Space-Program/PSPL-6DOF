# PSPL 6DOF
Six Degrees-of-Freedom model for PSP Liquids. This code is meant to be versatile and to be able to model a variety of rockets, from small scale solids to large scale liquid systems. The primary use case of the code is to run the modeling for the PSP Liquids launch vehicles, but is more generally capable.

This code is intended to be dynamic and model a variety of different rocket characteristics in a general way. This code is also intended to model sensors and other characteristics neccesary for GNC of the rocket.

## Getting Started

An example script ```MainRK4``` is included in the repo as an example of a simulation. Simply running this simulation will output the results of the simulation.

### Creating a custom simulation

The main simulation is comprised of a few main components. All simulations start by defining a ```Rocket``` object, which contains information on the  
structure of the rocket system. A default rocket is created upon calling ```Rocket```, but custom parameters can be added upon initialization of the object.

Any simulation also needs to define an ```Environment```, which contains data on the launch location, as well as weather, elevation, local gravity, and launch date. Most of this data comes for free with the specification of the launch location! A default location and time are automatically set for this object.