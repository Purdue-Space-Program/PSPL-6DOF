# PSP Liquids 6-DOF  Guide

This document details the use and the implementation of the 6-DoF model for PSP Liquids. This code is meant to be versatile and to be able to model a variety of rockets, from small scale solids to large scale liquid systems. The primary use case of the code is to run the modeling for the PSP Liquids launch vehicles, but is more generally capable. This represents over 4000 lines of code and several years of combined effort.

This code is intended to be dynamic and model a variety of different rocket characteristics in a general way. This code is also intended to model sensors and other characteristics neccesary for GNC of the rocket.

This document will not dive deeply into the methodology and mathematics behind the model. For that, please refer to our 'Flight Dynamics Bible' book series, which covers these topics in great depth:
- [Volume I](https://www.overleaf.com/read/cgjkwxwwzxmc#1818cf)
- [Volume II](https://www.overleaf.com/project/67c9db5f55304a0119638f83)
- [Volume III](https://www.overleaf.com/read/ttdvwvvznjgb#da02d8)

## Getting Started

### Requirements:

To get started with this tool, please clone this repository to your local machine to get access to the tool.

 The GUI is designed for operation with MATLAB 2025a or 2025b. Compatibility with other versions of MATLAB is not guaranteed. Python 3.13 with the following packages is also required:
- openmeteo_requests
- pandas
- requests_cache
- retry_requests
- numpy

The required packages may be installed with:

```
pip install openmeteo-requests
pip install requests-cache retry-requests numpy pandas
```

More information about the install and use of the OpenMeteo packages can be found at [OpenMeteo](https://open-meteo.com/en/docs).

### Using the GUI


The main user interaction with the 6-DoF in the MATLAB app. This can be accessed via the `RocketGUI.mlapp` file within MATLAB. This file is found in the `TheSixDoF` folder.


 Upon opening the `RocketGUI.mlapp`, the user must first select the root directory. The `TheSixDoF` folder should be selected as the root directory. Any other folder will throw an error for the user.

 After loading this, the user can begin designing a rocket or load one of the pre-made rocket samples, such as `CMS.mat`. This rocket object defines all of the properties and components that are associated with your rocket.

 When creating a rocket from scratch, it is important to input all of the parameters on the main page. Especially important is the RASAero data. This data should be exported from RASAero II "Aero Plots". At the current time, our models cannot directly compute aerodynamics based off the input geometry (this is a top priority for us). Save the rocket file before performing any further actions.

 If working from a loaded object, make sure to save the rocket after changing any of the parameters.

 ### Plot Interaction
 The plot inside the GUI is specially designed for operation with rocket objects. For this reason, the standard zoom and pan tools in the toolbar should not be used. Simply click the plot to drag it around. In 3D, use the top toolbar to adjust the orientation and size of the model.

 ### Adding Components
 Adding components to the rocket is a straighforward process. In the components tab, any number of components can be added. Each component has specific properties which must be input. For a component, *every* property must have an input for the component to be generated. All inputs are in standard SI units (showing units is WIP).

 For the tanks component, the selection of 'Fuel' or 'Ox' must be input exactly as 'Fuel' or 'Ox'. Changing this to a drop down is a work in progress.

 For each component, the distance is given as X,Y,Z. X is the distance from the nose (positive towards back), and Y and Z are distances from the centerline.

 In the bottom of the window, a component browser shows all the components of the rocket. Double clicking on any component allows the user to delete that component from the rocket. No two components may share the same name.

 ## Running a simulation
 After creating a rocket object, the user may run a simulation in the simulation tab. The user can input the location and time of launch they wish to simulate. Weather is retrieved from OpenMeteo and a confirmation window appears upon success of this operation.

 Then, the user may simulate the launch, which runs the simulation until apogee and displays a number of plots and outputs associated with the vehicle.

 More options for the simulation of the vehicle launch are coming in the future.


If you have any further questions, feel free to reach out to [Hudson Reynolds](mailto:reyno250@purdue.edu) or refer to any of our books.
