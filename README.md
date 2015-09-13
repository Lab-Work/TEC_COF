# Traffic Estimation and Control (TEC) toolbox based on Convex Optimizatoin Framework (COF).
Yanning Li, Sep 13 2015, emlynlyn@gmail.com

### This toolbox supports:
1. traffic estimation on a network (connection, merge, diverge), given the boudnary conditions. 
2. traffic control on a network with controllable ramps.

### Features:
1. The toolbox supports uneven discretization of both time and space. 
2. The toobox has full support of internal and density conditions. 
3. It has a iterative grid updating component to resolve the discretization error.
4. It has the built-in entropic objective function component, which gives the entropic solution for each step. When adding control objectives, the user need to check the theory to make sure the new control objective function satisfies the entropic conditions.
5. Compared to previous version, this version is class based; uses dictionaries (MATLAB struct) insteaded of using global indexing; uses MATLAB index (starts from 1) instead of original index (based on java starts from 0); has much more intuitive and precise function and variable names.
6. The user only need to write the main scipt which creates objects of classes for different estimation and control scenarios.

### Structure of the toolbox:
-script: This script defines the network and application. Run as main function.
-initNetwork: the class that contains all necessary topology information and data for a network.
-optProgram: the class that builds an optimization program and output solutions.
-postSolution: this class visualize the solution. It also iteratively update the discretization grid to resolve the discretization error.

### Installation:
1. Install CPLEX.
2. Add the following path to MATLAB:
  - (CPLEX installation folder)/IBM/ILOG/CPLEX_Stuidio1262/cplex/matlab
  - (this toolbox folder)
  - (this toolbox folder)/COF_toolbox/matlab
  - (this toolbox folder)/UCBerkeley_LWR_solver/matlab_package/source

### Run the examples:
  Simply run the main script in MATLAB.
Three exampels are provided, respectively
- testToolboxForEstimation.m: the estimation of traffic density on a merge.
- testToolboxForMPC.m: the model predictive control in a work zone. Need AIMSUN to run the control simulation.
