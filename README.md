# Traffic Estimation and Control (TEC) toolbox based on Convex Optimizatoin Framework (COF).
Yanning Li, Feb 18th 2016, yli171@illinois.edu

##1) Overview:
This toolbox provides a numerical scheme for computing the traffic evolution on a network containing an unsignalized merge or diverge using a convex optimization. Given the measurement data at the network boundaries (the endpoints of links not connected to a junction), this toolbox computes the how traffic flows interact at the junction and accordingly computes the traffic density on all links.

The toolbox can also be used for optimal traffic control by constructing an optimization program. An optimal on-ramp metering control example is presented in this toolbox and discussed in the associated article.

##2) License

This software is licensed under the *University of Illinois/NCSA Open Source License*:

**Copyright (c) 2013 The Board of Trustees of the University of Illinois. All rights reserved**

**Developed by: Department of Civil and Environmental Engineering University of Illinois at Urbana-Champaign**

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal with the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimers in the documentation and/or other materials provided with the distribution. Neither the names of the Department of Civil and Environmental Engineering, the University of Illinois at Urbana-Champaign, nor the names of its contributors may be used to endorse or promote products derived from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.

##3) Installation:
1. Download and Install [IBM CPLEX optimization studio](http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud) which is free for academic use. 
2. Add the following paths to MATLAB environment:
  - (CPLEX installation folder)/IBM/ILOG/CPLEX_Stuidio1262/cplex/matlab
  - [(this toolbox folder)](https://github.com/Lab-Work/TEC_COF)
  - [(this toolbox folder)/COF_toolbox/matlab](https://github.com/Lab-Work/TEC_COF/tree/master/COF_toolbox/matlab)
  - [(this toolbox folder)/UCBerkeley_LWR_solver/matlab_package/source](https://github.com/Lab-Work/TEC_COF/tree/master/UCBerkeley_LWR_solver/matlab_package/source)

##4) Run the examples:
Two examples are included in the toolbox, respectively a merge solver and an optimal on-ramp metering controller.

### A merge solver:
Run the script [merge_solver.m](https://github.com/Lab-Work/TEC_COF/blob/master/merge_solver.m). Refer to the published html [merge_solver.html](https://github.com/Lab-Work/TEC_COF/blob/master/html/merge_solver.html) for a more readable description. 

### Optimal on-ramp metering control
An optimal on-ramp metering controller example is provided in this toolbox. To run this example, you need to use [AIMSUN](https://www.aimsun.com/wp/) to simulate a microscopic traffic environment.
- Install [AIMSUN](https://www.aimsun.com/wp/).
- Download the microscopic model we built in AIMSUN and the python scripts we developed for AIMSUN to interact with MATLAB from [link](). We refer to the documentation in [link]() for the configuration of AIMSUN for running this example.
- Copy the [AIMSUN_MATLAB_COM foler](https://github.com/Lab-Work/TEC_COF/tree/master/AIMSUN_MATLAB_COM) to **E: Drive**. You may copy to other path, but you need to modify paths in the [runMPC_workzone.m script](https://github.com/Lab-Work/TEC_COF/blob/master/runMPC_workzone.m), the [com_config file.txt](https://github.com/Lab-Work/TEC_COF/blob/master/AIMSUN_MATLAB_COM/COM_CONFIG.txt), and the python scripts for AIMSUN.
- Run the script [runMPC_workzone.m script](https://github.com/Lab-Work/TEC_COF/blob/master/runMPC_workzone.m) and start simulation in AIMSUN. 

##5) Folders:

### [\AIMSUN_MATLAB_COM](https://github.com/Lab-Work/TEC_COF/tree/master/AIMSUN_MATLAB_COM)
This folder includes the files needed for running the optimal on-ramp metering control example. 
- [COM_CONFIG.txt](https://github.com/Lab-Work/TEC_COF/blob/master/AIMSUN_MATLAB_COM/COM_CONFIG.txt) This file defines the path files shared by MALTAB and AIMSUN for communication.
- [historical_data.txt](https://github.com/Lab-Work/TEC_COF/blob/master/AIMSUN_MATLAB_COM/historical_data.txt) This file contains the historical data for the model predictive control. It is same as the demand configured in AIMSUN.
- [\file_format](https://github.com/Lab-Work/TEC_COF/tree/master/AIMSUN_MATLAB_COM/file_format) This folder contains the example files used in communication between AIMSUN and MATLAB.

### [\COF_toolbox](https://github.com/Lab-Work/TEC_COF/tree/master/COF_toolbox)
This folder contains the code developed for the toolbox. Specifically, the following classes are developed for the solver and the optimal controller. Type *doc classname* in MATLAB command window for the documentation. 
- [initNetwork.m](https://github.com/Lab-Work/TEC_COF/blob/master/COF_toolbox/matlab/initNetwork.m) This class constructs a network object containing the network information, initial and boudnary conditions. 
- [optProgram.m](https://github.com/Lab-Work/TEC_COF/blob/master/COF_toolbox/matlab/optProgram.m) This class builds a convex program including the constraints and the objective function, and uses CPLEX to solve the convex program. 
- [postSolution.m](https://github.com/Lab-Work/TEC_COF/blob/master/COF_toolbox/matlab/postSolution.m) This class contains functions for visualizing the solution and iteratively updating the boudnary grid to compute the exact solution. 
- [rampController.m](https://github.com/Lab-Work/TEC_COF/blob/master/COF_toolbox/matlab/rampController.m) This class defines the model predictive control framework, which is used to construct an optimal on-ramp metering controller.
- Other classes should not be used directly. 

### [\UCBerkeley_LWR_solver](https://github.com/Lab-Work/TEC_COF/tree/master/UCBerkeley_LWR_solver) 
This folder contains the code developed by Berkeley for computing the solution of LWR PDE on a single link.

### [\html](https://github.com/Lab-Work/TEC_COF/tree/master/html) 
This folder contains the more readable documentation of the two examples scripts. 










