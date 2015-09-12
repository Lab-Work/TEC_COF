Update in MPC version
Improvements:
1. The new version of the code uses matlab struct (similar to dictionary) instead of matrix indexing. The readability of the code is much better and there is no need for the initEnv.m file for setting the global indexing variables.
2. Rewrote the dynamic grid algorithm, which is now cleaner.
3. Re-derived the theory. The entropic conditions are now documented in pdf. 
4. Added junction types onrampjunc, and offrampjunc, where onramp and offramps are controllable.
4. All array are defined as column vectors. 
5. Re-name functions and variables to make the code more readable and more precise, such as entropicSolutionAtJunc, sendingFuncAtLink, findBoundaryIntersection.


Details:
1. initNetwork. We no longer assume even discretization of the space. The input to setInitialCon is now a struct with fields: IC, X_grid_cum.
2. changed the INDEX_UP and INDEX_DOWN to struct (linkStr).upstream; (linkStr).downstream
3. Added a samplePoint function which only compute the vehicle id at points.
4. Added extractDensity function, which can be used to extract the density estimate at certain time.
5. In postSolution. If t_interval is shorter than dt_tol, return the middle point as the approximated found intersection point.
6. Added several functions for extracting the traffic density at any time, including searchFuncSlope, searchShocksOnLine, extractDensity.
7. Added utility functions such as groupSameElement to class.
8. Debugged the current control toolbox is working for estimation.
9. Changed the variable names in setIneqConstraints to more intuitive names.
10. Added an errors structure which supports different levels of error for each type of measurement.


Todo:
- The code was originally developed in Java in which the index begins with 0. Modified the index for internal and density conditions, and now index starts from 1.
- enable internal and density conditions to the framework. The density condition is helpful for setting the objective as minimizing the L1 error.
- update the dv_index structure.
- When converting min(s1, s2) to linear constraints by adding bools variables, we used C=500000 which we thought was large enough. Now we used inf directly.
- investigate AIMSUN ramp metering scheme and define the signal format.
- auto logging with time stamps
- unit conversion. Easier to read the density and flow.


Issues:
- We have not defined an offramp actuator yet. We use a green red cycle ramp meter for onramp control.
- Need to properly define a meaningful control objective.
- In find slope functions, if the interval is too small, the functions just return NaN which will introduce an error. Fix by computing an approximated slope.




Goal:
2. add density and internal conditions
3. Write a python script for to test if the toolbox works for control
4. Remote back and work on the API.




