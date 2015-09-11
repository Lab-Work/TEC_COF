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



Todo:
- properly define a meaningful control objective.
- enable internal conditions to the framework; change the variable names in the code
- enable density conditions to the framework; change the variable names in the code
- Offramp actuator not defined yet.
- investigate AIMSUN ramp metering scheme and define the signal format.
- auto logging with time stamps
- unit conversion. Easier to read the density and flow.
- In find slope functions, if the interval is too small, the functions just return NaN which will introduce an error. Fix by computing an approximated slope.


Goal today:
2. add density and internal conditions
3. Write a python script for to test if the toolbox works for control
4. Remote back and work on the API.




