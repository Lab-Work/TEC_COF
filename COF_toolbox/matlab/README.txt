Yanning Li, Sep 01, 2015

This toolbox is based on the network optimizer toolbox.
The network optimizer toolbox contains code for optimal traffic control on a network, assuming all intersections are equipped with actuators (e.g. traffic lights).

The original toolbox contains the following code:
1. The core code for setting the constraints defined in classes.
2. The code for integrating the satellite image and road topology. 
3. Several matlab scripts that to be executed consecutively to set up the network, define the optimization program, solve the program, and plot the result.


Modification:
1. The core code has been modified which can handle uneven discretization of the space to remove the discretization error. 
2. We no longer need the road topology code.
3. The matlab scripts will be too complicated in a large problem. Hence, those code are completely rewritten in classes. 
4. Overall, this new toolbox contains completely new code which enables traffic estimation involving multiple sensors on the stretch of road a a highway intersection.