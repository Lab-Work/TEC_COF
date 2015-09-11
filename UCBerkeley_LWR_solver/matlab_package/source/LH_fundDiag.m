classdef LH_fundDiag
%LH_fundDiag is an abstract class for any fundamental diagram. Any
%fundamental diagram should be a subclass of this class.

   properties
       kappa;%maximum density
   end

   methods (Abstract)
       q = flow(k)
       %fundamental diagram
       
       v = wspeed(k)
       %derivative of the fundamental diagram
       
       r = R(v)
       %Convex transform of the fundamental diagram
       
       [k1,k2] = densities(v,g)
       %Calculates the possible densities if a flow g is observed by an
       %observer at speed v, ie k1,k2 are solutions of flow(k) = k*v + g
       %Remark that for v = 0, it inverses the fundamental diagram.
       
       k = density(v)
       %Opposite of the derivative of R, or equivalently inverse of wspeed.
       %Needed only to compute exact densities. v may be a vector.
   end
end