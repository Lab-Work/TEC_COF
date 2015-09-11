classdef LH_Tfd < LH_fundDiag
%LH_TFD Triangular fundamental diagram class

   properties
       v_f;%free flow speed
       w;%congestion wave speed
       k_c;%critical density
   end

   methods
       function obj = LH_Tfd(ffspeed,swspeed,kmax)
           obj.v_f = ffspeed;
           obj.w = swspeed;
           obj.kappa = kmax;
           obj.k_c = -swspeed*kmax/(ffspeed-swspeed);
       end
       
       function q = flow(obj,k)
           if(k < obj.k_c)
               q = k*obj.v_f;
           else
               q = obj.w*(k-obj.kappa);
           end
       end
       
       function v = wspeed(obj,k)
           if(k < obj.k_c)
               v = obj.v_f;
           else
               v = obj.w;
           end
       end
       
       function r = R(obj,v)
           r=obj.k_c * (obj.v_f-v);
       end
       
       function k = density(obj,v)
           k = NaN*ones(size(v));
           velIsPossible = v >= obj.w & v <= obj.v_f;
           k(velIsPossible) = obj.k_c;
       end
       
       function res = densities(obj,v,g)
           k1 = g/(obj.v_f-v);
           k2 = (obj.w*obj.kappa+g)/(obj.w-v);
           res = [k1 k2];
       end
           
   end
end 
