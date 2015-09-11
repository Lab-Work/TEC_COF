classdef LH_Greenshields < LH_fundDiag
%LH_GREENSHIELDS Greenshields fundamental diagram class
%   Detailed explanation goes here

   properties
       v_f;%free flow speed
   end

   methods
       function obj = LH_Greenshields(vf,kmax)
           obj.v_f = vf;
           obj.kappa = kmax;
       end
       
       function q = flow(obj,k)
           q = obj.v_f*k*(1-k/obj.kappa);
       end
       
       function v = wspeed(obj,k)
           v = obj.v_f*(1-2*k/obj.kappa);
       end
       
       function r = R(obj,v)
           r = obj.kappa * (obj.v_f-v).^2 / (4 * obj.v_f);
       end
       
       function k = density(obj,v)
           k = NaN*ones(size(v));
           velIsPossible = v >= -obj.v_f & v <= obj.v_f;
           k(velIsPossible) = obj.kappa * (obj.v_f-v(velIsPossible)) / (2 * obj.v_f);
       end
       
       function res = densities(obj,v,g)
           d = sqrt((1 - v/obj.v_f)^2 - 4*g/(obj.kappa * obj.v_f));
           k1 = (1 - v/obj.v_f - d) * obj.kappa / 2;
           k2 = (1 - v/obj.v_f + d) * obj.kappa / 2;
           res = [k1 k2];
       end
           
   end
end 
