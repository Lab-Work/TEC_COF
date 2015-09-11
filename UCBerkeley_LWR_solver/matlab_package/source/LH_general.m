classdef LH_general < handle
%LH_general represents a particular HJ traffic problem.

    properties
        valConditions=[];
        funddiag;
        xmin;
        xmax;
    end

    methods
        
        % 2 lignes epxlications sur le fait qu il met FD etc. dans funddiag
        % etc.
        % This defines a problem environment, ie an object that contains 
        % a fundamental diagram FD
        % a minimum space value xmin (where upstream conditions apply)
        % a maximum space value xmin (where downstream conditions apply)
        function obj=LH_general(FD,xup,xdown)
            obj.funddiag=FD;
            obj.xmin=xup;
            obj.xmax=xdown;
        end
        
        % Sets initial conditions corresponding to the given initial
        % densities: for all x in [x(i),x(i+1)], density(0,x) = kx(i)
        % x is a vector of size [1 n+1] and kx of size [1 n].
        function setIniDens(obj,x,kx)
            if(size(x,1)~=1 || size(kx,1) ~= 1 || size(x,2)-1~=size(kx,2))
                error('check dimensions of x and kx');
            end
            N=obj.getStartingCount(0,x(1));
            lengths=diff(x);
            xbeg=x(1:(size(x,2)-1));
            xend=x(2:size(x,2));
            nCars=-kx.*lengths;
            totCars=N+[0 cumsum(nCars(1:end-1))];
            z = zeros(size(kx'));
            obj.valConditions=[obj.valConditions;z z xbeg' xend' nCars' totCars'];
        end
        
        % Sets upstream condition corresponding to the given upstream
        % flows: for all t in [t(i),t(i+1)], flow(t,xmin) = qt(i)   
        % t is a vector of size [1 n+1] and qt of size [1 n].
        function setUsFlows(obj,t,qt)
            if(size(t,1)~=1 || size(qt,1) ~= 1 || size(t,2)-1~=size(qt,2))
                error('check dimensions of x and kx');
            end
            N=obj.getStartingCount(t(1),obj.xmin);
            timeIntervals=diff(t);
            tbeg=t(1:(size(t,2)-1));
            tend=t(2:size(t,2));
            nCars=qt.*timeIntervals;
            totCars=N+[0 cumsum(nCars(1:end-1))];
            o = ones(size(qt'))*obj.xmin;
            obj.valConditions=[obj.valConditions;tbeg' tend' o o nCars' totCars' ];
        end
        
        % Sets upstream condition corresponding to the given downstream
        % flows: for all t in [t(i),t(i+1)], flow(t,xmax) = qt(i)   
        % t is a vector of size [1 n+1] and qt of size [1 n].  
        function setDsFlows(obj,t,qt)
            if(size(t,1)~=1 || size(qt,1) ~= 1 || size(t,2)-1~=size(qt,2))
                error('check dimensions of x and kx');
            end
            N=obj.getStartingCount(t(1),obj.xmax);
            timeIntervals=diff(t);
            tbeg=t(1:(size(t,2)-1));
            tend=t(2:size(t,2));
            nCars=qt.*timeIntervals;
            totCars=N+[0 cumsum(nCars(1:end-1))];
            o = ones(size(qt'))*obj.xmax;
            obj.valConditions=[obj.valConditions;tbeg' tend' o o nCars' totCars' ];            
        end
  
        % Assigns a value to the beginning of a value condition condition: to
        % enforce continuity, that value should be the solution of the
        % problem at that point. If no solution can be found, output=0;
        function N=getStartingCount(obj,t,x)
            N=obj.quickSol(t,x);
            if(isnan(N))
                N=0;
                if(~isempty(obj.valConditions))
                    disp('warning: could not determine the correct offset for the value condition.');
                end
            end
        end
      
        %explicit Solution:
        %Calculates the count function and the active value condition
        %at points (t,x) where t and x may be matrices.
        function res = explSol(obj,t,x)
            if(size(t) ~= size(x))
                error 'space and time matrices have different dimensions';
            end
            N=NaN * ones(size(t));
            i=0;
            index = NaN*ones(size(t));
            for myCond = obj.valConditions'
                i=i+1;
                tmp = LH_partSol(obj.funddiag,t,x,myCond);
                isLower = tmp == min(tmp,N); %Equivalent to tmp <= N, but weird condition to avoid problems when N = NaN
                %Modified specifically for the network control===========
                if (sum(sum(isLower,1)~=0)>2)   
                index(isLower) = i;
                end
                %=======================================================
                N(isLower) = tmp(isLower);
            end
            res = {N, index};
        end

        function k = density(obj,t,x,index)
            %Calculates the density and the count function at points
            %(t,x) on val Conditions #index. Needs a fundamental diagram with the Rprime method.
            k=NaN * ones(size(t));
            if(size(t) ~= size(x))
                error 'space and time matrices have different dimensions';
            end
            if(size(index) ~= size(t))
                error 'wrong index matrix';
            end
            i=0;
            for myCond = obj.valConditions'
                i=i+1;
                k(index==i) = LH_partDens(obj.funddiag,t(index==i),x(index==i),myCond);
            end


            
        end

        function N = quickSol(obj,t,x)
            %Quick solution: you get only the count function.
            N=NaN*ones(size(t));
            for myCond = obj.valConditions'
                N = min(N,LH_partSol(obj.funddiag,t,x,myCond));
            end
        end

        
        %From now on, old functions:
        function addFirstIntCond(obj,t0,t1,x0,x1,g)
            obj.valConditions = [t0 t1 x0 x1 g*(t1-t0) 0];
        end

        function addIntCond(obj,t0,t1,x0,x1,g)
            tmp = [t0 t1 x0 x1 g*(t1-t0) obj.getStartingCount(t0,x0)];
            obj.valConditions = [obj.valConditions;tmp];
        end
        
        function addIntCondL(obj,t0,t1,x0,x1,g,L)
            tmp = [t0 t1 x0 x1 g*(t1-t0) L];
            obj.valConditions = [obj.valConditions;tmp];
        end

        function addFirstIniCond(obj,x0,x1,k)
            %addIniCond adds the first initial condition (offset 0)
            tmp = [0 0 x0 x1 -k*(x1-x0) 0];
            obj.valConditions = tmp;
        end

        function addIniCond(obj,x0,x1,k)
            %addIniCond adds an initial condition and computes its offset.
            tmp = [0 0 x0 x1 -k*(x1-x0) obj.quickSol(0,x0)];
            obj.valConditions = [obj.valConditions;tmp];
        end

    end
end