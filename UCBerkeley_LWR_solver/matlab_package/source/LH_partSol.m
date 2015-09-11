function N = LH_partSol(fd,t,x,myCond)
%Computes a solution component in (t,x) for the fund. diag. fd, under the
%value condition myCond. TODO: modify so that t and x can be (1,n) vectors.

    if(size(t)~=size(x))
        error('t and x have different sizes');
    end

    N=NaN*ones(size(t));
    t0=myCond(1);
    t1=myCond(2);%We assume t0 <= t1
    x0=myCond(3);
    x1=myCond(4);
    m =myCond(5);
    n0=myCond(6);
    
    if(t1==t0)
        k1 = 0;
        k2 = -m/(x1-x0);
    else
        dens = fd.densities((x1-x0)/(t1-t0),m/(t1-t0));
        k1=dens(1);
        k2=dens(2);
    end
    
    isAfterBeginning = t>= t0;
    isAfterEnd = t >= t1;

    virtVel0 = N;
    virtVel1 = N;
    
    velObserver = (x1-x0) / (t1-t0);
    
    virtVel0(isAfterBeginning) =...
    (x(isAfterBeginning)-x0)./(t(isAfterBeginning)-t0);  % original
    %virtVel0(isAfterBeginning) =...
    %(x(isAfterBeginning)-x0+sign(sum(x(isAfterBeginning)-x0))*10e-5)./(t(isAfterBeginning)-t0);
    virtVel1(isAfterEnd) = (x(isAfterEnd)-x1)./(t(isAfterEnd)-t1);
    
    isInActionZone = (min(virtVel0,virtVel1) <= fd.wspeed(0) ...
        & max(virtVel0,virtVel1) >= fd.wspeed(fd.kappa));
    
    isInEndDomain = virtVel1 <= fd.wspeed(k1) & virtVel1 >= fd.wspeed(k2);
    N(isInEndDomain) = n0 + m ...
        + (t(isInEndDomain)-t1).*fd.R(virtVel1(isInEndDomain));

    isInCharacDomain = virtVel0 >= fd.wspeed(k2) ...
        & (virtVel0 <= fd.wspeed(k1) | (virtVel0 <= velObserver & virtVel1 < fd.wspeed(0)))... %the second case is for initial boundary conditions
        & ~isInEndDomain;
    
    isInUpCharacDomain = virtVel0 > velObserver & isInCharacDomain;
    N(isInUpCharacDomain) = n0 + (x0-x(isInUpCharacDomain))*k1 ...
        + (t(isInUpCharacDomain)-t0)*fd.flow(k1);
    
    isInDownCharacDomain = virtVel0 <= velObserver & isInCharacDomain;
    N(isInDownCharacDomain) = n0 + (x0-x(isInDownCharacDomain))*k2 ...
        + (t(isInDownCharacDomain)-t0)*fd.flow(k2);
    
    isInBegDomain = isInActionZone & ~isInEndDomain & ~isInCharacDomain;
    N(isInBegDomain) = n0 + ...
        (t(isInBegDomain)-t0).*fd.R(virtVel0(isInBegDomain));
    
    isBegPoint = t==t0 & x==x0;   %3 m
    N(isBegPoint) = n0;
    
    
    %Modified to avoid numerical problem===================================
    isEndPoint = t==t1 & abs(x-x1)<0.1;   %3 m
    N(isEndPoint)=n0+m;
    %======================================================================
    
end
