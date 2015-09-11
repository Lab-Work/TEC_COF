function k = LH_partDens(fd,t,x,myCond)
%Computes the density in (t,x) for the fund. diag. fd, under the
%value condition myCond.

    if(size(t)~=size(x))
        error('t and x have different sizes');
    end

    k=NaN*ones(size(t));
    t0=myCond(1);
    t1=myCond(2);
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

    virtVel0 = k;
    virtVel1 = k;
    
    velObserver = (x1-x0) / (t1-t0);
    %Modified to deal with the first element issue
    
    %virtVel0(isAfterBeginning) =...
    %(x(isAfterBeginning)-x0)./(t(isAfterBeginning)-t0);   %original
    virtVel0(isAfterBeginning) = (x(isAfterBeginning)-x0+sign(sum(x(isAfterBeginning)-x0))*10e-15)./(t(isAfterBeginning)-t0);
    virtVel1(isAfterEnd) = (x(isAfterEnd)-x1)./(t(isAfterEnd)-t1);
    
    %Modified to solve the boundary problem====================
    isBegPoint = abs(t-t0)<10e-5 & abs(x-x0)<10e-5;
    virtVel0(isBegPoint) = 10e-15*sign(sum(x(isAfterBeginning)-x0));
    isEndPoint = abs(t-t1)<10e-5 & abs(x-x1)<10e-5;
    virtVel0(isEndPoint) = 10e-15*sign(sum(x(isAfterBeginning)-x0));
    %==========================================================
    
    isInActionZone = (min(virtVel0,virtVel1) <= fd.wspeed(0) ...
        & max(virtVel0,virtVel1) >= fd.wspeed(fd.kappa));
    
    isInEndDomain = virtVel1 <= fd.wspeed(k1) & virtVel1 >= fd.wspeed(k2);
    k(isInEndDomain) = fd.density(virtVel1(isInEndDomain));

    isInCharacDomain = virtVel0 >= fd.wspeed(k2) ...
        & (virtVel0 <= fd.wspeed(k1) | (virtVel0 <= velObserver & virtVel1 < fd.wspeed(0)))... %the second case is for initial boundary conditions
        & ~isInEndDomain;
    
    %isInUpCharacDomain = virtVel0 >= velObserver & isInCharacDomain;
    isInUpCharacDomain = virtVel0 > velObserver & isInCharacDomain;     %original
    k(isInUpCharacDomain) = k1; 
    
    %isInDownCharacDomain = virtVel0 < velObserver & isInCharacDomain;  
    isInDownCharacDomain = virtVel0 <= velObserver & isInCharacDomain;  %original
    k(isInDownCharacDomain) = k2;  %original
    %k(isInDownCharacDomain) = k1;
    
    isInBegDomain = isInActionZone & ~isInEndDomain & ~isInCharacDomain;
    k(isInBegDomain) = fd.density(virtVel0(isInBegDomain));
    
%Modified to solve the boundary problem====================
%     isBegPoint = abs(t-t0)<10e-5 & abs(x-x0)<10e-5;
%     k(isBegPoint) = k1;
%     %k(isBegPoint) = k2;         %original
%     
%     isEndPoint = abs(t-t1)<10e-5 & abs(x-x1)<10e-5;
%     k(isEndPoint)=k2;  %original
%     %k(isBegPoint) = k1;
%==========================================================
end
