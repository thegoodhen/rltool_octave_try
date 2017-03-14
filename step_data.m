function pkv=peakVal(y,~,~)
    pkv=max(y);
end


function ev=endVal(y,~,~)
    ev=y(end);
end

function overshoot=OS(y,~,~)
    overshoot=((peakVal(y,0,0)-endVal(y,0,0))/(endVal(y,0,0)))*100;
end

function rT=raiseTime(y,t,~,from,to)
    endValue=endVal(y,t,0);
    fromValue=from*endValue;
    toValue=to*endValue;
    greaterThanFrom=y>fromValue;
    greaterThanTo=y>toValue;
    raiseStart=find(greaterThanFrom,1);
    raiseEnd=find(greaterThanTo,1);
    rT=t(raiseEnd)-t(raiseStart);
end

function sT=settlingTime(y,t,~,range)
    eV=endVal(y,t,0);
    tol=abs(eV*range);
    yFlipped=flipud(y);
    badVec=abs(yFlipped-eV)>tol;
    firstBadIndex=find(badVec,1);
    sT=t(size(t,1)-firstBadIndex);
    
end
