
warning ("off", "Octave:broadcast");

function tfOut=rlTest(gain)
    %ommitting any check for the sake of speed 

    %reg=global regG;
    %sys=sysG;
    %H=HG;
    %F=FG;

  reg=tf([1 1],[2 1]);   
  sys=tf([2 1],[1 1 2]);
  F=1;
  H=1;

    G=reg*sys
    tfout=F*((G*gain)/(1+G*H*gain))
    %tfout(1,:)=F*((G*gain(1,:))/(1+G*H*gain(1,:)));
end


function [realPart, imagPart]=rlPlot(OLNum, OLDen, H, F,startGain,endGain, stepGain)
    %%TODO: handle this as matrix operations
tic


%sysNum=cell2mat(sys.num);
%sysDen=cell2mat(sys.den);
%regNum=cell2mat(reg.num);
%regDen=cell2mat(reg.den);

%we calculate the total transfer of the regulator together with the regulated system (open loop)
%a=conv(regDen, sysDen);
%b=conv(regNum, sysNum);


a=OLDen;
b=OLNum;


%completeSys=sys*reg;
%a=cell2mat(completeSys.den)
%b=cell2mat(completeSys.num)

%this will not change unless the settings are altered; vector with the steps for the root locus plot 
r=linspace(startGain,endGain,stepGain);
%r=0:1:20;
r=r';


%these will not be changed unless a zero or a pole is added to the regulator
aS=size(a,2)
bS=size(b,2)
rS=size(r,1)




%could be done better I think
aBig=zeros(rS,aS)+a;
bBig=zeros(rS, bS)+b;

if(bS>aS)
    aBig=[zeros(bS-aS,rS) aBig];
end

if(aS>bS)
    bBig=[zeros(rS,aS-bS) bBig];
end

%we calculate the charasteristic polynomials of the CL-transfers with all the possible gains here (in a matrix form)
denPolys=getFeedbackDen(aBig,bBig,r);

hold on;

foundRoots=rootsMatrix(denPolys);
foundRoots;

%return;

%%Some points of the root locus will be too far from oneanother; we will try to maintain a constant distance by finding those point couples and then filling the gaps

%maxDistance=0.5;%maximum distance between 2 points in the chart allowed
%maxDistanceVect=zeros(1,size(foundRoots,2));
%maxDistanceVect(:)=maxDistance;
%r

rScale=1;
iScale=4;

    distanceMat=(foundRoots(1:end-1,:)-foundRoots(2:end,:));%we subtract the matrix of all roots from itself with one row offset. This basically means we get the distances between each 2 consecutive roots
    distanceMatR=abs(real(distanceMat));%then we split the results into real and imaginary parts
    distanceMatI=abs(imag(distanceMat));
    avgDistanceR=mean(max(distanceMatR,[],2));
    avgDistanceI=mean(max(distanceMatI,[],2));
    avgDistance=avgDistanceR;%should be*2, but it works;
    iScale=avgDistanceR/avgDistanceI;
    maxDistance=avgDistance*1.5;

    distanceMat=distanceMatR.*rScale+distanceMatI.*iScale;%Weighted sum to give x-distance a different weight from y-distance
    foundDistanceMaxVect=max(distanceMat,[],2);%ǧet a column vector of maximum distances for each row (that is, find the biggest one from the distances  the poles have traveled between two values of r)
    maxDistanceVect=repmat(maxDistance,size(foundRoots,1)-1,1);%we just make a row vector and fill it with the maximum traveled distance allowed
    isDistanceBadVect=foundDistanceMaxVect>maxDistanceVect;%then we compare the worst traveled distance and see if it is too much
    %size(foundRoots)

for it2=0:2%%improvement iterations
    it2


    distanceMat=(foundRoots(1:end-1,:)-foundRoots(2:end,:));%we subtract the matrix of all roots from itself with one row offset. This basically means we get the distances between each 2 consecutive roots
    distanceMatR=abs(real(distanceMat));%then we split the results into real and imaginary parts
    distanceMatI=abs(imag(distanceMat));

    distanceMat=distanceMatR.*rScale+distanceMatI.*iScale;%Weighted sum to give x-distance a different weight from y-distance
    foundDistanceMaxVect=max(distanceMat,[],2);%ǧet a column vector of maximum distances for each row (that is, find the biggest one from the distances  the poles have traveled between two values of r)
    maxDistanceVect=repmat(maxDistance,size(foundRoots,1)-1,1);%we just make a row vector and fill it with the maximum traveled distance allowed
    isDistanceBadVect=foundDistanceMaxVect>maxDistanceVect;%then we compare the worst traveled distance and see if it is too much
    %size(foundRoots)

it=2;
it_bv=1;%iterator over the isDistanceBadVect
while it<=size(foundRoots,1)
    %size(foundRoots,1);
    %rootVecCurrent=foundRoots(it,:);
    %rootVecPrev=foundRoots(it-1,:);
    %distanceVect=(rootVecCurrent-rootVecPrev);
    %distanceVectI=imag(distanceVect);
    %distanceVectR=real(distanceVect);

    %distanceVect=distanceVectR*rScale+distanceVectI*iScale;%using a manhattan distance for fun (and speed)


    %isDistanceBadVect=(distanceVect)>maxDistanceVect;
    %isDistanceBad=sum(isDistanceBadVect);%when we sum it, it is only zero if all elements were zero, so it's like OR-ing it


    if(isDistanceBadVect(it_bv)!=0)%distance is too big!
        maxMeasuredDistance=foundDistanceMaxVect(it_bv);%max(distanceVect);
        %123456
        refiningStepsCount=ceil(maxMeasuredDistance/maxDistance)+1;%%+1 to be sure
        kCurrent=r(it);
        kPrev=r(it-1);
        kDiff=abs(kCurrent-kPrev);
        stepSize=kDiff/refiningStepsCount;
        kFine=(kPrev+stepSize:stepSize:kCurrent-stepSize)';
        if(size(kFine,1)==0)
            kFine=(kCurrent+kPrev)/2;
        end
        a;
        b;
        [aBig bBig]=getaBigbBig(a,b,size(kFine,1));

        denPolys=getFeedbackDen(aBig,bBig,kFine);
        currentRoots=rootsMatrix(denPolys);
        currentRoots;
        kFine;
        foundRoots=[foundRoots(1:it-1,:);currentRoots;foundRoots(it:end,:)];
        r=[r(1:it-1,:);kFine;r(it:end,:)];
        it+=size(currentRoots,1);
    end
    it++;
    it_bv++;
end

end

hold off;
%return;

    %for i=0:0.1:10
        %theTf=rlTest(sys, reg, H, F,i);
        %theRoots=roots(cell2mat(theTf.den));
        %h=plot(theRoots);
        %set (h, "linestyle", "none"); 
        %set (h, "marker", "*"); 
    %end
    %hold off;
    realPart=real(foundRoots);
    imagPart=imag(foundRoots);
    toc
end

function foundRoots=rootsMatrix(polyMatrix)

foundRoots=zeros(size(polyMatrix,1),size(polyMatrix,2)-1);

%%here we find all the roots of the characteristic polynomials for all the gains
for it=1:size(polyMatrix,1)
    it;
    %(roots(denPolys(i,:)))
    foundRoots(it,:)=roots(polyMatrix(it,:))+0.001*i;%0.001i=ugly hack to fix the plotting
    %h=plot(foundRoots);
    %set (h, "linestyle", "none"); 
    %set (h, "marker", "*"); 
end
end

%this function calculates the characteristic polynomials of the closed loop transfer, described by the denominator of the total OL transfer function a, the nominator b and gain k. Both a and b are either row vectors or
%matrices. In the latter case, each row represents one transfer function num/den. The parameter k is either a scalar or a column vector.
function den=getFeedbackDen(a, b, k)
%whoopsie, I kinda forgot to implement other structures than reg->sys in feedback loop
den=a+b.*k;
end

%following function takes a vector aj
function [aBig bBig]=getaBigbBig(a, b, rowCount)
aS=size(a,2);
bS=size(b,2);
rowCount;



aBig=zeros(rowCount,aS)+a;
bBig=zeros(rowCount, bS)+b;

if(bS>aS)
    aBig=[zeros(bS-aS,rowCount) aBig];
end

if(aS>bS)
    bBig=[zeros(rowCount,aS-bS) bBig];
end

end
