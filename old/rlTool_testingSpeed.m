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


function tfs=rlPlot(sys, reg, H, F)
%%TODO: handle this as matrix operations
tic


sysNum=cell2mat(sys.num);
sysDen=cell2mat(sys.den);
regNum=cell2mat(reg.num);
regDen=cell2mat(reg.den);

a=conv(regDen, sysDen);
b=conv(regNum, sysNum);


%completeSys=sys*reg;
%a=cell2mat(completeSys.den)
%b=cell2mat(completeSys.num)

%this will not change unless the settings are altered; vector with the steps for the root locus plot (TODO: make the step size dynamic)
r=0:0.5:10;
r=r';


%these will not be changed unless a zero or a pole is added to the regulator
aS=size(a,2)
bS=size(b,2)
rS=size(r,1)


    for it2=0:0.1:10
        regNum=regNum+[0 0.1];
        a=conv(regDen, sysDen);
        b=conv(regNum, sysNum);

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

        %%here we find all the roots of the characteristic polynomials for all the gains
        for it=1:size(r,1)
            it;
            (roots(denPolys(it,:)));
            %h=plot(roots(denPolys(it,:))+0.001*i);%0.001i=ugly hack to fix the plotting
            %set (h, "linestyle", "none"); 
            %set (h, "marker", "*"); 
        end
        hold off;

        %for i=0:0.1:10
        %theTf=rlTest(sys, reg, H, F,i);
        %theRoots=roots(cell2mat(theTf.den));
        %h=plot(theRoots);
        %set (h, "linestyle", "none"); 
        %set (h, "marker", "*"); 
        %end
        %hold off;
    end
        toc
end

function den=getFeedbackDen(a, b, k)
    den=a+b.*k;
end
