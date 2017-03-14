source rlTool_adaptiveStep_vectorised.m;%also kinda ugly, this file just contains function definitions
pkg load control;

function commsTest()
sys=tf([2 1],[1 1 2]);
reg=tf([1 1],[2 1]);
%A=eye(100);
x=0:0.1:2*pi;
y=sin(x)';
inc=0.01;
%tic;
m=0;

initTxt=0;

startVal=0;
endVal=10;
stepSize=1;



sysNum=cell2mat(sys.num);
sysDen=cell2mat(sys.den);
regNum=cell2mat(reg.num);
regDen=cell2mat(reg.den);

theGain=1;%the current set gain of the root locus

%we calculate the total transfer of the regulator together with the regulated system (open loop)
OLDen=conv(regDen, sysDen);
OLNum=conv(regNum, sysNum);
OLNumRespectingGain=conv(regNum, sysNum)*theGain;


regPoles=roots(regDen);
regZeros=roots(regNum);

    %save("/home/thegoodhen/workspace/pipeTest/OOJI","-ascii", "initTxt");


while(1)
    disp('loop begin');
if(m>1 || m<0)
    inc=-inc;
end;
    m+=inc;
    y2=y.*m;
    fname="/home/thegoodhen/sketchbook/floatTest/JOOI";
    %[theText]=textread(fname,"%s");
    disp('Gonna read now! :3');
    [output theText]=system("cat /home/thegoodhen/workspace/pipeTest/JOOI")%temporary workaround
    %kokodak=load(fname,"-ascii")
    %disp(theText);
    if(strcmp(theText,"P\n"))%Root locus-plot
        disp("root locus plot request")

        sendAck();
        paramVec=loadVectorOfDoubles()
        startGain=paramVec(1);
        endGain=paramVec(2);
        stepCountGain=paramVec(3);
        [rP iP]=rlPlot(OLNum, OLDen, 1,1,startGain,endGain,stepCountGain);
        %disp("kokodak");
        [rowCount columnCount]=size(rP);
        disp('Gonna write now');

        sendPlotMatrix(rowCount,columnCount, rP, iP);
    disp('written, wee!');
    end;


    if(strcmp(theText, "G\n"))
        disp("gain change request")
        sendAck();
        paramVec=loadVectorOfDoubles()
        startGain=paramVec(1);
        theGain=startGain;
        endGain=startGain+1;
        stepCountGain=2;
        OLNumRespectingGain=conv(regNum, sysNum)*theGain;
        [rP iP]=rlPlot(OLNum, OLDen, 1,1,startGain,endGain,stepCountGain);
        [rowCount columnCount]=size(rP);
        disp('Gonna write now');
        sendPlotMatrix(rowCount,columnCount, rP, iP);
        disp('written, wee!');
    end;


    if(strcmp(theText,"PO\n"))%Pole change request
        disp('Pole change req');

        %regPoles(1)+=0.1;%TODO
        sendAck();
           
        regPoles=loadVectorOfDoubles();

        OLDen=recalcZP(sysDen,regPoles);
        sendAck();
    end;

    if(strcmp(theText,"ZE\n"))%Zero change request
        disp('Zero change request');
        sendAck();
        regZeroes=loadVectorOfDoubles();
        OLNum=recalcZP(sysNum,regZeroes);
        sendAck();
    end;

    if(strcmp(theText,"S\n"))%Step response request
        %OLNumRespectingGain=conv(regNum, sysNum)*theGain;
        bq=OLNumRespectingGain;
        ap=OLDen;


        
        aS=size(ap,2)
        bS=size(bq,2)

        %TODO: refactor this into a function
        aBig=zeros(1,aS)+ap;
        bBig=zeros(1, bS)+bq;

        if(bS>aS)
            aBig=[zeros(bS-aS,1) aBig];
        end

        if(aS>bS)
            bBig=[zeros(1,aS-bS) bBig];
        end



        sys=tf(bBig,aBig+bBig);%this might be kind of computational intensive, therefore: TODO: see if it can be simplified, test how much time this takes (about 0.001s, not as bad as I thought);
        sendAck();
        paramVec=loadVectorOfDoubles()
        startTime=paramVec(1);
        endTime=paramVec(2);
        stepCountTime=paramVec(3);
        [y t x]=step(sys);
        tP=t;
        yP=y;
        %[xP yP]=[x y];
        [rowCount columnCount]=size(tP);
        disp('Gonna write now');
        sendPlotMatrix(rowCount,columnCount, tP, yP);
        disp('written, wee!');
    end

    %save("/home/thegoodhen/sketchbook/floatTest/testFifo","-ascii", "y2");
    %toc;
end;
end;

function sendPlotMatrix(rowCount, columnCount,rP,iP)


    rowCountVec=zeros(1,size(rP,2));
    rowCountVec(1)=rowCount;
    columnCountVec=zeros(1,size(rP,2));
    columnCountVec(1)=columnCount;


    saveMatrix=[rowCountVec;columnCountVec;rP;iP];

    save("/home/thegoodhen/workspace/pipeTest/OOJI","-ascii", "saveMatrix");
end

function vec=loadVectorOfDoubles()
        [output theText]=system("cat /home/thegoodhen/workspace/pipeTest/JOOI")%now for the actual data
        indivStrings=strsplit(deblank(theText));
        vec=str2double(indivStrings)
end

function sendAck()
        ACK=1;%some reply
        disp('sending ack');
        save("/home/thegoodhen/workspace/pipeTest/OOJI","-ascii", "ACK");
        disp('sent');
end

function OLDen=recalcZP(sysDen, regPoles)%multiplies the ponylomial given as the first argument with a set of linear ones, described by their roots given as a vector in the 2nd argument
    OLDen=sysDen;
    size(regPoles)
    for i=1:size(regPoles,2)
        regPoles(i);
        OLDen=conv(OLDen,[1 -regPoles(i)]);
    end;
end;
