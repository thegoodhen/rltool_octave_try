function writeTest()
%A=eye(100);
x=0:0.1:2*pi;
y=sin(x)';
inc=0.01;
tic;
m=0;

while(1)
if(m>1 || m<0)
    inc=-inc;
end;
    m+=inc;
    y2=y.*m;
    save("/home/thegoodhen/sketchbook/floatTest/testFifo","-ascii", "y2");
    toc;
end;
end;
