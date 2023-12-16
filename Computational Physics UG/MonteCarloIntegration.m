%% Proj 7 Ethan Keller

Rad1=input('Enter radius for circle 1: ');
Rad2=input('Enter radius for circle 2: ');
d=input('Enter seperation distance: ');
L=2*max(Rad1,Rad2);

GraphOverlappingCircles(Rad1,Rad2,d);

disp('Area of object via Monte Carlo Integration: ');
disp('     ');
pError=input('Enter an acceptable % error: ');


area=TotalArea(Rad1,Rad2,d,pError);


dmax=Rad1+Rad2;
% dd=dmax/5;
% for d_vis=dd:dd:dmax
%     GraphOverlappingCircles(Rad1,Rad2,d_vis);
% end

dd=dmax/50;
x_plot=zeros(50,1);
for k=1:50
    d_gr=dd*k;
    x_plot(k)=d_gr;
end
y_plot=zeros(50,25);
figure;
hold on
for j=1:25
    for k=1:50
         d_gr=dd*k;
         y_plot(k,j)=TotalArea(Rad1,Rad2,d_gr,pError);
    end
    plot(x_plot,y_plot(:,j),'.');
end

ymean=mean(y_plot,2);
plot(x_plot,ymean,'r','LineWidth',2);

function [A,nPts] = TotalArea(Rad1,Rad2,d,pError)
[xmin,xmax,ymin,ymax]=GetBoxCoord(Rad1,Rad2,d);

Abox=(xmax-xmin)*(ymax-ymin);
%disp(Abox);

R=min(Rad1,Rad2);
Amin=pi*R*R;
aveAreaPerPoint=pError*Amin/100;
nPts=50*ceil(Abox/aveAreaPerPoint);
%disp(nPts);

rvec=rand(nPts,2);
rvec(:,1)=xmin+(xmax-xmin)*rvec(:,1);
rvec(:,2)=ymin+(ymax-ymin)*rvec(:,2);
%hold on
%plot(rvec(:,1),rvec(:,2),'.');

amtIn=0;
Rad1Sq=Rad1*Rad1;
Rad2Sq=Rad2*Rad2;
    for k=1:nPts
        q1Sq=rvec(k,1)*rvec(k,1)+rvec(k,2)*rvec(k,2);
        x=rvec(k,1)-d;
        q2Sq=x*x+rvec(k,2)*rvec(k,2);
        if q1Sq<Rad1Sq
            amtIn=amtIn+1;
        elseif q2Sq<Rad2Sq
            amtIn=amtIn+1;
        end
        %[q1Sq,q2Sq,Rad1Sq,Rad2Sq,rvec(k,1:end)]
    end
%disp(amtIn);
frac=amtIn/nPts;
A=frac*Abox;
end
     

function [xmin,xmax,ymin,ymax]=GetBoxCoord(Rad1,Rad2,d)
xmin=-Rad1;
xmax=d+Rad2;
ymax=max(Rad1,Rad2);
ymin=-ymax;
end


function GraphOverlappingCircles(Rad1,Rad2,d)
N=200;
t=(2*pi/N)*(1:N);
figure;
hold on
axis equal
fill(Rad1*cos(t),Rad1*sin(t),'k');
fill(d+Rad2*cos(t),Rad2*sin(t),'k');
[xmin,xmax,ymin,ymax]=GetBoxCoord(Rad1,Rad2,d);
x=zeros(5,1);
y=zeros(5,1);

xmin=xmin-.01;
xmax=xmax+.01;
ymin=ymin-.01;
ymax=ymax+.01;

x(1)=xmin;
x(2)=xmin;
x(3)=xmax;
x(4)=xmax;
x(5)=xmin;
y(1)=ymin;
y(2)=ymax;
y(3)=ymax;
y(4)=ymin;
y(5)=ymin;

plot(x,y,'r');

end




