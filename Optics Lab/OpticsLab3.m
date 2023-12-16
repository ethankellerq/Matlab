%%                                                initialize plot settings
set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
width = 6;                                              % Width in inches
height = 4;                                             % Height in inches
alw = 1.5;                                              % AxesLineWidth 
fsz = 14;                                               % Fontsize 
lw = 1.5;                                               % LineWidth 
msz = 8;                                                % MarkerSize 
set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);   
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw);
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition',[defpos(1) defpos(2) width*100,height*100]); 
set(0,'defaultFigurePosition',[400, 50, width*100, height*110]); 


% 
% %n=length(PositionCh12M1);
% posVec=.06:.001:.14;
% m=length(posVec);
% avgVec=zeros(1,m);
% %sum=0;
% % for index1 = 1:n-1
% % 
% % if PositionCh12M1(index1,1) == PositionCh12M1(index1+1,1)
% % sum=sum+LightIntensityChAMax1(index1,1);
% % else
% % sumVec(1,index1)=sum;    
% % end
% % 
% % end
% 
% bound1=find(abs(PositionCh12M1)==.06,1);
% bound2=find(abs(PositionCh12M1)==.14,1);
% PositionNew=abs(PositionCh12M1(bound1:bound2,1)).';
% IntensityNew=LightIntensityChAMax1(bound1:bound2,1).';
% NormInt=(IntensityNew-min(IntensityNew))./(max(IntensityNew)-min(IntensityNew));
% 
% for index2=1:m
% y0=0;
% x0=find(abs(PositionNew-posVec(1,index2))<.0001);
% k=length(x0);
% x0min=min(x0);
% x0max=max(x0);
% avgVec(1,index2)=sum(NormInt(1,x0min:x0max))/k;
% end


%% Raw Data
%{
figure(1)
hold on
scatter(Positionsingle_red_160um,single_red_160um,'red');
scatter(Positionsingle_red_80um,single_red_80um,'green');
scatter(Positionsingle_red_40um,single_red_40um,'blue');
xlabel("Position (m)");
ylabel("Intensity Reading");

figure(2)
hold on
scatter(Positionsingle_green_160um,single_green_160um,'red');
scatter(Positionsingle_green_80um,single_green_80um,'green');
scatter(Positionsingle_green_40um,single_green_40um,'blue');
xlabel("Position (m)");
ylabel("Intensity Reading");

figure(3)
hold on
scatter(Positiondouble_red_40um_d125_lightsoff,double_red_40um_d125_lightsoff,'red');
scatter(Positiondouble_red_40um_d250,double_red_40um_d250,'green');
scatter(Positiondouble_red_40um_d500,double_red_40um_d500,'blue');
xlabel("Position (m)");
ylabel("Intensity Reading");

figure(4)
hold on
scatter(Positiondouble_green_40um_d125,double_green_40um_d125,'red');
scatter(Positiondouble_green_40um_d250,double_green_40um_d250,'green');
scatter(Positiondouble_green_40um_d500,double_green_40um_d500,'blue');
xlabel("Position (m)");
ylabel("Intensity Reading");
%}

% figure(2)
% plot(posVec,avgVec);
% 
% figure(3)
% scatter(PositionNew,IntensityNew);
% 
% figure(4)
% scatter(PositionNew,NormInt);

%% Processed Data

NormInt1=(single_red_160um-min(single_red_160um))./(max(single_red_160um)-min(single_red_160um));
NormInt2=(single_red_80um-min(single_red_80um))./(max(single_red_80um)-min(single_red_80um));
NormInt3=(single_red_40um-min(single_red_40um))./(max(single_red_40um)-min(single_red_40um));

NormInt4=(single_green_160um-min(single_green_160um))./(max(single_green_160um)-min(single_green_160um));
NormInt5=(single_green_80um-min(single_green_80um))./(max(single_green_80um)-min(single_green_80um));
NormInt6=(single_green_40um-min(single_green_40um))./(max(single_green_40um)-min(single_green_40um));

NormInt7=(double_red_40um_d125_lightsoff-min(double_red_40um_d125_lightsoff))./(max(double_red_40um_d125_lightsoff)-min(double_red_40um_d125_lightsoff));
NormInt8=(double_red_40um_d250-min(double_red_40um_d250))./(max(double_red_40um_d250)-min(double_red_40um_d250));
NormInt9=(double_red_40um_d500-min(double_red_40um_d500))./(max(double_red_40um_d500)-min(double_red_40um_d500));

NormInt10=(double_green_40um_d125-min(double_green_40um_d125))./(max(double_green_40um_d125)-min(double_green_40um_d125));
NormInt11=(double_green_40um_d250-min(double_green_40um_d250))./(max(double_green_40um_d250)-min(double_green_40um_d250));
NormInt12=(double_green_40um_d500-min(double_green_40um_d500))./(max(double_green_40um_d500)-min(double_green_40um_d500));

% figure(5)
% hold on
% scatter(Positionsingle_red_160um,NormInt1,'red');
% scatter(Positionsingle_red_80um,NormInt2,'green');
% scatter(Positionsingle_red_40um,NormInt3,'blue');
% xlabel("Position (m)");
% ylabel("Intensity Reading");
% 
% figure(6)
% hold on
% scatter(Positionsingle_green_160um,NormInt4,'red');
% scatter(Positionsingle_green_80um,NormInt5,'green');
% scatter(Positionsingle_green_40um,NormInt6,'blue');
% xlabel("Position (m)");
% ylabel("Intensity Reading");

% figure(7)
% hold on
% scatter(Positiondouble_red_40um_d125_lightsoff,NormInt7,'red');
% scatter(Positiondouble_red_40um_d250,NormInt8,'green');
% scatter(Positiondouble_red_40um_d500,NormInt9,'blue');
% xlabel("Position (m)");
% ylabel("Intensity Reading");
% 
% figure(8)
% hold on
% scatter(Positiondouble_green_40um_d125,NormInt10,'red');
% scatter(Positiondouble_green_40um_d250,NormInt11,'green');
% scatter(Positiondouble_green_40um_d500,NormInt12,'blue');
% xlabel("Position (m)");
% ylabel("Intensity Reading");


% figure(9)
% hold on
% plot(Positionsingle_red_160um,NormInt1,'red');
% plot(Positionsingle_red_80um,NormInt2,'green');
% plot(Positionsingle_red_40um,NormInt3,'blue');
% xlabel("Position (m)");
% ylabel("Intensity Reading");

slitsize=[40,80,160];
separation=[125,250,500];

single_red_w=[0.0095
0.0054
0.0031
];
single_green_w=[0.0086
0.0045
0.003
];

double_red_w=[0.00155
0.0009
0.0004
];

double_green_w=[0.0013
0.0007
0.00023
];

slitsize_array=40:1:160;
y1=(-4.9821e-05).*slitsize_array +.0107;
y2=(-4.26785714285714e-05).*slitsize_array + 0.00935000000000000;

separation_array=120:1:500;
y3=-2.91428571428571e-06.*separation_array+0.00180000000000000;
y4=-2.71428571428571e-06.*separation_array+0.00153500000000000;

% figure(10)
% hold on
% scatter(slitsize,single_red_w,'red');
% scatter(slitsize,single_green_w,'green');
% plot(slitsize_array,y1,'red');
% plot(slitsize_array,y2,'green');
% xlabel("Slit Size (μm)");
% ylabel("0th Peak Width (m)");

figure(11)
hold on
scatter(separation,double_red_w,'red');
scatter(separation,double_green_w,'green');
plot(separation_array,y3,'red');
plot(separation_array,y4,'green');
xlabel("Separation (μm)");
ylabel("0th Peak Width (m)");
