%% Code to create graphics for Lab

%
lambda1=w4;
trans1=GrellowGivenTransSpectrum10I1000S0B;
[n,~]=size(trans1);
for i=1:n
    if trans1(n,1) < 0
        trans1(n,1)=0;
    end
end
translog1=log10(trans1);

%

[leftb1,~]=find(lambda1==403.85);
[rightb1,~]=find(lambda1==440.16);
[leftb2,~]=find(lambda1==473.52);
[rightb2,~]=find(lambda1==496.68);

lambdaR1=lambda1(leftb1:rightb1,1);     % set up new subvectors
lambdaR2=lambda1(leftb2:rightb2,1);
translogR1=translog1(leftb1:rightb1,1);
translogR2=translog1(leftb2:rightb2,1);

[P1,S1]=polyfit(lambdaR1,translogR1,1);   % this tells you slope and intercept 
[P2,S2]=polyfit(lambdaR2,translogR2,1);   % in that order
[yfit1,delta1]=polyval(P1,lambda1,S1);          
[yfit2,delta2]=polyval(P2,lambda1,S2);          
%}

figure(1)
scatter(lambda1,trans1,20);
ylim([0,100]);
%xlim([450,900]);
title("Transmission vs Wavelength for Green Given Color Filter");
ylabel("Transmission (%)");
xlabel("Wavelength (nm)");
%xline(710);
%xline(760);


figure(2)
%{
set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
width = 6;                                                                 % Width in inches
height = 4;                                                                % Height in inches
alw = 1.5;                                                                 % AxesLineWidth 
fsz = 14;                                                                  % Fontsize 
lw = 1.5;                                                                  % LineWidth 
msz = 8;                                                                   % MarkerSize 
set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);   
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw);
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
set(0,'defaultFigurePosition', [400, 50, width*100, height*110]); 
%}
scatter(lambda1,translog1,20);
title("log of Transmission vs Wavelength for Green Given Color Filter");
ylabel("log(Transmission)");
xlabel("Wavelength (nm)");
% xline(350);
% xline(500);
% xline(570);
% xline(600);
xlim([350,700]);
ylim([1.5,2]);

%   

hold on
%plotting linear fits for spectral cutoffs
plot(lambda1,yfit1,lambda1,yfit2,'LineWidth',1.5);
%plotting 95% prediction interval for both spectral cutoffs
plot(lambda1,yfit1+2*delta1,'m--',lambda1,yfit1-2*delta1,'m--');
plot(lambda1,yfit2+2*delta2,'m--',lambda1,yfit2-2*delta2,'m--');

%legend('Data','Linear Fit 1','Linear Fit 2','95% Prediction Interval')

%}
%}

%% Plotting Only

% plot(w1,IncandSourceSpectrum10I100S0BAt40,w1,IncandSourceSpectrum10I100S0BAt50,w1,IncandSourceSpectrum10I100S0BAt55);
% title("Intensity vs Wavelength for Incandescent Spectrum");
% ylabel("Intensity (#)");
% xlabel("Wavelength (nm)");

% plot(wavelength,DarkSpectrum10I100S0B,wavelength,DarkSpectrum10I100S1B,wavelength,DarkSpectrum10I100S10B,wavelength,DarkSpectrum10I100S100B);
% title("Intensity vs Wavelength for Dark Spectrum");
% ylabel("Intensity (#)");
% xlabel("Wavelength (nm)");

% plot(wavelength,DarkSpectrum10I1S0B,wavelength,DarkSpectrum100I1S0B,wavelength,DarkSpectrum1000I1S0B);
% title("Intensity vs Wavelength for Dark Spectrum");
% ylabel("Intensity (#)");
% xlabel("Wavelength (nm)");

% plot(wavelength,DarkSpectrum10I100S0B,wavelength,DarkSpectrum10I100S0B_PointedComputer,wavelength,DarkSpectrum10I100S0B_PointedDoor);
% title("Intensity vs Wavelength for Dark Spectrum");
% ylabel("Intensity (#)");
% xlabel("Wavelength (nm)");

% plot(wavelength,DarkSpectrum10I100S0B,wavelength,DarkSpectrum10I100S0B_PointedComputer,wavelength,DarkSpectrum10I100S0B_PointedDoor);
% title("Intensity vs Wavelength for Dark Spectrum");
% ylabel("Intensity (#)");
% xlabel("Wavelength (nm)");
%{
hydrogenminusdark=hydrogensource-Darkspec;
figure(10);
scatter(wavelength,hydrogenminusdark);
title("Intensity vs Wavelength for Hydrogen Spectrum");
ylabel("Intensity (#)");
xlabel("Wavelength (nm)");

[leftb1,~]=find(wavelength==653.45);       % decide on boundary for fittings                                                                       
[rightb1,~]=find(wavelength==659.38); 
lambdaR1=wavelength(leftb1:rightb1,1);
intR1=hydrogenminusdark(leftb1:rightb1,1);
Gfit=fit(lambdaR1,intR1,'gauss1');


n=10001;
x=linspace(653.45,659.38,n);
dx=x(2)-x(1);
g1=(4.769e+04)*exp(-((x-656.5)/0.8794).^2);
CA=zeros(1,n);                          % Preallocation
%CB=zeros(1,n);
MLA=zeros(1,n);
%MLB=zeros(1,n);


for i = 1:n-1                           % CA and CB are directly from the                            
    CA(i)=(g1(i+1)+g1(i))*dx/2;         % the equation for a trapezoid.
    MLA(i+1)=MLA(i)+CA(i);              % MLA and MLB allow me to sum 
    %CB(i)=(NB(i+1)+NB(i))*dt/2;         % CA and CB over i which produces
    %MLB(i+1)=MLB(i)+CB(i);              % the approximation of the integral
end  


sigma=MLA(end)/((4.848e+4)*(2*pi)^(0.5));
figure(11);
plot(Gfit,lambdaR1,intR1);
title("Intensity vs Wavelength for Hydrogen Spectrum");
ylabel("Intensity (#)");
xlabel("Wavelength (nm)");
%}

%
set(0,'DefaultFigureColor','white')
fig.InvertHardcopy = 'off';
width = 6;                                                                 % Width in inches
height = 4;                                                                % Height in inches
alw = 1.5;                                                                 % AxesLineWidth 
fsz = 14;                                                                  % Fontsize 
lw = 1.5;                                                                  % LineWidth 
msz = 8;                                                                   % MarkerSize 
set(0,'defaultAxesFontSize',fsz); 
set(0,'defaultLineLineWidth',lw);   
set(0,'defaultLineMarkerSize',msz); 
set(0,'defaultAxesLineWidth',alw);
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]); 
set(0,'defaultFigurePosition', [400, 50, width*100, height*110]); 
%}
