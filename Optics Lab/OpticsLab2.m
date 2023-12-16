unblockgreenint=.48;
unblockredint=.83;
greenlasergreenfilterint=.462;
greenlaserbluefilterint=.125;
redlaserredfilterint=.48;

%% Setup 1 Only Polarizers

% Normalize wrt MAX
normRedIntSetup1=IntensitypolarizedLightOfRedLaserSetup1./unblockredint;
normGreenIntSetup1=IntensitypolarizedLightOfGreenLaserSetup1./unblockgreenint;

% Create sin^2 and cos^2 for fitting of Setup 1
x1=linspace(0,350,360);
y1=0.4233.*(cosd(x1)).^2 +.06;
y2=max(normRedIntSetup1).*(sind(x1+10)).^2;

% do not forget to mention 10 degree phase shift on red and increasing
% amplitude of green due to 'warming up' of laser and also vertical
% shifting of the green data (suggesting partial polarization)

figure(1)
scatter(Angle,normGreenIntSetup1,'green');
hold on
plot(x1,y1,'green');
scatter(Angle, normRedIntSetup1,'red');
plot(x1,y2,'red');
ylabel('Normalized Intensity (%)')
xlabel('Angle (Degrees)')

% Creating polar plots for the same data
Radians=Angle.*(pi/180);
Radians1=x1.*(pi/180);

figure(3)
polarscatter(Radians, normRedIntSetup1,'red');
hold on
polarplot(Radians1, y2,'red');
polarscatter(Radians,normGreenIntSetup1,'green')
polarplot(Radians1,y1,'green');


%% Setup 2 Only Polarizers

% normalize wrt Max
normRedIntSetup2=IntensitypolarizedLightOfRedLaserP190P20Setup2./unblockredint;
normGreenIntSetup2=IntensitypolarizedLightOfGreenLaserP10P20Setup2./unblockgreenint;

% Create sin^2 and cos^2 for fitting of Setup 2

y3=0.34.*(cosd(x1)).^2;
y4=max(normRedIntSetup2).*(sind(x1)).^2;

% interesting, phase shift from red is now gone with two polarizers
% same story for green with increasing amplitude due to laser 'warming up'

figure(2)
scatter(Angle,normGreenIntSetup2,'green');
hold on
plot(x1,y3,'green');
scatter(Angle, normRedIntSetup2,'red');
plot(x1,y4,'red');
ylabel('Normalized Intensity (%)')
xlabel('Angle (Degrees)')

% Creating polar plots for the same data

figure(4)
polarscatter(Radians, normRedIntSetup2,'red');
hold on
polarplot(Radians1, y4,'red');
polarscatter(Radians,normGreenIntSetup2,'green')
polarplot(Radians1,y3,'green');

%% Setup 2 with color filters

%normalize wrt max
blueFilterNormInt=IntensityblueFilterP10p20./unblockgreenint;
greenFilterNormInt=IntensityGreenFilterP10P20./unblockgreenint;
redFilterNormInt=IntensityRedFilterP190P20./unblockredint;

% Create sin^2 and cos^2 for fitting of Setup 2 w/ filters

y5=max(greenFilterNormInt).*(cosd(x1)).^2;
y6=max(redFilterNormInt).*(sind(x1)).^2;
y7=max(blueFilterNormInt).*(cosd(x1)).^2;

figure(5)
scatter(Angle,blueFilterNormInt,'blue');
hold on
plot(x1,y7,'blue');
scatter(Angle,greenFilterNormInt,'green');
plot(x1,y5,'green');
scatter(Angle,redFilterNormInt,'red');
plot(x1,y6,'red');
ylabel('Normalized Intensity (%)')
xlabel('Angle (Degrees)')

% Creating polar plots for the same data

figure(6)
polarscatter(Radians, redFilterNormInt,'red');
hold on
polarplot(Radians1, y6,'red');
polarscatter(Radians,greenFilterNormInt,'green')
polarplot(Radians1,y5,'green');
polarscatter(Radians,blueFilterNormInt,'blue')
polarplot(Radians1,y7,'blue');

%% Data analysis

% Insertion losses for all graphs
insertLossRed1=1-max(normRedIntSetup1);
insertLossRed2=1-max(normRedIntSetup2);
insertLossGreen1=1-max(normGreenIntSetup1);
insertLossGreen2=1-max(normGreenIntSetup2);
insertLossRedFilt=1-max(redFilterNormInt);
insertLossGreenFilt=1-max(greenFilterNormInt);
insertLossBlueFilt=1-max(blueFilterNormInt);

% Extinction ratio for all graphs
exRatioRed1=10*log10(max(normRedIntSetup1)/min(normRedIntSetup1));
exRatioRed2=10*log10(max(normRedIntSetup2)/min(normRedIntSetup2));
exRatioGreen1=10*log10(max(normGreenIntSetup1)/min(normGreenIntSetup1));
exRatioGreen2=10*log10(max(normGreenIntSetup2)/min(normGreenIntSetup2));
exRatioRedFilt=10*log10(max(redFilterNormInt)/min(redFilterNormInt));
exRatioGreenFilt=10*log10(max(greenFilterNormInt)/min(greenFilterNormInt));
exRatioBlueFilt=10*log10(max(blueFilterNormInt)/min(blueFilterNormInt));

% Degree of linear polarization;

dopRed1=max(normRedIntSetup1)/unblockredint;
dopRed2=max(normRedIntSetup2)/unblockredint;
dopGreen1=max(normGreenIntSetup1)/unblockgreenint;
dopGreen2=max(normGreenIntSetup2)/unblockgreenint;
dopRedFilt=max(redFilterNormInt)/unblockredint;
dopGreenFilt=max(greenFilterNormInt)/unblockgreenint;
dopBlueFilt=max(blueFilterNormInt)/unblockgreenint;

















%% Publication Graphic Settings

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