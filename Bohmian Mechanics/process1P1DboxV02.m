% This script analyzes 1P1D Qtrajectories generated by qTraj1P1DboxV02
% Starter code for Ethan, Pavan, and others
%
% Part I: We set the boundary conditions based on a particle in the box in
%         thermal equilibrium at a specified temperature. The box size is
% fixed. The evolution of the system is based on the Schrodinger equation
% assuming that the system is closed. This means the parameters of the
% wavefunction are fixed and are constant. To characterize the dynamics of
% a closed system, version 02 of the code starts all particle trajectories
% at the same location. However, the intial phase angles are randomized to
% different values. The idea is to compare how different the trajectories
% are from one another based on differences in phase that is initially 
% placed on the expansion coefficients of the energy eigenfunctions.
%
% Part II: We acknowledge the system is not closed. As an open quantum
%          system that is in thermal equilibrium with its environment at a
% constant temperature, energy will fluctuate within the system, and thus
% energy is not conserved. All particle trajectories will start at the
% same location. In addition, initially the same randomized phase will be
% placed on the expansion coeficients of the energy eigenfunctions. Due to
% the system being open, the solutions to the time dependent Schrodinger 
% equation are dynamically altered by the environment as a stochastic 
% process. As an open system it is acknowledged that the length of the 1D
% box will fluctuate. A model assuming a constant length is unrealistic. 
% As such, there are three modes of diffusive processes being modeled for
% the open system, listed as: (1) the length of the box fluctuates since 
% it is impossible to keep the walls perfectly static. (2) The magnitude 
% of the expansion coefficients fluctuate because temperature is causing
% energy exchange, and (3) randomization in phase on the expansion 
% coefficients occurs due to the influence of the environment. The user
% will define the diffusion coefficients for each of these distinct 
% diffusive processes. In this way, and number of diffusion processes can
% be turned off by setting the diffusion coefficient to zero. This will
% allow us to view all particle trajectories to start identically and as
% time passes, the trajectories will diverge. Lyapunov exponents could be
% calculated using these trajectories.  
% ------------------------------------------------------------------------
clear all;
close all;
%%                                                           set constants
kb = 8.6183294*10^-5;                         % Boltzmann Constant in ev/K
hbar = 4.14/(2*pi);                                 % dimension is (eV*fs)
AngChar = char(197);                        % gives the Angstrom character
% -------------------------------------------- convert AMU to evfs^2/Ang^2
% JtoeV = 6.242*10^18;                                            J --> eV
% Regc = 3*10^8;                            speed of light in units of m/s
% Specalc = 3*10^3;                         speed of light in units ang/fs
%AMUtoSU1 = (mass * AMUtokg ) * Regc^2;                   AMU -> E units J
%AMUtoSU2 = ( AMUtoSU1 * JtoeV ) / Specalc^2;        E/c^2 -> evfs^2/Ang^2
%      massAMU*1.6605*(10^-27)*( (Regc^2)/Specialc^2 )*6.242*10^18
%      massAMU*1.6605*6.242*(10^-9)*10^10 = mass*1.6605*6.242*10
% AMUtoKg = 1.6605*10^-27;                                    % AMU --> kg
%%                                                              user input
%baseName = input('Enter base file name: ','s');
baseName = '1P1D300Lx100s12345RK4dt0_037394m0_00054858Ps0s4321Dp0_01Da0_01DL0open';
fileNameLog = [baseName,'Log.txt'];
fid=fopen(fileNameLog);
lines = cell(0,1);
while ~feof(fid)
    line=fgetl(fid);
    lines{end+1} = line;
end
fclose(fid);

mass_match = regexp(lines{6}, 'mass= (\d+\.\d+)', 'tokens');
    if ~isempty(mass_match)
        massAMU = str2double(mass_match{1}{1});
    end
ax_match = regexp(lines{9}, 'ax= (\d+)', 'tokens');
    if ~isempty(ax_match)
        ax = str2double(ax_match{1}{1});
    end
Emax_match = regexp(lines{20}, 'Emax = (\d+\.\d+)', 'tokens');
    if ~isempty(Emax_match)
        Emax = str2double(Emax_match{1}{1});
    end
nSamples_match = regexp(lines{21}, 'nSamples= (\d+)', 'tokens');
    if ~isempty(nSamples_match)
        nSamples = str2double(nSamples_match{1}{1});
    end
runTime_match = regexp(lines{12}, 'Trun= (\d+\.\d+)', 'tokens');
    if ~isempty(runTime_match)
        runTime = str2double(runTime_match{1}{1});
    end
dtObsv_match = regexp(lines{13}, 'dtObsv= (\d+\.\d+)', 'tokens');
    if ~isempty(dtObsv_match)
        dtObsv = str2double(dtObsv_match{1}{1});
    end
dtRK4_match = regexp(lines{14}, 'dtRK4= (\d+\.\d+)', 'tokens');
    if ~isempty(dtRK4_match)
        dtRK4 = str2double(dtRK4_match{1}{1});
    end
% read log file to pull any information that might be needed. FIX ME
% For example, need to read mass and Emax directly from log file. FIX ME
%ax = input('Enter the length of the box (Angstroms): ');
%ax=100;
%massAMU = input('Enter the particle mass (AMU): ');
%massAMU = 0.000548580;
mass = massAMU*103.6484;                        % in units of evfs^2/Ang^2
%Emax = input('Enter maximum energy state: Emax (eV): ');
%Emax = 0.24115;
%nSamples = input('Enter the number of trajectories generated: ');
%nSamples = 10;
%runTime = input('Enter the total run time: (fs) ');
%runTime = 14957.768497518;
%dtObsv = input('Enter observation time resolution => dtObsv (fs): ');
%dtObsv = 0.149577685;
%dtRK4 = input('Enter the dif-eq time resolution given by dtRK4 (fs): ');
%dtRK4 = 0.037394421;
xShift = 0.05*ax;
xLow = 0 - xShift;
xBig = ax + xShift;
%%                                       read results to output data file
fileNameP_x = [baseName,'P_x.txt'];
fileNameV_x = [baseName,'V_x.txt'];
fileNameQpe = [baseName,'Qpe.txt'];
posX = readmatrix(fileNameP_x);
velX = readmatrix(fileNameV_x);
q_PE = readmatrix(fileNameQpe);
[mTest,npts] = size(posX);
  if( mTest ~= nSamples )
  error('mTest is not equal to nSamples');
  end
tArray = dtObsv*(0:npts-1);
  if( abs(tArray(end) - runTime) > 0.5*dtObsv )
  error('tArray is not accurate');
  end
%%                                              replace potential outliers
qUmax = 100*abs(Emax);
qUmin = -qUmax;
Lmin = (q_PE < qUmin );
q_PE(Lmin) = qUmin;
Lmax = (q_PE > qUmax );
q_PE(Lmax) = qUmax;
%%                                              publication quality
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
%%                                                        plot information
% ---------------------------------- plot the last 10% of the trajectories
nFinal = length(tArray);
nStart = ceil(0.8*nFinal);
%%                                                   plot position vs time
nFig = 0;
% --------------------------------------------------------- plot X vs time
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
  for s=1:nSamples
  plot(tArray(nStart:nFinal),posX(s,nStart:nFinal));
  end
ylim([xLow,xBig]);
yline(ax);
yline(0);
xlabel('t  (fs)');
ylabel(['X(t)  (',AngChar,')']);
title(['nS= ',num2str(nSamples),'  dtRK4= ',num2str(dtRK4), ...
       ' fs  mass= ',num2str(massAMU),' AMU']); 
%%                                              plot each position vs time
% --------------------------------------------------------- plot X vs time
mFig = 99;
  for s=1:nSamples
  mFig = mFig + 1;
  figure(mFig);
  clf;
  plot(tArray,posX(s,:));
  ylim([xLow,xBig]);
  yline(ax);
  yline(0);
  xlabel('t  (fs)');
  ylabel(['X(t)  (',AngChar,')']);
  title(['traj= ',num2str(s),'  dtRK4= ',num2str(dtRK4), ...
         ' fs  mass= ',num2str(massAMU),' AMU']);
  end 
%%                                                   plot velocity vs time
% -------------------------------------------------------- plot Vx vs time
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
  for s=1:nSamples
  plot(tArray(nStart:nFinal),velX(s,nStart:nFinal));
  xlabel('t  (fs)');
  ylabel(['Vx(t)  (',AngChar,'/fs)']);
  title(['trajectory: ',num2str(s),'  dtRK4= ',num2str(dtRK4), ...
         ' fs  mass= ',num2str(massAMU),' AMU']);
  end
%%                                                     plot energy vs time
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
  for s=1:nSamples
  Earray = 0.5*mass*velX(s,nStart:nFinal).^2 + q_PE(s,nStart:nFinal);
  plot(tArray(nStart:nFinal),Earray);
  xlabel('t  (fs)');
  ylabel('E(t)  (eV)');
  title(['trajectory: ',num2str(s),'  dtRK4= ',num2str(dtRK4), ...
       ' fs  mass= ',num2str(massAMU),' AMU']);
  aveEnergy = mean(Earray);
  disp(['trajectory: ',num2str(s),'   aveEnergy= ',num2str(aveEnergy)]);
  end
%%                                               plot velocity vs position
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
  for s=1:nSamples
  plot(posX(s,nStart:nFinal),velX(s,nStart:nFinal), 'LineWidth', 0.75);
  xlabel(['X(t)  (',AngChar,')']);
  ylabel(['Vx(t)  (',AngChar,'/fs)']);
  title(['trajectory: ',num2str(s),'  dtRK4= ',num2str(dtRK4), ...
       ' fs  mass= ',num2str(massAMU),' AMU']); 
  end
%%                                      plot quantum potential vs position
nFig = nFig + 1;
figure(nFig);
clf;
hold on;
  for s=1:nSamples
  qPEarray = q_PE(s,:);
  plot(posX(s,nStart:nFinal),qPEarray(nStart:nFinal));
  xlabel(['X(t)  (',AngChar,')']);
  ylabel('qU(x)  (eV)'); 
  title(['trajectory: ',num2str(s),'  dtRK4= ',num2str(dtRK4), ...
       ' fs  mass= ',num2str(massAMU),' AMU']);
  end
%% Velocity movie
% To get this to work I need to get psi for all space stop at each dt.
%{
runTime = 10000; %input('What is the run time for the experiment (fs)? ');
dt      = 1; %input('input the value of dt: ');
to      = 0; %input('Enter the intial time, to (fs): ');
tFinal = to + runTime;
tArray = to:dt:tFinal;
% --------------------------
xArray4 = logspace(-5,-3,40);
% --------------------------
dx = 0.001;
b3 = 0.01*b;
xL = dx/2;
xArray3 = xL:dx:b3;
% --------------------------
dx = 0.01;
b2 = 0.1*b;
xL = b3;
xArray2 = xL:dx:b2;
% --------------------------
dx = 0.1;
b1 = 0.5*b;
xL = b2;
xArray1 = xL:dx:b1;
% --------------------------
xArray5 = b - xArray1;
xArray6 = b - xArray2;
xArray7 = b - xArray3;
xArray8 = b - xArray4;
xArray = unique([xArray1,xArray2,xArray3,xArray4, ...
                 xArray5,xArray6,xArray7,xArray8]);
% ------------------------------------------------------------------------
nPts = length(xArray);
velx = zeros(1,nPts);
% -------------------------------------------------------- estimate limits
vmin = 1.0e50;
vmax = -vmin;
dt2 = (tFinal - to)/20;
for t=to:dt2:tFinal
    for j=1:nPts
    x = xArray(j);
    vTest = GetVel(qSystem,x,t);
      if( vTest < vmin )
      vmin = vTest;
      end
      if( vTest > vmax )
      vmax = vTest;
      end
    end
end
% cut peaks by half;
vmin = 0.5*vmin;
vmax = 0.5*vmax;
% ----------------------------------------------------------- plot "movie"
figure(2)
clf;
  for t=tArray
    for j=1:nPts
    x = xArray(j);
    velx(j) = GetVel(qSystem,x,t);
    end
  plot(xArray,velx,'k','LineWidth',1.8);
  xlim([0,b]);
  %ylim([vmin,vmax]);
  ylim([-5,5]);
  pause(0.05);
  end
%}