clear all
close all
%%                                                        define constants
hbar = 4.14/(2*pi);                                 % dimension is (eV*fs)
c = 3000;                                  % speed of light in Angstrom/fs
eEnergy = 510990.6;                        % rest energy of electron in eV
eMass = eEnergy/c^2;                         % in units of (fs^2)*eV/Ang^2
%%                                                         user directives
mF = 75;                             % number of points spread over sigmaF
aFx = 20;                % aFx*sigmaF = amplitude for field in x-direction
aFy = 20;                % aFy*sigmaF = amplitude for field in y-direction
fShift = 0.75;              % 0 => 0 degree shift   1 => 2*pi degree shift
% fShift = 0.50 gives stable short arcs (back and forth motion on arcs)
%               total area for (x-direction, y-direction) = (2,1)
% fShift = 0.75 gives total area for (x-dir, y-dir) = (1.006, 1.006)
%               Produces ~chatic loop motion that is very localized
% fShift = pi/5
%   fShift = -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1 are interesting
mass = 1*eMass;             % mass of particle in units of (fs^2)*eV/Ang^2
ks = 0.0001;                                  % spring constant (ev/Ang^2)
spaceBuffer = 6;
% --------------------------------------------------- select option number
% 1: check 1D wave motion formulas
% 2: check 2D elliptical motion using 2 orthogonal separable wavefunctions
% 3: check numerical method to integrate equations of motion for particle
% 4: check if reverse elliptical orbits can be created (opposite of #2)
% 5: without superposition, view simutaneous forward & reverse motions
% 6: trajectories for 2D interference of two wavepackets pulsating
% ------------------------------------------------------------------------
opt = 6;
nLaps = 25;                  % used to extend number of periods to observe
% ------------------------------------------ initialize all cases to false
checkTestWaveMotion = false;      %1
checkEllipticalMotion = false;    %2
checkNumericalMethod = false;     %3
checkReversedMotions = false;     %4
checkForwardReversed = false;     %5
exact2DInterference = false;      %6            the case we wanted to see!
% ------------------------------------------------------------ select case
  switch opt
      case 1
      checkTestWaveMotion = true;
      case 2
      checkEllipticalMotion = true;
      case 3
      checkNumericalMethod = true;
      case 4
      checkReversedMotions = true;
      case 5
      checkForwardReversed = true;
      case 6
      exact2DInterference = true;
      otherwise
      error('option is unavailable');
  end
% ----------------------------------------------------------- set up cases
  if( checkEllipticalMotion || checkNumericalMethod )
  nWaves = 1;
  waveData = cell(1,nWaves);
    for nW=1:nWaves
    waveData{nW} = struct;
    end
% ------------------------------------------- set initial field conditions
  waveData{1}.a_xFo = aFx;            % amplitude for field in x-direction
  waveData{1}.a_yFo = aFy;            % amplitude for field in y-direction
  waveData{1}.f_xFo = 0;               % fraction of period in x-direction
  waveData{1}.f_yFo = 0.25;            % fraction of period in y-direction
  end
% ------------------------------------------------------------------------
  if( checkReversedMotions )
  nWaves = 1;
  waveData = cell(1,nWaves);
    for nW=1:nWaves
    waveData{nW} = struct;
    end
% ------------------------------------------- set initial field conditions
  waveData{1}.a_xFo = aFx;            % amplitude for field in x-direction
  waveData{1}.a_yFo = aFy;            % amplitude for field in y-direction
  waveData{1}.f_xFo = -0.5;            % fraction of period in x-direction
  waveData{1}.f_yFo = -0.75;           % fraction of period in y-direction
  checkEllipticalMotion = true;       % => only changed initial conditions
  end
% ------------------------------------------------------------------------
  if( checkForwardReversed || exact2DInterference )
  nWaves = 2;
  waveData = cell(1,nWaves);
    for nW=1:nWaves
    waveData{nW} = struct;
    end
% --------------------------------------- set forward motion wavefunctions
  waveData{1}.a_xFo = aFx;            % amplitude for field in x-direction
  waveData{1}.a_yFo = aFy;            % amplitude for field in y-direction
  waveData{1}.f_xFo = 0;               % fraction of period in x-direction
  waveData{1}.f_yFo = 0.25;            % fraction of period in y-direction
% --------------------------------------- set reverse motion wavefunctions
  waveData{2}.a_xFo = aFx;            % amplitude for field in x-direction
  waveData{2}.a_yFo = aFy;            % amplitude for field in y-direction
  waveData{2}.f_xFo = -0.50 + fShift;  % fraction of period in x-direction
  waveData{2}.f_yFo = -0.75 + fShift;  % fraction of period in y-direction
  end
%%                                                  dependent consequences
omega=sqrt(ks/mass);              % angular frequency of field disturbance
Tperiod = 2*pi/omega;                        % period of field disturbance
sigmaF=sqrt(hbar/(mass*omega));        % sigmaF => std for spread in field
const=(mass*omega)/hbar;
c1 = sqrt(const);                                   % spatial scale factor
%%                                                    scaling of dt and ds
dtFlow = Tperiod/1000;
%ds = sigmaF/333;                                     % spatial resolution
ds = sigmaF/mF;
% -------------------------------------- compare to numerical calculations
fixConst1 = 0.1;                           % can be anything less than 0.5
b = hbar./(2*mass);
dt = ds*ds*fixConst1./b;        % for calculations of Schrodinger equation
Nrepeat = ceil(dtFlow/dt);
dt = dtFlow/Nrepeat;
% ----------------------------------------------------- record information
disp('   ');
disp('------------------------------------------------------ parameters');
   if( (opt == 6) || (opt == 5) )
   disp(['  fShift = ',num2str(fShift)]);
   end
disp(['    mass = ',num2str(mass/eMass),' eMass']);
disp(['    Kcpw = ',num2str(ks),' ev/Ang^2']);
disp(['   omega = ',num2str(omega),' 1/fs']);
disp(['  period = ',num2str(Tperiod),' fs']);
disp(['  sigmaF = ',num2str(sigmaF),' Ang']);
disp(['     aFx = ',num2str(aFx),' Ang']);
disp(['     aFy = ',num2str(aFy),' Ang']);
disp(['flow: dt = ',num2str(dtFlow),' fs']);       % plotting trajectories
disp(['grid: dt = ',num2str(dt),' fs']);      % integration & evolving psi
disp(['grid: ds = ',num2str(ds),' Ang']);                   % evolving psi
disp('-----------------------------------------------------------------');
%%                                         set default graphics parameters
% REMARKS: These set functions are to make the plots ~publication quality.
%          Using the default set functions for axes, font size and figure
%          color are for convenience only.
%fig.InvertHardcopy = 'off';         
set(0,'DefaultFigureColor','white')
set(0,'defaultAxesFontSize',14); 
set(0,'defaultLineLineWidth',1.4);   
set(0,'defaultLineMarkerSize',7); 
set(0,'defaultAxesLineWidth',1.4);
%set(0,'defaultGeographicAxes','on');
%%                                               test wave function motion
if( checkTestWaveMotion )
amplitude = 25;                                      % multiples of sigmaF                    
fractionT = 0.5;                                      % fraction of period
Atwf = amplitude*sigmaF;                          % amplitude of oscilaton
xMaximum = Atwf + 10*sigmaF;             % add 10 sigma more at boundaries
xMinimum = -xMaximum;
xtwf = xMinimum:ds:xMaximum;
Ftwf = fractionT;                                   % fraction of 1 period
Nx = length(xtwf);
% ---------------------------------------------------- dependent variables
t0 = Ftwf*Tperiod;
tArray = 0:dtFlow:(Tperiod - dtFlow/4);
ampTWF = c1*(const/pi)^(-1/4);
% ------------------------------------------- work with test wave function
Xtwf=c1*xtwf;
X0=c1*Atwf;
T0=omega*t0;
jMax = length(tArray);
% -------------------------------------------- augment summary information
disp(['    Nx = ',num2str(Nx)]);
disp(['  jMax = ',num2str(jMax)]);
totalProb = zeros(1,jMax);
k = 0;
  for j=1:5:jMax
  t = tArray(j);
  T = t*omega;
  k = k + 1;
% ------------------------------------------------- set unitless variables
  XoT1 = X0*cos(T-T0);
  K1 = -1*X0*sin(T-T0);
  Phi1 = (-1/4)*(X0^2)*sin(2*(T-T0));
  dsTWF = Xtwf - XoT1;
  imTWF = Xtwf*K1 - Phi1 - 0.5*(T - T0);
  psiTWF = ampTWF*exp( -0.5*dsTWF.^2 + 1i*imTWF );
  pdf = abs(psiTWF).^2;
  totalProb(k) = sum(pdf)*ds;
  end
kMax = k;
temp = mean( totalProb(1:kMax) );
totalProb = totalProb/temp;
temp = 1/sqrt(temp);
ampTWF = temp*ampTWF;
disp(['scale factor for normalization constant = ',num2str(ampTWF)]);
% ----------------------------------------------------------- plot results
figure(1);
clf;
plot( (1:5:jMax)*dtFlow, totalProb(1:kMax), 'k-' );
ylim( [0.9,1.1] );
% --------------------------------------------------------------
figure(2);
clf;
dsTWF = Xtwf - XoT1;
imTWF = Xtwf*K1 - Phi1 - 0.5*(T - T0);
psiTWF = ampTWF*exp( -0.5*dsTWF.^2 + 1i*imTWF );
pdf = abs(psiTWF).^2;
plot( xtwf, pdf, 'k-' );
xlabel('x coordinate');
ylabel('| \psi |^2');
% --------------------------------------------------------------
figure(3);
clf;
dsTWF = Xtwf - XoT1;
imTWF = Xtwf*K1 - Phi1 - 0.5*(T - T0);
psiTWF = ampTWF*exp( -0.5*dsTWF.^2 + 1i*imTWF );
plot( xtwf, real(psiTWF), 'k-' );
xlabel('x coordinate');
ylabel('RE[ \psi ]');
% --------------------------------------------------------------
figure(4);
clf;
dsTWF = Xtwf - XoT1;
imTWF = Xtwf*K1 - Phi1 - 0.5*(T - T0);
psiTWF = ampTWF*exp( -0.5*dsTWF.^2 + 1i*imTWF );
plot( xtwf, imag(psiTWF), 'k-' );
xlabel('x coordinate');
ylabel('IM[ \psi ]');
% ----------------------------------- plot real part as a function of time
figure(5);
clf;
  for j=1:jMax
  t = tArray(j);
  T = t*omega;
% ------------------------------------------------- set unitless variables
  XoT1 = X0*cos(T-T0);
  K1 = -1*X0*sin(T-T0);
  Phi1 = (-1/4)*(X0^2)*sin(2*(T-T0));
  dsTWF = Xtwf - XoT1;
  imTWF = Xtwf*K1 - Phi1 - 0.5*(T - T0);
  psiTWF = ampTWF*exp( -0.5*dsTWF.^2 + 1i*imTWF );
  figure(5);
  clf;
  hold on
  plot( Xtwf/c1, real(psiTWF), 'k-' );
  plot( [Atwf,Atwf], [-0.25,0.25], 'r--' );
  plot( [-Atwf,-Atwf], [-0.25,0.25], 'r--' );
  title('real part of wavefunction');
  pause(0.05)
  end
end
%%                                                 check elliptical motion
   if( checkEllipticalMotion )
   ampF = c1*(const/pi)^(-1/4);
   xFo = waveData{1}.a_xFo*sigmaF;
   yFo = waveData{1}.a_yFo*sigmaF;
   Txo = waveData{1}.f_xFo*Tperiod;
   Tyo = waveData{1}.f_yFo*Tperiod;
   xFmax = xFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   yFmax = yFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   Fmax = max(xFmax,yFmax);
   xFmin = -xFmax;
   yFmin = -yFmax;
   xF = xFmin:ds:xFmax;
   yF = yFmin:ds:yFmax;
   Nx = length(xF);
   Ny = length(yF);
   tArray = 0:dtFlow:(Tperiod - dtFlow/4);
   Nt = length(tArray);
% -------------------------------------------- augment summary information
   disp(['    Nx = ',num2str(Nx)]);
   disp(['    Ny = ',num2str(Ny)]);
   disp(['    Nt = ',num2str(Nt)]);
% ------------------------------------------ work elliptical wave function
   Xs = c1*xF;
   Xs0 = c1*xFo;
   Tx0 = omega*Txo;
% ------------------- 
   Ys = c1*yF;
   Ys0 = c1*yFo;
   Ty0 = omega*Tyo;     
% ----------------------------------------------- normalize wave functions
   totalProb = zeros(1,Nt);
   k = 0;
     for j=1:5:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSx = Xs0*cos(T-Tx0);
     negSINx =-Xs0*sin(T-Tx0);
     x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
     dXs = Xs - posCOSx;
     imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
     psi = ampF*exp( -0.5*dXs.^2 + 1i*imXpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(temp);
   ampFx = temp*ampF;
   disp(['scale factor for X-normalization constant = ',num2str(temp)]);
   pdfXmax = max(pdf)*temp^2;
% ---------------------------------- normalize wavefunction in y-direction
   totalProb = 0*totalProb;
   k = 0;
     for j=1:5:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSy = Ys0*cos(T-Ty0);
     negSINy =-Ys0*sin(T-Ty0);
     y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
     dYs = Ys - posCOSy;
     imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
     psi = ampF*exp( -0.5*dYs.^2 + 1i*imYpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(temp);
   ampFy = temp*ampF;
   disp(['scale factor for Y-normalization constant = ',num2str(temp)]);
   pdfYmax = max(pdf)*temp^2;   
% --------------------------------------------------- calculate trajectory
   shift_xFo = min(0.5*xFo,1.5*sigmaF);
   shift_yFo = min(0.5*yFo,1.5*sigmaF);
   xPo1 = xFo;
   yPo1 = yFo;
   xPo2 = xFo + shift_xFo;
   yPo2 = yFo + shift_yFo;
   xPo3 = xFo - shift_xFo;
   yPo3 = yFo - shift_yFo;
% % % %    xPo1 = 1.00*xFo;
% % % %    yPo1 = 1.00*yFo;
% % % %    xPo2 = 1.50*xFo;
% % % %    yPo2 = 1.50*yFo;
% % % %    xPo3 = 0.50*xFo;
% % % %    yPo3 = 0.50*yFo;
   %xP1 = xFo*cos( omega*(tArray - Txo) ) + xPo1 - xFo;
   %yP1 = yFo*cos( omega*(tArray - Tyo) ) + yPo1 - yFo;
   xP2 = xFo*cos( omega*(tArray - Txo) ) + xPo2 - xFo;
   yP2 = yFo*cos( omega*(tArray - Tyo) ) + yPo2 - yFo;
   xP3 = xFo*cos( omega*(tArray - Txo) ) + xPo3 - xFo;
   yP3 = yFo*cos( omega*(tArray - Tyo) ) + yPo3 - yFo;
   xFpeak = xFo*cos( omega*(tArray - Txo) );
   yFpeak = yFo*cos( omega*(tArray - Tyo) );
% --------------------------------------------- plot |psi|^2 decomposition 
     for j=1:10:Nt
     figure(1);
     clf; 
     axis equal
     t = tArray(j);
     T = t*omega;
% -------------------------------------------- show trajectory of particle
     subplot(2,2,2);
     hold on
     plot(xFpeak,yFpeak,'k','linewidth',1.3);
     plot(xP2,yP2,'color',[0.6,0.6,0.6],'linewidth',1.5);
     plot(xP3,yP3,'color',[0.6,0.6,0.6],'linewidth',1.5);
     xlim( [-Fmax,Fmax] );
     ylim( [-Fmax,Fmax] );
     xP1t = xFo*cos( omega*(t - Txo) );
     yP1t = yFo*cos( omega*(t - Tyo) );
     xP2t = xFo*cos( omega*(t - Txo) ) + xPo2 - xFo;
     yP2t = yFo*cos( omega*(t - Tyo) ) + yPo2 - yFo;
     xP3t = xFo*cos( omega*(t - Txo) ) + xPo3 - xFo;
     yP3t = yFo*cos( omega*(t - Tyo) ) + yPo3 - yFo;
     scatter(xP1t,yP1t,60,'b','fill');
     scatter(xP2t,yP2t,60,'r','fill');
     scatter(xP3t,yP3t,60,'m','fill');
% -------------------------------------------------------- pdf x-direction
     posCOSx = Xs0*cos(T-Tx0);
     negSINx =-Xs0*sin(T-Tx0);
     x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
     dXs = Xs - posCOSx;
     imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
     psi = ampFx*exp( -0.5*dXs.^2 + 1i*imXpart );
     pdf = abs(psi).^2;
     subplot(2,2,4);
     hold on
     plot(xF,pdf,'k-' );
     plot( [xFo,xFo],  [0,1.1*pdfXmax],'r--' );
     plot( [-xFo,-xFo],[0,1.1*pdfXmax],'r--' );
     [~,indx] = sort( abs(xF - xP1t) );
     indx1 = indx(1);
     [~,indx] = sort( abs(xF - xP2t) );
     indx2 = indx(1);
     [~,indx] = sort( abs(xF - xP3t) );
     indx3 = indx(1);
     scatter(xP1t,pdf(indx1),60,'b','fill');
     scatter(xP2t,pdf(indx2),60,'r','fill');
     scatter(xP3t,pdf(indx3),60,'m','fill');
     xlim( [-Fmax,Fmax] );
% -------------------------------------------------------- pdf y-direction
     posCOSy = Ys0*cos(T-Ty0);
     negSINy =-Ys0*sin(T-Ty0);
     y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
     dYs = Ys - posCOSy;
     imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
     psi = ampFy*exp( -0.5*dYs.^2 + 1i*imYpart );
     pdf = abs(psi).^2;
     subplot(2,2,1);
     hold on
     plot(pdf,yF,'k-' );
     plot( [0,1.1*pdfYmax], [yFo,yFo],'r--' );
     plot( [0,1.1*pdfYmax], [-yFo,-yFo],'r--' );
     [~,indx] = sort( abs(yF - yP1t) );
     indx1 = indx(1);
     [~,indx] = sort( abs(yF - yP2t) );
     indx2 = indx(1);
     [~,indx] = sort( abs(yF - yP3t) );
     indx3 = indx(1);
     scatter(pdf(indx1),yP1t,60,'b','fill');
     scatter(pdf(indx2),yP2t,60,'r','fill');
     scatter(pdf(indx3),yP3t,60,'m','fill');
     ylim( [-Fmax,Fmax] );
     pause(0.05)
     end
   end
%%                                                  check numerical method
   if( checkNumericalMethod )
   ampF = c1*(const/pi)^(-1/4);
   xFo = waveData{1}.a_xFo*sigmaF;
   yFo = waveData{1}.a_yFo*sigmaF;
   Txo = waveData{1}.f_xFo*Tperiod;
   Tyo = waveData{1}.f_yFo*Tperiod;
   xFmax = xFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   yFmax = yFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   Fmax = max(xFmax,yFmax);
   xFmin = -xFmax;
   yFmin = -yFmax;
   xF = xFmin:ds:xFmax;
   yF = yFmin:ds:yFmax;
   Nx = length(xF);
   Ny = length(yF);
   tFinal = nLaps*Tperiod + 3*dtFlow;
   tArray = 0:dtFlow:tFinal;
   Nt = length(tArray);
% -------------------------------------------- augment summary information
   disp(['    Nx = ',num2str(Nx)]);
   disp(['    Ny = ',num2str(Ny)]);
   disp(['    Nt = ',num2str(Nt)]);
% ------------------------------------------ work elliptical wave function
   Xs = c1*xF;
   Xs0 = c1*xFo;
   Tx0 = omega*Txo;
% ------------------- 
   Ys = c1*yF;
   Ys0 = c1*yFo;
   Ty0 = omega*Tyo;     
% ----------------------------------------------- normalize wave functions
   totalProb = zeros(1,Nt);
   k = 0;
     for j=1:25:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSx = Xs0*cos(T-Tx0);
     negSINx =-Xs0*sin(T-Tx0);
     x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
     dXs = Xs - posCOSx;
     imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
     psi = ampF*exp( -0.5*dXs.^2 + 1i*imXpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(temp);
   ampFx = temp*ampF;
   disp(['scale factor for X-normalization constant = ',num2str(temp)]);
   pdfXmax = max(pdf)*temp^2;
% ---------------------------------- normalize wavefunction in y-direction
   totalProb = 0*totalProb;
   k = 0;
     for j=1:25:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSy = Ys0*cos(T-Ty0);
     negSINy =-Ys0*sin(T-Ty0);
     y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
     dYs = Ys - posCOSy;
     imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
     psi = ampF*exp( -0.5*dYs.^2 + 1i*imYpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(temp);
   ampFy = temp*ampF;
   disp(['scale factor for Y-normalization constant = ',num2str(temp)]);
   pdfYmax = max(pdf)*temp^2;   
% --------------------------------------------------- calculate trajectory
   dj = 2;                                       % any even number will do
   Nrepeat2 = round(Nrepeat*dj/2);
% ----------------------------------
   shift_xFo = min(0.5*xFo,1.5*sigmaF);
   shift_yFo = min(0.5*yFo,1.5*sigmaF);
   xPo1 = xFo;
   yPo1 = yFo;
   xPo2 = xFo + shift_xFo;
   yPo2 = yFo + shift_yFo;
   xPo3 = xFo - shift_xFo;
   yPo3 = yFo - shift_yFo;
% % % %    xPo1 = 1.00*xFo;
% % % %    yPo1 = 1.00*yFo;
% % % %    xPo2 = 1.50*xFo;
% % % %    yPo2 = 1.50*yFo;
% % % %    xPo3 = 0.50*xFo;
% % % %    yPo3 = 0.50*yFo;
   %xP1 = xFo*cos( omega*(tArray - Txo) ) + xPo1 - xFo;
   %yP1 = yFo*cos( omega*(tArray - Tyo) ) + yPo1 - yFo;
   xP2 = xFo*cos( omega*(tArray - Txo) ) + xPo2 - xFo;
   yP2 = yFo*cos( omega*(tArray - Tyo) ) + yPo2 - yFo;
   xP3 = xFo*cos( omega*(tArray - Txo) ) + xPo3 - xFo;
   yP3 = yFo*cos( omega*(tArray - Tyo) ) + yPo3 - yFo;
   xFpeak = xFo*cos( omega*(tArray - Txo) );
   yFpeak = yFo*cos( omega*(tArray - Tyo) );
% ------------------------------------------------------------------------
   xP1Locate = xFo*cos( omega*(-dj*dtFlow/2 - Txo) ) + xPo1 - xFo;
   yP1Locate = yFo*cos( omega*(-dj*dtFlow/2 - Tyo) ) + yPo1 - yFo;
   xP2Locate = xFo*cos( omega*(-dj*dtFlow/2 - Txo) ) + xPo2 - xFo;
   yP2Locate = yFo*cos( omega*(-dj*dtFlow/2 - Tyo) ) + yPo2 - yFo;
   xP3Locate = xFo*cos( omega*(-dj*dtFlow/2 - Txo) ) + xPo3 - xFo;
   yP3Locate = yFo*cos( omega*(-dj*dtFlow/2 - Tyo) ) + yPo3 - yFo;
% ------------------------------------------------------------------------
   estimatedP1x = zeros(1,Nt);
   estimatedP1y = zeros(1,Nt);
   estimatedP2x = zeros(1,Nt);
   estimatedP2y = zeros(1,Nt);
   estimatedP3x = zeros(1,Nt);
   estimatedP3y = zeros(1,Nt);
   exactP1x = zeros(1,Nt);
   exactP1y = zeros(1,Nt);
   exactP2x = zeros(1,Nt);
   exactP2y = zeros(1,Nt);
   exactP3x = zeros(1,Nt);
   exactP3y = zeros(1,Nt);
% ------------------------- apply Bohmain procedure for quantum trajectory
   dShiftNeg = -3*ds;
   dShiftPos = 3.1*ds;
   iSpan = 1:7;
   xFit = ds*(iSpan - 1);
   k = 0;
     for j=1:dj:Nt
     t = tArray(j);
     T = t*omega;     
% {     
% ---------------------------------------------------- determine Xs and Ys 
       for nR=1:Nrepeat2
       % -------------------------------------------------------
       x1 = (xP1Locate + dShiftNeg):ds:(xP1Locate + dShiftPos);
       y1 = (yP1Locate + dShiftNeg):ds:(yP1Locate + dShiftPos);
       % -------------------------------------------------------
       x2 = (xP2Locate + dShiftNeg):ds:(xP2Locate + dShiftPos);
       y2 = (yP2Locate + dShiftNeg):ds:(yP2Locate + dShiftPos);
       % -------------------------------------------------------
       x3 = (xP3Locate + dShiftNeg):ds:(xP3Locate + dShiftPos);
       y3 = (yP3Locate + dShiftNeg):ds:(yP3Locate + dShiftPos);
       % ------------------------------------------------------ scale x,y
       Xs1 = c1*x1;
       Xs2 = c1*x2;
       Xs3 = c1*x3;
       Ys1 = c1*y1;
       Ys2 = c1*y2;
       Ys3 = c1*y3;
% ------------------------------------------------ calculate wavefunctions
       posCOSx = Xs0*cos(T-Tx0);
       negSINx =-Xs0*sin(T-Tx0);
       x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
% -------------------------------------------------------
       dXs1 = Xs1 - posCOSx;
       imXpart = Xs1*negSINx - x_phase - 0.5*(T - Tx0);
       psi1x = ampFx*exp( -0.5*dXs1.^2 + 1i*imXpart );
% -------------------------------------------------------
       dXs2 = Xs2 - posCOSx;
       imXpart = Xs2*negSINx - x_phase - 0.5*(T - Tx0);
       psi2x = ampFx*exp( -0.5*dXs2.^2 + 1i*imXpart );
% -------------------------------------------------------
       dXs3 = Xs3 - posCOSx;
       imXpart = Xs3*negSINx - x_phase - 0.5*(T - Tx0);
       psi3x = ampFx*exp( -0.5*dXs3.^2 + 1i*imXpart );
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       posCOSy = Ys0*cos(T-Ty0);
       negSINy =-Ys0*sin(T-Ty0);
       y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
% -------------------------------------------------------
       dYs1 = Ys1 - posCOSy;
       imYpart = Ys1*negSINy - y_phase - 0.5*(T - Ty0);
       psi1y = ampFy*exp( -0.5*dYs1.^2 + 1i*imYpart );
% -------------------------------------------------------
       dYs2 = Ys2 - posCOSy;
       imYpart = Ys2*negSINy - y_phase - 0.5*(T - Ty0);
       psi2y = ampFy*exp( -0.5*dYs2.^2 + 1i*imYpart );
% -------------------------------------------------------
       dYs3 = Ys3 - posCOSy;
       imYpart = Ys3*negSINy - y_phase - 0.5*(T - Ty0);
       psi3y = ampFy*exp( -0.5*dYs3.^2 + 1i*imYpart );
% ============================================================= get action
       qPhase = atan2( imag(psi1x), real(psi1x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel1 = hbar*( aCoef(2) + 2*aCoef(1)*xP1Locate )/mass;
       vel = -omega*xFo*sin( omega*(t - Txo) );
         if( abs(vel - vel1) > 1.0e-9 )
         disp( [vel1,vel] );
         end
       xP1Locate = xP1Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi2x), real(psi2x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP2Locate )/mass;
       xP2Locate = xP2Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi3x), real(psi3x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP3Locate )/mass;
       xP3Locate = xP3Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi1y), real(psi1y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP1Locate )/mass;
       yP1Locate = yP1Locate + vel*dt;
% ----------------------------------------------------------       
       qPhase = atan2( imag(psi2y), real(psi2y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP2Locate )/mass;
       yP2Locate = yP2Locate + vel*dt;
% ----------------------------------------------------------       
       qPhase = atan2( imag(psi3y), real(psi3y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP3Locate )/mass;
       yP3Locate = yP3Locate + vel*dt;
       end    
%}   
% --------------------------------------- calculate trajectory of particle
     xP1t = xFo*cos( omega*(t - Txo) );
     yP1t = yFo*cos( omega*(t - Tyo) );
     xP2t = xFo*cos( omega*(t - Txo) ) + xPo2 - xFo;
     yP2t = yFo*cos( omega*(t - Tyo) ) + yPo2 - yFo;
     xP3t = xFo*cos( omega*(t - Txo) ) + xPo3 - xFo;
     yP3t = yFo*cos( omega*(t - Tyo) ) + yPo3 - yFo;
% ----------------------------------------------------- record information
     k = k + 1;
     exactP1x(k) = xP1t;
     exactP1y(k) = yP1t;
     %-------------------
     exactP2x(k) = xP2t;
     exactP2y(k) = yP2t;
     %-------------------
     exactP3x(k) = xP3t;
     exactP3y(k) = yP3t;
     %----------------------------
     estimatedP1x(k) = xP1Locate;
     estimatedP1y(k) = yP1Locate;
     %----------------------------
     estimatedP2x(k) = xP2Locate;
     estimatedP2y(k) = yP2Locate;
     %----------------------------
     estimatedP3x(k) = xP3Locate;
     estimatedP3y(k) = yP3Locate;     
% ---------------------------------------------------- determine Xs and Ys 
       for nR=1:Nrepeat2
       % -------------------------------------------------------
       x1 = (xP1Locate + dShiftNeg):ds:(xP1Locate + dShiftPos);
       y1 = (yP1Locate + dShiftNeg):ds:(yP1Locate + dShiftPos);
       % -------------------------------------------------------
       x2 = (xP2Locate + dShiftNeg):ds:(xP2Locate + dShiftPos);
       y2 = (yP2Locate + dShiftNeg):ds:(yP2Locate + dShiftPos);
       % -------------------------------------------------------
       x3 = (xP3Locate + dShiftNeg):ds:(xP3Locate + dShiftPos);
       y3 = (yP3Locate + dShiftNeg):ds:(yP3Locate + dShiftPos);
       % ------------------------------------------------------ scale x,y
       Xs1 = c1*x1;
       Xs2 = c1*x2;
       Xs3 = c1*x3;
       Ys1 = c1*y1;
       Ys2 = c1*y2;
       Ys3 = c1*y3;
% ------------------------------------------------ calculate wavefunctions
       posCOSx = Xs0*cos(T-Tx0);
       negSINx =-Xs0*sin(T-Tx0);
       x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
% -------------------------------------------------------
       dXs1 = Xs1 - posCOSx;
       imXpart = Xs1*negSINx - x_phase - 0.5*(T - Tx0);
       psi1x = ampFx*exp( -0.5*dXs1.^2 + 1i*imXpart );
% -------------------------------------------------------
       dXs2 = Xs2 - posCOSx;
       imXpart = Xs2*negSINx - x_phase - 0.5*(T - Tx0);
       psi2x = ampFx*exp( -0.5*dXs2.^2 + 1i*imXpart );
% -------------------------------------------------------
       dXs3 = Xs3 - posCOSx;
       imXpart = Xs3*negSINx - x_phase - 0.5*(T - Tx0);
       psi3x = ampFx*exp( -0.5*dXs3.^2 + 1i*imXpart );
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       posCOSy = Ys0*cos(T-Ty0);
       negSINy =-Ys0*sin(T-Ty0);
       y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
% -------------------------------------------------------
       dYs1 = Ys1 - posCOSy;
       imYpart = Ys1*negSINy - y_phase - 0.5*(T - Ty0);
       psi1y = ampFy*exp( -0.5*dYs1.^2 + 1i*imYpart );
% -------------------------------------------------------
       dYs2 = Ys2 - posCOSy;
       imYpart = Ys2*negSINy - y_phase - 0.5*(T - Ty0);
       psi2y = ampFy*exp( -0.5*dYs2.^2 + 1i*imYpart );
% -------------------------------------------------------
       dYs3 = Ys3 - posCOSy;
       imYpart = Ys3*negSINy - y_phase - 0.5*(T - Ty0);
       psi3y = ampFy*exp( -0.5*dYs3.^2 + 1i*imYpart );
% ============================================================= get action
       qPhase = atan2( imag(psi1x), real(psi1x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel1 = hbar*( aCoef(2) + 2*aCoef(1)*xP1Locate )/mass;
       vel = -omega*xFo*sin( omega*(t - Txo) );
         if( abs(vel - vel1) > 1.0e-9 )
         disp( [vel1,vel] );
         end
       xP1Locate = xP1Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi2x), real(psi2x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP2Locate )/mass;
       xP2Locate = xP2Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi3x), real(psi3x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP3Locate )/mass;
       xP3Locate = xP3Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi1y), real(psi1y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP1Locate )/mass;
       yP1Locate = yP1Locate + vel*dt;
% ----------------------------------------------------------       
       qPhase = atan2( imag(psi2y), real(psi2y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP2Locate )/mass;
       yP2Locate = yP2Locate + vel*dt;
% ----------------------------------------------------------       
       qPhase = atan2( imag(psi3y), real(psi3y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP3Locate )/mass;
       yP3Locate = yP3Locate + vel*dt;
       end        
                
     
%{     
% -------------------------------------------------------- pdf x-direction
     posCOSx = Xs0*cos(T-Tx0);
     negSINx =-Xs0*sin(T-Tx0);
     x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
     dXs = Xs - posCOSx;
     imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
     psi = ampFx*exp( -0.5*dXs.^2 + 1i*imXpart );
     pdf = abs(psi).^2;
     subplot(2,2,4);
     hold on
     plot(xF,pdf,'k-' );
     plot( [xFo,xFo],  [0,1.1*pdfXmax],'r--' );
     plot( [-xFo,-xFo],[0,1.1*pdfXmax],'r--' );
     [~,indx] = sort( abs(xF - xP1t) );
     indx1 = indx(1);
     [~,indx] = sort( abs(xF - xP2t) );
     indx2 = indx(1);
     [~,indx] = sort( abs(xF - xP3t) );
     indx3 = indx(1);
     scatter(xP1t,pdf(indx1),60,'b','fill');
     scatter(xP2t,pdf(indx2),60,'r','fill');
     scatter(xP3t,pdf(indx3),60,'m','fill');
     xlim( [-Fmax,Fmax] );
% -------------------------------------------------------- pdf y-direction
     posCOSy = Ys0*cos(T-Ty0);
     negSINy =-Ys0*sin(T-Ty0);
     y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
     dYs = Ys - posCOSy;
     imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
     psi = ampFy*exp( -0.5*dYs.^2 + 1i*imYpart );
     pdf = abs(psi).^2;
     subplot(2,2,1);
     hold on
     plot(pdf,yF,'k-' );
     plot( [0,1.1*pdfYmax], [yFo,yFo],'r--' );
     plot( [0,1.1*pdfYmax], [-yFo,-yFo],'r--' );
     [~,indx] = sort( abs(yF - yP1t) );
     indx1 = indx(1);
     [~,indx] = sort( abs(yF - yP2t) );
     indx2 = indx(1);
     [~,indx] = sort( abs(yF - yP3t) );
     indx3 = indx(1);
     scatter(pdf(indx1),yP1t,60,'b','fill');
     scatter(pdf(indx2),yP2t,60,'r','fill');
     scatter(pdf(indx3),yP3t,60,'m','fill');
     ylim( [-Fmax,Fmax] );
%} 
     
     %pause(0.02)
     %disp(k);
     end
% ----------------------------------------------------------- plot results
   kMax = k;
   figure(1);
   clf; 
   hold on
   %-------------------------------------------------------
   plot(estimatedP1x(1:kMax),estimatedP1y(1:kMax),'b');
   plot(estimatedP2x(1:kMax),estimatedP2y(1:kMax),'r');
   plot(estimatedP3x(1:kMax),estimatedP3y(1:kMax),'m');
   %-------------------------------------------------------
   plot(xFpeak,yFpeak,'k','linewidth',1.3);
   plot(xP2,yP2,'color',[0.6,0.6,0.6],'linewidth',1.5);
   plot(xP3,yP3,'color',[0.6,0.6,0.6],'linewidth',1.5);
   xlim( [-Fmax,Fmax] );
   ylim( [-Fmax,Fmax] );
   title(['exact vs. numerical trajectory comparison  #Laps = ', ...
          num2str(nLaps)]);
   xlabel('x position  (Ang)');
   ylabel('y position  (Ang)');
   legend('est1','est2','est3','exct1','exct2','exct3', ...
          'location','southeast');
   % --------------------------------------------------
   dError1x = estimatedP1x(1:kMax) - exactP1x(1:kMax);
   dError2x = estimatedP2x(1:kMax) - exactP2x(1:kMax);
   dError3x = estimatedP3x(1:kMax) - exactP3x(1:kMax);
   % --------------------------------------------------
   dError1y = estimatedP1y(1:kMax) - exactP1y(1:kMax);
   dError2y = estimatedP2y(1:kMax) - exactP2y(1:kMax);
   dError3y = estimatedP3y(1:kMax) - exactP3y(1:kMax);
   % --------------------------------------------------
   figure(2);
   clf;
   hold on;
   plot( 1:kMax, dError1x );
   plot( 1:kMax, dError1y );
   plot( 1:kMax, dError2x );
   plot( 1:kMax, dError2y );
   plot( 1:kMax, dError3x );
   plot( 1:kMax, dError3y ); 
   title('Errors: Raw differences in position');
   ylabel('estimated - exact  (Ang)');
   xlabel(['time steps  dt=',num2str(dtFlow),' fs']);
   legend('E1x','E1y','E2x','E2y','E3x','E3y','location','southeast');
   end
%%                   check simultaneous forward/reverse elliptical motions
% REMARK: The forwrd & reverse wavefunctions are not linearly superimposed
   if( checkForwardReversed )
   ampF = c1*(const/pi)^(-1/4);
%  look at waveData{1}
% --------------------- by symmetry forward & reverse setups are identical
   xFo = waveData{1}.a_xFo*sigmaF;
   yFo = waveData{1}.a_yFo*sigmaF;
   Txo = waveData{1}.f_xFo*Tperiod;
   Tyo = waveData{1}.f_yFo*Tperiod;
   xFmax = xFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   yFmax = yFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   Fmax = max(xFmax,yFmax);
   xFmin = -xFmax;
   yFmin = -yFmax;
   xF = xFmin:ds:xFmax;
   yF = yFmin:ds:yFmax;
   Nx = length(xF);
   Ny = length(yF);
   tArray = 0:dtFlow:(Tperiod - dtFlow/4);
   Nt = length(tArray);
% -------------------------------------------- augment summary information
   disp(['    Nx = ',num2str(Nx)]);
   disp(['    Ny = ',num2str(Ny)]);
   disp(['    Nt = ',num2str(Nt)]);
% ------------------------------------------ work elliptical wave function
   Xs = c1*xF;
   Xs0 = c1*xFo;
   Tx0 = omega*Txo;
% ------------------- 
   Ys = c1*yF;
   Ys0 = c1*yFo;
   Ty0 = omega*Tyo;     
% ----------------------------------------------- normalize wave functions
   totalProb = zeros(1,Nt);
   k = 0;
     for j=1:5:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSx = Xs0*cos(T-Tx0);
     negSINx =-Xs0*sin(T-Tx0);
     x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
     dXs = Xs - posCOSx;
     imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
     psi = ampF*exp( -0.5*dXs.^2 + 1i*imXpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(temp);
   ampFx = temp*ampF;
   disp(['scale factor for X-normalization constant = ',num2str(temp)]);
   pdfXmax = max(pdf)*temp^2;
% ---------------------------------- normalize wavefunction in y-direction
   totalProb = 0*totalProb;
   k = 0;
     for j=1:5:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSy = Ys0*cos(T-Ty0);
     negSINy =-Ys0*sin(T-Ty0);
     y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
     dYs = Ys - posCOSy;
     imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
     psi = ampF*exp( -0.5*dYs.^2 + 1i*imYpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(temp);
   ampFy = temp*ampF;
   disp(['scale factor for Y-normalization constant = ',num2str(temp)]);
   pdfYmax = max(pdf)*temp^2; 
% ----------------------------- setup for forward and reverse trajectories
   Txo = zeros(1,2);
   Tyo = zeros(1,2);
   Txo(1) = waveData{1}.f_xFo*Tperiod;
   Tyo(1) = waveData{1}.f_yFo*Tperiod;
   Txo(2) = waveData{2}.f_xFo*Tperiod;
   Tyo(2) = waveData{2}.f_yFo*Tperiod;
% --------------------------------------------------- calculate trajectory
   shift_xFo = min(0.5*xFo,1.5*sigmaF);
   shift_yFo = min(0.5*yFo,1.5*sigmaF);
   xPo1 = xFo;
   yPo1 = yFo;
   xPo2 = xFo + shift_xFo;
   yPo2 = yFo + shift_yFo;
   xPo3 = xFo - shift_xFo;
   yPo3 = yFo - shift_yFo;
% % % %    xPo1 = 1.00*xFo;
% % % %    yPo1 = 1.00*yFo;
% % % %    xPo2 = 1.50*xFo;
% % % %    yPo2 = 1.50*yFo;
% % % %    xPo3 = 0.50*xFo;
% % % %    yPo3 = 0.50*yFo;
   %xP1 = xFo*cos( omega*(tArray - Txo) ) + xPo1 - xFo;
   %yP1 = yFo*cos( omega*(tArray - Tyo) ) + yPo1 - yFo;
   xP2 = xFo*cos( omega*( tArray - Txo(1) ) ) + xPo2 - xFo;
   yP2 = yFo*cos( omega*( tArray - Tyo(1) ) ) + yPo2 - yFo;
   xP3 = xFo*cos( omega*( tArray - Txo(1) ) ) + xPo3 - xFo;
   yP3 = yFo*cos( omega*( tArray - Tyo(1) ) ) + yPo3 - yFo;
   xFpeak = xFo*cos( omega*( tArray - Txo(1) ) );
   yFpeak = yFo*cos( omega*( tArray - Tyo(1) ) );
% --------------------------------------------- plot |psi|^2 decomposition 
     for j=1:10:Nt
     figure(1);
     clf; 
     axis equal
     t = tArray(j);
     T = t*omega;
% -------------------------------------------- show trajectory of particle
     subplot(2,2,2);
     hold on
     plot(xFpeak,yFpeak,'k','linewidth',1.3);
     plot(xP2,yP2,'color',[0.6,0.6,0.6],'linewidth',1.5);
     plot(xP3,yP3,'color',[0.6,0.6,0.6],'linewidth',1.5);
     xlim( [-Fmax,Fmax] );
     ylim( [-Fmax,Fmax] );
       for jj=1:2
       xP1t = xFo*cos( omega*( t - Txo(jj) ) );
       yP1t = yFo*cos( omega*( t - Tyo(jj) ) );
       xP2t = xFo*cos( omega*( t - Txo(jj) ) ) + xPo2 - xFo;
       yP2t = yFo*cos( omega*( t - Tyo(jj) ) ) + yPo2 - yFo;
       xP3t = xFo*cos( omega*( t - Txo(jj) ) ) + xPo3 - xFo;
       yP3t = yFo*cos( omega*( t - Tyo(jj) ) ) + yPo3 - yFo;
         if( jj == 1 ) % => forward
         scatter(xP1t,yP1t,60,'b','fill');
         scatter(xP2t,yP2t,60,'r','fill');
         scatter(xP3t,yP3t,60,'m','fill');
         else          % => reverse
         scatter(xP1t,yP1t,60,'b','linewidth',2);
         scatter(xP2t,yP2t,60,'r','linewidth',2);
         scatter(xP3t,yP3t,60,'m','linewidth',2);
         end
       end
% -------------------------------------------------------- pdf x-direction
     subplot(2,2,4);
     hold on
       for jj=2:-1:1
       xP1t = xFo*cos( omega*( t - Txo(jj) ) );
       xP2t = xFo*cos( omega*( t - Txo(jj) ) ) + xPo2 - xFo;
       xP3t = xFo*cos( omega*( t - Txo(jj) ) ) + xPo3 - xFo;
       Tx0 = omega*Txo(jj);
       posCOSx = Xs0*cos(T-Tx0);
       negSINx =-Xs0*sin(T-Tx0);
       x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
       dXs = Xs - posCOSx;
       imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
       psi = ampFx*exp( -0.5*dXs.^2 + 1i*imXpart );
       pdf = abs(psi).^2;
         if( jj == 1 )     % => forward
         plot(xF,pdf,'k-' );
         else              % => reverse
         L = ( pdf > 0.00005 );
         temp_xF = xF(L);
         temp_pdf = pdf(L);
         mmm = length(temp_pdf);
         sss = 1:4:mmm;
         fill(temp_xF(sss),temp_pdf(sss),'y-');
         end
       
       plot( [xFo,xFo],  [0,1.1*pdfXmax],'r--' );
       plot( [-xFo,-xFo],[0,1.1*pdfXmax],'r--' );
       [~,indx] = sort( abs(xF - xP1t) );
       indx1 = indx(1);
       [~,indx] = sort( abs(xF - xP2t) );
       indx2 = indx(1);
       [~,indx] = sort( abs(xF - xP3t) );
       indx3 = indx(1);
         if( jj == 1 )    % => forward
         scatter(xP1t,pdf(indx1),60,'b','fill');
         scatter(xP2t,pdf(indx2),60,'r','fill');
         scatter(xP3t,pdf(indx3),60,'m','fill');
         end
       xlim( [-Fmax,Fmax] );
       end
% -------------------------------------------------------- pdf y-direction
     subplot(2,2,1);
     hold on
       for jj=2:-1:1
       yP1t = yFo*cos( omega*( t - Tyo(jj) ) );
       yP2t = yFo*cos( omega*( t - Tyo(jj) ) ) + yPo2 - yFo;
       yP3t = yFo*cos( omega*( t - Tyo(jj) ) ) + yPo3 - yFo;
       Ty0 = omega*Tyo(jj);
       posCOSy = Ys0*cos(T-Ty0);
       negSINy =-Ys0*sin(T-Ty0);
       y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
       dYs = Ys - posCOSy;
       imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
       psi = ampFy*exp( -0.5*dYs.^2 + 1i*imYpart );
       pdf = abs(psi).^2;
         if( jj == 1 )     % => forward
         plot(pdf,yF,'k-' );
         else              % => reverse
         L = ( pdf > 0.00005 );
         temp_yF = yF(L);
         temp_pdf = pdf(L);
         mmm = length(temp_pdf);
         sss = 1:4:mmm;
         fill(temp_pdf(sss),temp_yF(sss),'y' );
         end
       plot( [0,1.1*pdfYmax], [yFo,yFo],'r--' );
       plot( [0,1.1*pdfYmax], [-yFo,-yFo],'r--' );
       [~,indx] = sort( abs(yF - yP1t) );
       indx1 = indx(1);
       [~,indx] = sort( abs(yF - yP2t) );
       indx2 = indx(1);
       [~,indx] = sort( abs(yF - yP3t) );
       indx3 = indx(1);
         if( jj == 1 )    % => forward
         scatter(pdf(indx1),yP1t,60,'b','fill');
         scatter(pdf(indx2),yP2t,60,'r','fill');
         scatter(pdf(indx3),yP3t,60,'m','fill');
         end
       ylim( [-Fmax,Fmax] );
       end
     pause(0.15)
     end
   end    
%%                     follow trajectores for linear supposition
% REMARK: linear superposition of countr-clckwise & clockwse wavefunctions
   if( exact2DInterference )
   dj = 2;                                       % any even number will do
   ampF = c1*(const/pi)^(-1/4);
%  look at waveData{1}
% --------------------- by symmetry forward & reverse setups are identical
   xFo = waveData{1}.a_xFo*sigmaF;
   yFo = waveData{1}.a_yFo*sigmaF;
   Txo = waveData{1}.f_xFo*Tperiod;
   Tyo = waveData{1}.f_yFo*Tperiod;
   xFmax = xFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   yFmax = yFo + spaceBuffer*sigmaF;       % add spaceBuffer at boundaries
   Fmax = max(xFmax,yFmax);
   xFmin = -xFmax;
   yFmin = -yFmax;
   xF = xFmin:ds:xFmax;
   yF = yFmin:ds:yFmax;
   Nx = length(xF);
   Ny = length(yF);
   tFinal = nLaps*Tperiod + 2*dj*dtFlow;
   tArray = 0:dtFlow:tFinal;
   Nt = length(tArray);
% -------------------------------------------- augment summary information
   disp(['    Nx = ',num2str(Nx)]);
   disp(['    Ny = ',num2str(Ny)]);
   disp(['    Nt = ',num2str(Nt)]);
   xMIN = min(xF);
   yMIN = min(yF);
% ------------------------------------------ work elliptical wave function
   Xs = c1*xF;
   Xs0 = c1*xFo;
   Tx0 = omega*Txo;
% ------------------- 
   Ys = c1*yF;
   Ys0 = c1*yFo;
   Ty0 = omega*Tyo;     
% ----------------------------------------------- normalize wave functions
   totalProb = zeros(1,Nt);
   k = 0;
     for j=1:25:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSx = Xs0*cos(T-Tx0);
     negSINx =-Xs0*sin(T-Tx0);
     x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
     dXs = Xs - posCOSx;
     imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
     psi = ampF*exp( -0.5*dXs.^2 + 1i*imXpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(2*temp);           % must decrease amplitude by 1/sqrt(2)
   ampFx = temp*ampF;
   disp(['scale factor for X-normalization constant = ',num2str(temp)]);
   pdfXmax = max(pdf)*temp^2;
% ---------------------------------- normalize wavefunction in y-direction
   totalProb = 0*totalProb;
   k = 0;
     for j=1:25:Nt
     t = tArray(j);
     T = t*omega;
     k = k + 1;
% ------------------------------------------------- set unitless variables
     posCOSy = Ys0*cos(T-Ty0);
     negSINy =-Ys0*sin(T-Ty0);
     y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
     dYs = Ys - posCOSy;
     imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
     psi = ampF*exp( -0.5*dYs.^2 + 1i*imYpart );
     pdf = abs(psi).^2;
     totalProb(k) = sum(pdf)*ds;
     end
   kMax = k;
   temp = mean( totalProb(1:kMax) );
   totalProb = totalProb/temp;
   temp = 1/sqrt(2*temp);           % must decrease amplitude by 1/sqrt(2)
   ampFy = temp*ampF;
   disp(['scale factor for Y-normalization constant = ',num2str(temp)]);
   pdfYmax = max(pdf)*temp^2; 
% ----------------------------- setup for forward and reverse trajectories
   Txo = zeros(1,2);
   Tyo = zeros(1,2);
   Txo(1) = waveData{1}.f_xFo*Tperiod;
   Tyo(1) = waveData{1}.f_yFo*Tperiod;
   Txo(2) = waveData{2}.f_xFo*Tperiod;
   Tyo(2) = waveData{2}.f_yFo*Tperiod;
% --------------------------------------------------- calculate trajectory
   shift_xFo = min(0.5*xFo,1.5*sigmaF);
   shift_yFo = min(0.5*yFo,1.5*sigmaF);
   xPo1 = xFo;
   yPo1 = yFo;
   xPo2 = xFo + shift_xFo;
   yPo2 = yFo + shift_yFo;
   xPo3 = xFo - shift_xFo;
   yPo3 = yFo - shift_yFo;
% % % % %    xPo1 = 1.00*xFo;
% % % % %    yPo1 = 1.00*yFo;
% % % % %    xPo2 = 1.50*xFo;
% % % % %    yPo2 = 1.50*yFo;
% % % % %    xPo3 = 0.50*xFo;
% % % % %    yPo3 = 0.50*yFo;
% ------------------------------------------------------------------------
   xP1Locate = xFo*cos( omega*( -dj*dtFlow/2 - Txo(1) ) ) + xPo1 - xFo;
   yP1Locate = yFo*cos( omega*( -dj*dtFlow/2 - Tyo(1) ) ) + yPo1 - yFo;
   xP2Locate = xFo*cos( omega*( -dj*dtFlow/2 - Txo(1) ) ) + xPo2 - xFo;
   yP2Locate = yFo*cos( omega*( -dj*dtFlow/2 - Tyo(1) ) ) + yPo2 - yFo;
   xP3Locate = xFo*cos( omega*( -dj*dtFlow/2 - Txo(1) ) ) + xPo3 - xFo;
   yP3Locate = yFo*cos( omega*( -dj*dtFlow/2 - Tyo(1) ) ) + yPo3 - yFo;
                                                 % ^ -> no ghost particles
% ------------------------------------------------------------------------
   recordP1x = zeros(1,Nt);
   recordP1y = zeros(1,Nt);
   recordP2x = zeros(1,Nt);
   recordP2y = zeros(1,Nt);
   recordP3x = zeros(1,Nt);
   recordP3y = zeros(1,Nt); 
% ------------------------- apply Bohmain procedure for quantum trajectory
   dShiftNeg = -3*ds;
   dShiftPos = 3.1*ds;
   iSpan = 1:7;
   xFit = ds*(iSpan - 1);
   mSpan = length(iSpan);
   k = 0;
     for j=1:dj:Nt
     figure(1);
     clf; 
     axis equal
     t = tArray(j);
     T = t*omega; 
% ------------------------------------------ integrate equations of motion
       for nR=1:Nrepeat
       % -------------------------------------------------------
       x1 = (xP1Locate + dShiftNeg):ds:(xP1Locate + dShiftPos);
       y1 = (yP1Locate + dShiftNeg):ds:(yP1Locate + dShiftPos);
       % -------------------------------------------------------
       x2 = (xP2Locate + dShiftNeg):ds:(xP2Locate + dShiftPos);
       y2 = (yP2Locate + dShiftNeg):ds:(yP2Locate + dShiftPos);
       % -------------------------------------------------------
       x3 = (xP3Locate + dShiftNeg):ds:(xP3Locate + dShiftPos);
       y3 = (yP3Locate + dShiftNeg):ds:(yP3Locate + dShiftPos);
       % ------------------------------------------------------- scale x,y
       Xs1 = c1*x1;
       Xs2 = c1*x2;
       Xs3 = c1*x3;
       Ys1 = c1*y1;
       Ys2 = c1*y2;
       Ys3 = c1*y3;
% ------------------------------------------------ calculate wavefunctions
       psi1x = zeros(1,mSpan);
       psi2x = zeros(1,mSpan);
       psi3x = zeros(1,mSpan);
       psi1y = zeros(1,mSpan);
       psi2y = zeros(1,mSpan);
       psi3y = zeros(1,mSpan);
         for jj=1:2
         Tx0 = omega*Txo(jj);
         Ty0 = omega*Tyo(jj);
         posCOSx = Xs0*cos(T-Tx0);
         negSINx =-Xs0*sin(T-Tx0);
         x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
% --------------------------------------------------------------
         dXs1 = Xs1 - posCOSx;
         imXpart = Xs1*negSINx - x_phase - 0.5*(T - Tx0);
         psi1x = psi1x + ampFx*exp( -0.5*dXs1.^2 + 1i*imXpart );
% --------------------------------------------------------------
         dXs2 = Xs2 - posCOSx;
         imXpart = Xs2*negSINx - x_phase - 0.5*(T - Tx0);
         psi2x = psi2x + ampFx*exp( -0.5*dXs2.^2 + 1i*imXpart );
% --------------------------------------------------------------
         dXs3 = Xs3 - posCOSx;
         imXpart = Xs3*negSINx - x_phase - 0.5*(T - Tx0);
         psi3x = psi3x + ampFx*exp( -0.5*dXs3.^2 + 1i*imXpart );
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         posCOSy = Ys0*cos(T-Ty0);
         negSINy =-Ys0*sin(T-Ty0);
         y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
% --------------------------------------------------------------
         dYs1 = Ys1 - posCOSy;
         imYpart = Ys1*negSINy - y_phase - 0.5*(T - Ty0);
         psi1y = psi1y + ampFy*exp( -0.5*dYs1.^2 + 1i*imYpart );
% --------------------------------------------------------------
         dYs2 = Ys2 - posCOSy;
         imYpart = Ys2*negSINy - y_phase - 0.5*(T - Ty0);
         psi2y = psi2y + ampFy*exp( -0.5*dYs2.^2 + 1i*imYpart );
% --------------------------------------------------------------
         dYs3 = Ys3 - posCOSy;
         imYpart = Ys3*negSINy - y_phase - 0.5*(T - Ty0);
         psi3y = psi3y + ampFy*exp( -0.5*dYs3.^2 + 1i*imYpart );
         end
% ============================================================= get action
       qPhase = atan2( imag(psi1x), real(psi1x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP1Locate )/mass;
       xP1Locate = xP1Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi2x), real(psi2x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP2Locate )/mass;
       xP2Locate = xP2Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi3x), real(psi3x) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*xP3Locate )/mass;
       xP3Locate = xP3Locate + vel*dt;
% ----------------------------------------------------------
       qPhase = atan2( imag(psi1y), real(psi1y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP1Locate )/mass;
       yP1Locate = yP1Locate + vel*dt;
% ----------------------------------------------------------       
       qPhase = atan2( imag(psi2y), real(psi2y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP2Locate )/mass;
       yP2Locate = yP2Locate + vel*dt;
% ----------------------------------------------------------       
       qPhase = atan2( imag(psi3y), real(psi3y) );
       qPhase = unwrap(qPhase);                              % track phase
% ------------------------------------------------------- smoothing method
       aCoef = polyfit(xFit,qPhase,2);
       vel = hbar*( aCoef(2) + 2*aCoef(1)*yP3Locate )/mass;
       yP3Locate = yP3Locate + vel*dt;
       end 
% ------------------------------------------- record particle trajectories
     k = k + 1;
     recordP1x(k) = xP1Locate;
     recordP1y(k) = yP1Locate;
     recordP2x(k) = xP2Locate;
     recordP2y(k) = yP2Locate;
     recordP3x(k) = xP3Locate;
     recordP3y(k) = yP3Locate;
% -------------------------------------------------------- pdf x-direction
     psi = zeros( size(Xs) );
       for jj=1:2
       Tx0 = omega*Txo(jj);
       posCOSx = Xs0*cos(T-Tx0);
       negSINx =-Xs0*sin(T-Tx0);
       x_phase = (-1/4)*(Xs0^2)*sin(2*(T-Tx0));
       dXs = Xs - posCOSx;
       imXpart = Xs*negSINx - x_phase - 0.5*(T - Tx0);
       psi = psi + ampFx*exp( -0.5*dXs.^2 + 1i*imXpart );
       end     
% ----------------------------------------------------------- plot |psi|^2     
     pdf = abs(psi).^2;
     netProb = sum(pdf)*ds;
     subplot(2,2,4);
     hold on
     plot(xF,pdf,'k-' );
     plot( [xFo,xFo],  [0,1.1*pdfXmax],'r--' );
     plot( [-xFo,-xFo],[0,1.1*pdfXmax],'r--' );
     xlim( [-Fmax,Fmax] );
     title(['total area = ',num2str(netProb)]);
% -------------------------------------------------------- pdf y-direction
     psi = zeros( size(Ys) ); 
       for jj=1:2
       Ty0 = omega*Tyo(jj);
       posCOSy = Ys0*cos(T-Ty0);
       negSINy =-Ys0*sin(T-Ty0);
       y_phase = (-1/4)*(Ys0^2)*sin(2*(T-Ty0));
       dYs = Ys - posCOSy;
       imYpart = Ys*negSINy - y_phase - 0.5*(T - Ty0);
       psi = psi + ampFy*exp( -0.5*dYs.^2 + 1i*imYpart );
       end
% ----------------------------------------------------------- plot |psi|^2     
     pdf = abs(psi).^2;
     netProb = sum(pdf)*ds;
     subplot(2,2,1);
     hold on
     plot(pdf,yF,'k-' );
     %disp( [pdf; yF] );
     plot( [0,1.1*pdfYmax], [yFo,yFo],'r--' );
     plot( [0,1.1*pdfYmax], [-yFo,-yFo],'r--' );
     ylim( [-Fmax,Fmax] );
     title(['total area = ',num2str(netProb)]);
% --------------------------------------------- plot particle trajectories
     subplot(2,2,2);
     hold on
     plot( recordP1x(1:k), recordP1y(1:k),'b','linewidth',1.6);
     plot( recordP2x(1:k), recordP2y(1:k),'r','linewidth',1.6);
     plot( recordP3x(1:k), recordP3y(1:k),'m','linewidth',1.6);
     xlim( [-Fmax,Fmax] );
     ylim( [-Fmax,Fmax] );
     mLap = floor( j*dtFlow/Tperiod );
     title(['trajectory  #Laps = ', ...
            num2str(mLap)]); 
     xlabel('x-coordinate');
     ylabel('y-coordinate');
     %pause
     end
   end    
