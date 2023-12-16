% This script simulates one particle in a 1D box in two parts.
%
% Part I: We set the boundary conditions based on a particle in the box in
%         thermal equilibrium at a specified temperature. The box size is
% fixed. The evolution of the system is based on the Schrodinger equation
% assuming that the system is closed. This means the parameters of the
% wavefunction are fixed and are constant. To characterize the dynamics of
% a closed system, version 02 of the code allows the user to start each 
% trajectory at randomized position or to have all particle trajectories
% start at the same location. In both cases, the intial phase angles are]
% randomized to different values. Even for the case that all trajectories 
% start at the same location, the trajectories will diverge due to having
% different starting phases. In this way, it is possible to study how the
% trajectories diverge from one another based on differences in phase that
% is initially placed on the expansion coefficients of the energy 
% eigenfunctions. However, in other cases, if one is trying to see as many
% different phase space orbits as possible, it makes sense to start the 
% trajectories at randomized locations. 
%
% Part II: We acknowledge the system is not closed. As an open quantum
%          system that is in thermal equilibrium with its environment at a
% constant temperature, energy will fluctuate within the system, and thus
% energy is not conserved. All particle trajectories will start at the
% same location. In addition, the same initial randomized phase will be
% placed on the expansion coeficients of the energy eigenfunctions. Due to
% the system being open, the solutions to the time dependent Schrodinger 
% equation are dynamically altered by the environment as a stochastic 
% process. As an open system it is acknowledged that the length of the 1D
% box will fluctuate. A model assuming a constant length is unrealistic. 
% As such, there are three modes of diffusive processes being modeled for
% the open system, listed as: (1) the length of the box fluctuates around
% a well-defined mean length (given by ax) since it is impossible to keep
% the walls of the box perfectly static. (2) The magnitude of expansion 
% coefficients fluctuate because temperature is causing energy exchange, 
% which affects the weighting of energy states, and (3) randomization of 
% phase on the expansion coefficients occurs due to the influence of the
% environment. The user will define the diffusion coefficients for each of
% these distinct diffusive processes. A specific diffusion processes can
% be turned off by setting the respective diffusion coefficient to zero.
% This will particle trajectories to start identically and as time passes,
% watch the trajectories diverge. Lyapunov exponents could be calculated
% using these sycronized trajectories. 
%
% INPUT
% Many model and user defined parameters.
%
% PROCESS
% Calculation of the Bohmian equations of motion based on the momentum
% operator expression (or current density). This is in contrast to a 
% direct F = ma approach. 
%
% OUTPUT
% A log file. 
% A trajectory file. 
%
% A Bohmian analysis program needs to be used on the output files that are
% generated from this program. 
%
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
disp('   ');
%disp('mass of electron in Atomic Mass Units (AMU) = 0.000548579909 AMU');
%disp('mass of  proton  in Atomic Mass Units (AMU) = 1.00727647 AMU');
%disp('mass of neutron  in Atomic Mass Units (AMU) = 1.008665 AMU');
%massAMU = input('what mass will be used? (AMU): ');
massAMU = 0.000548579909;
mass = massAMU*103.6484;                        % in units of evfs^2/Ang^2
  if( mass < 1.0e-25 )
  error('mass cannot be negative');
  end
%ax = input('Enter box size (in Ang): ');
ax = 100;
  if( ax < 1 )
  error('Box size must be greater than 1 Angstrom');
  end
%oneStartNum=input('Enter 0,1 for {many,single} initial positions: ');
oneStartNum = 1;
%iParity =   input('Enter 0,1,2 for {no,odd,even} parity symmetry: ');
iParity =  0;
%systemType= input('Enter {-1,1} to model a {closed, open} system: ');
systemType = 1;
%nSamples =  input('Enter the number of trajectories to calculate: ');
nSamples = 10;
%T =         input('Enter the temperature of the box  (in Kelvin): ');
T = 300;
%phaseAmpDeg=input('Enter initial random phase angle spread (deg): ');
phaseAmpDeg = 0;
%seed1 =     input('Enter SEED1 to setup an initial quantum state: ');
seed1 = 12345;
  if( systemType > 0 )  % => open system
  disp(['-------------------------------------------------------', ...
        '------------ open systems'])
  %seed2 = input('Enter SEED2 to simulate the diffusive process: ');
  seed2 = 4321;
  %Dp =    input('Enter Dp in (degrees^2)/ps    for phase angle: ');
  Dp = 0.01;
  %Da =    input('Enter Da in (percent^2)/ps  on prob-amplitude: ');
  Da = 0.01;
  %DL =    input('Enter DL in (percent^2)/ps on box side Length: ');
  DL = 0;
  else
  Dp = 0;
  Da = 0;
  DL = 0;
  end
% ---------------------------------- when desired: enforce parity symmetry
%  no-symmetry        odd-parity         even-parity
% nx0=1; dnx=1;      nx0=1; dnx=2;      nx0=2; dnx=2;
  switch iParity
    case 0
    nx0=1; 
    dnx=1;
    case 1
    nx0=1; 
    dnx=2;
    case 2
    nx0=2; 
    dnx=2;
    otherwise
    error('parity type is not specified!');
  end
%%                                              set type of quantum system
  if( systemType > 0 ) % => open system
  openSystem = true;
  else
  openSystem = false;
  end
%%                                                             start setup
rng(seed1); 
d = hbar/mass;
dd =((hbar^2 * pi^2)/(2*mass));
d1 = (hbar*pi^2)/(2*mass);           % d1 = dd/hbar = (hbar*pi^2)/(2*mass)
d2 = (hbar^2)/(2*mass);
vRMS1d = sqrt((kb*T/mass));                % in 1D: classical RMS velocity
vmax = 10*vRMS1d;
% --------------------------------------------- estimate maximum quantum #
imax = round((mass*vmax*ax)/(pi*hbar)); 
imax = max(imax,3);
xShift = 0.05*ax;
xLow = 0 - xShift;
xBig = ax + xShift;
%%                                                           sanity checks
  if( (nSamples < 1) || (nSamples > 1000))
  disp('nSamples must be at least 1 and not greater than 1000');
  end
  if( vmax > 300 )
  disp('lower the temperature: Relativistic effects are appearing.')
  end
  if( oneStartNum == 1 )
  oneStart = true;            % => a single starting position
  else
  oneStart = false;           % => many initial positions based on |psi|^2
  end
%%                                                          get Cnm and Qx
kT = kb*T;
qModes = nx0:dnx:imax;                    % list of quantum modes (qModes)
Estates = dd*(qModes/ax).^2;
Eground = dd*(nx0/ax)^2;
dEstates = Estates - Eground;
sWeights = exp( -dEstates/kT );                      % statistical weights
% ------------------------------------------------------------------------
totalWeight = sum(sWeights);
prob = sWeights/totalWeight;
Cprob = cumsum(prob);
L = ( Cprob < 0.9999995 );
nmodes = max(sum(L),3);                      % must track at least 3 modes
qModes = qModes(1:nmodes);
Estates = dd*(qModes/ax).^2;
omega = Estates/hbar;
rn = sqrt( prob(1:nmodes) );                          % note: rn = abs(Cn)
theta = 2*pi*rand(1,nmodes);             % => fix the initial random phase
%Cn = rn.*exp(1i*phi);                                       % apply phase
Qx = (pi*qModes)/ax;
% ---------------------------------------------------------- saniety check
KEclassical = 0.5*mass*vRMS1d^2;
r2total = sum( rn.^2 );
  if( abs(r2total - 1) > 0.000001 )
  %disp('Cn coeficients were renormalized!');
  rn = rn/sqrt(r2total);
  r2total = sum( rn.^2 );
    if( abs(r2total - 1) > 0.000001 )
    disp('Cn coeficients cannot be normalized!')
    end
  end
aveE0 = d2*sum( (rn.*Qx).^2 );
percentError = 100*( KEclassical - aveE0)/aveE0;
% -------------------------------- prevent negative diffusion coefficients
  if( Dp < 1.0-35 )   % for phase
  Dp = 0;
  end
  if( Da < 1.0-35 )   % for amplitude
  Da = 0;
  end
  if( DL < 1.0e-35 )  % for side of box
  DL = 0;
  end
% ----------------------------------------------------------- give warning
  if( nmodes > 5000 )
  disp(['nmodes = ',num2str(nmodes)]);
  disp('Warning: Number of modes is > 5000. compute time will be long.')
  end
Emax = max(Estates);
%%                                            scale the time of simulation
dEmove = aveE0 - Eground;
vx = sqrt(2*dEmove/mass) + 1.0e-10;
Tx = min(ax/vx,10000000);   % time scale for particle to cross side of box
mObs = 1000;
dtObsv = Tx/mObs;           % this program sets this very small as a check
%Nsweeps = input('Enter the number of sweeps to be observed: ');
Nsweeps = 100;
          disp(['  dt for observations = ',num2str(dtObsv),' fs ']);
nObservations = Nsweeps*mObs;
%dtRK4 = input('Enter the dt (fs) for solving dif-eq using RK4: ');
dtRK4 = dtObsv/4;
runTime = Nsweeps*Tx;
nRK4 = ceil(dtObsv/dtRK4);
dtRK4 = dtObsv/nRK4;
%%       define diffusion coefficients for each process for an open system
  if( systemType > 0 )
% ---------------------------------------------------- set diffusion rates
% Dp = sig0^2/ps => R2 = Dp*(1ps) = sig0^2 = (sig^2/dtRK4)*1000 
  sigDp = sqrt(dtRK4*Dp/1000);            % => sqrt( (dtRK4*sig0^2)/1000 )
  sigDa = sqrt(dtRK4*Da/1000); 
  sigDL = sqrt(dtRK4*DL/1000);
  else
  sigDp = 0;
  sigDa = 0;
  sigDL = 0;
  end
%%                                   create output log and data file names
strMass = num2str(massAMU);
mmm = length(strMass);
  for mmmm=1:mmm
     if( strMass(mmmm) == '.' )
     strMass(mmmm) = '_';
     end
  end
% -------------------------------
strDtRK4 = num2str(dtRK4);
mmm = length(strDtRK4);
  for mmmm=1:mmm
     if( strDtRK4(mmmm) == '.' )
     strDtRK4(mmmm) = '_';
     end
  end
% -------------------------------
  if( openSystem )
  strDp = num2str(Dp);
  mmm = length(strDp);
    for mmmm=1:mmm
      if( strDp(mmmm) == '.' )
      strDp(mmmm) = '_';
      end
    end
  %------------------------------
  strDa = num2str(Da);
  mmm = length(strDa);
    for mmmm=1:mmm
      if( strDa(mmmm) == '.' )
      strDa(mmmm) = '_';
      end
    end
  %------------------------------
  strPs = num2str(phaseAmpDeg);
  mmm = length(strPs);
    for mmmm=1:mmm
      if( strPs(mmmm) == '.' )
      strPs(mmmm) = '_';
      end
    end
  %------------------------------
  strDL = num2str(DL);
  mmm = length(strDL);
    for mmmm=1:mmm
      if( strDL(mmmm) == '.' )
      strDL(mmmm) = '_';
      end
    end
  baseName = ['1P1D',num2str(T),'Lx',num2str(ax),'s',num2str(seed1), ...
             'RK4dt',strDtRK4,'m',strMass,'Ps',num2str(phaseAmpDeg),'s',num2str(seed2),'Dp', ...
             strDp,'Da',strDa,'DL',strDL,'open'];
  else
  baseName = ['1P1D',num2str(T),'Lx',num2str(ax),'s',num2str(seed1), ...
             'RK4dt',strDtRK4,'m',num2str(massAMU),'Ps',num2str(phaseAmpDeg),'closed'];
  end 
%%                           minimum check on output folder and file names
folderName = '1P1D';
currentFolder = pwd;
testEnd = currentFolder(end-3:end);
  if( strcmp(testEnd,folderName) )         % => already in the data folder
  baseFileName = baseName;  %=> no pre-appending folder name is neccessary
  else                                 % => in directory above data folder
% ------------------------------------------------------- create directory
    if( ~isfolder(folderName) )  % folder does not exist: create a new one
    mkdir(folderName);
    end
    if( isunix )                        % now folder exists (empty or not)
    folderPrefix = '1P1D/';
    else
    folderPrefix = '1P1D\'; 
    end
  baseFileName = [folderPrefix,baseName];        % add folder path to name
  end
fileNameP_x = [baseFileName,'P_x.txt'];
fileNameV_x = [baseFileName,'V_x.txt'];
fileNameQpe = [baseFileName,'Qpe.txt'];
flagError = -1;
  if( isfile(fileNameP_x) )
  flagError = 1;
  end
  if( isfile(fileNameV_x) )
  flagError = 1;
  end
  if( isfile(fileNameQpe) )
  flagError = 1;
  end
  if( flagError > 0 )
  error('This simulation has been done previously. No auto-overwrite');
  end
%%                                       build system quantum wavefunction
bb = sqrt(2/ax);
qSystem = struct;
qSystem.nmodes = nmodes;
qSystem.qModes = qModes;
qSystem.d = d;
qSystem.dd = dd;
qSystem.d2 = d2;
qSystem.ax = ax;
qSystem.rn = rn;  % note: rn = abs(Cn) where Cn are expansion coefficients
qSystem.Qx = Qx;
qSystem.bb = bb;
qSystem.theta = theta;           % initial phase on expansion coefficients
qSystem.omega = omega;                         % frequencies for each mode
qSystem.openSystem = openSystem;       % true/false toggle for open system
qSystem.sigDp = sigDp;             % STD for accumulative Gaussian process
qSystem.sigDa = sigDa;               % STD for stationary Gaussian process
qSystem.sigDL = sigDL;               % STD for stationary Gaussian process
%%                          generate initial random phases for each sample
phaseAmp = phaseAmpDeg*pi/180;                % convert degrees to radians
phaseShift = cell(1,nSamples);
dPhi = zeros(1,nmodes);                            % simply to preallocate 
  if( systemType > 0 )  % => open system
    for s=1:nSamples
    phaseShift{s} = zeros(1,nmodes);
    end
  else                  % systemType < 0 => closed system
    for s=1:nSamples
    phaseShift{s} = phaseAmp*randn(1,nmodes);           % units in radians
    end
  end
%%                                        calculate equilibrium properties
nGrid = 50000;
dx = ax/nGrid;
x0 = dx/2;
xArray = x0:dx:ax;
iPlaces = length(xArray);
rhoArray = zeros(1,iPlaces);
xLocate = zeros(1,iPlaces);
aveVx = 0;
aveVx2 = 0;
aveQPE = 0;
tTryArray = 0:10:100;
reducedProb = 0;
  for tTry=tTryArray
    parfor i=1:iPlaces
    x = xArray(i);
    rhoArray(i) = getProbDensity1d(qSystem,x,tTry);
    xLocate(i) = x;
    end
% ------------------------------------------- calculate average properties
  prob = rhoArray*dx;
  totalProb = sum(prob);   % should be 1 if grid resolution is high enough
  prob = prob/totalProb;
% ------------------------------------- extract initial starting positions
    if( tTry == 0 )
    save_xLocate = xLocate;
    save_prob = prob;
    end
% ---------------------------------------------------------- sample places
  nOutliers = 0;
    parfor i=1:iPlaces
    x = xLocate(i);
    kx = getVelocity1d(dPhi,qSystem,x,tTry);
    quantumPE = getQuantumPE1P1D(dPhi,qSystem,x,tTry);
      if( (abs(quantumPE) < 100*Emax) && (kx < 2000) )  % removes outliers
      aveVx  = aveVx  + prob(i)*kx;
      aveVx2 = aveVx2 + prob(i)*kx*kx;
      aveQPE = aveQPE + prob(i)*quantumPE;
      reducedProb = reducedProb + prob(i);
      else                                            % ignore: an outlier
      nOutliers = nOutliers + 1;
      %disp(quantumPE);
      end
    end
  % disp(aveVx)
  end
aveVx2 = aveVx2/reducedProb;
meanKE = 0.5*mass*aveVx2;
aveQPE = aveQPE/reducedProb;
aveBohmEnergy = meanKE + aveQPE;   % note avePE = 0 classical contribution
percentError2 = 100*(aveBohmEnergy - aveE0)/aveE0;
%%                                             set up time of observations
tArray = 0:dtObsv:runTime;
runTime = tArray(end);
nPts = length(tArray);
%%                                          place test particles in system
  if( oneStart )
  xo = ax/sqrt(2); 
  xoArray = xo*ones(1,nSamples);                       % only one location
  else                            % let us randomize particle location too
  [~,indx] = sort(save_prob,'descend');
  mTemp = ceil(iPlaces/2);
  dmS = max(floor(mTemp/nSamples),1);
  mSmax = nSamples*dmS;
  xoArray = save_xLocate(indx(1:dmS:mSmax));
  end
%%                                                    summarize parameters
  if( systemType > 0 )
  systemTypeName = 'open system';
  else
  systemTypeName = 'closed system';
  end
disp('   ');
disp(['             base file name= ',baseName]);
disp(['  initial conditions: SEED1= ',num2str(seed1)]);
disp(['                       mass= ',num2str(mass),' evfs^2/Ang^2']);
disp(['                       mass= ',num2str(massAMU,'%.9f'),' AMU']);
disp(['                          T= ',num2str(T),' Kelvin']);
disp(['           classical vRMS1d= ',num2str(vRMS1d),' Ang/fs']);
disp(['                         ax= ',num2str(ax),' Ang']);
disp(['               sweep period= ',num2str(Tx),' fs']);
disp(['                    Nsweeps= ',num2str(Nsweeps)]);
disp(['                       Trun= ',num2str(runTime,'%.9f'),' fs']);
disp(['                     dtObsv= ',num2str(dtObsv,'%.9f'),' fs  ']);
disp(['                      dtRK4= ',num2str(dtRK4,'%.9f'),' fs']);
disp(['        dtObsv subdivisions= ',num2str(nRK4), ...
                                    '  => # of loops per observation']);
disp(['    total # of observations= ',num2str(nObservations)]);
if( oneStart )
disp(['the single initial position= ',num2str(xo),' Ang']);
else
disp([' use many initial positions: ','randomize over |psi|^2']);
end
disp(['  parity symmetry type: nx0= ',num2str(nx0), ...
                                    '  dnx= ',num2str(dnx)]);
disp(['        number of Cnm terms= ',num2str(nmodes)]);
disp(['maximum energy state: Emax = ',num2str(Emax),' eV']);
disp(['                   nSamples= ',num2str(nSamples)]);
disp('----------------------------------------------------------------');
disp(['            simulation type= ',systemTypeName]);
if( openSystem )
disp([' diffusive processes: SEED2= ',num2str(seed2)]);
disp([' initial random phase stdev= ',num2str(phaseAmpDeg),'  deg']);
disp(['  phase angle diffusion: Dp= ',num2str(Dp,'%.9f'), ...
                                                 '  (degree^2)/ps']);
disp(['   amplitude modulation: Da= ',num2str(Da,'%.9f'), ...
                                                 '  (percnt^2)/ps']);
disp(['   ax-length modulation: DL= ',num2str(DL,'%.9f'), ...
                                                 '  (percnt^2)/ps']);
end
disp(['                KEclassical= ',num2str(KEclassical,'%.9f'),' eV']);
disp(['     Copenhagen mean energy= ',num2str(aveE0,'%.9f'),' eV']);
disp(['                mean BohmKE= ',num2str(meanKE,'%.9f'),' eV']);
disp(['             mean quantumPE= ',num2str(aveQPE,'%.9f'),' eV']);
disp(['              aveBohmEnergy= ',num2str(aveBohmEnergy,'%.9f'), ...
                                                            ' eV']);
disp(['  number of energy outliers= ',num2str(nOutliers)]);
disp(['Copenhagen/classical %error= ',num2str(percentError)]);
disp(['  Copenhagen/Bohmian %error= ',num2str(percentError2)]);
%error('stop here for now');
%%                                                   set open system seed2
  if( openSystem )
  rng(seed2);                 % reset random numbers for diffusive process
  end
%%                             solve Bohmian equations of motion using RK4
cpu0 = cputime;
h  = dtRK4;
posX = zeros(nSamples,nPts);
velX = zeros(nSamples,nPts);
q_PE = zeros(nSamples,nPts);
% ------------------------------------------------ parallelize the samples
to = tArray(1);                     % Note: to never changes => a constant
  parfor s=1:nSamples              % use for or parfor for parallelization
  xo = xoArray(s);
  dPhi = phaseShift{s};
  vxo = getVelocity1d(dPhi,qSystem,xo,to);
  [pX,vX,qU] = propagateMotion1d(dPhi,qSystem,xo,vxo,tArray,h);
  posX(s,:) = pX;
  velX(s,:) = vX;
  q_PE(s,:) = qU;
  end
cpu1 = cputime;
timeCPU = cpu1 - cpu0;
disp(['                CPU seconds= ',num2str(timeCPU)]);
%%                                       write results to output data file
writematrix(posX,fileNameP_x);
writematrix(velX,fileNameV_x);
writematrix(q_PE,fileNameQpe);
%%                                        write results to output log file
fileNameLog = [baseFileName,'Log.txt'];
fid = fopen(fileNameLog,'w');
fprintf(fid,'%s \n',baseName);                  % 1st line is the basename
msg ='simulation of quantum trajectories for confined particles';
fprintf(fid,'%s \n',msg);
msg = ['------------------------------------------------------------', ...
       '----- summary'];
fprintf(fid,'%s \n',msg);
msg = ['  initial conditions: SEED1= ',num2str(seed1)];
fprintf(fid,'%s \n',msg);
msg = ['                       mass= ',num2str(mass),' evfs^2/Ang^2'];
fprintf(fid,'%s \n',msg);
msg = ['                       mass= ',num2str(massAMU,'%.9f'),' AMU'];
fprintf(fid,'%s \n',msg);
msg = ['                          T= ',num2str(T),' Kelvin'];
fprintf(fid,'%s \n',msg);
msg = ['           classical vRMS1d= ',num2str(vRMS1d),' Ang/fs'];
fprintf(fid,'%s \n',msg);
msg = ['                         ax= ',num2str(ax),' Ang'];
fprintf(fid,'%s \n',msg);
msg = ['               sweep period= ',num2str(Tx),' fs'];
fprintf(fid,'%s \n',msg);
msg = ['                    Nsweeps= ',num2str(Nsweeps)];
fprintf(fid,'%s \n',msg);
msg = ['                       Trun= ',num2str(runTime,'%.9f'),' fs'];
fprintf(fid,'%s \n',msg);
msg = ['                     dtObsv= ',num2str(dtObsv,'%.9f'),' fs  '];
fprintf(fid,'%s \n',msg);
msg = ['                      dtRK4= ',num2str(dtRK4,'%.9f'),' fs'];
fprintf(fid,'%s \n',msg);
msg = ['        dtObsv subdivisions= ',num2str(nRK4), ...
                                     '  => # of loops per observation'];
fprintf(fid,'%s \n',msg);
msg = ['    total # of observations= ',num2str(nObservations)];
fprintf(fid,'%s \n',msg);
if( oneStart )
msg = ['the single initial position= ',num2str(xo,'%.9f'),' Ang'];
fprintf(fid,'%s \n',msg);
else
msg = [' use many initial positions: ','randomize over |psi|^2'];
fprintf(fid,'%s \n',msg);
end
msg = ['  parity symmetry type: nx0= ',num2str(nx0), ...
                                     '  dnx= ',num2str(dnx)];
fprintf(fid,'%s \n',msg);
msg = ['        number of Cnm terms= ',num2str(nmodes)];
fprintf(fid,'%s \n',msg);
msg = ['maximum energy state: Emax = ',num2str(Emax),' eV'];
fprintf(fid,'%s \n',msg);
msg = ['                   nSamples= ',num2str(nSamples)];
fprintf(fid,'%s \n',msg);
msg = '----------------------------------------------------------------';
fprintf(fid,'%s \n',msg);
msg = ['            simulation type= ',systemTypeName];
fprintf(fid,'%s \n',msg);
if( openSystem )
msg = [' diffusive processes: SEED2= ',num2str(seed2)];
fprintf(fid,'%s \n',msg);
msg = [' initial random phase stdev= ',num2str(phaseAmpDeg),'  deg'];
fprintf(fid,'%s \n',msg);
msg = ['  phase angle diffusion: Dp= ',num2str(Dp,'%.9f'), ...
                                                  '  (degree^2)/ps'];
fprintf(fid,'%s \n',msg);
msg = ['   amplitude modulation: Da= ',num2str(Da,'%.9f'), ...
                                                  '  (percnt^2)/ps'];
fprintf(fid,'%s \n',msg);
msg = ['   ax-length modulation: DL= ',num2str(DL,'%.9f'), ...
                                                  '  (percnt^2)/ps'];
fprintf(fid,'%s \n',msg);
end
msg = ['                KEclassical= ',num2str(KEclassical,'%.9f'),' eV'];
fprintf(fid,'%s \n',msg);
msg = ['     Copenhagen mean energy= ',num2str(aveE0,'%.9f'),' eV'];
fprintf(fid,'%s \n',msg);
msg = ['                mean BohmKE= ',num2str(meanKE,'%.9f'),' eV'];
fprintf(fid,'%s \n',msg);
msg = ['             mean quantumPE= ',num2str(aveQPE,'%.9f'),' eV'];
fprintf(fid,'%s \n',msg);
msg = ['              aveBohmEnergy= ',num2str(aveBohmEnergy,'%.9f'), ...
                                                             ' eV'];
fprintf(fid,'%s \n',msg);
msg = ['  number of energy outliers= ',num2str(nOutliers)];
fprintf(fid,'%s \n',msg);
msg = ['Copenhagen/classical %error= ',num2str(percentError)];
fprintf(fid,'%s \n',msg);
msg = ['  Copenhagen/Bohmian %error= ',num2str(percentError2)];
fprintf(fid,'%s \n',msg);
fclose(fid);



%{
%%                                              replace potential outliers
qUmax = 100*abs(Emax);
qUmin = -qUmax;
Lmin = (q_PE < qUmin );
q_PE(Lmin) = qUmin;
Lmax = (q_PE > qUmax );
q_PE(Lmax) = qUmax;
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
  plot(posX(s,nStart:nFinal),velX(s,nStart:nFinal));
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
%}  
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





%%                                                       list of functions

function rho = getProbDensity1d(qSystem,x,t)
% Probability density at (x,t) for a particle in a 1d-box is calculated
bb = qSystem.bb;
rn = qSystem.rn;
theta = qSystem.theta;                       % initial phase for each mode
omega = qSystem.omega;                         % frequencies for each mode
Cn = rn.*exp( -1i*theta );                % initial expansion coefficients
Qx = qSystem.Qx;
psi = bb*sum( Cn.* exp( -1i*t*omega ).* sin(Qx*x) );
rho = abs(psi)^2;
end

function vx = getVelocity1d(dPhi,qSystem,x,t)
% get vx given x,t and quantum system details for a 1d-box
% ------------------------------------------- unpackage needed information
% bb  = system.bb;     % not needed because |psi|^2 need not be normalized
d     = qSystem.d;                                             % hbar/mass
rn    = qSystem.rn;                  % amplitude of expansion coefficients
theta = qSystem.theta;                       % initial phase for each mode
omega = qSystem.omega;                         % frequencies for each mode
Cn    = rn.*exp( -1i*theta );             % initial expansion coefficients
Qx    = qSystem.Qx;
% --------------------------------------------------------- do calculation
CnExpOmega = Cn.* exp( -1i*(t*omega + dPhi) );
psi           = sum( CnExpOmega.*sin(Qx*x) );                        % *bb
partial_Psi_x = sum( CnExpOmega.*cos(Qx*x).*Qx );                    % *bb
conj_partial_Psi_x = conj( partial_Psi_x );
JBohmx             = d* imag( psi* conj_partial_Psi_x );
rho                = abs(psi)^2;
vx                 = JBohmx/rho;
end

function quantumPE = getQuantumPE1P1D(dPhi,qSystem,x0,t)
% get quantum PE given x0,t and quantum system details for a 1d-box
% ------------------------------------------- unpackage needed information
% bb  = system.bb;     % not needed because |psi|^2 need not be normalized
d     = qSystem.d;                                             % hbar/mass
d2    = qSystem.d2;                                    % (hbar^2)/(2*mass)
rn    = qSystem.rn;                  % amplitude of expansion coefficients
theta = qSystem.theta;                       % initial phase for each mode
omega = qSystem.omega;                         % frequencies for each mode
Cn    = rn.*exp( -1i*theta );             % initial expansion coefficients
Qx    = qSystem.Qx;
% --------------------------------------------------------- do calculation
dx = 0.00001;
dxSQ = dx^2;
CnExpOmega = Cn.* exp( -1i*(t*omega + dPhi) );
xN = x0 - dx;
xP = x0 + dx;
psiP = sum( CnExpOmega.*sin(Qx*xP) );                                % *bb
psi0 = sum( CnExpOmega.*sin(Qx*x0) );                                % *bb
psiN = sum( CnExpOmega.*sin(Qx*xN) );                                % *bb
Rpos = abs(psiP);
Rmid = abs(psi0);
Rneg = abs(psiN);
secondDerivativeR = (Rpos - 2*Rmid + Rneg)/dxSQ;
quantumPE = -d2*secondDerivativeR/Rmid;
%disp( [d2,secondDerivativeR, Rmid]);
%pause
end

function [pX,vX,qU] = propagateMotion1d(dPhi,qSystem,xo,vxo,tArray,h)
% performs RK4 method using adaptive time scale near box walls
% In the case of open system, a quasistatic approximation is used by 
% augmenting a diffusive perturbation to the energy eigenfunction 
% expansion that describes particle motion. 
dtShift = 0.01*h;
h2 = h/2;
h6 = h/6;
ax    = qSystem.ax;
dAbs = 0.01*ax;
dAbs10 = dAbs/10;
kFix = 0;
tEnd = tArray(end);
% ------------------------------------------------------------------------
phiPts = length(dPhi);
nPts = length(tArray);
pX = zeros(1,nPts);
vX = pX;
qU = pX;                                        % quantum potential energy
pX(1) = xo;
vX(1) = vxo;
t0 = tArray(1);
qU(1) = getQuantumPE1P1D(dPhi,qSystem,xo,t0);
% ----------------------------------------- augment diffusive perturbation
openSystem = qSystem.openSystem;
  if( openSystem )
  sigDp = qSystem.sigDp;
  sigDa = 5*qSystem.sigDa;   % 5 is needed due to internal running average
  sigDL = 5*qSystem.sigDL;   % 5 is needed due to internal running average
  fDa = 0;  fDL = 0;         % running average fDa= 0.923*fDa + 0.077*dfDa
    for iii=1:10000
    dfDa = sigDa*randn;
    fDa = 0.923*fDa + 0.077*dfDa; 
    dfDL = sigDL*randn;
    fDL = 0.923*fDL + 0.077*dfDL; 
    end
  rn = qSystem.rn;
  ax = qSystem.ax;
  qModes = qSystem.qModes;
  hbar = 4.14/(2*pi);                               % dimension is (eV*fs)
  dd = qSystem.dd;
  %else
  % dPhi is fixed (nothing to do)
  end
% ------------------------------------- set diffusion rate for phase angle
    for nt=2:nPts
    t0  = tArray(nt-1);
    x0  = pX(nt-1);
    vx0 = vX(nt-1);
    tf  = tArray(nt) - dtShift;
% ------------------------------------------------------- implementing RK4
    t = t0;
      while( t < tf )
      dt = h;                                                % presumption
      k0x = vx0;
      x1 = x0 + k0x*h2;
      k1x = getVelocity1d(dPhi,qSystem,x1,t+h2);
      x2 = x1 + k1x*h2;
      k2x = getVelocity1d(dPhi,qSystem,x2,t+h2);
      x3 = x2 + k2x*h;
      k3x = getVelocity1d(dPhi,qSystem,x3,t+h);
% -------------------------------------------------------- try normal step
      xx = x0 + (k0x + 2*k1x + 2*k2x + k3x)*h6;
      dxAbs = abs(xx - x0);
% ------------------------------------------------------- check boundaries
        while( xx > ax || xx < 0 || dxAbs > dAbs )
        dt = dt/2;
          if( dt < 0.0001 )               % => try to shake off divergence
          dt = 0.0001;
          xx = x0 + dAbs10*randn;            % perform a small random step
            while( xx > ax || xx < 0 || dxAbs > dAbs )
            xx = x0 + dAbs10*randn;
            end
          %disp(['BC: random step made: t= ',num2str(t), ...
          %     '  xx= ',num2str(xx)]);                  % can comment out
          kFix = kFix + 1;     
            if( kFix > 10 )
            error('velocity field is too erratic to integrate');
            end
          else
          kFix = 0;
          hB = dt;
          h2B = dt/2;
          h6B = hB/6;
          x1 = x0 + k0x*h2B;
          k1x = getVelocity1d(dPhi,qSystem,x1,t+h2B);
          x2 = x1 + k1x*h2B;
          k2x = getVelocity1d(dPhi,qSystem,x2,t+h2B);
          x3 = x2 + k2x*hB;
          k3x = getVelocity1d(dPhi,qSystem,x3,t+hB);
% ----------------------------------------------------- try shortened step
          xx = x0 + (k0x + 2*k1x + 2*k2x + k3x)*h6B;
          end
        %disp(['BC: ',num2str([t/tEnd,dt,xx/ax])]);      % can comment out
        dxAbs = abs(xx - x0);
        end
      t = t + dt;
      x0 = xx;
        if( openSystem )
        dPhi = dPhi + sigDp*randn(1,phiPts);
        dfDa = sigDa*randn;
        fDa = 0.923*fDa + 0.077*dfDa; 
        rnNew = rn.*(1 + fDa);
        probArray = rnNew.^2;
        totalProb = sum(probArray);
        rnNew = sqrt(probArray/totalProb);
        dfDL = sigDL*randn;
        fDL = 0.923*fDL + 0.077*dfDL;
        axNew = ax*(1 + fDL);
        bbNew = sqrt(2/axNew);
        QxNew = (pi*qModes)/axNew;
        Estates = dd*(qModes/axNew).^2;
        omegaNew = Estates/hbar; 
        qSystem.ax = axNew;
        qSystem.bb = bbNew;
        qSystem.rn = rnNew;
        qSystem.Qx = QxNew;
        qSystem.omega = omegaNew;
        end
      vx0 = getVelocity1d(dPhi,qSystem,xx,t);
      end
% ------------------------------------------------------------ record @ nt
% Remark: last point at tArray(nt) is the initial point at tArray(nt-1)
    pX(nt) = x0;
    vX(nt) = vx0;
    tArray(nt) = t;
    qU(nt) = getQuantumPE1P1D(dPhi,qSystem,x0,t);
    end
end







