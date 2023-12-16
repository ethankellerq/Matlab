  % Numerical Verific ation of Schrodinger Eq in 1D
clear all
 
tic
% Set up units
% Basic units are femtosecond, angstrom, and electron-volt
%%                                               set fundamental constants
%                                mass of proton = 1836.152*(electron mass)
hbar = 4.14/(2*pi);                                 % dimension is (eV*fs)
c = 3000;                                  % speed of light in Angstrom/fs
eEnergy = 510990.6;                        % rest energy of electron in eV
eMass = eEnergy/c^2;                         % in units of (fs^2)*eV/Ang^2

%%                                  Set ranges for parameters and variables
%   
% sigma=10*rand(1);                       % sigma ranges from 0-10 Angstrom
% phi=2*pi*rand(1);     
% x_0=10*rand(1);
% v_0=10*rand(1);
% a=1 *rand(1);
%x=1*rand(1);
%t=1*rand(1);

sigma=10;
phi=0/sigma^2;
%phi=0;
x_0=0;
v_0=10;
a=0;
m=1*eMass;

numPulses=2;
numParticles=2;
PAR0=zeros(numPulses,6,numParticles);

% Creating parameter 3D array where each page represents a new particle

for r=1:numParticles
    
for w=1:numPulses
 
 PAR0(w,1,r)=sigma; %+10*rand(1);
 PAR0(w,2,r)=phi; %1*rand(1);
 PAR0(w,3,r)=x_0; %+10*rand(1);
 PAR0(w,4,r)=v_0; %+10*rand(1);
 PAR0(w,5,r)=a; %+1*rand(1);
 PAR0(w,6,r)=m;    

end
end

PAR0(1,4,2)=11;
PAR0(2,4,2)=9;
%PAR0(1,4,1)=10;
%PAR0(2,4,1)=-10;

% PAR0(1,1,2)=1;
% PAR0(2,1,2)=1;


%error("stop")
% PAR0(1,1)=sigma;
% PAR0(1,2)=phi;
% PAR0(1,3)=x_0;
% PAR0(1,4)=v_0;
% PAR0(1,5)=a;
% PAR0(1,6)=m;
% 
% PAR0(2,1)=sigma;
% PAR0(2,2)=phi;
% PAR0(2,3)=x_0;
% PAR0(2,4)=v_0;
% PAR0(2,5)=a;
% PAR0(2,6)=m;
% 
% PAR0(3,1)=sigma;
% PAR0(3,2)=phi;
% PAR0(3,3)=x_0;
% PAR0(3,4)=v_0;
% PAR0(3,5)=a;
% PAR0(3,6)=m;

%% Setting up Schrodinger's Eq

%{
k=21;

MLA_f=zeros(1,k);

for j=1:k
t=j-1;

n=100001;

lambda=sqrt((sigma^2)*(1-2*hbar.*t*phi/m).^2+(hbar.*t/(2*m*sigma)).^2);

x_i=-5*lambda + (x_0+v_0*t+0.5*a*t^2);
x_f=5*lambda + (x_0+v_0*t+0.5*a*t^2);
x=linspace(x_i,x_f,n);
dx=x(2)-x(1);

% t_i=0;
% t_f=100;
% t=linspace(t_i,t_f,n);
% dt=t(2)-t(1);

const1=1/(4*sigma^2) + 1i*phi;

const2=(2*pi*sigma^2)^(-1/4);

gamma=sqrt(1 + 2*1i*const1.*t*hbar/m);
%D1_gamma_t=1i*hbar*const1./(m.*gamma);

alpha=(x-x_0-v_0.*t-(1/2)*a.*t.^2);
%D1_alpha_t=-1*(v_0+a.*t);

beta= x-x_0-v_0.*t/2;
delta= (x-x_0).*t-0.5*v_0.*t.^2-(1/6)*a.*t.^3;

psi=const2*(1./gamma).*exp( (-const1*alpha.^2)/gamma.^2+ ...
                1i*m*v_0.*beta/hbar + ...
                1i*m*x_0*a.*t/hbar + 1i*m*a.*delta/hbar ); 

end            
% D1_h_x=-2*const1.*alpha./gamma.^2 + 1i*m*v_0/hbar + 1i*m*a.*t/hbar;
% 
% D1_h_t=(2*const1*alpha.*(D1_gamma_t.*alpha - gamma.*D1_alpha_t))/gamma.^3 +...
%     (-1i*m*v_0^2)/(2*hbar) + 1i*m*a*x_0/hbar + 1i*m*a.*alpha/hbar;
% 
% D1_psi_x=psi.*D1_h_x;
% D2_psi_x=psi.*((-2*const1./gamma.^2)+(D1_h_x).^2);
% 
% D1_psi_t=psi.*((-1*D1_gamma_t./gamma)+D1_h_t);
% 
% wave_eq_lhs=1i*hbar.*D1_psi_t;
% wave_eq_rhs=-1*((D2_psi_x.*hbar^2)/(2*m))-(m*a*x.*psi);
% % split this up into KE ^ and               PE^
% KE=(D2_psi_x.*hbar^2)/(2*m);
% PE=(m*a*x.*psi);
% con1= abs(wave_eq_lhs) - abs(wave_eq_rhs);

%}

%% Numerical integration
% Trapezoidal Integration Technique 

% psi_star= conj(psi);
% iG_prob= psi_star .* psi;
% iG_prob = abs(psi).^2;
% iG_energy_1= psi_star .* wave_eq_lhs;
% iG_energy_2= psi_star .* wave_eq_rhs;
% iG_momentum= psi_star .* -1i*hbar.*D1_psi_x;

% Integration by Trapeziodal Approximation
%{
CA=zeros(1,n);                          % Preallocation
%CB=zeros(1,n);
MLA=zeros(1,n);
%MLB=zeros(1,n);


for i = 1:n-1                           % CA and CB are directly from the                            
    CA(i)=(iG_momentum(i+1)+iG_momentum(i))*dx/2;         % the equation for a trapezoid.
    MLA(i+1)=MLA(i)+CA(i);              % MLA and MLB allow me to sum 
    %CB(i)=(NB(i+1)+NB(i))*dt/2;         % CA and CB over i which produces
    %MLB(i+1)=MLB(i)+CB(i);              % the approximation of the integral
end

MLA_f(j)=MLA(end);

psi_real=real(psi);

figure(1)
plot(x,psi_real);
xlabel('Position (Angstrom)');
ylabel('Psi');
title('Real component of Psi as func of x');

end

% Graphing only every 100 points, Changing DDt effects smoothness of
% graph but does not effect accuracy which is determined by n. 
%DDx=1000;
%MLAs=MLA(1:DDx:end);
%%MLBs=MLB(1:DDt:end);
%m=length(MLAs);
%T=linspace(x_i,x_f,m);


% figure(1)
% plot(T,MLAs);
% xlabel('Position (Angstrom)');
% ylabel('Total area under the curve');
% title('Trapezoidal Integration');

%}


%% Plotting

%psi_real=real(psi);
%psi_imag=imag(psi);

% Plots for position as indep
%{
figure(2)
plot(x,psi_real);
xlabel('Position (Angstrom)');
ylabel('Psi');
title('Real component of Psi as func of x');

figure(3)
plot(x,psi_imag);
xlabel('Position (Angstrom)');
ylabel('Psi');
title('Imaginary component of Psi as func of x');

figure(4)
plot(x,iG_prob);
xlabel('Position (Angstrom)');
ylabel('Psi Modulus Squared');
title('Psi Modulus Squared as func of x');
%}

% Plots for time as indep
%{
figure(2)
plot(t,psi_real);
xlabel('Time (femtosecond)');
ylabel('Psi');
title('Real component of Psi as func of t');

figure(3)
plot(t,psi_imag);
xlabel('Time (femtosecond)');
ylabel('Psi');
title('Imaginary component of Psi as func of t');

figure(4)
plot(t,iG_prob);
xlabel('Time (femtosecond)');
ylabel('Psi Modulus Squared');
title('Psi Modulus Squared as func of t');
%}

%% Bohmian Trajectories

%{

q_p=100000;
q=zeros(1,q_p);
t=zeros(1,q_p);
accel=zeros(1,q_p);
velocity=zeros(1,(q_p+1));
dt=.001;
q(1)=0*sigma;
%RK4 Algorithm
for z=1:q_p
    
temp_x_1=q(z);
temp_t_1=t(z);
K1=Bohm_Velocity_Func(temp_x_1,temp_t_1,PAR0);

temp_x_2=temp_x_1+K1*(dt/2);
temp_t_2=t(z)+dt/2;
K2=Bohm_Velocity_Func(temp_x_2,temp_t_2,PAR0);

temp_x_3=temp_x_1+K2*(dt/2);
temp_t_3=t(z)+dt/2;
K3=Bohm_Velocity_Func(temp_x_3,temp_t_3,PAR0);

temp_x_4=temp_x_1+K3*dt;
temp_t_4=t(z)+dt;
K4=Bohm_Velocity_Func(temp_x_4,temp_t_4,PAR0);

avg_Slope=(K1+2*K2+2*K3+K4)*(1/6);
 
velocity(z)=avg_Slope;

if velocity(z) > 300
error('Error: Relativistic effects');
end

% if z>2 && z<(q_p-1)
% accel(z)=(velocity(z)-velocity(z-1))/(dt);
% end

q(z+1)=q(z)+avg_Slope*dt;
t(z+1)=t(z)+dt;
end

hold on
figure(5)
plot(t,q);
xlabel('Time (femtosecond)');
ylabel('Position (Angstrom)');
title('Position as func of time');

figure(6)
plot(t,velocity);
xlabel('Time (femtosecond)');
ylabel('Velocity (Angstrom/fs)');
title('Velocity as func of time');
% 
% figure(7)
% plot(t,accel);
% xlabel('Time (femtosecond)');
% ylabel(' (Angstrom)');
% title('Position as func of time');

%}


%% Superposition of gaussians 

%Q=q(1:1000:end);

%
k=2;
linP=1000;

for j=1:k

t1=(j-1)*100;
psi_super_array=zeros(numParticles,linP);
prod_state=ones(1,linP);
x_array=zeros(numParticles,linP);

for u=1:numParticles

packet_pos=(PAR0(:,3,u)+PAR0(:,4,u)*t1+0.5*PAR0(:,5,u)*t1^2);
    
lambda=sqrt((PAR0(:,1,u).^2).*(1-2*hbar*t1.*PAR0(:,2,u)./PAR0(:,6,u)).^2 + ...
    (hbar*t1./(2*PAR0(:,6,u).*PAR0(:,1,u))).^2);

x_i=min(packet_pos)-5*max(lambda);
x_f=max(packet_pos)+5*max(lambda);
x_array(u,:)=linspace(x_i,x_f,linP);

for i=1:numPulses

const1=1/(4*PAR0(i,1,u)^2) + 1i*PAR0(i,2,u);
const2=(2*pi*PAR0(i,1,u)^2)^(-1/4);

gamma=sqrt(1 + 2*1i*const1.*t1*hbar/PAR0(i,6,u));

alpha=(x_array(u,:)-PAR0(i,3,u)-PAR0(i,4,u).*t1-(1/2)*PAR0(i,5,u).*t1.^2);

beta= x_array(u,:)-PAR0(i,3,u)-PAR0(i,4,u).*t1/2;

delta= (x_array(u,:)-PAR0(i,3,u)).*t1-0.5*PAR0(i,4,u).*t1.^2-(1/6)*PAR0(i,5,u).*t1.^3;

psi=const2*(1./gamma).*exp( (-const1*alpha.^2)/gamma.^2+ ...
                1i*PAR0(i,6,u)*PAR0(i, 4,u).*beta/hbar + ...
                1i*PAR0(i,6,u)*PAR0(i,3,u)*PAR0(i,5,u).*t1/hbar + ...
                1i*PAR0(i,6,u)*PAR0(i,5,u).*delta/hbar ); 

psi_super_array(u,:)=psi+psi_super_array(u,:);

end            

%prod_state=prod_state.*psi_super_array(l,:);
prod_state=zeros(1000);

for n=1:linP
for m=1:linP
prod_state(n,m)=psi_super_array(1,n) * psi_super_array(2,m); 
end
end

end

psi_product_prob=abs(prod_state).^2;

figure(10)
clf
x_max=max(max(x_array,[],2));
x_min=min(min(x_array,[],2));
x_plot1=linspace(x_min,x_max,linP);

%h=surf(x_plot,x_plot,psi_product_prob);
% h=heatmap(x_plot1,x_plot1,psi_product_prob);
% h.XDisplayData = [];
% h.YDisplayData = [];
%h=imagesc(x_plot1,x_plot1,psi_product_prob);
h=imagesc(x_array(1,:),x_array(2,:),psi_product_prob);

set(gca,'YDir','normal') 

%set(h,'LineStyle','none')
colormap('turbo');

pause(1.5);
%plot(x,psi_product_prob);
%xline(Q(j))
%xlabel('Position (Angstrom)');
%ylabel('Psi Modulus Squared');
title('Psi Modulus Squared as func of x1 and x2');

end

%}
toc
%%

function [D1_q_t] = Bohm_Velocity_Func(x,t,PAR0)
% Function to Calculate Psi, Psi*, and D1_psi_x to give 
% initialize variables
hbar = 4.14/(2*pi);
[numD,~]=size(PAR0);
psi_total=0;
D1_psi_x_total=0;

for i=1:numD

const1=1/(4*PAR0(i,1)^2) + 1i*PAR0(i,2);
const2=(2*pi*PAR0(i,1)^2)^(-1/4);

gamma=sqrt(1 + 2*1i*const1.*t*hbar/PAR0(i,6));

alpha=(x-PAR0(i,3)-PAR0(i,4).*t-(1/2)*PAR0(i,5).*t.^2);

beta= x-PAR0(i,3)-PAR0(i,4).*t/2;

delta= (x-PAR0(i,3)).*t-0.5*PAR0(i,4).*t.^2-(1/6)*PAR0(i,5).*t.^3;

psi=const2*(1./gamma).*exp( (-const1*alpha.^2)/gamma.^2+ ...
                1i*PAR0(i,6)*PAR0(i, 4).*beta/hbar + ...
                1i*PAR0(i,6)*PAR0(i,3)*PAR0(i,5).*t/hbar + ...
                1i*PAR0(i,6)*PAR0(i,5).*delta/hbar ); 

D1_h_x= -2*const1.*alpha./gamma.^2 + 1i*PAR0(i,6)*PAR0(i,4)/hbar + ...
    1i*PAR0(i,6)*PAR0(i,5).*t/hbar;

D1_psi_x=psi.*D1_h_x;

psi_total=psi_total+psi;

D1_psi_x_total=D1_psi_x_total+D1_psi_x;

end

psi_total_star=conj(psi_total);
prob_Total=abs(psi_total)^2;
% Getting velocity thru prob density

j=(hbar/PAR0(1,6))*imag(psi_total_star*D1_psi_x_total);
D1_q_t=(j/prob_Total);

end
