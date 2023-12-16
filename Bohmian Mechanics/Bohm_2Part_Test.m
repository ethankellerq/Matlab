%% Testing 2d entanglement
clear All
%%                                               set fundamental constants
%                                mass of proton = 1836.152*(electron mass)
hbar = 4.14/(2*pi);                                 % dimension is (eV*fs)
c = 3000;                                  % speed of light in Angstrom/fs
eEnergy = 510990.6;                        % rest energy of electron in eV
eMass = eEnergy/c^2;                         % in units of (fs^2)*eV/Ang^2

%%                                  Set ranges for parameters and variables

sigma=10;
phi=0/sigma^2;
x_0=0;
v_0=0;
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

PAR0(1,4,2)=-5;
PAR0(2,4,2)=5;
PAR0(1,4,1)=0;
PAR0(2,4,1)=5;

%

q_p=100000;
q1=zeros(1,q_p);
q2=zeros(1,q_p);
t=zeros(1,q_p); 
accel=zeros(1,q_p);
velocity1=zeros(1,(q_p));
velocity2=zeros(1,(q_p));
dt_0=.0001;
dt=dt_0;
q1(1)=1*sigma;
q2(1)=2*sigma;

%% RK4 Algorithm

for z=1:q_p

temp_x1_1=q1(z);
temp_x2_1=q2(z);
temp_t_1=t(z);
[K1_1,K1_2]=Bohm_Velocity_Func(temp_x1_1,temp_x2_1,temp_t_1,PAR0);

temp_x1_2=temp_x1_1+K1_1*(dt/2);
temp_x2_2=temp_x2_1+K1_2*(dt/2);
temp_t_2=t(z)+dt/2;
[K2_1,K2_2]=Bohm_Velocity_Func(temp_x1_2,temp_x2_2,temp_t_2,PAR0);

temp_x1_3=temp_x1_1+K2_1*(dt/2);
temp_x2_3=temp_x2_1+K2_2*(dt/2);
temp_t_3=t(z)+dt/2;
[K3_1,K3_2]=Bohm_Velocity_Func(temp_x1_3,temp_x2_3,temp_t_3,PAR0);

temp_x1_4=temp_x1_1+K3_1*dt;
temp_x2_4=temp_x2_1+K3_2*dt;
temp_t_4=t(z)+dt;
[K4_1,K4_2]=Bohm_Velocity_Func(temp_x1_4,temp_x2_4,temp_t_4,PAR0);

avg_Slope1=(K1_1+2*K2_1+2*K3_1+K4_1)*(1/6);
avg_Slope2=(K1_2+2*K2_2+2*K3_2+K4_2)*(1/6);

velocity1(z)=avg_Slope1;
velocity2(z)=avg_Slope2;
 
if (velocity1(z)||velocity2(z)) > 3000
error('Error: Cannot be faster than speed of light');
end 

% Adaptive time
%{
if  z>1 && abs(velocity1(z-1)-velocity1(z)) > 0.001
dt=dt_0;
elseif z>1 && abs(velocity2(z-1)-velocity2(z)) > 0.001
dt=dt_0;
elseif z>1 && abs(velocity1(z-1)-velocity1(z)) < .000001 ...
        && abs(velocity2(z-1)-velocity2(z)) < .000001
dt=2*dt;
end

if dt>1
dt=1;
end

%}  

q1(z+1)=q1(z)+avg_Slope1*dt;
q2(z+1)=q2(z)+avg_Slope2*dt;
t(z+1)=t(z)+dt;
end

%}

figure(1)
plot(t,q1);
xlabel('Time (femtosecond)');
%xlim([0 10]);
ylabel('Position of Particle1 (Angstrom)');
%ylim([0 100]);
title('Position as func of time');

figure(2)
plot(t,q2);
xlabel('Time (femtosecond)');
%xlim([0 10]);
ylabel('Position of Particle2 (Angstrom)');
%ylim([0 100]);
title('Position as func of time');

figure(3)
tv=t(1:q_p);
plot(tv,velocity1)
hold on 
plot(tv,velocity2)
xlabel('Time (femtosecond)');
ylabel('Velocity of Particles (Angstrom/fs)');

%% Creating 3D PDF 
% {
Q1=q1(1:10000:end);
Q2=q2(1:10000:end);

k=11;
linP=1000;

for j=1:k

t1=(j-1)*1;
psi_super_array1=zeros(numParticles,linP);
psi_super_array2=zeros(numParticles,linP);
D1_psi_x_array=zeros(numParticles,linP);
prod_state=ones(1,linP);
x_array=zeros(numParticles,linP);

packet_pos=(PAR0(:,3,:)+PAR0(:,4,:)*t1+0.5*PAR0(:,5,:)*t1^2);

lambda=sqrt((PAR0(:,1,:).^2).*(1-2*hbar*t1.*PAR0(:,2,:)./PAR0(:,6,:)).^2 + ...
    (hbar*t1./(2*PAR0(: ,6,:).*PAR0(:,1,:))).^2);

x_f1=max(packet_pos(:,:,1))+5*max(lambda(:,:,1));
x_i1=-1*x_f1;

x_f2=max(packet_pos(:,:,2))+5*max(lambda(:,:,2));
x_i2=-1*x_f2;

x_array(1,:)=linspace(x_i1,x_f1,linP);
x_array(2,:)=linspace(x_i2,x_f2,linP);

for u=1:numParticles

for i=1:numPulses

const1=1/(4*PAR0(i,1,u)^2) + 1i*PAR0(i,2,u);
const2=(2*pi*PAR0(i,1,u)^2)^(-1/4);

gamma=sqrt(1 + 2*1i*const1.*t1*hbar/PAR0(i,6,u));

alpha=(x_array-PAR0(i,3,u)-PAR0(i,4,u).*t1-(1/2)*PAR0(i,5,u).*t1.^2);

beta= x_array-PAR0(i,3,u)-PAR0(i,4,u).*t1/2;

delta= (x_array-PAR0(i,3,u)).*t1-0.5*PAR0(i,4,u).*t1.^2-(1/6)*PAR0(i,5,u).*t1.^3;

psi=const2*(1./gamma).*exp( (-const1*alpha.^2)/gamma.^2+ ...
                1i*PAR0(i,6,u)*PAR0(i,4,u).*beta/hbar + ...
                1i*PAR0(i,6,u)*PAR0(i,3,u)*PAR0(i,5,u).*t1/hbar + ...
                1i*PAR0(i,6,u)*PAR0(i,5,u).*delta/hbar ); 
if u==1
psi_super_array1(1,:)=psi(1,:)+psi_super_array1(1,:);
psi_super_array1(2,:)=psi(2,:)+psi_super_array1(2,:);
elseif u==2
psi_super_array2(1,:)=psi(1,:)+psi_super_array2(1,:);
psi_super_array2(2,:)=psi(2,:)+psi_super_array2(2,:);
end
        
end            

prod_state=zeros(1000);

for n=1:linP
for m=1:linP
prod_state(n,m)=psi_super_array1(1,m).*psi_super_array2(2,n)-...
   psi_super_array1(2,n).*psi_super_array2(1,m);
end
end

end

psi_product_prob=abs(prod_state).^2;

figure(10)
clf
h=imagesc(x_array(1,:),x_array(2,:),psi_product_prob);
set(gca,'YDir','normal') 
xlabel('x1')
ylabel('x2')
colormap('turbo');
xline(packet_pos(1,1,1),'w','LineWidth',2.5);
xline(packet_pos(2,1,1),'w','LineWidth',2.5);
yline(packet_pos(1,1,2),'w','LineWidth',2.5);
yline(packet_pos(2,1,2),'w','LineWidth',2.5);

xline(Q1(j),'g','LineWidth',2.5);
yline(Q2(j),'g','LineWidth',2.5);


% pause(1.5);
% plot(x,psi_product_prob);
% xline(Q(j))
% xlabel('Position (Angstrom)');
% ylabel('Psi Modulus Squared');
% title('Psi Modulus Squared as func of x');

end

%}

%% Function
function [D1_q1_t , D1_q2_t] = Bohm_Velocity_Func(x1,x2,t,PAR0)
% Function to Calculate Psi, Psi*, and D1_psi_x to give 
% initialize variables
hbar = 4.14/(2*pi);
[numPulses,~,numParticles]=size(PAR0);
psi_super_array=zeros(numParticles,2);
D1_psi_x_array=zeros(numParticles,2);
x_array=[x1 x2];

for u=1:numParticles
   
for i=1:numPulses
    
const1=1/(4*PAR0(i,1,u)^2) + 1i*PAR0(i,2,u);
const2=(2*pi*PAR0(i,1,u)^2)^(-1/4);

gamma=sqrt(1 + 2*1i*const1.*t*hbar/PAR0(i,6,u));

alpha=(x_array-PAR0(i,3,u)-PAR0(i,4,u).*t-(1/2)*PAR0(i,5,u).*t.^2);

beta= x_array-PAR0(i,3,u)-PAR0(i,4,u).*t/2;

delta= (x_array-PAR0(i,3,u)).*t-0.5*PAR0(i,4,u).*t.^2-(1/6)*PAR0(i,5,u).*t.^3;

psi=const2*(1./gamma).*exp( (-const1*alpha.^2)/gamma.^2+ ...
                1i*PAR0(i,6,u)*PAR0(i,4,u).*beta/hbar + ...
                1i*PAR0(i,6,u)*PAR0(i,3,u)*PAR0(i,5,u).*t/hbar + ...
                1i*PAR0(i,6,u)*PAR0(i,5,u).*delta/hbar ); 

D1_h_x= -2*const1.*alpha./gamma.^2 + 1i*PAR0(i,6,u)*PAR0(i,4,u)/hbar + ...
    1i*PAR0(i,6,u)*PAR0(i,5,u).*t/hbar;

D1_psi_x=psi.*D1_h_x;

psi_super_array(u,:)=psi+psi_super_array(u,:);

D1_psi_x_array(u,:)=D1_psi_x+D1_psi_x_array(u,:);
            
end

end
% Getting velocity thru prob density
%disp(psi_super_array)
psi_total_ent=psi_super_array(1,1)*psi_super_array(2,2)-...
    psi_super_array(1,2)*psi_super_array(2,1);
%disp(psi_total_ent);

D1_psi_total_x1=D1_psi_x_array(1,1)*psi_super_array(2,2)-...
    D1_psi_x_array(2,1)*psi_super_array(1,2);

D1_psi_total_x2=D1_psi_x_array(2,2)*psi_super_array(1,1)-...
    D1_psi_x_array(1,2)*psi_super_array(2,1);

psi_total_conj=conj(psi_total_ent);

%prob_total=abs(psi_total_ent)^2;
epsilon=abs(psi_total_ent);

j1=(hbar/PAR0(1,6,1))*imag(psi_total_conj*D1_psi_total_x1);
j2=(hbar/PAR0(1,6,2))*imag(psi_total_conj*D1_psi_total_x2);
D1_q1_t=(1/epsilon)*(j1/epsilon);
D1_q2_t=(1/epsilon)*(j2/epsilon);
% disp(D1_q1_t)
% disp(D1_q2_t)

end
