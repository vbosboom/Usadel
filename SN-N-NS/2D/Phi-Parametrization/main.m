%General execution script for solving the Usadel equations in the
%Phi-parametrization on a 2D SN-N-NS bridge
clear all

%add general functions folder to path
mydir = pwd;
sep = ['\','/'];
sepindx = contains(mydir,sep(2));
idcs = strfind(mydir,sep(sepindx+1));
GenPath = mydir(1:idcs(end-2)-1);
GenEnd = strcat(sep(sepindx+1),'General');
addpath(strcat(GenPath,GenEnd));

%computational parameters
minvert = 400; %minimum number of nodes in computational domain
n_mats = 20; %number of Matsubara frequencies used in self-consistency calculation
maxit = 400;    %maximum number of iterations for self-consistency calculation
tol = 10^(-5); %tolerance for convergence

%declare the physical parameters during our simulation
gamma_B = 1.5; %interface parameter for both boundaries  gamma_{B}
gamma = 1; %proximit effect parameter
ksi = 1; %coherence length inside superconductors relative to normal metal
T=0.5; %temperature relative to T_c
phase = pi/2; %phase difference over the junction

%spatial dimensions of the junction
s=1;
l = 5;P
dn = 0.2;
ds = 5;

%2D Gaussian quadrature points and weights in barycentric coordinates
weights = [0.109951743655322,0.109951743655322,0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011];
points = [0.816847572980459,0.091576213509771,0.091576213509771;0.091576213509771,0.816847572980459,0.091576213509771;0.091576213509771,0.091576213509771,0.816847572980459...
    ;0.108103018168070,0.445948490915965,0.445948490915965;0.445948490915965,0.108103018168070,0.445948490915965;0.44594849091596,0.44594849091596,0.108103018168070];

%1D Gaussian quadrature points and weights
weights1D = [128/225,(322+13*sqrt(70))/900,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900,(322-13*sqrt(70))/900];
points1D = [1/2,1/2;1-1/6*sqrt(5-2*sqrt(10/7))-1/2,1/6*sqrt(5-2*sqrt(10/7))+1/2;1+1/6*sqrt(5-2*sqrt(10/7))-1/2,-1/6*sqrt(5-2*sqrt(10/7))+1/2;1-1/6*sqrt(5+2*sqrt(10/7))-1/2,1/6*sqrt(5+2*sqrt(10/7))+1/2;1+1/6*sqrt(5+2*sqrt(10/7))-1/2,-1/6*sqrt(5+2*sqrt(10/7))+1/2];

%Generate geometry information
GI = GenerateGeometry(s,l,dn,ds,minvert,false);

%Solve the Usadel equations self-consistently
[phis1,phis2,phin,deltas1,deltas2,Difference] = SolveUsadel(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,phase,ksi,n_mats,maxit,tol,false);

%Calculate the current density inside the junction
[Currxs1,Currxs2,Currxn,Currys1,Currys2,Curryn,indx] = CalculateCurrent(GI,n_mats,phis1,phis2,phin,T,gamma_B,gamma,ksi);

%Plot the current density throughout the junction
figure
hold on
fill3(GI.x(GI.elmats1)',GI.y(GI.elmats1)',kron(ones(1,size(GI.elmats1,2)),Currxs1)',kron(ones(1,size(GI.elmats1,2)),Currxs1)');
fill3(GI.x(GI.elmats2)',GI.y(GI.elmats2)',kron(ones(1,size(GI.elmats2,2)),Currxs2)',kron(ones(1,size(GI.elmats2,2)),Currxs2)');
fill3(GI.x(GI.elmatn)',GI.y(GI.elmatn)',kron(ones(1,size(GI.elmatn,2)),Currxn)',kron(ones(1,size(GI.elmatn,2)),Currxn)');
xlabel('x/\xi_N')
ylabel('y/\xi_N')
title('Current density (eI_cR_N/2\piT_c)')
colorbar