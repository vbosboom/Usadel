%General execution script for solving the Usadel equations in the
%Phi-parametrization on a 2D SF-F-FS bridge
clear all

%add general functions folder to path
mydir = pwd;
idcs = strfind(mydir,'\');
GenPath = mydir(1:idcs(end-2)-1);
addpath(strcat(GenPath,'\General'));

%computational parameters
minvert=400; %minimum number of nodes in computational domain
n_mats = 20; %number of Matsubara frequencies used in self-consistency calculation
maxit = 400;    %maximum number of iterations for self-consistency calculation
tol = 10^(-5); %tolerance for convergence

%declare the physical parameters during our simulation
gamma_B = 1.5; %interface parameter for both boundaries  gamma_{B}
gamma = 1; %proximit effect parameter
ksi = 1; %coherence length inside superconductors relative to normal metal
T=0.5; %temperature relative to T_c
phase = pi/2; %phase difference over the junction
H = 0.3; %Exchange field

%spatial dimensions of the junction
s=1;
l = 5;
dn = 0.2;
ds = 5;

%parameters for the Gaussian quadrature calculation
weights = [0.109951743655322,0.109951743655322,0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011];
points = [0.816847572980459,0.091576213509771,0.091576213509771;0.091576213509771,0.816847572980459,0.091576213509771;0.091576213509771,0.091576213509771,0.816847572980459...
    ;0.108103018168070,0.445948490915965,0.445948490915965;0.445948490915965,0.108103018168070,0.445948490915965;0.44594849091596,0.44594849091596,0.108103018168070];

%1D parameters for the Gaussian quadrature
weights1D = [128/225,(322+13*sqrt(70))/900,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900,(322-13*sqrt(70))/900];
points1D = [1/2,1/2;1-1/6*sqrt(5-2*sqrt(10/7))-1/2,1/6*sqrt(5-2*sqrt(10/7))+1/2;1+1/6*sqrt(5-2*sqrt(10/7))-1/2,-1/6*sqrt(5-2*sqrt(10/7))+1/2;1-1/6*sqrt(5+2*sqrt(10/7))-1/2,1/6*sqrt(5+2*sqrt(10/7))+1/2;1+1/6*sqrt(5+2*sqrt(10/7))-1/2,-1/6*sqrt(5+2*sqrt(10/7))+1/2];

%create the geometry info structure
GI = GenerateGeometry(s,l,dn,ds,minvert,false);

%Solve the Usadel equations self-consistently
[phis1P,phis1CM,phis2P,phis2CM,phifP,phifCM,deltas1,deltas2,Difference] = SolveUsadel(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,H,phase,ksi,n_mats,maxit,tol,false);

%Calculate the current inside the junction
[Currxs1,Currxs2,Currxf,Currys1,Currys2,Curryf,indx] = CalculateCurrent(GI,n_mats,phis1P,phis1CM,phis2P,phis2CM,phifP,phifCM,T,H,gamma,ksi);