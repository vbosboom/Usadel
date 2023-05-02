%General execution script for solving the Usadel equations in the
%Phi-parametrization on a 1D SN-N-NS bridge
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
n=500; %number of nodes in the S and L layers
n_mats = 20; %number of Matsubara frequencies used in self-consistency calculation
maxit = 400;    %maximum number of iterations for self-consistency calculation
tol = 10^(-5); %tolerance for convergence

%physical parameters
T = 0.5; %temperature; normalized by T_c
phase = 0.3; %phase difference
gamma = 0.4; %Interface parameter

%bridge dimensions
S = 1; %center length;
L = 5;  %total strip length
nS = n; %number of vertices in center
nL = n; %number of vertices outside of center

%Guassian quadrature points and weights
weights = [128/225;(322+13*sqrt(70))/900;(322+13*sqrt(70))/900;(322-13*sqrt(70))/900;(322-13*sqrt(70))/900];
points = [0;1/3*sqrt(5-2*sqrt(10/7));-1/3*sqrt(5-2*sqrt(10/7));1/3*sqrt(5+2*sqrt(10/7));-1/3*sqrt(5+2*sqrt(10/7))];

%Generate Geometry information
GI = GenerateGeometry1D(S,L,nS,nL,points);

%solve the Usadel equations
phin = SolveUsadel(GI,weights,T,gamma,phase,n_mats,maxit,tol);

%calculate the supercurrent trough the metal layer
Curr = CalculateCurrent(GI,n_mats,phin,T);

%plot the current density in the metal layer
plot(GI.x,Curr)
xlabel('x/\xi_N')
ylabel('eI_cR_N/2\piT_c')
title('Current density')