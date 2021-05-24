%General execution script for solving the Usadel equations in the
%Theta,Chi-parametrization on a 1D SN-N-NS bridge
clear all

%add general functions folder to path
mydir = pwd;
idcs = strfind(mydir,'\');
GenPath = mydir(1:idcs(end-2)-1);
addpath(strcat(GenPath,'\General'));

%computational parameters
n=500; %number of nodes in the S and L layers
maxit = 400;    %maximum number of iterations for self-consistency calculation
tol = 10^(-5); %tolerance for convergence

%physical parameters
T = 0.5; %temperature; normalized by T_c
phase = 0.3; %phase difference
gamma = 0.4; %Interface parameter
E = -0.6; %Energy, negative value easier to implement theta_S BC

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

%Solve the Usadel equations for energy E
[chi,theta,success] = RealEnergy(GI,points,weights,T,phase,E,gamma,maxit,tol,false);