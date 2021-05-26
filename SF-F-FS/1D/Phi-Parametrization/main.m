%General execution script for solving the Usadel equations in the
%Phi-parametrization on a 1D SF-F-FS bridge
clear all

%add general functions folder to path
mydir = pwd;
idcs = strfind(mydir,'\');
GenPath = mydir(1:idcs(end-2)-1);
addpath(strcat(GenPath,'\General'));

%computational parameters
n=500; %number of nodes in the S and L layers
n_mats = 20; %number of Matsubara frequencies used in self-consistency calculation
maxit = 400;    %maximum number of iterations for self-consistency calculation
tol = 10^(-5); %tolerance for convergence

%physical parameters
T = 0.5; %temperature; normalized by T_c
phase = 0.3; %phase difference
gamma = 0.4; %Interface parameter
H = 0.1; %Exchange field

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
[phiP,phiCM] = SolveUsadel(GI,weights,T,gamma,H,phase,n_mats,maxit,tol,false);

%calculate the supercurrent trough the junction
Curr = CalculateCurrent(GI,n_mats,phiP,phiCM,T,H);