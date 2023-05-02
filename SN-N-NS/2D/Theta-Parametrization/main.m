%General execution script for solving the Usadel equations in the
%Theta-parametrization on a 2D SN-N-NS bridge
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
minvert=400; %minimum number of nodes in computational domain
maxit = 400;    %maximum number of iterations for self-consistency calculation
tol = 10^(-5); %tolerance for convergence

%declare the physical parameters during our simulation
gamma_B = 1.5; %interface parameter for both boundaries  gamma_{B}
gamma = 1; %proximit effect parameter
ksi = 1; %coherence length inside superconductors relative to normal metal
T=0.5; %temperature relative to T_c
phase = 0.3; %phase difference over the junction
E = -0.6; %Energy, negative is easier for bulk BC on thetaS

%spatial dimensions of the junction
s=1;
l = 5;
dn = 0.2;
ds = 5;

%2D Gaussian quadrature points and weights
weights = [0.109951743655322,0.109951743655322,0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011];
points = [0.816847572980459,0.091576213509771,0.091576213509771;0.091576213509771,0.816847572980459,0.091576213509771;0.091576213509771,0.091576213509771,0.816847572980459...
    ;0.108103018168070,0.445948490915965,0.445948490915965;0.445948490915965,0.108103018168070,0.445948490915965;0.44594849091596,0.44594849091596,0.108103018168070];

%1D Gaussian quadrature points and weights
weights1D = [128/225,(322+13*sqrt(70))/900,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900,(322-13*sqrt(70))/900];
points1D = [1/2,1/2;1-1/6*sqrt(5-2*sqrt(10/7))-1/2,1/6*sqrt(5-2*sqrt(10/7))+1/2;1+1/6*sqrt(5-2*sqrt(10/7))-1/2,-1/6*sqrt(5-2*sqrt(10/7))+1/2;1-1/6*sqrt(5+2*sqrt(10/7))-1/2,1/6*sqrt(5+2*sqrt(10/7))+1/2;1+1/6*sqrt(5+2*sqrt(10/7))-1/2,-1/6*sqrt(5+2*sqrt(10/7))+1/2];

%Define the geometry info structure
GI = GenerateGeometry(s,l,dn,ds,minvert,false);

%Obtain pair potential info
if isfile('DeltaInfo.mat')
     %load pair potential info from external file if it exists
     load('DeltaInfo.mat')
else
     %calculate pair potential info if it does not exist
     PhiPath = strcat(mydir(1:idcs(end)-1),sep(sepindx+1),'Phi-Parametrization');
     addpath(PhiPath);
     rmpath(mydir);
     cd(PhiPath);
     [~,~,~,deltas1,deltas2]=SolveUsadel(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,phase,ksi,20,maxit,tol,false);
     addpath(mydir);
     rmpath(PhiPath);
     cd(mydir);
     save('DeltaInfo.mat','deltas1','deltas2');
end

%solve the Usadel equation at energy E
[thetas1,chis1,thetas2,chis2,thetan,chin,success] = RealEnergy(GI,points,weights,points1D,weights1D,T,gamma_B,gamma,phase,ksi,E,deltas1,deltas2,maxit,tol,false);

%Plot the density of states throughout the junction
figure
hold on
trisurf(GI.convs1(GI.elmats1),GI.xs1,GI.ys1,real(cos(thetas1)))
trisurf(GI.convs2(GI.elmats2),GI.xs2,GI.ys2,real(cos(thetas2)))
trisurf(GI.convn(GI.elmatn),GI.xn,GI.yn,real(cos(thetan)))
xlabel('x')
ylabel('y')
title('DOS (N(E)/N_0)')
colorbar