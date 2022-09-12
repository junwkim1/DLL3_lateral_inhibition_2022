% This program is adapted from Lateral Inhibition Tutorial (2014) 
% by Formosa-Jordan and Sprinzak. For more details please see
% Formosa-Jordan, P., Sprinzak, D. Modeling Notch signaling: a practical tutorial. 
% Methods Mol Biol. 2014;1187:285-310. doi: 10.1007/978-1-4939-1139-4_22.

function paramsearch_LI

% This code plots the log(Hmax), log(Hmin), and log(Hmax/Hmin)
params.Y=1;
params.ke=0.1;
params.kdr=1;
params.kf=0.1;
params.Ys=1;
params.p=3; %Cooperativity value p
params.krs=300000;
params.Yr=1;
params.betaD3=0;%20
params.betaH=1000000;
params.m=1; %Cooperativity value m
params.sigma=0.1;
params.kc=0.1; 
params.kt=1; 
params.P=6;
params.Q=6;
k=params.Q*params.P;
% variable parameters
betaD=logspace(-1,5,40); % creates a series of betaD from 0.1 to 100000
betaN=logspace(-1,5,40); % creates a series of betaN from 0.1 to 100000
ind=0; h=waitbar(0,'% of progress'); % generates a waitbar
for i=1:length(betaD)
    for j=1:length(betaN)
        params.betaD=betaD(i);
        params.betaN=betaN(j);
        ind=ind+1; waitbar(ind/(length(betaD)*length(betaN))) 
        [yout,tout] = transcis_multicell_LI(params); % calling the LI solver
        % finding max and min values of H
        Hmax(i,j)=max(yout(end,k+1:2*k));
        Hmin(i,j)=abs(min(min(yout(end,k+1:2*k))));
    end
end

close(h)

figure(1)
imagesc(log10(betaN),log10(betaD),log10(Hmin));
set(gca,'YDir','normal')
xlabel('log(\beta_N)','fontsize',14);
ylabel('log(\beta_D)','fontsize',14);
title('log(H_{min})','fontsize',14)
colorbar

figure(2)
imagesc(log10(betaN),log10(betaD),log10(Hmax./Hmin));
set(gca,'YDir','normal')
xlabel('log(\beta_N)','fontsize',14);
ylabel('log(\beta_d3)','fontsize',14);
title('log(H_{max}/H_{min})','fontsize',14)
colorbar

figure(3)
imagesc(log10(betaN),log10(betaD),log10(Hmax));
set(gca,'YDir','normal')
xlabel('log(\beta_d)','fontsize',14);
ylabel('log(\beta_d3)','fontsize',14);
title('log(H_{max})','fontsize',14)
colorbar


function [yout,tout,params,F] = transcis_multicell_LI(params)

% DLL3 is included to transcis_multicell_LI from Formosa-Jordan and Sprinzak 2014, 
% which simulates trans-annihilation with cis-inactivation in a hexagonal lattice. 
% The structure params contains the model parameters of the system. 
% TOUT is a vector containing the time points of the solution between 0 and Tmax. 
% YOUT is a matrix containing the numerical solution for each variable for each time point. 
% Each row in YOUT is a vector of the size of TOUT. 

Tmax=50; tspan=[0 Tmax]; % set time for simulation

if(nargin < 1)
    params=defaultparams; % get the default parameters if none provided
end

P=params.P;  % number of cells per column
Q=params.Q;  % number of columns - MUST BE EVEN
k=P*Q; % number of cells

% get the connectivity matrix
params.connectivity=getconnectivityM(P,Q);

% setting the initial conditions + noise
y0=getIC(params,k);

% run simulation with lateral inhibition
[tout,yout] = ode15s(@li,tspan,y0,[],params);

function dy = li(t,y,params) 

betaD=params.betaD;
betaN=params.betaN;
betaH=params.betaH;
betaD3=params.betaD3;
m=params.m;
M=params.connectivity;
k=length(M);
kc=params.kc;
kt=params.kt;
Y=params.Y;
ke=params.ke;
kdr=params.kdr;
kf=params.kf;
Ys=params.Ys;
p=params.p;
krs=params.krs;
Yr=params.Yr;

D = y(1:k);         % levels of Delta in cells 1 to k
R = y(k+1:2*k);     % levels of Hes1 (Repressor in Formosa-Jordan and Sprinzak, 2014) in cells 1 to k
N = y(2*k+1:3*k);   % levels of Notch in cells 1 to k
D3 = y(3*k+1:4*k);  % levels of DLL3 in cells 1 to k

Dneighbor=M*y(1:k);       % Delta level in the neighboring cells
Nneighbor=M*y(2*k+1:3*k); % Notch level in the neighboring cells

dN = betaN - Y.*N - N.*Dneighbor./kt - N.*D3./ke - N.*D./kc;
dD = betaD*kdr./(kdr + R.^m) - Y.*D - Nneighbor.*D./kt- N.*D./kc;
dD3 = betaD3 - 1*Y.*D3 - N.*D3./ke;
dR = betaH.*(N.*Dneighbor./(kt*Ys)).^p./(krs + (N.*Dneighbor./(kt*Ys)).^p) - Yr.*R;
dy = [dD;dR;dN;dD3];

function params=defaultparams

params.Y=1;
params.ke=2;
params.kdr=1;
params.kf=.1;
params.Ys=1;
params.p=3;
params.krs=300000;
params.Yr=1;
params.betaD=10;%
params.betaN=10; %
params.betaH=1000000;%
params.betaD3=350;
params.m=1;%
params.sigma=0.1;
params.kc=2; %
params.kt=1; %
params.P=12;
params.Q=12;


function M=getconnectivityM(P,Q)

k=P*Q; % number of cells
M=zeros(k,k); % connectivity matrix
w=1/6; % weight for interactions

% calculating the connectivity matrix
for s=1:k
    kneighbor=findneighborhex(s,P,Q); % finds the neighbors of cell s
    for r=1:6
        M(s,kneighbor(r))=w;
    end
end

function y0=getIC(params,k)

U=rand(k,1) - 1/2; % a uniform random distribution
epsilon=1e-5;   % multiplicative factor of Delta initial condition
D0=epsilon*params.betaD.*(1 + params.sigma*U); % initial Delta levels 
R0=zeros(k,1);  % initial repressor levels
N0=params.betaN.*(1 + params.sigma*U);  % initial Notch levels
D30=epsilon*params.betaD3.*(1 + params.sigma*U); % initial DLL3 levels
y0=[D0;R0;N0;D30];  % vector of initial conditions

function out = findneighborhex(ind,P,Q)
[p,q] = ind2pq(ind,P);

%above and below:
out(1) = pq2ind(mod(p,P)+1,q,P);
out(2) = pq2ind(mod(p-2,P)+1,q,P);

%left side:
qleft = mod(q-2,Q)+1;
qright = mod(q,Q)+1;

if q/2~=round(q/2),
    pup = p;
    pdown = mod(p-2,P)+1;
else 
    pup = mod(p,P)+1;
    pdown = p;
end;
out(3) = pq2ind(pup,qleft,P);
out(4) = pq2ind(pdown,qleft,P);
out(5) = pq2ind(pup,qright,P);
out(6) = pq2ind(pdown,qright,P);

function ind=pq2ind(p,q, P)
ind = p + (q-1)*P;

function [p,q]=ind2pq(ind, P)
q = 1+floor((ind-1)/P);
p = ind - (q-1)*P;
