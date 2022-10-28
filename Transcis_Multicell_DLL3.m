% This program is adapted from Lateral Inhibition Tutorial (2014) 
% by Formosa-Jordan and Sprinzak. For more details please see
% Formosa-Jordan, P., Sprinzak, D. Modeling Notch signaling: a practical tutorial. 
% Methods Mol Biol. 2014;1187:285-310. doi: 10.1007/978-1-4939-1139-4_22.

function [yout,tout,params,F] = transcis_multicell_LI(params)

% DLL3 is included to transcis_multicell_LI from Formosa-Jordan and Sprinzak 2014, 
% which simulates trans-annihilation with cis-inactivation in a hexagonal lattice. 
% The structure params contains the model parameters of the system. 
% TOUT is a vector containing the time points of the solution between 0 and Tmax. 
% YOUT is a matrix containing the numerical solution for each variable for each time point. 
% Each row in YOUT is a vector of the size of TOUT. F is a movie showing the simulation. 

Tmax=200; tspan=[0 Tmax]; % set time for simulation

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

% show time traces of two cells with lateral inhibition

plot2cells(tout,yout,k)

% show lattice simulation

F=movielattice(tout,yout,P,Q,k);

function dy = li(t,y,params) 

betaD=params.betaD;% Parameters for DLL1
betaN=params.betaN;% Parameters for notch receptor
betaR=params.betaR;% Parameters for HES1
betaD3=params.betaD3;% Parameters for DLL3
m=params.m;% Cooperativity parameter m
M=params.connectivity;% Connectivity matrix used to identify neighboring cells
k=length(M);
kc=params.kc;%
kt=params.kt;%
Y=params.Y;
ke=params.ke;
kdr=params.kdr;
kf=params.kf;
Ys=params.Ys;
p=params.p;% Cooperativity parameter p
krs=params.krs;
Yr=params.Yr;

D = y(1:k);         % levels of DLL1 (representative classical notch ligant) in cells 1 to k
R = y(k+1:2*k);     % levels of HES1 (Repressor in Formosa-Jordan and Sprinzak, 2014) in cells 1 to k
N = y(2*k+1:3*k);   % levels of notch receptors in cells 1 to k
D3 = y(3*k+1:4*k);  % levels of DLL3 in cells 1 to k

Dneighbor=M*y(1:k);       % DLL1 level in the neighboring cells
Nneighbor=M*y(2*k+1:3*k); % Notch receptor level in the neighboring cells

% Differential equations for notch receptor, DLL1, DLL3, and HES1 
dN = betaN - Y.*N - N.*Dneighbor./kt - N.*D3./ke- N.*D./kc;
dD = betaD*kdr./(kdr + R.^m) - Y.*D - Nneighbor.*D./kt - N.*D./kc;
dD3 = betaD3 - 1*Y.*D3- N.*D3./ke ;
dR = betaR.*(N.*Dneighbor./(kt*Ys)).^p./(krs + (N.*Dneighbor./(kt*Ys)).^p) - Yr.*R;
dy = [dD;dR;dN;dD3];

function params=defaultparams

params.Y=1;
params.ke=.1;
params.kdr=1;
params.kf=.1;
params.Ys=1;
params.p=3; % Cooperativity value p
params.krs=300000;
params.Yr=1;
params.betaD=15; % Delta production rate 
params.betaN=20; % Notch production rate
params.betaR=1000000;
params.betaD3=25; % DLL3 production rate
params.m=1; % Cooperativity value m
params.sigma=0.1;
params.kc=0.1; %
params.kt=1; %
params.P=6;
params.Q=6;


function M=getconnectivityM(P,Q)

k=P*Q; % number of cells
M=zeros(k,k); % connectivity matrix
w=1/6; % weight for interactions

% calculating the connectivity matrix
for s=1:k
    kneighbor=findneighborhex(s,P,Q); % finds the neighbors of cells
    for r=1:6
        M(s,kneighbor(r))=w;
    end
end

function y0=getIC(params,k)

U=rand(k,1) - 1/2; % a uniform random distribution
epsilon=1e-5;   % multiplicative factor of Delta initial condition
D0=epsilon*params.betaD.*(1 + params.sigma*U); % initial DLL1 levels 
R0=zeros(k,1);  % initial HES1 levels
N0=params.betaN.*(1 + params.sigma*U);  % initial notch receptor levels 
D30=epsilon*params.betaD3.*(1 + params.sigma*U); % initial DLL3 levels
y0=[D0;R0;N0;D30];  % vector of initial conditions

function plot2cells(tout,yout,k) % Plotting time trace of cells

figure(1)
clf
for i=1:k
    subplot(1,1,1)
    semilogy(tout,yout(:,k+i),'Color',[0.4660 0.6740 0.1880],'linewidth',0.2)
    hold on
    ylim([1e-5,1e2])
    title(['cell #',num2str(i)])
    xlabel('time [a.u]'); ylabel('concentration [a.u]')
    legend('Hes1')
    set(gca,'FontSize',15)
end

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

function plotHexagon(p0,q0,c)

% This function plots a hexagon centered at hex lattice coordinates p,q

s32 = sqrt(3)/4;

q = q0*3/4;
p = p0*2*s32;

if q0/2 == round(q0/2),
   p = p+s32;
end;

x(1) = q-.5; x(2) = q-.25; x(3) = q+.25; x(4) = q+.5; x(5) = q+.25; x(6) = q-.25;

y(1) = p ; y(2) = p+s32; y(3) = p+s32; y(4) = p; y(5) = p-s32; y(6) = p-s32;

c=min(c,ones(1,3));

patch(x,y,c,'linewidth',2);

function F=movielattice(tout,yout,P,Q,k) % This function shows a time course movie of cell lattice

figure(4)
sy1 = sort(yout(end,1+k:2*k));
Norm = sy1(round(length(sy1)*0.95)); % find the Delta level in the high Delta cells
frameind=0;
for tind = 1:5:length(tout),   % shows every 5th frame
    clf;
    for i = 1:P,
        for j = 1:Q,
            ind = pq2ind(i,j,P);
            mycolor = min([yout(tind,ind+k)/ Norm,1]); % defined the normalized color of cell
            plotHexagon(i,j,[1-mycolor,1-mycolor,1]);
        end;
    end;
    axis image; axis off; box off;
    
    frameind=frameind+1;
    F(frameind) = getframe; % generates a movie variable
end;
movie2avi(F,'movielattice'); % save movie in avi format
