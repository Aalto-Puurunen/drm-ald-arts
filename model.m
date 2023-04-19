% General ALD growth model based on Yanguas-Gil et al., Chemical Vapor
% Deposition 18, 46-52 (2012). See Arts et al., J. Phys. Chem. C 123, 44, 
% 27030–27035 (2019) for further details.

clear; clc; close all;

%% ==== Definitions =======================================================
% n:     number density of gas-phase reactant species (m-3)
% t:     dosing time (s)
% z:     distance into high-AR structure (m)
% L:     lenth of high-AR structure (m)
% S/V:   surface to volume ratio of high-AR structure (m-1)

% H:     gap height of semi-infinite trech (PillarHall structure) (m) 
%        In this geometry, S/V =approx (2*width*L)/(H*width*L) = 2/H
% AR:    Aspect ratio (AR = L/H) of semi-infinite trench (-)

% s0:    initial sticking probability
% r:     surface recombination probability

% v_th:  mean thermal velocity of reactant molecules (m/s) 
%        Calculated as v_th = sqrt(8*k_b*T/pi*M)  

% D:     diffusion coefficient (m2/s)
%        Calculated as D=(2/3)*v_th*H [KNUDSEN APPROXIMATION USED]


% A0:    average surface area per absorption site (m2)
% theta: fraction of available adsorption sites (-)
%        [SO HERE THE SURFACE COVERAGE EQUALS 1-theta]              (!)

% // NORMALIZED, DIMENSIONLESS PARAMETERS //

% x = n/n0:         normalized reactant density (-)
% eta = z/L:        scaled distance into high-AR structure (-)
% tau = tD/L^2:     dimensionless dosing time (-)

% gamma: Excess number: number of gas-phase reactant species present in
%        high-AR struc-ture per adsorption site (-) 
%        Calculated as gamma = V*n0*s0/S  (<< 1 under typical conditions)

% eta:   Ratio between surface collision rate and diffusion rate of
%        gas-phase reactant molecules (-)
%        Calculated as eta = (1/4)*L^2*(S/V)*v_th/D

% alpha = eta*s0:   Ratio between adsorption rate and diffusion rate of 
%                   gas-phase reactant molecules (-)
% nu = eta*r:       Ratio between recombination rate and diffusion rate of 
%                   gas-phase reactant molecules (-)

% For molecular diffusion in a semi-infinite trench:
% alpha = (3/4)*(L/H)^2*s0 = (3/4)*AR^2*s0
% nu    = (3/4)*(L/H)^2*r = (3/4)*AR^2*r

%% ==== CALCULATION ======================================================= 
% // pdepe (see documentation): ODE solver 2nd order partial diff. eqs. //
% // Used functions are defined at bottom of script //////////////////////

% // INPUT VALUES
s0 = 1E-2;     % SET VALUE Initial sticking probability (-)
r = 0;         % SET VALUE Surface recombination probability (-)
t_max = 0.1;     % SET VALUE Dosing time (or pulse length) (s) 
T = 523.15;    % SET VALUE Temperature (K) 
pA0 = 100;      % SET VALUE Partial pressure of reactant at z=0 (Pa) 
MA = 0.1;      % SET VALUE Molar mass of reactant(kg/mol) 
q = 4E18;      % SET VALUE Adsorption capacity (m-2) 
H = 500E-9;    % SET VALUE Cavity height (m)
AR = 600;      % SET VALUE Aspect ratio (-)

% // Not used in Knudsen approximation
dA = 0.6E-9;   % NOT USED HERE Diameter of reactant (m) 
pI=50;         % NOT USED HERE Partial pressure of inert gas (Pa)
MI = 0.028;    % NOT USED HERE Molar mass of inert gas (kg/mol)
dI = 0.373E-9; % NOT USED HERE Diameter of inert gas (m)
Pd = 0.01;     % NOT USED HERE Desorption probability per unit time (s-1) 

% // Constants
R = 8.31446;   % Gas constant J/(mol*K); 
NA = 6.022E23; % Avogadro's number (mol-1);   

% // Calculated using set values
L = H*AR;                       % << calculated using set values
v_th = sqrt((8*R*T)/(pi*MA));   % << calculated using set values 
n0 = pA0*NA/(R*T);              % << calculated using set values   
A0 = 1/q;                       % << calculated using set values   
gamma = n0*A0*H/2;              % << calculated using set values 
D = (2/3)*v_th*H;               % << calculated using set values       (!)
tau_max = t_max*D/L^2;          % << calculated using set values
alpha = (3/4)*s0*AR^2;          % << calculated using set values 
nu = (3/4)*r*AR^2;              % << calculated using set values 
n_rec = 1;                      % Fixed at 1 (first order recombination)
setGlobal(AR,s0,gamma,alpha,nu,n_rec)

% // Number of spatial and temporal points in numerical solution
N_tau = 400;                    % SET VALUE
zN = L/200;                     % SET VALUE Distance between points (m)
%zN = 0.25E-6;                  % Used by Jihong (m)
N_AR = floor(L/zN);             % << calculated using set values
m = 0;                          % Fixed at 0 for slab geometry
eta = linspace(0,1,N_AR);   
tau = linspace(0,tau_max,N_tau);   

% // Calculate solution //////////////////////////////////////////////////
sol = pdepe(m,@pdefun,@icfun,@bcfun,eta,tau);   
x = sol(:,:,1);
theta = sol(:,:,2);
Coverage = 1-theta;
pA = x*n0*R*T/NA;

%% ==== EXTRACTED s0 ======================================================

for iz = 2:N_AR-1
   if Coverage(N_tau,iz)>(1/2) && Coverage(N_tau,iz+1)<(1/2)
       Slope_HTPD=(Coverage(N_tau,iz+1)-Coverage(N_tau,iz-1))/...
                  (AR*(eta(iz+1)-eta(iz-1)));
   end
end
s0_out = 13.9*Slope_HTPD^2;

disp(['s0_in = ',num2str(s0)]);
disp(['s0_out = ',num2str(s0_out)]);
disp(['s0_out / s0_in = ',num2str(s0_out/s0)]);

%% ==== PLOTS ============================================================= 

% // Surface coverage profile reached at t_max
figure (1);
plot(eta*AR,Coverage(N_tau,:),'k','LineWidth',3);
ylabel('Coverage'); ylim([0 1]);xlim([0 1000]);
xlabel('Distance / Cavity Height')
set(gca,'FontSize',24);
set(gca,'LineWidth',3);

% Reactant density profile reached at t_max
figure (2);
plot(eta*AR,pA(N_tau,:),'k','LineWidth',3);
ylabel('Partial pressure (Pa)'); %ylim([0 15]);xlim([0 1000]);
xlabel('Distance / Cavity Height');
set(gca,'FontSize',24);
set(gca,'LineWidth',3);

%% ==== SAVE OUTPUT =======================================================
% s0 = 1E-2;     % SET VALUE Initial sticking probability (-)
% r = 0;         % SET VALUE Surface recombination probability (-)
% t_max = 2;     % SET VALUE Dosing time (or pulse length) (s) 
% T = 523.15;    % SET VALUE Temperature (K) 
% pA0 = 10;      % SET VALUE Partial pressure of reactant at z=0 (Pa) 
% MA = 0.1;      % SET VALUE Molar mass of reactant(kg/mol) 
% q = 4E18;      % SET VALUE Adsorption capacity (m-2) 
% H = 500E-9;    % SET VALUE Cavity height (m)
% AR = 1000;     % SET VALUE Aspect ratio (-)
clear X
X(:,1) = {'s0 (-)','r (-)','t_dose (s)','T (K)','pA0 (Pa)','MA (kg/mol)','q (m-2)','H (m)','AR=L/H (-)'};
X(:,2) = {num2str(s0),num2str(r),num2str(t_max),num2str(T),num2str(pA0),num2str(MA),num2str(q),num2str(H),num2str(AR)};
X(1,4) = {'Distance (m)'};
X(1,5) = {'Distance / Cavity Height'};
X(1,6) = {'Coverage'};
X(1,7) = {'Partial pressure (Pa)'};
for i = 1:length(eta)
    X(i+1,4)={num2str(eta(i)*AR*H)};
    X(i+1,5)={num2str(eta(i)*AR)};
    X(i+1,6)={num2str(Coverage(N_tau,i))};
    X(i+1,7)={num2str(pA(N_tau,i))};
end
xlswrite('Output_Arts.xlsx',X)

%% ==== FUNCTIONS =========================================================
% --------------------------------------------------------------
% // Components of the 2nd order partial differential equation
function [c,f,s] = pdefun(x,t,u,DuDx)          
para = getGlobal;
gamma = para(3);
alpha = para(4); 
nu = para(5);
n_rec = para(6);
c = [1; 1]; 
f = [1; 0].*DuDx; 
s = [-alpha*u(1)*u(2)-nu*u(1)^n_rec; -alpha*gamma*u(1)*u(2)];  
end
% --------------------------------------------------------------
% // Initial conditions
function u0 = icfun(x)                        
u0 = [0; 1];
end
% --------------------------------------------------------------
% // Boundary conditions
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
para = getGlobal;
AR = para(1);
alpha = para(4);         
pl = [ul(1)-1; 0]; 
ql = [0; 1]; 
pr = [(alpha/(2*AR))*ur(1)*ur(2); 0]; 
qr = [1; 1]; 
end
% --------------------------------------------------------------
% Set global variables to use in function 
function setGlobal(val_AR,val_s0,val_gamma,val_alpha,val_nu,val_n_rec)
global x_para
x_para = [val_AR,val_s0,val_gamma,val_alpha,val_nu,val_n_rec];
end
% --------------------------------------------------------------
% Get global variables (used in pdefun and bcfun)
function val_para = getGlobal
global x_para
val_para = x_para;
end
% ----------------------------------------------------------