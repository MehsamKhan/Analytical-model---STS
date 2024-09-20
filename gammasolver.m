%% Code for determining equivalent shear stiffness parameter from withdrawal stiffness
%  
% Author : Mehsam Khan
% email: mehsam@u.northwestern.edu
% Important: Symbolic MATLAB Toolbox required for this code to run

clc
clear all

syms gamm

% Young's modulus of screw (MPa)
Es=208200; % Input parameter

 % Inner core diameter of screw (mm)
d=5; % Input parameter

As=pi*d*d/4;

% Young's modulus of wood in longitduinal direction of screw (MPa)
Ew=530.64; % Input parameter

% Effective area of wood (mm^2)
Aw=16032; % Input parameter (will depend on the withdrawal test setup,
% refer to Figure 3 and Equation 5 in the manuscript)

% Effective penetration length of screw (mm)
leff=72; % Input parameter

% Withdrawal stiffness from screw withdrawal test (N/mm)
kw=17.41*1000; % Input parameter 
% [Multiply by pi*(outer diameter of screw)*(effective length) if unit is in force/length^3]

beta=1/(Es*As)+1/(Ew*Aw);

% Inital Guess for root finding algorithm
iguess=1; % Input parameter (Keep 1 as default)

gamma_e=vpasolve(kw==(pi*d*leff*gamm*tanh(sqrt(pi*d*gamm*beta*leff*leff)))/...
    (sqrt(pi*d*gamm*beta*leff*leff)),gamm,iguess);

fprintf('Equivalent shear stiffness parameter (MPa/mm): %f\n ', gamma_e)
