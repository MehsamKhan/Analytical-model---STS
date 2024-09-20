%% Code for determining maximum axial stress vs effective penetration length
% Author : Mehsam Khan
% email: mehsam@u.northwestern.edu

clc
clear all

% Effective penetration length range (mm)
leff=0:1:500; % Input parameter [starting number (keep as 0): step size (≥ 1):
% effective penetration length (specify the screw
% effective penetration length)]

for ii=1:length(leff)
    % Inner core diameter of screw (mm)
    d=5; % Input parameter

    % Outer core diameter of screw (mm)
    do=8; % Input parameter

    % Shear stiffness parameter (Can be solved using gammasolver.m)
    gam=20.21; % Input parameter

    As=pi/4*(d)^2;

    % Young's modulus of screw (MPa)
    Es=208.2*1000; % Input parameter

    % Expected axial load (N)
    P=5*1000; % Input parameter
    
    % Effective area of wood (mm^2)
    Aweff1=6*15*do^2-(pi/4)*do^2; % For the connection condition described ...
    % in section 5 of manuscript
    Aweff2=pi/2*((leff(ii)/6+do/2)^2-(do/2)^2);
    
    % Effective swelling coefficient of wood in longitduincal direction of screw
    alpha=0.0029; % Input parameter

    % Moisture content change percentage (in whole numbers)
    MC=6; % Input parameter

    % Young's modulus of wood in longitduinal direction of screw (MPa)
    Ew=620; % Input parameter

    ks=sqrt((4/(d*Es)+(pi*d)/(Aweff2*Ew))*gam);
    beta=1/(As*Es)+1/(Aweff1*Ew);
    omega=sqrt(pi*d*gam*beta*leff(ii)*leff(ii));
    x1=0:.1:leff(ii);
    x2=0:.1:leff(ii)/2;
    x3=leff(ii)/2:.1:leff(ii);
    fffx=(4*P)/(pi*d*d)*sinh(omega-(omega*x1/leff(ii)))/sinh(omega);
    fffm1=(4*alpha*MC*gam)/(ks*ks*d)*(1-exp(-ks*x2));
    fffm2=(4*alpha*MC*gam)/(ks*ks*d)*(1-exp(ks*(x3-leff(ii))));
    k=length(x2);
    m=length(x1);
    fffm_total=zeros(1,length(fffx));
    fffm_total=fffm1;
    for i=1:k-1
        fffm_total(1,i+k)=fffm2(i+1);
    end
    ffinal=fffx+fffm_total;
    plot(x1,ffinal);
    hold on;
    plot(x1,fffx);
    plot(x1,fffm_total);
    stress=transpose(ffinal);
    maxstress_values(ii)=max(stress);
end

title('Plot of stress Distribution for different effective penetration lengths')
xlabel('Length along screw-axis (mm)')
ylabel('Stress (MPa)')

figure()
plot(leff,maxstress_values)
title('Maximum axial stress  at different penetration lengths')
xlabel('L_e_f_f (mm)')
ylabel('Maximum Stress (MPa)')
