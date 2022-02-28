clear;
clc;
%% ============ Computational Methods: 1-D Heat Equation =============== %%
% The code is iterated to solve for 1-D heat equation using explicit and
% implicit method. For explicit methodology, second order central
% difference in space and forward Euler in time is considered. 
%
% Case: T(x,t=0) = 6*sin(pi*x/L)
%       T(x=0,t) = 0
%       T(x=1,t) = 0
%       alpha    = 1
%
% Tabulated Variables:
%
% alpha     = Diffusivity Constant (m^2/sec)
% beta      = Itervative Constant
% dt        = Time Step
% dx        = Space Step
% i,j       = Domain Coordinates
% L         = Space Domain (m) 
% R_rms     = Root mean square residuals
% Td        = Time Domain (sec)
% tnstep    = Number of Divisions in time domain
% xnstep    = Number of divisions in space domain
%
% ======================================================================= %

% Defining the input variables:
dt      = 0.003;
dx      = 0.1;
alpha   = 1;
beta    = alpha*dt/((dx)^2);
Td      = 1;
L       = 1;

% Spaced steps:
xnstep   = (L/dx) + 1;
tnstep   = round((Td/dt)) + 1;
T        = zeros(xnstep,tnstep);
x        = linspace(0,1,xnstep);
T(1,:)   = 0;
T(:,1)   = 6*sin((pi*x)/L);
T(xnstep,:)  = 0;

% Iteration procedure:
for j = 2:tnstep
    for i = 2:xnstep-1
    Tn(i,j) = T(i,j-1)+ beta*(T(i+1,j-1)-2*T(i,j-1)+T(i-1,j-1));
    T(i,j)  = Tn(i,j);
    end
end

Temp    = T';

% Estimating the residuals of the iteration
for j = 2:tnstep
    for i = 2:xnstep-1
    R(i,j)   = Tn(i,j) - T(i,j-1)+ beta*(T(i+1,j-1)-2*T(i,j-1)+T(i-1,j-1));
    R_rms(j) = sqrt(sum((R(i,j))^2)/(xnstep-1)); 
    end
end

time_step   = [0:dt:1];

figure(1)
hold on
contourf(Temp,40,':','Linewidth',1)
colormap turbo
colorbar
set(gca,'XTick',[], 'YTick', [])
xlabel('Length Axis (L = 1 m); \Deltax = 0.1 m')
ylabel('Time Axis (T = 1 sec); \Deltat = 0.007 sec')
hold off

figure(2)
hold on
plot(time_step,R_rms,'--','Linewidth',3)
xlabel('Time Axis (T = 1 sec)')
ylabel('Residuals: Root Mean Square (R_{RMS})')
grid on
hold off