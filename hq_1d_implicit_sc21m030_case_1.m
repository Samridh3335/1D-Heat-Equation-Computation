clear;
clc;
%% ============ Computational Methods: 1-D Heat Equation =============== %%
% The code is iterated to solve for 1-D heat equation using explicit and
% implicit method. For explicit methodology, second order central
% difference in space and forward Euler in time is considered. For implicit
% scheme, Thomas algorithm is used.
%
% Case: T(x,t=0) = 20
%       T(x=0,t) = 20
%       T(x=1,t) = 100
%       alpha    = 1
%
% Tabulated Variables:
%
% alpha     = Diffusivity Constant (m^2/sec)
% beta      = Itervative Constant
% dt        = Time Step
% deltax    = Space Step
% i,j       = Domain Coordinates
% L         = Space Domain (m) 
% Rrms      = Root mean square residuals
% Temp      = Temperature variation between space and time (K)
% t         = Time Domain (sec)
%
% ======================================================================= %

% Defining the input variables:
alpha   = 1;
beta    = 0.8;
deltaX  = 0.25;
t       = 1;
L       = 1;

% Defining Iteration Scheme:
T0  = 20*ones(1+round((alpha*t)/(beta*((deltaX)^2))),1);
T1  = 100*ones(1+round((alpha*t)/(beta*((deltaX)^2))),1);
T11 = 20*ones(1+round(L/deltaX),1);
n1  = length(T0);
n2  = length(T11);
T(:,1)  = T11;
b(1) = 1; b(n2) = 1; a(n2) = 0; a(1) = 0; c(1) = 0; c(n2) = 0;
d1(n2) = 100; d1(1) = 20;

for i = 1:((n2)-2)
 b(i+1)  = ((2*beta)+1);
 a(i+1)  = -beta;
 c(i+1)  = -beta;
 d1(i+1) = T11(i);
end

d(:,1) = d1';
for i = 2:n1
    d(:,i) = thomasAlgo(a',b',c',d(:,i-1));
end
T       = d;
Temp    = T';


% Solving for residuals:
for i=2:n2-1
    for j=1:n1-1
        R(i,j)=((T(i,j+1)-T(i,j))/((alpha*((deltaX)^2))/beta))-((((T(i+1,j))-(2*(T(i,j)))+(T(i-1,j)))/((deltaX)^2))*(alpha)); 
        Rrms(j)=sqrt((sum((R(i,j))^2))/((n2)-1));
    end
end
t   = linspace(0,1,n1-1);

figure(1)
hold on
contourf(Temp,40,':','Linewidth',1)
colormap turbo
colorbar
set(gca,'XTick',[], 'YTick', [])
xlabel('Length Axis (L = 1 m); \Deltax = 0.25 m')
ylabel('Time Axis (T = 1 sec); \Deltat = 0.05 sec')
hold off

figure(2)
hold on
plot(t,Rrms,'--','Linewidth',3)
xlabel('Time Axis (T = 1 sec)')
ylabel('Residuals: Root Mean Square (R_{RMS})')
grid on
hold off


% Defining Thomas Algorithm
function T = thomasAlgo(a,b,c,d)
N       = length(b);
beta    = b(1);
gamma(1)= d(1)/beta(1);
for i = 2:1:N
    beta(i)     = b(i)-a(i)*c(i-1)/beta(i-1);
    gamma(i)    = (d(i)-a(i)*gamma(i-1))/beta(i);
end
T(N) = gamma(N);
for j = N-1:-1:1
    T(j) = gamma(j)-c(j)*T(j+1)/beta(j);
end
end

