clc
close all
clear

% Element Points

P = [0, 0;
    2, 0;
    0, 2];

% Steifigkeitsmatrix

E = 210e6;
v = 0.3;
G = E/(1+v)/2;
C = 2*G*[(1 + v/(1-2*v)), v/(1-2*v), 0
    v/(1-2*v), (1 + v/(1-2*v)), 0
    0, 0, 1];   
load gauss_lambda.mat
load gauss_tau.mat
order = 7;
lambda = lambda(order,1:order);
tau = tau(order,1:order);

%K = integrateElement(P, C, lambda, tau);

f = @(x) ones(size(x));

x1 = P(1,1);
x2 = P(2,1);
x3 = P(3,1);
y1 = P(1,2);
y2 = P(2,2);
y3 = P(3,2);

w = lambda;
% detJ berechnen
dJ = abs((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));

I  = integrate(f, 0, 1, tau, w);
        % äußeres Integral
        f = @(eta) I.*(1-eta);
        A = dJ*integrate(f, 0, 1, tau, w);
