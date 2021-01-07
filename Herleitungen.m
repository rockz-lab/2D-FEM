clc
close all
clear

%% zur Herleitung von Ansatzfkt, etc.


%% Lineares Dreieckselement
A = [1, 0, 0
     1, 1, 0
     1, 0, 1];
 
invA = inv(A);

P1 = [1 0 0];
P2 = [0 1 0];
P3 = [0 0 1];

Nlin1 = A\P1';
Nlin2 = A\P2';
Nlin3 = A\P3';
x1 = 1
x2 = 2
x3 = 3
y1 = 0
y2 = 0
y3 = 4;
test = [0, 1, 0; 0, 0, 1]*invA*[[x1; x2; x3], [y1; y2; y3]];
test2 = (invA*[[x1; x2; x3], [y1; y2; y3]])' * [0, 1, 0; 0, 0, 1]';
inv(test)

% Ansatzfunktionen x, y,
P = [0, 0;
     2, 1;
     3, 1];
 M = [[1;1;1], P];
 alpha1 = M\P1';
 alpha2 = M\P2';
 alpha3 = M\P3';

%% quadratisches Dreieckselement

% a + b*xi + c*eta + d*xi*eta + e*xi^2 + d*eta^2
B = [1, 0, 0, 0, 0, 0;
     1, 1, 1, 0, 0, 0;
     1, 0, 0, 1, 1, 0;
     1, 0.5, 0.25, 0, 0, 0;
     1, 0.5, 0.25, 0.5, 0.25, 0.25;
     1, 0, 0, 0.5, 0.25, 0];
 
P1 = [1 0 0 0 0 0];
P2 = [0 1 0 0 0 0]; 
P3 = [0 0 1 0 0 0];
Pm1 = [0 0 0 1 0 0];
Pm2 = [0 0 0 0 1 0];
Pm3 = [0 0 0 0 0 1];
 
% Koeffizienten für die Ansatzfunktionen:

Nquad1 = B\P1';
Nquad2 = B\P2';
Nquad3 = B\P3';
Nquad4 = B\Pm1';
Nquad5 = B\Pm2';
Nquad6 = B\Pm3';



%% C-Matrix
%eps = sym(grad(u)) = 1/2(grad(u) + grad(u)T)
% 2D: u = ux ex + uy ey
% ebener Verzerrungszustand
% grad u = [dux/dx, dux/dy
%           duy/dx, duy/dy];
% eps = [dux/dx, 1/2*(dux/dy+ duy/dx)
%        1/2*(dux/dy+ duy/dx), duy/dy
% eps = [dux/dx, duy/dy, 1/2(dux/dy + duy/dx)]
% sigma = 2G(eps + v/(1-2*v)sp(eps)*I)
% sigma = [2G*dux/dx + 2*G*v/(1-2v)*(dux/dx + duy/dy), 2G*duy/dy + 2*G*v/(1-2v)*(dux/dx + duy/dy), 2G*1/2(dux/dy + duy/dx)]

E = 210e6;
v = 0.3;
G = E/(1+v)/2;
C = 2*G*[(1 + v/(1-2*v)), v/(1-2*v), 0
    v/(1-2*v), (1 + v/(1-2*v)), 0
    0, 0, 1];
C/E*(1+v)*(1-2*v)/(1-v)





