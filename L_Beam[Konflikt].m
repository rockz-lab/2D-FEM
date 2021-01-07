clc
%close all
clear

%% Geometrieparameter

h = 10;
w = 10;
r = 3;
th = 2;
tv = 2; 

laux1 = w-th-r;
laux2 = h-tv-r;

x_m = th + r;
y_m = tv + r;


%% Material


E = 210e3;
v = 0.3;

Model2D = 1;    % 1: Ebener Verzerrungszustand
                % 2: Ebener Spannungszustand

%% Vernetzung:

% Randlinien erstellen

hmesh = 0.1;

% Linien ausgehend vom Ursprung zeichnen

n1 = round(w/hmesh)+1;
n2 = round(tv/hmesh)+1;
n3 = round(laux1/hmesh)+1;
n4 = round(laux2/hmesh)+1;
n5 = round(th/hmesh)+1;
n6 = round(h/hmesh)+1;

L1 = [linspace(0, w, n1)', zeros(n1, 1)];
L2 = [zeros(n2, 1) + w, linspace(0, tv, n2)'];
L3 = [linspace(w, w-laux1, n3)', zeros(n3, 1) + tv];
L4 = [zeros(n4, 1) + th, linspace(tv+r, h, n4)'];
L5 = [linspace(th, 0, n5)', zeros(n5, 1) + h];
L6 = [zeros(n6, 1), linspace(h, 0, n6)'];

nt = round(2*pi*r/hmesh/4)+1;
t = linspace(-pi/2, -pi, nt);
P_Kreis = [r*cos(t)' + x_m, r*sin(t)' + y_m];

Bnodes = [L1; L2; L3; P_Kreis; L4; L5; L6];
nges = n1+n2+n3+n4+n5+n6+nt;
Bedges = [(1:nges-1)', (2:nges)'; nges, 1];
        
[nodes, etri, elements,tnum] = refine2(Bnodes, Bedges, [], [], hmesh);


iB = ismember(nodes, Bnodes, 'rows');

plotmeshing(nodes, elements, Bedges, Bnodes);


Bfun = 'linear';


% ebener Verzerrungszustand

tic

if Model2D == 1
    % Ebener Verzerrungszustand
    C = E/(1+v)/(1-2*v)*[1-v, v, 0
             v, 1-v, 0
             0, 0, (1-2*v)/2]; 
    
elseif Model2D == 2
    % Ebenber Spannungszustand
    C = E/(1-v^2)*[1, v, 0
                 v, 1, 0
                 0, 0, (1-v)/2]; 
end
    

load gauss_lambda.mat
load gauss_tau.mat
order = 1;
lambda = lambda(order,1:order);
tau = tau(order,1:order);

Nn = length(nodes)*2;
Ne = length(elements);

switch Bfun
    case 'quadratic'
        
        [nodesQuad, elementsQuad] = createMidsideNodes(nodes, elements);
        Kges = assemble(nodesQuad, elementsQuad, C, lambda, tau, Bfun);
    case 'linear'

        Kges = assemble(nodes, elements, C, lambda, tau, Bfun);
end
toc
%% Randbedingungen festlegen

u = zeros(Nn, 1);
f = zeros(Nn, 1);

% Verschiebungsrandbedingungen


idxU = find(nodes(:,1) == 0);       % alles an Linker Kante feste Einspannung

u(idxU*2) = 0;          % y
u(idxU*2-1) = 0;        % x
% Kraftrandbedingungen

idxF1 = find(nodes(:,1) == w);
%idxF2 = find(nodes(:,2) == h);

%f(idxF1*2) = -5e1;          % y
f(idxF1*2) = -0.5e1;          % y
%f(idxF1*2-1) = -1e2;        % x

%% Lösen

tic
idxU = [idxU*2-1; idxU*2];
[Kred, P] = ReduceSparseMatrix(Kges, idxU);
    
rhs = P*(f - Kges*u);


%ured = bicgstab(Kred, rhs, 1e-03, 500);
ured = Kred\rhs;

u = P'*ured + u;

toc


%% Plotten & Post

% Verformte Geometrie
scale = 1e1;

nodesDef = nodes + scale*[u(1:2:Nn), u(2:2:Nn)];
%Randknoten:
BnodesDef = nodesDef(iB,:);

% Spannungen

sigma = evaluateStress(nodes, elements, u, C);

% Gesamtverschiebung

utot = sqrt(u(1:2:Nn).^2 + u(2:2:Nn).^2);

plot_rst(nodes, elements, Bedges, Bnodes, etri, nodesDef, BnodesDef, sigma, utot)



