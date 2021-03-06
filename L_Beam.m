clc
close all
clear

%% Geometrieparameter

h = 10;
w = 12;
r = 6;
th = 3;
tv = 3;
thickness = 10;

laux1 = w-th-r;
laux2 = h-tv-r;

x_m = th + r;
y_m = tv + r;


%% Material


E = 210e3;
v = 0.3;

Model2D = 2;    % 1: Ebener Verzerrungszustand
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


Bfun = 'quadratic';


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
order = 4;
lambda = lambda(order,1:order);
tau = tau(order,1:order);

Nn = length(nodes)*2;
Ne = length(elements);

switch Bfun
    case 'quadratic'
        
        [nodesQuad, elementsQuad] = createMidsideNodes(nodes, elements);
        Kges = assemble(nodesQuad, elementsQuad, C, thickness, lambda, tau, Bfun);
        
        Nn = length(nodesQuad)*2;
        Ne = length(elementsQuad);
    case 'linear'

        Kges = assemble(nodes, elements, C, thickness, lambda, tau, Bfun);
        
        Nn = length(nodes)*2;
        Ne = length(elements);
end


toc
%% Randbedingungen festlegen

u = zeros(Nn, 1);
f = zeros(Nn, 1);

switch Bfun
    
    case 'linear'
        idxU = find(nodes(:,1) == 0);       % alles an Linker Kante feste Einspannung
        idxF1 = find(nodes(:,1) == w);
        idxF2 = find(nodes(:,2) == h);
        idxF3 = find(nodes(:,1) == w & nodes(:,2) == h);

        u(idxU*2) = 0;          % y
        u(idxU*2-1) = 0;        % x
    case 'quadratic'
        idxU = find(nodesQuad(:,1) == 0);       % alles an Linker Kante feste Einspannung
        idxF1 = find(nodesQuad(:,1) == w);
        idxF2 = find(nodesQuad(:,2) == h);
        idxF3 = find(nodesQuad(:,1) == w & nodesQuad(:,2) == h);
        idxF4 = find(nodesQuad(:,2) == tv);
        
        u(idxU*2) = 0;          % y
        u(idxU*2-1) = 0;        % x
end

% Kraftrandbedingungen

%f(idxF1*2-1) = 10e3/length(idxF1);        % x
%f(idxF2*2-1) = -1e3/length(idxF2);
%f(idxF1*2) = +2e2/length(idxF1);        % y
f(idxF4*2) = -2e2/length(idxF4);        % y

%% L�sen

tic
idxU = [idxU*2-1; idxU*2];
[Kred, P] = ReduceSparseMatrix(Kges, idxU);
    
rhs = P*(f - Kges*u);

%L = ichol(Kred);
%ured = bicgstabl(Kred, rhs, 1e-03, 500, L, L');
ured = Kred\rhs;

u = P'*ured + u;

toc


%% Plotten & Post

% Verformte Geometrie
scale = 1e2;

switch Bfun
    
    case 'linear'
        nodesDef = nodes + scale*[u(1:2:Nn), u(2:2:Nn)];
        %Randknoten:
        BnodesDef = nodesDef(iB,:);

        % Spannungen

        sigma = evaluateStress(nodes, elements, u, C, Bfun);

        % Gesamtverschiebung

        utot = sqrt(u(1:2:Nn).^2 + u(2:2:Nn).^2);
        disp(['maximale Verschiebung:', num2str(max(utot)), ' mm'])
        plot_rst(nodes, elements, Bedges, Bnodes, etri, nodesDef, BnodesDef, sigma, utot)
        
    case 'quadratic'
        nCorner = length(nodes);
        nEdge = length(nodesQuad) - nCorner;

        nodesDef = nodesQuad + scale*[u(1:2:Nn), u(2:2:Nn)];
        %Randknoten:
        BnodesDef = nodesDef(iB,:);

        % Spannungen

        sigma = evaluateStress(nodesQuad, elementsQuad, u, C, Bfun, nEdge);

        % Gesamtverschiebung

        utot = sqrt(u(1:2:Nn).^2 + u(2:2:Nn).^2);
        disp(['maximale Verschiebung:', num2str(max(utot)), ' mm'])
        %plot_rst(nodes, elements, Bedges, Bnodes, etri, nodesDef, BnodesDef, sigma, utot)
        plot_rst_Quad(nodesQuad, elementsQuad, elements, Bedges, Bnodes, etri, nodesDef, BnodesDef, utot)

end





% Van Mieses Spannungen
if Model2D == 1
    % Ebener Verzerrungszustand (plain strain)
    sigmavM = sqrt((sigma(:, 1).^2 + sigma(:, 2).^2)*(v^2 - v + 1) - sigma(:, 1).*sigma(:,2)*(2*v^2 - 2*v - 1) + 3*sigma(:,3).^2); 
    
elseif Model2D == 2
    % Ebener Spannungszustand (plain stress)
    sigmavM = sqrt(sigma(:, 1).^2 + sigma(:, 2).^2 - sigma(:, 1).*sigma(:,2) + 3*sigma(:,3).^2);
end

switch Bfun
    
    case 'linear'
        figure
        set(gcf, 'Position', get(0, 'Screensize'));
        patch('faces',elements(:,1:3),'vertices',nodes, ...
            'FaceVertexCData',sigmavM,'FaceColor','flat', ...
            'facealpha', 1, ...
            'edgecolor',[.1,.1,.1],...
            'edgealpha', 0.1);
        hold on; axis image off;
        colorbar('southoutside'); 
        colormap jet
        title('GEH');
        ax = gca;
        ax.Clipping = 'off';
        
    case 'quadratic'
        
        figure
        set(gcf, 'Position', get(0, 'Screensize'));
        
        hold on; axis image off;
        colorbar('southoutside'); 
        colormap(jet(200))
        title('GEH');
        ax = gca;
        ax.Clipping = 'off';
        patch('faces',elementsQuad(:,1:3),'vertices',nodesQuad, ...
            'FaceVertexCData',sigmavM,'FaceColor','interp', ...
            'edgecolor','none',...
            'facealpha', 1);
        patch('faces',elements(:,1:3),'vertices',nodes, ...
            'facecolor','none', ...
            'facealpha', 1, ...
            'edgecolor',[.1,.1,.1],...
            'edgealpha', 0.1);
end