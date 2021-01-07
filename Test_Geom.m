clc
close all
clear

%% Modell Platte mit Loch

% Rechteck
h = 3;
l = 4;
thickness = 20;
% Loch
r = 1;
x_m = 2;
y_m = h/2;

E = 210e3;
v = 0.3;

Model2D = 2;            % 1: Plain Strain 2: Plain Stress
%% Netzerstellung

% Vernetzungsparameter
h1 = 0.1;
% Randpunkte festlegen:

% Rechteck
n_hor = round(l/h1)+1;
n_vert = round(h/h1)+1;
K1 = [linspace(0, l, n_hor)', zeros(n_hor, 1)];
K2 = [zeros(n_vert, 1), linspace(0, h, n_vert)'];
K3 = [linspace(0, l, n_hor)', zeros(n_hor, 1)+h];
K4 = [zeros(n_vert, 1)+l, linspace(0, h, n_vert)'];

% Kreis
h2 = 0.1;
nt = round(2*pi*r/h2)+1;
t = linspace(0, 2*pi, nt);
P_Kreis = [r*cos(t)' + x_m, r*sin(t)' + y_m];

% Kanten und Eckpunkte definieren
n_out = 2*n_hor+2*n_vert;
Bnodes = [K1; K4; flipud(K3); flipud(K2); P_Kreis(1:end-1,:)];
Bedges = [(1:1:n_out-1)', (2:1:n_out)'; n_out, 1;
     (n_out+1:1:nt+n_out-2)', (n_out+2:1:nt+n_out-1)';
     nt+n_out-1, n_out+1];
% Bnodes = [K1; K4; flipud(K3); flipud(K2)];
% Bedges = [(1:1:n_out-1)', (2:1:n_out)'; n_out, 1];

hfun = h1;
% Venetzung erstellen
[nodes, etri, elements,tnum] = refine2(Bnodes, Bedges, [], [], hfun);
%[nodes, etri, elements,tnum] = smooth2(nodes, etri, elements,tnum);
% indices von Randknoten finden

iB = ismember(nodes, Bnodes, 'rows');

plotmeshing(nodes, elements, Bedges, Bnodes);

%% Assemblierung

order = 4;
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

lambda = lambda(order,1:order);
tau = tau(order,1:order);



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
        idxU2 = find(nodes(:,1) == 0);       % alles an Linker Kante feste Einspannung
        idxF1 = find(nodes(:,1) == l);
        idxF2 = find(nodes(:,2) == h);
        idxF3 = find(nodes(:,1) == l & nodes(:,2) == h);

        u(idxU2*2) = 0;          % y
        u(idxU2*2-1) = 0;        % x
    case 'quadratic'
        idxU1 = find(nodesQuad(:,1) == 0);       % alles an Linker Kante feste Einspannung
        idxU2 = find(nodesQuad(:,1) == l);       % alles an Linker Kante feste Einspannung
        idxF1 = find(nodesQuad(:,1) == l);
        idxF2 = find(nodesQuad(:,2) == h);
        idxF3 = find(nodesQuad(:,1) == l & nodesQuad(:,2) == h);

        u(idxU1*2) = 0;          % y
        u(idxU1*2-1) = 0;        % x
        
        u(idxU2*2) = 0;          % y
        u(idxU2*2-1) = 0;        % x
end

% Kraftrandbedingungen

%f(idxF1*2-1) = 10e3/length(idxF1);        % x
%f(idxF2*2-1) = -1e3/length(idxF2);
f(idxF2*1) = -1e3/length(idxF2);        % y

%% Lösen

tic
idxU = [idxU1*2-1; idxU1*2; idxU2*2-1; idxU2*2];
[Kred, P] = ReduceSparseMatrix(Kges, idxU);
    
rhs = P*(f - Kges*u);

% L = ichol(Kred);
% ured = pcg(Kred, rhs, 1e-03, 500, L, L');
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