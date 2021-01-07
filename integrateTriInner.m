function I = integrateTriInner(a, b, tau, w)

% Stützstellen
tau = tau*(b-a) + a;
% Gewichte
w = w*(b-a);
% Summe bilden
I = sum(w);
