function I = integrateTriOuter(a, b, tau, w)

% St�tzstellen
tau = tau*(b-a) + a;
% Gewichte
w = w*(b-a);
% Summe bilden
I = w*(1-tau)';
