function I = integrateTriInner(a, b, tau, w)

% St�tzstellen
tau = tau*(b-a) + a;
% Gewichte
w = w*(b-a);
% Summe bilden
I = sum(w);
