function I = integrate(f, a, b, tau, w)

% St�tzstellen
tau = tau*(b-a) + a;
% Gewichte
w = w*(b-a);
% Summe bilden
I = w*f(tau)';

