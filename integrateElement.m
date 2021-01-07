function K = integrateElement(P, C, t, tau, w, Bfun)

x1 = P(1,1);
x2 = P(2,1);
x3 = P(3,1);
y1 = P(1,2);
y2 = P(2,2);
y3 = P(3,2);
        
        
switch Bfun
    
    case 'linear'

        

        % detJ berechnen
      
        dJ = abs((x1-x3)*(y2-y3) - (x2-x3)*(y1-y3));
        % B-Vektor

        % Koeffizienten der Ansatzfunkion des Elements berechnen
%          M = [[1;1;1], P];
% 
%          alpha1 = M\[1; 0; 0];
%          alpha2 = M\[0; 1; 0];
%          alpha3 = M\[0; 0; 1];

        % Bi = (dNi/dx, 0; 0, dNi/y; dNi/dy, dNi/dx)
%         B = [alpha1(2), 0, alpha2(2), 0, alpha3(2), 0;
%              0, alpha1(3), 0, alpha2(3), 0, alpha3(3);
%              alpha1(3) , alpha1(2), alpha2(3) , alpha2(2), alpha3(3) , alpha3(2)];
        
        B = 1/dJ*[y2-y3, 0, y3-y1, 0, y1-y2, 0;
                0, x3-x2, 0, x1-x3, 0, x2-x1;
                x3-x2, y2-y3, x1-x3, y3-y1, x2-x1, y1-y2];

        F = B'*C*B;  
        K = zeros(6);


        for i = 1:6
            for j = 1:6
        %         f = @(xi) F(i,j)*ones(size(xi));
        %         % inneres Integral
        %         I = integrate(f, 0, 1, tau, w);
        %         % ‰uﬂeres Integral
        %         f = @(eta) I.*(1-eta);

                %I = integrateTriInner(0, 1, tau, w);
                K(i,j) =  t*1/2*F(i,j)*dJ;
            end
        end
        
    case 'quadratic'
        
        % evaluation at midside Points
        xi = 0;
        eta = 0.5;
        [F1, detJ] = integrandQuadTri (P, xi, eta, C);
        xi = 0.5;
        eta = 0;
        [F2, ~] = integrandQuadTri (P, xi, eta, C);
        xi = 0.5;
        eta = 0.5;
        [F3, ~] = integrandQuadTri (P, xi, eta, C);
        
        K = t*1/2*abs(detJ)/3*(F1 + F2 + F3);
end