function sigma = evaluateStress(nodes, elts, u, C, Bfun, varargin)

Nn = length(nodes);
Ne = length(elts);



switch Bfun
    
    case 'linear'
        
        sigma = zeros(Ne, 3);
        % Spannungen konstant für jedes Element
        for i = 1:Ne

            ind = elts(i,:);
            P = nodes(ind, :);

            M = [[1;1;1], P];

            alpha1 = M\[1; 0; 0];
            alpha2 = M\[0; 1; 0];
            alpha3 = M\[0; 0; 1];

            % Bi = (dNi/dx; dNi/y; 1/2*(dNi/dx+dNi/dy))
            B = [alpha1(2), 0, alpha2(2), 0, alpha3(2), 0;
                 0, alpha1(3), 0, alpha2(3), 0, alpha3(3);
                 alpha1(3) , alpha1(2), alpha2(3) , alpha2(2), alpha3(3) , alpha3(2)];

             sigma(i,:) = C*B*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2)]';

        end
        
    case 'quadratic'
        
        nEdge = varargin{1};
        nCorner = length(nodes) - nEdge;
        % berechnung der Spannungen an den Integrationspunkten
        % (Midside-Nodes)
        sigma = zeros(Nn, 3);
        sigma_accumulator = zeros(Ne, 3, 6);
        for i = 1:Ne
        
            ind = elts(i,:);
            P = nodes(ind, :);
            
            
            B1 = BpostQuad(P, 0, 0);
            B2 = BpostQuad(P, 1, 0);
            B3 = BpostQuad(P, 0, 1);
            

            B4 = BpostQuad(P, 0.5, 0);
            B5 = BpostQuad(P, 0.5, 0.5);
            B6 = BpostQuad(P, 0, 0.5);
            
            sigma_accumulator(i,:, 1) = C*B1(:, 1:12)*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2), ...
                                                        u(ind(4)*2-1), u(ind(4)*2), u(ind(5)*2-1), u(ind(5)*2), u(ind(6)*2-1), u(ind(6)*2)]';
            sigma_accumulator(i,:, 2) = C*B2(:, 1:12)*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2), ...
                                                        u(ind(4)*2-1), u(ind(4)*2), u(ind(5)*2-1), u(ind(5)*2), u(ind(6)*2-1), u(ind(6)*2)]';
            sigma_accumulator(i,:, 3) = C*B3(:, 1:12)*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2), ...
                                                        u(ind(4)*2-1), u(ind(4)*2), u(ind(5)*2-1), u(ind(5)*2), u(ind(6)*2-1), u(ind(6)*2)]';
            
            sigma_accumulator(i,:, 4) = C*B4(:, 1:12)*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2), ...
                                                        u(ind(4)*2-1), u(ind(4)*2), u(ind(5)*2-1), u(ind(5)*2), u(ind(6)*2-1), u(ind(6)*2)]';
            sigma_accumulator(i,:, 5) = C*B5(:, 1:12)*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2), ...
                                                        u(ind(4)*2-1), u(ind(4)*2), u(ind(5)*2-1), u(ind(5)*2), u(ind(6)*2-1), u(ind(6)*2)]';
            sigma_accumulator(i,:, 6) = C*B6(:, 1:12)*[u(ind(1)*2-1), u(ind(1)*2), u(ind(2)*2-1), u(ind(2)*2), u(ind(3)*2-1), u(ind(3)*2), ...
                                                        u(ind(4)*2-1), u(ind(4)*2), u(ind(5)*2-1), u(ind(5)*2), u(ind(6)*2-1), u(ind(6)*2)]';
            
             if any(sigma(ind(1), :), 'all') == 0
                sigma(ind(1), :) = sigma_accumulator(i,:, 1);
            else
                sigma(ind(1), :) = sigma(ind(1), :)/2 + sigma_accumulator(i,:, 1)/2;
            end
            if any(sigma(ind(2), :), 'all') == 0
                sigma(ind(2), :) = sigma_accumulator(i,:, 2);
            else
                sigma(ind(2), :) = sigma(ind(2), :)/2 + sigma_accumulator(i,:, 2)/2;
            end
            if any(sigma(ind(3), :), 'all') == 0
                sigma(ind(3), :) = sigma_accumulator(i,:, 3);
            else
                sigma(ind(3), :) = sigma(ind(3), :)/2 + sigma_accumulator(i,:, 3)/2;
            end
            if any(sigma(ind(4), :), 'all') == 0
                sigma(ind(4), :) = sigma_accumulator(i,:, 4);
            else
                sigma(ind(4), :) = sigma(ind(4), :)/2 + sigma_accumulator(i,:, 4)/2;
            end
            if any(sigma(ind(5), :), 'all') == 0
                sigma(ind(5), :) = sigma_accumulator(i,:, 5);
            else
                sigma(ind(5), :) = sigma(ind(5), :)/2 + sigma_accumulator(i,:, 5)/2;
            end
            if any(sigma(ind(6), :), 'all') == 0
                sigma(ind(6), :) = sigma_accumulator(i,:, 6);
            else
                sigma(ind(6), :) = sigma(ind(6), :)/2 + sigma_accumulator(i,:, 6)/2;
            end
            
        end
end