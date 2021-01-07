function Kges = assemble(Nodes, Elts, C, t, lambda, tau, Bfun)


Ne = length(Elts);
Nn = length(Nodes);

ii = [];
jj = [];
KK = [];

switch Bfun
    
    case 'linear'
        
        for i = 1:Ne

            ind = Elts(i, :);
            P = Nodes(ind,:);

            K = integrateElement(P, C, t, lambda, tau, Bfun);
            ii = [ii, repmat([ind(1)*2-1, ind(1)*2, ind(2)*2-1, ind(2)*2, ind(3)*2-1, ind(3)*2], 1, 6)];
            jj = [jj, ones(1,6)*(ind(1)*2-1), ones(1,6)*(ind(1)*2), ones(1,6)*(ind(2)*2-1)...
                , ones(1,6)*(ind(2)*2), ones(1,6)*(ind(3)*2-1), ones(1,6)*(ind(3)*2)];
        %     KK = [KK, K(:,1)', K(:,2)', K(:,3)', K(:,4)', K(:,5)', K(:,6)'];
            KK = [KK, K(1,:), K(2,:), K(3,:), K(4,:), K(5,:), K(6,:)];

        end
    case 'quadratic'
        
        for i = 1:Ne

            ind = Elts(i, :);
            P = Nodes(ind,:);

            K = integrateElement(P, C, t, lambda, tau, Bfun);
            ii = [ii, repmat([ind(1)*2-1, ind(1)*2, ind(2)*2-1, ind(2)*2, ind(3)*2-1, ind(3)*2, ind(4)*2-1, ind(4)*2, ind(5)*2-1, ind(5)*2, ind(6)*2-1, ind(6)*2], 1, 12)];
            jj = [jj, ones(1,12)*(ind(1)*2-1), ones(1,12)*(ind(1)*2), ones(1,12)*(ind(2)*2-1)...
                , ones(1,12)*(ind(2)*2), ones(1,12)*(ind(3)*2-1), ones(1,12)*(ind(3)*2), ones(1,12)*(ind(4)*2-1), ones(1,12)*(ind(4)*2)...
                , ones(1,12)*(ind(5)*2-1), ones(1,12)*(ind(5)*2), ones(1,12)*(ind(6)*2-1), ones(1,12)*(ind(6)*2)];
        %     KK = [KK, K(:,1)', K(:,2)', K(:,3)', K(:,4)', K(:,5)', K(:,6)'];
            KK = [KK, K(:)'];

        end
        
end

Kges = sparse(ii, jj, KK);