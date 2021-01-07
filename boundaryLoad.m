function [fx, fy] = boundaryLoad(points, edges, load, norm)

% calculates boundray Load vectors fx, fy for given line Load load(x)
% norm specifies the direction normal vector of each point of the loaded
% boundary line.


nedges = length(edges);

for i = 1:nedges
    ind = edges(i,:);
    P = points(ind, :);
    x1 = P(1,1);
    x2 = P(2,1);
    y1 = P(1,2);
    y2 = P(2,2);
    N = @(x) [(x-x2)/(x2-x1), (x1-x)/(x2-x1)];
    
end