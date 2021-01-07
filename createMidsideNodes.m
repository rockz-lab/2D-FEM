function [nodesQuad, elementsQuad] = createMidsideNodes(nodes, elements)

Ne = length(elements);
Nn = length(nodes);


elementsQuad = zeros(Ne, 6);



[temp, ia, ib] = unique([nodes(elements(:, 2), :)/2 + nodes(elements(:, 1), :)/2; ...
                                 nodes(elements(:, 3), :)/2 + nodes(elements(:, 2), :)/2; ...
                                 nodes(elements(:, 3), :)/2 + nodes(elements(:, 1), :)/2], 'sorted', 'rows');
n = length(temp);

nodesQuad = zeros(Nn+n, 2);
nodesQuad(1:Nn, :) = nodes;
nodesQuad(Nn+1:end, :) = temp;

elementsQuad(:, 1:3) = elements;
%ib = reshape(ib, [3, Ne])';
ib = reshape(ib, [Ne, 3]);
elementsQuad(:, 4:6) = ib+Nn;
