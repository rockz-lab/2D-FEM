function [ Kred, P ] = ReduceSparseMatrix( K, idxU )

N = length(K);

Kred = K;

% Zeilen/Spalten löschen
Kred(idxU, :) = [];
Kred(:, idxU) = [];

% Permutationsmatrix

ii = 1:N;
jj = 1:N;
P = sparse(ii, jj, 1);
 
P(idxU, :) = [];