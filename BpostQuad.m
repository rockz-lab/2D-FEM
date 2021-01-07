function B = BpostQuad(P, xi, eta)



x1 = P(1,1);
x2 = P(2,1);
x3 = P(3,1);
x4 = P(4,1);
x5 = P(5,1);
x6 = P(6,1);

y1 = P(1,2);
y2 = P(2,2);
y3 = P(3,2);
y4 = P(4,2);
y5 = P(5,2);
y6 = P(6,2);

a1 = 4*x4 - x2 - 3*x1;
a2 = 2*x1 + 2*x2 - 4*x4;
a3 = 4*x6 - x3 - 3*x1;
a4 = 2*x1 + 2*x3 - 4*x6;
a5 = 4*x1 - 4*x4 + 4*x5 - 4*x6;

b1 = 4*y4 - y2 - 3*y1;
b2 = 2*y1 + 2*y2 - 4*y4;
b3 = 4*y6 - y3 - 3*y1;
b4 = 2*y1 + 2*y3 - 4*y6;
b5 = 4*y1 - 4*y4 + 4*y5 - 4*y6;

J = [a1 + 2*a2*xi + eta*a5, a3 + 2*a4*eta + a5*xi;
b1 + 2*b2*xi + b5*eta, b3 + 2*b4*eta + b5*xi];

detJ = abs(det(J));

dN_d_xi = [-3+4*xi+4*eta, 4*xi-1, 0, 4-8*xi-4*eta, 4*eta, -4*eta];
dN_d_eta = [-3+4*eta+4*xi, 0, 4*eta-1, -4*xi, 4*xi, 4-8*eta-4*xi];

% make adjoint
adj_J = [J(2,2), - J(1,2);
        -J(2,1), J(1,1)];
    
dNdxdy = (1/detJ*adj_J)'*[dN_d_xi; dN_d_eta];

B = [dNdxdy(1, 1), 0, dNdxdy(1, 2), 0, dNdxdy(1, 3), 0, dNdxdy(1, 4), 0, dNdxdy(1, 5), 0, dNdxdy(1, 6), 0;...
0, dNdxdy(2, 1), 0, dNdxdy(2, 2), 0, dNdxdy(2, 3), 0, dNdxdy(2, 4), 0, dNdxdy(2, 5), 0, dNdxdy(2, 6);...
dNdxdy(2, 1), dNdxdy(1, 1), dNdxdy(2, 2), dNdxdy(1, 2), dNdxdy(2, 3), dNdxdy(1, 3), dNdxdy(2, 4), dNdxdy(1, 4), dNdxdy(2, 5), dNdxdy(1, 5), dNdxdy(2, 6), dNdxdy(1, 6)];