function [L,S] = generate_L_S(n,p,a,b)
% w ~ U(a,b)

A = generateRandomPlanarGraph(p);
L = full(A);
w = L2w(L);
shift = a*sign(w);
w = b*rand(size(w)).*w + shift;
L = w2L(w);
S = generate_DATA_matrix(L,n);
end