function g = grad_by_w_of_f(w,S,lambda,L,p)
n = size(S,1);
e_D = 10^(-10);
J = 1/n; 
M = L + J;
ind = (1:n).' == (1:n);
M(ind) = M(ind) + e_D;
M = inv(M);
if p
	dM = diag(M);
	P = (dM + dM' - 2*M).^2;
	z = gMCP(w,lambda);
	g = (((Lstar(S-M))+z)./(-L2w(P)));
else
    g = Lstar(S-M);
    g = g + gMCP(w,lambda);
end
end