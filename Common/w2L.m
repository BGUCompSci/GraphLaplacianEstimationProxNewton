function L = w2L(w) 
m = length(w);
n = 0.5 + ((1+8*m)^0.5)/2;
ind = (1:n).' >= (1:n)+1;
L = zeros(n);
L(ind) = -w;
L = L + L.';
ind = (1:n).' == (1:n);
L(ind) = -sum(L);
end