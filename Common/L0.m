function L = L0(n)
L = 2*eye(n);
L(1,1) = 1;
L(2,1) = -1;
L(n,n) = 1;
L(n-1,n) = -1;
for t = 2:n-1
    L(t-1,t) = -1;
    L(t+1,t) = -1;
end
end

