function w = w_init0(S,lambda)
tr = lambda;
Lhat = L0(size(S,1));
w = abs(wthresh(L2w(Lhat),'h',tr));
end