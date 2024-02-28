function w = w_init_data(S,tr)
 %eD = 0.001;
 %Lhat = inv(S + eye(size(S,1))*eD);
%Lhat = pinv(S);
%tr = lambda;
Lhat = L0(size(S,1));

w = abs(wthresh(L2w(Lhat),'h',tr));

end