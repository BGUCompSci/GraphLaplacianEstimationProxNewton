function [L,obj,error,fs,nnzs,T,method] = PGD(Lt,S,lambda) % with Palomar's porjection
tic
wt = L2w(Lt);
n = size(S,1);
w = w_init0(S,lambda);

err = 10^(-4);
maxIter = 10000;
maxIter_line = 20;
eta = 1.;
beta = 0.5;
tr = 10^(-6);
k = 0;

g_control = 1;
precond = 0; %preconditioning  

f = @(w,k) objective_f(w,S,lambda);
g = @(w,L,k) grad_by_w_of_f(w,S,lambda,L,precond);

[obj(k+1),L] = f(w,k);
nnzs(1) = 2*nnz(w) + n;
normLt = norm(Lt,'fro');
error(k+1) = norm(Lt-L,'fro')/normLt;
fs(k+1) = F_score(w,wt);
T(1) = 0;
%============================= main loop ==================================
while k == 0 || k < maxIter 
    wold = w;
    d = -grad_control(g(w,L,k),g_control);
    [w,obj(k+2),L] = update_w_proj(w,d,f,obj(k+1),maxIter_line,eta , beta,tr);
    k = k + 1;
    error(k+1) = norm(Lt-L,'fro')/normLt;
    fs(k+1) = F_score(w,wt);
    T(k+1) = toc;
    nnzs(k+1) = 2*nnz(w) + n;
    stop_cond = norm(w-wold)/norm(wold) < err && nnzs(k+1) == nnzs(k);
     if stop_cond
        break
     end
end
method = 'PGD';
end

