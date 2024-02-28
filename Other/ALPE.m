function [L,obj,error,fs,nnzs,T,method] = ALPE(Lt,S,lambda) % with Palomar's porjection
tic
wt = L2w(Lt);
T(1) = 0 ;
n = size(S,1);
eps = 10^(-1);
err = 10^(-4);
tr = 10^(-6);
maxIter = 10000;
maxIter_line = 20;
eta0 = 1;
beta = 0.25;
g_control = 1;
precond = 0;
q = 2;

w = w_init0(S,lambda);

k = 0;
normLt = norm(Lt,'fro');
L = w2L(w);
error(k+1) = norm(Lt-L,'fro')/normLt;
fs(k+1) = F_score(w,wt);
nnzs(k+1) =  2*nnz(w) + n;
%T(k+1) = toc;
k= k+1;

warmStart = 0;
wML = NGL_GD_warmStartALPE(S,0,w);

a = 1./((abs(wML) + eps).^q);%.*supp;
m = max(a);
a = a/m;
f = @(w) objective_ALPE(w,S,a,lambda);

error(k+1) = norm(Lt-L,'fro')/normLt;
fs(k+1) = F_score(w,wt);
if warmStart == 1  
    w = wML;
end
[obj(k+1),L] = f(w);
nnzs(k+1) =  2*nnz(w) + n;
T(k+1) = toc;
%============================= main loop ==================================

g = @(L) grad_by_w_of_f(L,S,a,lambda,precond);

while k == 0 || k < maxIter  % add stopping condition
    eta = eta0;
    wold = w;
    grad = grad_control(g(L),g_control);%.*supp;
    f_old = f(w);
    for p = 1:maxIter_line 
        d = -grad;         
        w_new = wthresh(max(w + eta*d,0),'h',tr);
        [f_new,Lnew] = f(w_new);
        if f_new <= f_old% + dot(grad,w-w_new) + (norm(w-w_new).^2)/(2*eta)
        	w = w_new;
            L = Lnew;
        	break
        end
        eta = beta*eta;
        if p == maxIter_line
            w = w_new;
            L = Lnew;
        end
    end
    k = k+1;
    obj(k+1) = f_new;
    error(k+1) = norm(Lt-L,'fro')/normLt;
    fs(k+1) = F_score(w,wt);
    nnzs(k+1) =  2*nnz(w) + n;
    T(k+1) = toc;
    stop_cond = norm(w-wold)/norm(wold) < err && nnzs(k+1) == nnzs(k);
     if stop_cond
        break
     end
end
method = 'ALPE';
end
%=============================== Functions ================================
function [f,L] = objective_ALPE(w,S,a,lambda)
    n = size(S,1);
    J = 1/n; %  ones(n)./n;
    L = w2L(w);
    
    t = trace( L*S ); % trace term

    [v,flag] = chol(L + J);
    if flag
       e = norm(L,1)*1e-14; 
       [v,~] = chol(L + J + e*eye(n));
    end
    u = diag(v);
    d = 2*sum(log(u)); %logdet term
    r = lambda*dot(w,a);
    f = t - d + r;
end
%--------------------------------------------------------------------------
% function g = grad_by_w_of_f(L,S,a,lambda)
% e_D = 10^(-10);
% n = size(S,1);
% J = 1/n;
% g = Lstar(S - inv(L + J + e_D*eye(n))) + lambda*a;
% end

function g = grad_by_w_of_f(L,S,a,lambda,p)
n = size(S,1);
e_D = 10^(-10);
J = 1/n; %ones(n)./n;
M = L + J;
ind = (1:n).' == (1:n);
M(ind) = M(ind) + e_D;
M = inv(M);
if p
	dM = diag(M);
	P = (dM + dM' - 2*M).^2;
	g = (((Lstar(S-M))+lambda*a)./(-L2w(P)));
else
    g = Lstar(S-M);
    g = g + lambda*a;
end
end


