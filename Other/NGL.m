function [L,obj,error,fs,nnzs,T,method] = NGL(Lt,S,lambda) % with Palomar's porjection
tic
wt = L2w(Lt);
n = size(S,1);
err = 10^(-4);
maxIter = 50;
k = 0;
%init:
w = w_init0(S,lambda);
z = gMCP(w,lambda);
f = @(w) objective_palomar(w,S,z);
L = w2L(w);
obj(k+1) = f(w);
nnzs(k+1) = 2*nnz(w) + n;
normLt = norm(Lt,'fro');
error(k+1) = norm(Lt-L,'fro')/normLt;
fs(k+1) = F_score(w,wt);
T(1) = 0;
%============================= main loop ==================================
while k == 0 || k < maxIter  % add stopping condition
    wold = w;
    
    z = gMCP(w,lambda);
    f = @(w) objective_palomar(w,S,z);
    [w,L] = update_w(w,S,z,f,L);  
    
    k = k + 1;
    nnzs(k+1) = 2*nnz(w) + n;
    %obj(k+1) = f(w);
    obj(k+1) = objective_f(w,S,lambda);
    error(k+1) = norm(Lt-L,'fro')/normLt;
    fs(k+1) = F_score(w,wt);
    T(k+1) = toc;
    stop_cond = norm(w-wold)/norm(wold) < err && nnzs(k+1) == nnzs(k);
     if stop_cond
        break
     end
end
method = 'NGL';
L = w2L(w);
end
%=============================== Functions ================================
function [f,L] = objective_palomar(w,S,z)
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
    r = dot(z,w);
    f = t - d + r;
end
%--------------------------------------------------------------------------
% function g = grad_by_w_of_f(L,S,z)
% e_D = 10^(-10);
% n = size(S,1);
% J = 1/n;
% g = Lstar(S - inv(L + J + e_D*eye(n))) + z;
% end

function g = grad_by_w_of_f(L,S,z,p)
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
	g = (((Lstar(S-M))+z)./(-L2w(P)));
else
    g = Lstar(S-M);
    g = g + z;
end
end

%--------------------------------------------------------------------------
function [w,L] = update_w(w,S,z,f,L)
err = 10^(-4);
tr = 10^(-6);
maxIter = 200;
maxIter_line = 10;
eta0 = 1;
beta = 0.5;
g_control = 1;
precond = 0;

g = @(L) grad_by_w_of_f(L,S,z,precond);

for t = 1:maxIter
    eta = eta0;
    wold = w;
    grad = grad_control(g(L),g_control);
    f_old = f(w);
    for p = 1:maxIter_line 
        d = -grad;         
        w_new = wthresh(max(w + eta*d,0),'h',tr);
        [f_new,Lnew] = f(w_new);
        if f_new <= f_old %+ dot(grad,w-w_new) + (norm(w-w_new).^2)/(2*eta)
        	w = w_new;
            L = Lnew;
        	break
        end
        eta = beta*eta;
    end
    stop_cond = p == maxIter_line || norm(w_new-wold)/norm(wold) < err;
    if stop_cond
        w = w_new;
        L = Lnew;
        break
    end
end
%t
end

