function [L,obj,error,fs,nnzs,T,method] = NGL_MCP_NewtonCG(Lt,S,lambda,activeSet,CG,epsilon_param,precond_newton) % with Palomar's porjection
if nargin < 4
    activeSet = 1;
end
if nargin < 5
    CG = 1;
end
if nargin < 6
    epsilon_param = 0.1;
end
if nargin < 7
    precond_newton = 1;
end


T(1) = 0;
wt = L2w(Lt);
n = size(S,1);         
tr = 10^(-6);              %thresholding
err = 10^(-4);  %-4        % Stopping criteria: stop if: ||w - wnew||/||w|| < err
maxIter = 100;             % number of newton iterations
warmStart = 1;
ploting_flag = 0;
%lineseaech main loop cont's:
maxIter_line = 7;   % number of iterations for linesearch
beta = 0.5;
k = 0;
normLt = norm(Lt,'fro');
f = @(w) objective_f(w,S,lambda);

tic
w0 = w_init0(S,lambda);
[obj(k+1),L]  = f(w0); %f(w0)
nnzs(k+1) = 2*nnz(w0) + n;
error(k+1) = norm(Lt-L,'fro')/normLt;
fs(k+1) = F_score(w0,wt);
k = k+1;

eta = 1;
if warmStart %w0
	w = NGL_GD_warmStart(S,lambda,w0);
else
    w = w0;
end
T(k+1) = toc;
[obj(k+1),L]  = f(w); %f(w0)
nnzs(k+1) = 2*nnz(w) + n;
error(k+1) = norm(Lt-L,'fro')/normLt;
fs(k+1) = F_score(w,wt);

%============================= main loop ==================================
while k == 0 || k < maxIter 
    d = ProxNewtonCG(w,L,S,lambda,activeSet,CG,epsilon_param,precond_newton); %find newton direction
    wold = w;
    [w,obj(k+2),L] = update_w(w,d,f,obj(k+1),maxIter_line,eta ,beta,tr);
    T(k+2) = toc;
    k = k + 1;    
    nnzs(k+1) = 2*nnz(w) + n;
    error(k+1) = norm(Lt-L,'fro')/normLt;
    fs(k+1) = F_score(w,wt);
    stop_cond = norm(w-wold)/norm(wold) < err && nnzs(k+1) == nnzs(k);
    if stop_cond
        break
    end
end
method = 'NewGLE';
if ploting_flag
    ploting(error,obj,fs,nnzs,T,method);
end
end
