function w = NGL_GD_warmStart(S,lambda,w0) % with Palomar's porjection
w = w0;
maxIter = 5; %100
maxIter_line = 10;
eta = .7;
beta = 0.3;
tr = 10^(-5);
g_control = 1;
precond = 1; %preconditioning  
k = 0;

f = @(w) objective_f(w,S,lambda);
g = @(w,L) grad_by_w_of_f(w,S,lambda,L,precond);

[fold,L] = f(w);
obj(k+1) = fold;

%============================= main loop ==================================
while k == 0 || k < maxIter  % add stopping condition
    d = -grad_control(g(w,L),g_control);
    [w,fnew,L] = update_w_proj(w,d,f,fold,maxIter_line,eta , beta,tr);
    k = k + 1;   
    obj(k+1) = fnew;
    fold = fnew;
end
end