function w = NGL_GD_warmStartALPE(S,lambda,w0) % with Palomar's porjection
w = w0;
maxIter = 45; %45
err = 10^(-4);
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
    w_old = w;
    d = -grad_control(g(w,L),g_control);
    [w,fnew,L] = update_w(w,d,f,fold,maxIter_line,eta , beta,tr);
    k = k + 1;   
    obj(k+1) = fnew;
    fold = fnew;
    if norm(w-w_old)/norm(w_old) < err
        break
    end
end
end

function [w,fnew,Lnew] = update_w(w,d,f,fval,maxIter_line,eta,beta,tr)
for p = 1:maxIter_line
    w_new = wthresh(max(w + eta*d,0),'h',tr);
    [fnew,Lnew] = f(w_new);
    if fnew < fval% + dot(-d,w-w_new) + (norm(w-w_new).^2)/(2*eta)  
    	w = w_new;
        break
    end
    eta = eta*beta;
    if p == maxIter_line
        break
    end
end
if fnew > fval
    disp('???')
end
end

