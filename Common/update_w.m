function [w,fnew,Lnew] = update_w(w,d,f,fval, maxIter_line, eta, beta,tr)
for p = 1:maxIter_line
    w_new = wthresh(w + eta*d,'h',tr);
    [fnew,Lnew] = f(w_new);
    if fnew <= fval
    	w = w_new;
    	break
    end
    eta = eta*beta;
    if p == maxIter_line
        break
    end
end
end