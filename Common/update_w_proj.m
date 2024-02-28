function [w,fnew,Lnew] = update_w_proj(w,d,f,fval,maxIter_line,eta,beta,tr)
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
    disp('THIS SHOULD NOT BE POSSIBLE ???')
end
end