function d = grad_control(d,flag)
if flag
    max_grad = 0.1;
    nd = norm(d,inf);
    if nd > max_grad
        d = (d/nd).*max_grad;
    end
end
end