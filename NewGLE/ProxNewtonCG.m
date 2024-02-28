function delta = ProxNewtonCG(w,L,S,lambda,activeSet,CG,epsilon_param,precond_newton)

% CG with pre-cond - QN
M = 1;
if activeSet == 1
    g0 = grad_by_w_of_f(w,S,lambda,L,0);
    M = w > 0;
    M = M | ( g0 < -lambda);
end
gamma = 1.01;
epsilon_vec = epsilon_param*(M & (w < gamma.*lambda));
maxIter_Newton = 10;
maxIter_linesearch = 20;
n = size(S,1);
eta = 1;
beta = 0.5;
e_D = 10^(-6);
J = 1/n;
g_control = 0;
err = 0.01;


Q = inv(L + J + e_D*eye(n));

h = @(delta) Objective_Newton(delta,w,S,Q,epsilon_vec,lambda);
if precond_newton==1
    dM = diag(Q);
    P = 1./(-L2w((dM + dM' - 2.*Q).^2) + 2*epsilon_vec.^2);
    grad = @(delta,DQ) (gradObjNewtonByDelta(delta,w,S,Q,epsilon_vec,lambda,DQ).*P);
else
    grad = @(delta,DQ) gradObjNewtonByDelta(delta,w,S,Q,epsilon_vec,lambda,DQ);
end

delta = zeros(size(w)); 
[hold,DQ] = h(delta);
for t = 1:maxIter_Newton
	delta_old = delta;
    if t == 1 || CG==0
        g = grad_control(grad(w,DQ),g_control).*M;
        d = -g;
    else
        g_old = g;
        g = grad_control(grad(w,DQ),g_control).*M;
        d_old = d;
        norm_dd_old = norm(d_old)^2;
        y = g - g_old;
        p = 1 + max(0,-dot(d_old,y)/norm_dd_old);
        a = y + p*d_old;
        dot_d_old_a = dot(d_old,a);
        eta_k = 1 - dot(g,d_old)/dot_d_old_a;
        %beta_PDY = max(min((norm(g)^2)/dot_d_old_a,1),0);
        beta_PDY = (norm(g)^2)/dot_d_old_a;
        d = -eta_k*g + beta_PDY*d_old;
    end
    eta = min(1.0,2*eta);    
    for p = 1:maxIter_linesearch
        z = max(delta + eta*d,-w);
        [hnew,DQ] = h(z);
        lines_cond =  hnew < hold ;      
        if lines_cond
            delta = z;
            break
        end
        eta = eta*beta;
    end
    if p == maxIter_linesearch  || norm(delta - delta_old)/norm(delta) < err
        break
    end
    hold = hnew;
end
%t
end

