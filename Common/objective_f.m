function [f,L] = objective_f(w,S,lambda)
    n = size(S,1);
    %J = ones(n)./n;
    J = 1.0/n;
    L = w2L(w);
    
    t = sum(sum(L.*S)); % trace( L*S ); % trace term

    [v,flag] = chol(L + J);
    if flag
       e = norm(L,1)*1e-10; 
       [v,~] = chol(L + J + e*eye(n));
    end
    u = diag(v);
    d = 2*sum(log(u)); %logdet term
    r = MCP(w,lambda);
    f = t - d + r;
end