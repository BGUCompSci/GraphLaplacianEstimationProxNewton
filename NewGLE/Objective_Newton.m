function [h,DQ] = Objective_Newton(delta,w,S,Q,epsilon_vec,lambda)
    ed = epsilon_vec.*delta;
    Delta = w2L(delta);
    DQ = Delta*Q;
    trDQDQ = sum(sum(DQ.*DQ)); % trace((DQ)'*DQ)
    trDeltaS = sum(sum(Delta.*S)); % trace( Delta*S) 
    t = trDeltaS - trace(DQ) + norm(ed)^2 + 0.5*trDQDQ; 
    r = MCP(w + delta,lambda);
    h = t + r;
end
