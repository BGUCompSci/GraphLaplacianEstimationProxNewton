function g = gradObjNewtonByDelta(delta,w,S,Q,epsilon_vec,lambda,DQ)
z = gMCP(delta + w,lambda);
g = Lstar(S - Q + Q'*DQ) + 2*epsilon_vec.*epsilon_vec.*delta + z;
end