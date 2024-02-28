function r = MCP(w,lambda)
% where w >= 0 
    gama = 1.01; 
    gama_lambda = gama*lambda;
    helf_gama2_lambda = gama_lambda*lambda/2;
    S = w <= gama_lambda;
    MCPw = (lambda*w-(w.^2)./(2*gama) - helf_gama2_lambda).*S + helf_gama2_lambda;
    %MCPw = (lambda*w-(w.^2)./(2*gama)).*S + (0.5*(lambda^2)*gama)*(1-S);
    r = sum(MCPw);
end