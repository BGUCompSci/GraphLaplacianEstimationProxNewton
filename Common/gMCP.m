function z = gMCP(w,lambda)
%output z (size of w and w>=0)
gama = 1.01;
S =  w <= gama*lambda;
z = (lambda - w/gama).*S;
end