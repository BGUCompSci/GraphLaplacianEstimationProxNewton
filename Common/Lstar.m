function w = Lstar(Y)
w = L2w(  Y + Y.' - diag(Y)  - diag(Y)');
end