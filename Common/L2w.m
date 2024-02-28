function w = L2w(L)
m = (1:size(L,1)).' >= (1:size(L,2))+1;
w = -L(m);
end