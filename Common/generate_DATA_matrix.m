function S = generate_DATA_matrix(L,k)
n = size(L,1);
X = mvnrnd(zeros(n,1),pinv(L),k);
X = X - mean(X')';
S = cov(X);
end

