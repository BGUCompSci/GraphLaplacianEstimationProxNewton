close all;
k=1;
Res = load('GraphLap_NewGLE_lam_0.25.mat');
obj = Res.Res{1}; error = Res.Res{2};fs = Res.Res{3};nnzs = Res.Res{4}; T = Res.Res{5}; method = Res.Res{6}; 
objBatch{k} = obj; error_batch{k} = error; methodBatch{k} = method; fs_batch{k} = fs; nnzs_batch{k} = nnzs; T_batch{k} = T;

k=k+1;
Res = load('GraphLap_PGD_lam_0.25.mat');
obj = Res.Res{1}; error = Res.Res{2};fs = Res.Res{3};nnzs = Res.Res{4}; T = Res.Res{5}; method = Res.Res{6}; 
objBatch{k} = obj; error_batch{k} = error; methodBatch{k} = method; fs_batch{k} = fs; nnzs_batch{k} = nnzs; T_batch{k} = T;

k=k+1;
Res = load('GraphLap_NGL_lam_0.15.mat');
obj = Res.Res{1}; error = Res.Res{2};fs = Res.Res{3};nnzs = Res.Res{4}; T = Res.Res{5}; method = Res.Res{6}; 
objBatch{k} = obj; error_batch{k} = error; methodBatch{k} = method; fs_batch{k} = fs; nnzs_batch{k} = nnzs; T_batch{k} = T;

k=k+1;
Res = load('GraphLap_ALPE_lam_0.3.mat');
obj = Res.Res{1}; error = Res.Res{2};fs = Res.Res{3};nnzs = Res.Res{4}; T = Res.Res{5}; method = Res.Res{6}; 
objBatch{k} = obj; error_batch{k} = error; methodBatch{k} = method; fs_batch{k} = fs; nnzs_batch{k} = nnzs; T_batch{k} = T;

plotingBatch(k,objBatch,error_batch,fs_batch,nnzs_batch,T_batch,methodBatch);


