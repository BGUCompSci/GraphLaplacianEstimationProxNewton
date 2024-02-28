addpath('./Common/');
addpath('./Other/');
addpath('./NewGLE/');
close all
clear all
n = 500;
p = 1000;
a = 0.5;
b = 2.0;
[Ltrue,S] = generate_L_S(n,p,a,b);
% diary myDiaryFile
disp('NewGLE ~~~~~~~~~~~~~~~~~~~~~~~~~~')
for lambda = [0.25,0.2]
    [L,obj,error,fs,nnzs,T,method] = NGL_MCP_NewtonCG(Ltrue,S,lambda);
    % Res = {L,obj,error,fs,nnzs,T,method,Ltrue,S,lambda,n,p,a,b};
    Res = {obj,error,fs,nnzs,T,method,lambda,n,p,a,b};
    save(['GraphLap_NewGLE_lam_',num2str(lambda),'.mat'],'Res')
    disp([num2str(lambda),', ',num2str(error(end)),', ',num2str(fs(end)),', ',num2str(T(end)),',',num2str(length(obj))])
end
disp('PGD  ~~~~~~~~~~~~~~~~~~~~~~')
for lambda = [0.25,0.2]
    [L,obj,error,fs,nnzs,T,method] = PGD(Ltrue,S,lambda);
    disp([num2str(lambda),', ',num2str(error(end)),', ',num2str(fs(end)),', ',num2str(T(end)),',',num2str(length(obj))])
%     %Res = {L,obj,error,fs,nnzs,T,method,Ltrue,S,lambda,n,p,a,b};
    Res = {obj,error,fs,nnzs,T,method,lambda,n,p,a,b};
    save(['GraphLap_PGD_lam_',num2str(lambda),'.mat'],'Res')
end
disp('ALPE ~~~~~~~~~~~~~~~~~~~~~~')
for lambda = [0.3]
     [L,obj,error,fs,nnzs,T,method] = ALPE(Ltrue,S,lambda);
%     %Res = {L,obj,error,fs,nnzs,T,method,Ltrue,S,lambda,n,p,a,b};
    Res = {obj,error,fs,nnzs,T,method,lambda,n,p,a,b};
    save(['GraphLap_ALPE_lam_',num2str(lambda),'.mat'],'Res')
    disp([num2str(lambda),', ',num2str(error(end)),', ',num2str(fs(end)),', ',num2str(T(end)),',',num2str(length(obj))])
end
disp('NGL ~~~~~~~~~~~~~~~~~~~~~~')
for lambda = [0.15]    
     [L,obj,error,fs,nnzs,T,method] =  NGL(Ltrue,S,lambda);
%     %Res = {L,obj,error,fs,nnzs,T,method,Ltrue,S,lambda,n,p,a,b};
    Res = {obj,error,fs,nnzs,T,method,lambda,n,p,a,b};
    save(['GraphLap_NGL_lam_',num2str(lambda),'.mat'],'Res')
    disp([num2str(lambda),', ',num2str(error(end)),', ',num2str(fs(end)),', ',num2str(T(end)),',',num2str(length(obj))])
end
% diary off

