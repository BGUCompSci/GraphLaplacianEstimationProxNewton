function [] = ploting(error,obj,fs,nnzs,time,method_title)
    k = length(error)-1;
    figure;
    semilogy(1:k+1,error);
    xlabel('Iterations');
    ylabel('Error');
    title(['Error -  Frobenius norm - ',method_title] );
    legend('Initialization Based on Data')
    grid on
    
    figure;
    m = (1+ 10^(-3))*abs(min(obj));
    plot(1:k+1,obj);
    xlabel('Iterations');
    ylabel('f(L)');
    title(['Objective Function - ',method_title])
    %legend('Initialization Based on Data')
    grid on

    figure;
    plot(1:k+1,fs);
    xlabel('Iterations');
    ylabel('FS(L)');
    title(['F-Score - ',method_title])
    grid on
    
    %figure;
    %plot(1:k+1,nnzs);
    %xlabel('Iterations');
    %ylabel('nnz');
    %title('nnz ' );
    %legend('Initialization Based on Data')
    %grid on
    
    %figure;
    %semilogy(time,error);
    %xlabel('time [sec]');
    %ylabel('Error Norm');
    %title('Error vs time' );
    %legend('Initialization Based on Data')
    %grid on

end