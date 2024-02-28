function [] = plotingBatch(bs,b_obj,b_error,b_fs,b_nnzs,b_time,b_method_title)
    

    figure;
    for j = 1:bs
        error = b_error{j};
        k = length(error)-1;
        semilogy(1:k+1,error,'-','LineWidth',1.5);
        hold on
    end
    grid on
    xlabel('Iterations');
    ylabel('Error');
    title('Error -  Frobenius norm');
    legend(b_method_title)


    figure;
    for j = 1:bs
        error = b_error{j};
        time = b_time{j};
        semilogy(time,error,'-','LineWidth',1.5);
        hold on
    end
    grid on
    xlabel('Time [sec]');
    ylabel('Error');
    title('Error -  Frobenius norm');
    legend(b_method_title)
    


    figure;
    for j = 1:bs
        obj = b_obj{j};
        time = b_time{j};
        plot(time,obj,'-','LineWidth',1.5);
        hold on
    end
    grid on
    xlabel('Time [sec]');
    ylabel('Objective: F(L)');
    title('Objective Value');
    legend(b_method_title)


    figure;
    for j = 1:bs
        obj = b_obj{j};
        time = b_time{j};
        semilogy(time,abs((obj - obj(end))/(obj(1) - obj(end))),'-','LineWidth',1.5);
        hold on
    end
    grid on
    xlabel('Time [sec]');
    ylabel('Relative Decrease');
    title('Relative Decrease in Objective Value');
    legend(b_method_title)
        
    
    figure;
    for j = 1:bs
        fs = b_fs{j};
        time = b_time{j};
        plot(time,fs,'-','LineWidth',1.5);
        hold on
    end
    grid on
    xlabel('Time [sec]');
    ylabel('F-Score');
    title('F - Scores');
    legend(b_method_title)

end