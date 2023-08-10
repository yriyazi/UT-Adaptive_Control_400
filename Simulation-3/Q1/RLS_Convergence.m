function  RLS_Convergence(number_X,X,x,number_G,G,g,t,Titlework,Plot_true,plot_counter)
figure()
subplot(2,1,1)
for i=1:number_X
    legend_names{i} = [ x '' num2str(i)-1 ''];
end
plot(t,X(:,:), 'LineWidth' , 2) ;
legend(legend_names)
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title( [x ' convergence']) ;
% xlim('auto')
% ylim([-0.5 1.5])
grid on


subplot(2,1,2)
for i=1:number_G
    legend_names{i} = [ g '' num2str(i)-1 ''];
end
plot(t ,G(:,:), 'LineWidth' , 2) ;
legend(legend_names)
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title( [g ' convergence Den']) ;
% xlim('auto')
% ylim([-1 2])
grid on
    if Plot_true==1
        print(gcf,[Titlework , num2str(plot_counter) ' RLS Convegence.png'],'-dpng','-r400');
    end
end