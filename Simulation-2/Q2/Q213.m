
run('BASIC.m')
%% generate Data

Question_mark='Q213';
[uc,t,Status,tfinal]=Datagen(1,T_s,1200);
Titlework=[Question_mark,Status]
% General Input+white Noise

y = lsim(sys_discret  ,uc ,t);
plot(t,uc ,'LineWidth',2) ;
xlabel('Time (sec)') ;
ylabel('Value') ;
title('ressponse of Question 2') ;
grid on
legend('Input' , 'ّOutPut') ;
print(gcf,[Titlework , num2str(plot_counter) ' Refrence .png'],'-dpng','-r400');
plot_counter=plot_counter+1;
%% Assumption

%choose number of parameters
Parameters_in_den=3
Parameters_in_num=3
degre_canselled_zero=1
A_m=poly(linspace(0.01,0.7,(Parameters_in_num)))
B_m=sum(A_m)*[1 zeros(1,(numel(A_m)-1-1))];


deg_A=Parameters_in_num+1;
deg_B=Parameters_in_den  ;


DC_gain_model=sum(A_m)/sum(B_m);              betaa=DC_gain_model
Nv=Parameters_in_num+Parameters_in_den;
deg_B_plus =degre_canselled_zero+1;
Deg_B_minus=deg_B-deg_B_plus;
Deg_A_o=0;
%A_o=[1 zeros(1,Deg_A_o)];
A_o=[1];
A_c_prim=conv(A_m,A_o);
u_cont=uc;
%%     RLS

theta(1:Nv,1:40) = ones(Nv,40) ;
P = 1e16*eye(Nv) ;
phi=[];
N = numel(y) ;
S=zeros(Nv,Parameters_in_num);
R=zeros(Nv,Parameters_in_den);
T=zeros(Nv,Parameters_in_num);
for i = (deg_A+deg_B):N
    if i<50
        y(i) =  -A(2:end)*y(i-1:-1:i-3)+B*(u_cont(i-1:-1:i-3)) ;
        phi(:,i) = [(y(i-1:-1:i-Parameters_in_den))' , (u_cont(i-1:-1:i-Parameters_in_num))']';
        K = P*phi(:,i)*(1+phi(:,i)'*P*phi(:,i))^(-1) ;
        P = (eye(Nv) - K*phi(:,i)')*P ;
        theta(:,i) = theta(:,i-1) + K*(y(i) - phi(:,i)'*theta(:,i-1));

        Aest=[1 -theta(1:Parameters_in_num ,end)'];
        Best=theta((Parameters_in_num+1):end,end)';
        B_minus=Best(1);
        B_plus=B/B_minus;
    else
        y(i) =  -Aest(2:end)*y(i-1:-1:i-Parameters_in_num)+Best*(u_cont(i-1:-1:i-Parameters_in_den)) ;

        phi(:,i) = [(y(i-1:-1:i-Parameters_in_den))' , (u_cont(i-1:-1:i-Parameters_in_num))']';
        K = P*phi(:,i)*(1+phi(:,i)'*P*phi(:,i))^(-1) ;
        P = (eye(Nv) - K*phi(:,i)')*P ;
        theta(:,i) = theta(:,i-1) + K*(y(i) - phi(:,i)'*theta(:,i-1));

        Aest=[1 -theta(1:Parameters_in_num ,end)'];
        Best=theta((Parameters_in_num+1):end,end)';
        B_minus=Best(1);
        B_plus=Best/B_minus;



        [R_prim,Si] = Diophantine(Aest,Best ,A_c_prim);
        S(i,:)=Si;
        T(i,:)=conv(sum(A_m),[1 zeros(1,(numel(A_m)-1-1))])/Best(1);
;
        R(i,:)=B_plus;

        A_c=conv(A_c_prim,B_plus);

        var1=conv(B,T(i,:));      narv1=numel(var1)      ;
        var2=A_c(2:end)     ;     narv2=numel(A_c(2:end));
        var3=conv(A,T(i,:)) ;     narv3=numel(var3)      ;

        y(i)     =var1*uc(i:-1:i-narv1+1)-var2*y(i-1:-1:i-narv2);
        u_cont(i)=var3*uc(i:-1:i-narv3+1)-var2*u_cont(i-1:-1:i-narv2);
    end
end
%% General Input v.s. Output

plot(t,uc,'b',t,y,'r--','LineWidth',1)
title('Refrence vs Controled output')
legend('input','plant out put')
print(gcf,[Titlework , num2str(plot_counter) ' Refrence vs Controled output.png'],'-dpng','-r400');
plot_counter=plot_counter+1;

plot(t,u_cont)
title('Control Signal')
print(gcf,[Titlework , num2str(plot_counter) ' Control Signal.png'],'-dpng','-r400');
plot_counter=plot_counter+1;

plot(t,gensig('square' , tfinal/10 , tfinal ,T_s),'b',t,y,'r',t,u_cont,'g--','LineWidth',0.25)
title('All Combined')
legend('input','plant out put','Controller Signal')
xlim([0 tfinal])
ylim([0 2])
print(gcf,[Titlework , num2str(plot_counter) ' All Combined.png'],'-dpng','-r600');
plot_counter=plot_counter+1;
%% 
% RLS Convergence of R UND S

subplot(2,1,1)
    for i=1:deg_A-1
        legend_names{i} = ['s' num2str(i)-1 ''];
    end
    plot(t ,S(:,:), 'LineWidth' , 2) ;
    legend(legend_names)
    xlabel('Time (sec)') ;
    ylabel('Parameters') ;
    title('S convergence') ;
grid on
subplot(2,1,2)
    for i=1:deg_A-1
        legend_names{i} = ['r' num2str(i)-1 ''];
    end
    plot(t ,R(:,:), 'LineWidth' , 2) ;
    legend(legend_names)
    xlabel('Time (sec)') ;
    ylabel('Parameters') ;
    title('R convergence Den') ;
    grid on
print(gcf,[Titlework , num2str(plot_counter) ' RLS Convegence.png'],'-dpng','-r400');