clc;
clear all;
close all;
%% generate data

tic
%%
run('Q310_Basic.m')
%%
tfinal=200;
t = 0:T_s_close:tfinal;
u = gensig('sine' , tfinal/20 , tfinal ,T_s_close);
Noise=(-0.2+(0.2+0.2)*rand(numel(t),1))/2;
u=u+Noise;
y = lsim(sys_dis_close ,u ,t);
plot(t,u ,t , y ,'LineWidth',2) ;
xlabel('Time (sec)') ;
ylabel('Magnitute') ;
grid on
legend('Input','Output') ;
%%
N = numel(y) ;
Parameters_in_den=4
Parameters_in_num=4
Nv=Parameters_in_num+Parameters_in_den;
theta(:,1:Nv) = zeros(Nv , Nv) ;
P = 1e12*eye(Nv) ;
phi=[];
eror(1:Nv,1:N)=zeros(Nv,N);
for i = (max(Parameters_in_num,Parameters_in_den)+1):N
    phi(:,i) = [[y(i-1:-1:i-Parameters_in_den)]' , [u(i-1:-1:i-Parameters_in_num)]']';
    K = P*phi(:,i)*(1+phi(:,i)'*P*phi(:,i))^(-1) ;
    P = (eye(Nv) - K*phi(:,i)')*P ;
    theta(:,i) = theta(:,i-1) + K*(y(i) - phi(:,i)'*theta(:,i-1));
    eror(:,i)=theta(:,i)-[-d1(2:end),c1(2:end)]';
end
%% 
% Bode

ident_dis = tf(theta((Parameters_in_num+1):end,end)' ,[1 -theta(1:Parameters_in_num ,end)'], T_s_close)
ident_analog = d2c(ident_dis)
bode(ident_analog ,'g*',sys_cont_close)
legend('model ','system')
%% 
% RLS Convergence

subplot(2,1,1)
plot(t , theta((Parameters_in_num+1):end,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('RLS convergence Num') ;
grid on
legend('a_1','a_2','a_3')
xlim([0 6])
ylim([-1 1])
%--------------------------------------------------------------
subplot(2,1,2)
plot(t , -theta(1:Parameters_in_num ,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('RLS convergence Den') ;
grid on
legend('b_1','b_2','b_3')
xlim([0 6])
ylim([-7  7])
%%
toc