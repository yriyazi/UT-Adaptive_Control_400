clc;
clear all;
close all;
%%
run("Q310_Basic.m")
%% generate data

tfinal=100;
t = 0:T_s_close:tfinal;

Priemss=primes(100);
u=zeros(numel(t),1);
for i=20:numel(Priemss)
    input_dummy=gensig('sine' , tfinal/Priemss(1,i) , tfinal ,T_s_close);
    u=u+input_dummy;
end

var_e=0.7;
e=sqrt(var_e);
Noise=(-e+(e+e)*rand(numel(t),1));
u=u+Noise;
y = lsim(sys_dis_close ,u ,t);
plot(t,u ,t , y ,'LineWidth',2) ;
xlabel('Time (sec)') ;
ylabel('Magnitute') ;
grid on
legend('Input','Output') ;
%% Kalman Filter

N = numel(y) ;
% %--------------------------%
%choose number of parameters
Parameters_in_den=3
Parameters_in_num=3
% %--------------------------%
Nv=Parameters_in_num+Parameters_in_den;
%%
p_KF(1:Nv,1:Nv,1:N)=zeros(Nv,Nv,N);
%intitial Conditions
theta_hat_KF(1:Nv,1:N)=ones(Nv,N);
K_KF(1:Nv,1:N)=zeros(Nv,N);
var_e=0.05;

p_KF(1:Nv,1:Nv,1)=1e5*eye(Nv);p_KF(1:Nv,1:Nv,2)=p_KF(1:Nv,1:Nv,1);p_KF(1:Nv,1:Nv,3)=p_KF(1:Nv,1:Nv,2);


for i=(max(Parameters_in_num,Parameters_in_den)+1):N
    phi_KF(:,i)         =[(y(i-1:-1:i-Parameters_in_den))',(u(i-1:-1:i-Parameters_in_num))']';
    K_KF(:,i)           =p_KF(:,:,i-1)*phi_KF(:,i)*(1+phi_KF(:,i)'*p_KF(:,:,i-1)*phi_KF(:,i))^(-1) ;  
    p_KF(:,:,i)         =p_KF(:,:,i-1)-p_KF(:,:,i-1)*phi_KF(:,i)*(1+phi_KF(:,i)'*p_KF(:,:,i-1)*phi_KF(:,i))^(-1)*phi_KF(:,i)'*p_KF(:,:,i-1)+var_e;
    theta_hat_KF(:,i)   =theta_hat_KF(:,i-1)+K_KF(:,i)*(y(i)-phi_KF(:,i)'*theta_hat_KF(:,i-1));
end

%% 
% Bode

ident_dis = tf(theta_hat_KF((Parameters_in_num+1):end,end)' ,[1 -theta_hat_KF(1:Parameters_in_num ,end)'], T_s_close)
ident_analog = d2c(ident_dis)
bode(ident_analog ,'g*',sys_cont_close )
legend('model ','system')
%% 
% KF Convergence

subplot(2,1,1)
plot(t , theta_hat_KF((Parameters_in_num+1):end,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('KF convergence Num') ;
grid on
legend('a_1','a_2','a_3')
xlim([0 6])
% ylim([-1 1])
%--------------------------------------------------------------
subplot(2,1,2)
plot(t , -theta_hat_KF(1:Parameters_in_num ,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('KF convergence Den') ;
grid on
legend('b_1','b_2','b_3')
xlim([0 6])
% ylim([-7  7])