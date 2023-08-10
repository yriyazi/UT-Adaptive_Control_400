clc;
clear all;
close all;
%%
run("Q31_Basic.m")
%% generate data

tfinal=1000;
t = 0:T_s_close:tfinal;
u = 3*gensig('sine' , tfinal/20 , tfinal ,T_s_close);
Noise=(-0.2+(0.2+0.2)*rand(numel(t),1));
u=u+Noise;
y = lsim(sys_dis_close ,u ,t);
plot(t,u ,t , y ,'LineWidth',2) ;
xlabel('Time (sec)') ;
ylabel('Magnitute') ;
grid on
legend('Input','Output') ;
%% 
% kalman Filter

N = numel(y) ;
%--------------------------%
%choose number of parameters
Parameters_in_den=4
Parameters_in_num=4
%--------------------------%
Nv=Parameters_in_num+Parameters_in_den;
theta(:,1:Nv) = zeros(Nv , Nv) ;
P=[Nv,Nv];
P= 1e16*eye(Nv); 
phi=[Nv,N];
phi(1:Nv,1:N) = zeros(Nv , N) ;

Landa=1;
Sigma_e=(-0.1+(0.1+0.1)*rand(Nv,1));

for i = (max(Parameters_in_num,Parameters_in_den)+1):N
    phi(:,i) = [[y(i-1:-1:i-Parameters_in_den)]' , [u(i-1:-1:i-Parameters_in_num)]']';
    K = P*phi(:,i)*(Landa+phi(:,i)'*P*phi(:,i))^(-1);
    P =P-(K*phi(:,i)'*P/Landa)*(Landa + phi(:,i)'*P*phi(:,i))^(-1) ;
    theta(:,i) = theta(:,i-1) + K*(y(i) - phi(:,i)'*theta(:,i-1))+Sigma_e.^2; 
    Sigma_e=-0.05+(0.05+0.05)*rand(Nv,1);
end
%% 
% Bode

ident_dis = tf(theta((Parameters_in_num+1):end,end)' ,[1 -theta(1:Parameters_in_num ,end)'], T_s_close)
ident_analog = d2c(ident_dis)
bode(ident_analog ,'g*',sys_cont_close )
legend('model ','system')
%% 
% KF Convergence

plot(t , theta(:,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('Kalman Filter convergence') ;
grid on