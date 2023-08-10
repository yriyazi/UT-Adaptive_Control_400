clc;
clear all;
close all;
%% generate data

run ("Basics.m")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfinal=100;
t = 0:T_s:tfinal;
Priemss=primes(100);
u=zeros(numel(t),1);
% for i=15:numel(Priemss)
%     input_dummy=gensig('sine' , tfinal/Priemss(1,i) , tfinal ,T_s);
%     u=u+input_dummy;
% end
u = gensig('pulse' , tfinal/10 , tfinal ,T_s);
% Noise=0.01*rand(numel(t),1)
% u=u+Noise
y = lsim(sysd  ,u ,t);
% plot(t,u ,t , y ,'LineWidth',2) ;
% xlabel('Time (sec)') ;
% ylabel('V - position') ;
% title('square ressponse of Question 2') ;
% grid on
% legend(['U' , 'position']) ;
%% PA estimation

N = numel(y) ;
%choose number of parameters
Parameters_in_den=4
Parameters_in_num=4
Nv=Parameters_in_num+Parameters_in_den
theta(:,1:30) = 30*ones(Nv,30);
phi=[];

alfa=0.005
gama=1

for i = (max(Parameters_in_num,Parameters_in_den)+1):N
    phi(:,i) = [(y(i-1:-1:i-Parameters_in_den))',(u(i-1:-1:i-Parameters_in_num))']';
    error=gama*phi(:,i)*(y(i) - phi(:,i)'*theta(:,i-1))/(alfa+phi(:,i)'*phi(:,i));
    theta(:,i) = theta(:,i-1) +error; 
    norm(error);
end
theta(:,i)
ident_dis = tf(theta(Parameters_in_num+1:end,end)' ,[1 -theta(1:Parameters_in_num ,end)'], T_s)
ident_analog = d2c(ident_dis)
bode(ident_analog ,sys )

plot(t , theta(:,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('RLS convergence') ;
grid on