clc;
clear all;
close all;
%% generate data

run ("Basics.m")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tfinal=200;
t = 0:T_s:tfinal;
u = gensig('sine' , tfinal/20 , tfinal ,T_s);%+gensig('sine' , tfinal/50 , tfinal ,T_s)+gensig('square' , tfinal/10 , tfinal ,T_s);

Noise=[1:numel(t)];
Noise(1:10) = ones(1,10)';
for i=4:1:numel(t)
    Noise(i)=+0.4*Noise(i-1)+0.75*Noise(i-2)-0.15*Noise(i-3);
end
Noise2=-0.2+(0.2+0.2)*rand(numel(t),1);
u=u+Noise'+Noise2;

y = lsim(sysd  ,u ,t);
plot(t,u ,t , y ,'LineWidth',2) ;
xlabel('Time (sec)') ;
ylabel('Value') ;
title('ressponse of Question 2') ;
grid on
legend('Input' , 'ّOutPut') ;
%% recursive least esquare estimation

N = numel(y) ;
%choose number of parameters
Parameters_in_den=4
Parameters_in_num=4
Parameters_in_noise=3
Nv=Parameters_in_num+Parameters_in_den+Parameters_in_noise;
P = 1e12*eye(Nv) ;
theta=[Nv,N]
theta(1:Nv,1:30) = 5*ones(Nv,30) ;
e=[numel(y),1]
e=zeros(numel(y),1)
phi=[];
Eror=zeros(1,N)
tic
for i = (max(Parameters_in_num,Parameters_in_den+Parameters_in_noise))+1:N
    phi(:,i) = [[y(i-1:-1:i-Parameters_in_den)]' , [u(i-1:-1:i-Parameters_in_num)]' , [e(i-1:-1:i-Parameters_in_noise)]']';
    K = P*phi(:,i)*(1+phi(:,i)'*P*phi(:,i))^(-1) ;
    P = (eye(Nv) - K*phi(:,i)')*P ;
    theta(:,i) = theta(:,i-1) + K*(y(i) - phi(:,i)'*theta(:,i-1)); 
    Eror(i)=(Eror(i-1)+(y(i)-phi(:,i)'*theta(:,i))^2);
end
%% 
% Bode

ident_dis = tf(theta((Parameters_in_num+1):(end-Parameters_in_noise),end)' ,[1 -theta(1:Parameters_in_num ,end)'], T_s)
ident_analog = d2c(ident_dis)
bode(ident_analog ,'g*',sys )
legend('model ','system')
%% 
% RLS Convergence

plot(t , theta(:,:) , 'LineWidth' , 2) ;
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('RLS convergence') ;
grid on
%%
plot(1:1:N,Eror)
xlabel('Iteration') ;
ylabel('error square') ;
title('Cost function \times 2') ;
%%
tfinal=1000;
T_s=T_s
t = 0:T_s:tfinal;
u = gensig('square' , tfinal/20 , tfinal ,T_s);
u = u+rand(numel(t),1);
y = lsim(sysd,u ,t);

y_model = lsim(ident_dis ,u ,t);

plot(t,y_model ,'b*',t , y ,'LineWidth',1.25) ;
xlabel('Time (sec)') ;
ylabel('Value') ;
title('ressponse of Question 2') ;
grid on
legend('Under parameter Model' , 'ّSystem') ;