clc;
clear all;
close all;
%%
run("Q310_Basic.m")
%% generate data

tfinal=200;
t = 0:T_s:tfinal;
u = gensig('sine' , tfinal/20 , tfinal ,T_s);

% t = 0:T_s_close:tfinal;
% u = gensig('sine' , tfinal/20 , tfinal ,T_s_close);

var_e=0.05;

e=sqrt(var_e);
Noise=(-e+(e+e)*rand(numel(t),1));
u=u+Noise;
y = lsim(sys_dis_open ,u ,t);


% smooth parameter variation 


for i=1:numel(t)
    paras(:,i)=[d(2:end),c]';
end



uu=round(numel(t)/3,0)
uuu=round(numel(t)*2/3,0)
for temp=uu:uuu
    paras(:,temp)=[d(2:end)+d(2:end)*.05*(sin(.5*(temp-uu))),(c+c*.05*(sin(.5*(temp-uu))-1))]';
    y(temp)=[-(y(temp-1:-1:temp-3))',(u(temp:-1:temp-3))']*paras(:,temp);
end

for temp=uuu:numel(y)
    paras(:,temp)=[d(2:end)+d(2:end)*.08*(sin(.5*(temp-uu))),(c+c*.08*(sin(.5*(temp-uu))-1))]';
    y(temp)=[-(y(temp-1:-1:temp-3))',(u(temp:-1:temp-3))']*paras(:,temp);
end

sys_dis = tf(c+c*.0001*(sin(.05*(temp-uu))-1)' ,[1 -d(2:end)+d(2:end)*.0001*(sin(.05*(temp-uu)))'], T_s)
ident_analog2 = d2c(sys_dis)
legend('Input','Output') ;

[c2,d2]=tfdata(sys_dis,'v')

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
theta_hat_KF(1:Nv,1:N)=zeros(Nv,N);
K_KF(1:Nv,1:N)=zeros(Nv,N);
var_e=0.05;

p_KF(1:Nv,1:Nv,1)=1e5*eye(Nv);p_KF(1:Nv,1:Nv,2)=p_KF(1:Nv,1:Nv,1);p_KF(1:Nv,1:Nv,3)=p_KF(1:Nv,1:Nv,2);


for i=(max(Parameters_in_num,Parameters_in_den)+1):N
    phi_KF(:,i)=[(y(i-1:-1:i-Parameters_in_den))',(u(i-1:-1:i-Parameters_in_num))']';
    K_KF(:,i) =p_KF(:,:,i-1)*phi_KF(:,i)*(1+phi_KF(:,i)'*p_KF(:,:,i-1)*phi_KF(:,i))^(-1) ;  
    p_KF(:,:,i)=p_KF(:,:,i-1)-p_KF(:,:,i-1)*phi_KF(:,i)*(1+phi_KF(:,i)'*p_KF(:,:,i-1)*phi_KF(:,i))^(-1)*phi_KF(:,i)'*p_KF(:,:,i-1)+var_e;
    theta_hat_KF(:,i)=theta_hat_KF(:,i-1)+K_KF(:,i)*(y(i)-phi_KF(:,i)'*theta_hat_KF(:,i-1));
end

%% 
% Bode

ident_dis = tf(theta_hat_KF((Parameters_in_num+1):end,end)' ,[1 -theta_hat_KF(1:Parameters_in_num ,end)'], T_s)
ident_analog = d2c(ident_dis)
bode(ident_analog ,'g*',sys_cont_open ,ident_analog2)
legend('model ','last system','Fisrst system')
%% 
% KF Convergence


subplot(2,1,1)
plot(t , theta_hat_KF((Parameters_in_num+1):end,:) , '*-') ;
hold on
plot(t ,paras(5:end,:) , 'LineWidth' , 1.5)
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('KF convergence Num') ;
grid on
legend('a_1','a_2','a_3')
% xlim([0 6])
% ylim([-1 1])
hold off
%--------------------------------------------------------------
subplot(2,1,2)

plot(t , -theta_hat_KF(1:Parameters_in_num ,:) , '*-') ;
hold on
plot(t ,paras(1:Parameters_in_num ,:) , 'LineWidth' , 1.5)
xlabel('Time (sec)') ;
ylabel('Parameters') ;
title('KF convergence Den') ;
grid on
legend('b_1','b_2','b_3')
% xlim([0 6])
hold off
% ylim([-7  7])