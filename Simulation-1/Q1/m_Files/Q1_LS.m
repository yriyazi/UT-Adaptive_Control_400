
clc
clear
run ("Basics.m");
% System identification

tic
tfinal=200;
t = 0:T_s:tfinal;
u = zeros(numel(t),1);
% General Input+white Noise

% u = gensig('sine' , tfinal , tfinal ,T_s);
% Noise=-0.2+(0.2+0.2)*rand(numel(t),1);
% u=u+Noise;
% 1.Pulse Input

% u(1:50,1)=ones(50,1) ;
% 2.Step Input

% u=ones(numel(t),1);
% % u(round(numel(t)/10,0):end,1)=1;
% 3.Sine Input

% u = gensig('sine' , tfinal/15 , tfinal ,T_s);
% 4.Ramp Input

for i=1:numel(t)
    u(i)= 1.2*i;%randi(1);
end    
%% 
% Out Put Generating

y = lsim(sysd  ,u ,t);
plot(t,u ,t , y ,'LineWidth',2) ;
xlabel('Time (sec)') ;
ylabel('Value') ;
title('ressponse of Question 1') ;
grid on
legend('Input' , 'ّOutPut') ;
xlim([0 200])
% ylim([-1.2 1.2])
%% LS Identification

N = numel(y) ;
Parameters_in_den=4
Parameters_in_num=4
Nv=Parameters_in_num+Parameters_in_den
phi=[];
for i=(max(Parameters_in_num,Parameters_in_den)+1):N
    phi(i,:) = [(y(i-1:-1:i-Parameters_in_den))' , (u(i-1:-1:i-Parameters_in_num))'];
end
theta_hat=((phi'*phi)^(-1))*(phi'*y)
% norm([theta_hat]-[d,c(2:end)]')
% norm(Y-phi*theta_hat)
sysdd=tf(theta_hat((Parameters_in_num+1):end,end)' ,[1 -theta_hat(1:Parameters_in_num ,end)'], T_s)
% ‌BODE

ident_analog = d2c(sysdd)
figure
bode(ident_analog ,'g*',sys )
legend('model','system')
% Ploting discret system and Least square Model

figure
plot(y)
hold on
plot(phi*theta_hat,'r--')
xlabel('Sample cber')
ylabel('Output')
legend('\fontsize{12} discret system','\fontsize{12} Ls Model');
grid on;
% Ploting discret system and Least square Model via step input

figure
step(sysd,0:T_s:100*T_s)
hold on
step(sysdd,0:T_s:100*T_s,'r+')
legend('\fontsize{12} discret system','\fontsize{12} Ls Model');
grid on;
xlabel('time','fontsize',12);
% ylabel('x2','fontsize',16);
toc