clear all
close all
clc
%%
run("BASIC_nonMinimum.m")
C=poly([0.6 0.7]);

Question_mark='Q15-MV nonadaptive';
[uc,t,Status,tfinal,Noise]=Datagen_Stochastic(2,T_s,500,0.01,C);
Titlework=[Question_mark,Status]
%%
plot(t,Noise,'b','LineWidth',1)
title('Colored Noise input')
legend('Noise')
print(gcf,[Titlework , num2str(plot_counter) ' Colored Noise input.png'],'-dpng','-r400');
plot_counter=plot_counter+1;
N = numel(t) ;

% NMP systems
z=roots(B);
if (abs(z)>1)
    z=1/z;
    B=poly(z)*B(1,1)*z^-1;
else
    B=B;
end

Ac=conv(C,B);
d0=1;
D0=zeros(d0,1);
D0(1)=1;
if D0==[]
    Ac=conv(C,B);
else
    Ac=conv(Ac,D0);
end
%%
[R , S ] = Diophantine(A , B , Ac)      
%%
y=zeros(1,N);
y(1,1:10)=0.0*ones(1,10);%initial condotion
u=zeros(1,N);
u(1,1:10)=0.0*ones(1,10);%initial condotion
var_y=zeros(1,N);   mean_y=zeros(1,N);  ACLS_y=zeros(1,N);
var_u=zeros(1,N);   mean_u=zeros(1,N);  ACLS_u=zeros(1,N);
for i=3:numel(t)

    y(i)= -A(2:3)*[y(i-1:-1:i-2)]'+ B*[uc(i-1:-1:i-2)]+Noise(i);
    u(i)=(-R(2)*u(i-1)-S*[y(i:-1:i-1)]')./R(1);

    var_u (i)=var(u)    ;mean_u(i)=mean(u)  ;ACLS_u(i)=u(i)^2+ACLS_u(i-1)   ;
    var_y (i)=var(y)    ;mean_y(i)=mean(y)  ;ACLS_y(i)=y(i)^2+ACLS_y(i-1)   ;
end
%% Output and control signal plotting
% General Input v.s. Output

figure
plot(t,Noise,'b-',t,y,'r-','LineWidth',1)
title('Refrence vs Controled output')
legend('input','output')
xlim('auto')
ylim([-0.5 1.5])
print(gcf,[Titlework , num2str(plot_counter) ' Refrence vs Controled output.png'],'-dpng','-r400');
plot_counter=plot_counter+1;

plot(t,u)
title('Control Signal')
% xlim('auto')
% ylim([-1.5 1.5])
print(gcf,[Titlework , num2str(plot_counter) ' Control Signal.png'],'-dpng','-r400');
plot_counter=plot_counter+1;
%% Variance plotting

figure
subplot(211)
plot(t,var_y)
hold on
title('Output Variance')

subplot(212)
hold on
title('Variance of Controller Signal')
plot(t,var_u)
print(gcf,[Titlework , num2str(plot_counter) ' Variance.png'],'-dpng','-r400');
plot_counter=plot_counter+1;
%% Mean plotting

figure
subplot(211)
plot(t,mean_y)
hold on
title('Output Mean')

subplot(212)
hold on
plot(t,mean_u)
title('Mean of Controller Signal')
print(gcf,[Titlework , num2str(plot_counter) ' mean.png'],'-dpng','-r400');
plot_counter=plot_counter+1;
%% Accumulated loss plotting

figure
subplot(211)
hold on
plot(t,ACLS_y)
title('Output Accumulated Loss')

subplot(212)
hold on
title(' Accumulated Loss of Controller Signal')
plot(t,ACLS_u)
print(gcf,[Titlework , num2str(plot_counter) ' Accumulated Loss.png'],'-dpng','-r400');
plot_counter=plot_counter+1;