
clc
clear 
close

run('BASIC.m')

Titlework='Q11'
%% Desigered System

%run('am.mlx');
%A_m=den_discret_desierd%A=den_discret
A_m=poly([0.80 0.65 0.69]);


betaa=sum(A_m)/sum(B);
B_m=B*betaa;
B_plus=1;
%% MDPP with no zero canselation

A_o=[1 0 0];
A_c=conv(A_m,A_o)
[R_prim , S] = Diophantine(A , B , A_c)
T=conv(betaa,A_o)
R=conv(R_prim,B_plus)
tfinal=200;

t = 0:T_s:tfinal;
uc=gensig('square' , tfinal/3 , tfinal ,T_s);

u=zeros(numel(t),1);
y=zeros(numel(t),1);

for i=10:numel(t)
    var1=conv(B,T)  ;       narv1=numel(var1)      ;
    var2=A_c(2:end) ;       narv2=numel(A_c(2:end));
    var3=conv(A,T)  ;       narv3=numel(var3)      ;

    y(i)=var1*uc(i:-1:i-narv1+1)-var2*uc(i-1:-1:i-narv2);
    u(i)=var3*uc(i:-1:i-narv3+1)-var2*uc(i-1:-1:i-narv2);
end

%%
plot(t,uc,'b',t,y,'r--','LineWidth',1)
title('Refrence vs Controled output')
print(gcf,[Titlework ' Refrence vs Controled output.png'],'-dpng','-r400');

plot(t,u)
title('Control Signal')
print(gcf,[Titlework ' Control Signal.png'],'-dpng','-r400');