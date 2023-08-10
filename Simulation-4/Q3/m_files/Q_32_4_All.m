close all
clear all
clc
%%
run('BASIC.m')
n = length(A)-1;
m = length(B)-1;
d0=4;
%%
Question_mark='Q32_4';
for d=[2,4,5,8,10]
    n=length(A)-1;
%%
    T_s=0.3;
    tfinal=100;
    freq=5;
    variance=0.01;
    C=0;
    dist=1
    [uc,t,Status,tfinal,Noix]=Datagen(0,T_s,tfinal,freq,variance,C);
    if dist==1
        Status=['-disturbed'];
    end
    Titlework=[Question_mark,Status];
    uc(1:40,1)=ones(40,1);
%%
    N=numel(t);
    var_y=zeros(1,N);  mean_y=zeros(1,N); ACLS_y=zeros(1,N);
    var_u=zeros(1,N);  mean_u=zeros(1,N); ACLS_u=zeros(1,N);
    u=zeros(1,N);
    y=zeros(1,N);
    Normerror=0;
    Error=0;
    Q = poly(zeros(1,d+n-1));
    [F,G] = dab(A,1,Q);

    Q = poly(zeros(1,n-1));
    BF = conv(B,F);
    [R,R_bar] = dab(Q,1,BF);

    dg=d+length(B);
    theta_hat=zeros((n+dg),numel(t));
    theta_hat=ones((n+dg),40);
    p_per=10^8*eye(n+dg);
    for i=max(d0,d)+n+d0:numel(t)-d
        if and(i>floor(N/3),dist==1)
            y(i)=-A(2:end)*y(i-1:-1:i-n)'+B*u(i-d0:-1:i-length(B)+1-d0)'+.1;
        else
            y(i)=-A(2:end)*y(i-1:-1:i-n)'+B*u(i-d0:-1:i-length(B)+1-d0)';
        end
        phi=[-y(i-1:-1:i-n),u(i:-1:i-d-length(B)+1)]';
        k=(p_per*phi)*(((1+phi'*p_per*phi))^-1);
        theta_hat(:,i)=theta_hat(:,i-1)+k*(y(i)-phi'*theta_hat(:,i-1));
        p=(eye(n+dg)-k*phi')*p_per;
        p_per=p;
        A_hat=[1,theta_hat(1:2,i)'];
        B_hat=(theta_hat(end-1:end,i))';
        [alpha,beta]=d_step_maker(A_hat,B_hat,d);

        u(i) = (uc(i+d) -G*y(i:-1:i-length(G)+1)'-R_bar*u(i-1:-1:i-length(R_bar))')./sum(R);
        var_u (i)=var(u)    ;mean_u(i)=mean(u)  ;ACLS_u(i)=u(i)^2+ACLS_u(i-1)   ;
        var_y(i)=var(y)     ;mean_y(i)=mean(y)  ;ACLS_y(i)=y(i)^2+ACLS_y(i-1)   ;
        Error=(y(i)-uc(i))^2+Error;
        Normerror=norm(Error,2)+Normerror;
        if abs(u(i))>3
            u(i)=sign(u(i))*3;
        end
    end
%%
    figure(1);
    subplot(2,1,1);
    plot(t,y,t,uc','--r','linewidth',2);

    title(['delay= ',int2str(d) ,' ' Status],'Color','b','fontsize',14)
    legend('output','reference');

    subplot(212);
    plot(t,u,'linewidth',2);
    legend('control signal');
    print(gcf,[Titlework , ' d=' num2str(d) ' Result.png'],'-dpng','-r400');
%%
    disp(['mean of u=',num2str(mean_u(end-d)),'   mean of y=',num2str(mean_y(end-d))])
    disp(['Error of u=',num2str(Error),'   NormError of y=',num2str(Normerror)])
    disp(['--------------------------------------'])
end