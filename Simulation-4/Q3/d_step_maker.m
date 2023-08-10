function [alpha,beta] =d_step_maker(a1,a2,d)
a1=-a1(2:end);

a10=a1; da1=length(a1); a1=[a1,zeros(1,d-1)];
a20=a2; da2=length(a2); a2=[a2,zeros(1,d-1)];

for i = 1 : d-1
    a1(:,i+1:i+da1)=a1(:,i+1:i+da1)+a1(i)*a10;
    a2(:,i+1:i+da2)=a2(:,i+1:i+da2)+a1(i)*a20;
    a1(i)=0;
end
alpha=a1(end-1:end);
beta=a2;
end


    
    
    