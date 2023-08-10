function [uc,t,Status,tfinal,Noix]=Datagen_Stochastic(num_noise,T_s,tfinal,variance,C)
t = 0:T_s:tfinal;
uc = gensig('square' , tfinal/10 , tfinal ,T_s);
if num_noise==0
    Status=['-No NOISE-']

elseif num_noise==1
    Noix = sqrt(variance).*randn(numel(t),1);
    uc=uc+Noix;
    Status=['-white NOISE-']

elseif num_noise==2
    e = sqrt(variance).*randn(numel(t),1);
    %e = e-mean(e);
    Noix=zeros(numel(t),1);
    for i = numel(C):numel(t)
        Noix(i) = C*e(i:-1:i+1-(numel(C)));
    end
    uc=uc+(Noix);
    Status=['-Cololerd NOISE-']
end



