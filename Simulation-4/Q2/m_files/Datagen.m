function [uc,t,Status,tfinal,Noix]=Datagen(num_noise,T_s,tfinal,Freq,variance,C)
t = 0:T_s:tfinal;
uc = gensig('square' , tfinal/Freq , tfinal ,T_s);
if num_noise==0
    Status=['-No NOISE-'];
    Noix=0;
elseif num_noise==1
    Noix = sqrt(variance).*randn(numel(t),1);
    uc=uc+Noix;
    Status=['-white NOISE-'];

elseif num_noise==2
    e = sqrt(variance).*randn(numel(t),1);
    %e = e-mean(e);
    Noix=zeros(numel(t),1);
    for i = numel(C):numel(t)
        Noix(i) = C*e(i:-1:i+1-(numel(C)));
    end
    uc=uc+(Noix);
    Status=['-Cololerd NOISE-'];
end



