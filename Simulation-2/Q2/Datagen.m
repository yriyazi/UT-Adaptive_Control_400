function [uc,t,Status,tfinal]=Datagen(Noise,T_s,tfinal)
t = 0:T_s:tfinal;
uc = gensig('square' , tfinal/10 , tfinal ,T_s);
if Noise==0
    Status=['_No NOISE_']
else 
    c1 = 0.008; % Colored Noise
    ve = 0.01;
    e = sqrt(ve).*randn(numel(t),1);
    e = e-mean(e);
    Noise(1) = e(1);
    for i = 2:numel(t)
        Noise(i) = e(i)+c1*e(i-1);
    end
    uc=uc+(Noise/3)';
    Status=['_white NOISE_']
end



