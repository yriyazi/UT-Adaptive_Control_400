s=tf('s');
num=[0 1 3 7];
den=[1 8 7 7];
kc=7.79967495151787+182.26982236127/s;
sys_cont_open=tf(num,den);
fb_openloop = bandwidth(sys_cont_open);
T_s=0.05*2*pi/fb_openloop;%seconds
sys_dis_open = c2d(sys_cont_open, T_s,'zoh')
% kcdd = c2d(kc, Ts,'zoh')
[c,d]=tfdata(sys_dis_open,'v')
sys_cont_close=feedback(kc*sys_cont_open,1)
fb1 = bandwidth(sys_cont_close);
T_s_close=0.05*2*pi/fb1;%seconds
sys_dis_close = c2d(sys_cont_close, T_s_close,'zoh')
[c1,d1]=tfdata(sys_dis_close,'v')