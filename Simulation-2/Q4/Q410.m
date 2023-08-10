A=poly([-0.5 ,-3]);
B  =3;

omega_n=6;
zeta   =0.67;
A_m    =[1 2*zeta*omega_n omega_n^2];
B_m    =[omega_n^2];
sys=tf(B_m,A_m)
%step(sys,5)
%% Continous STR

P_0=1e12*eye(4);
theta_0=[2,0,0.5]

phi=[];
alpha=.3;
a_o=[1 10];
A_c=conv(a_o,A_m);
T=.0001;
%%
sim('Q410_sim_V2.slx')