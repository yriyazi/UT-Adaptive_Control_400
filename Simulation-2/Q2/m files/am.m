omega_n=2.5372
zeta=0.5912
g1=tf([1 2*zeta*omega_n omega_n^2],1);
p1=-15;
g2=zpk(p1,[],1);
z1=-10;
z2=-11
k2=-p1/(z1*z2)
g3=zpk([z1 z2],[],(omega_n^2)*k2);
sys=tf(g3/(g1*g2))
sys_discret=c2d(sys,0.3)
[num_discret,den_discret]=tfdata(sys_discret);
num_discret=cell2mat(num_discret);
num_discret=num_discret(2:end)      ;%B=num_discret
den_discret=cell2mat(den_discret)   %A=den_discret