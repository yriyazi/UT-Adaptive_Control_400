clc;
clear all;
close all;
warning off
plot_counter=1;

sys_discret=zpk([(1/0.74)],[0.15 0.82],1);
T_s=0.01;

[num_discret,den_discret]=tfdata(sys_discret);
num_discret=cell2mat(num_discret);
num_discret=num_discret(2:end)      ;B=num_discret;
den_discret=cell2mat(den_discret)   ;A=den_discret;