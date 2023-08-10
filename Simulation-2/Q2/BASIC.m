clc;
clear all;
close all;
warning off
plot_counter=1;

sys_cont=zpk([-1.6 -3.4],[-4/0.7 -1.6/0.3 1.1],0.1);
BW=bandwidth(sys_cont);
% to signal can be reconstructable we must have 2*Bandwich
Discret_ratio=15; % >2
T_s=2*pi/(BW*Discret_ratio);
sys_discret=c2d(sys_cont,T_s,'zoh')
[num_discret,den_discret]=tfdata(sys_discret);
num_discret=cell2mat(num_discret);
num_discret=num_discret(2:end)      ;B=num_discret;
den_discret=cell2mat(den_discret)   ;A=den_discret;