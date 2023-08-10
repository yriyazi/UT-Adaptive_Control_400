clc;
clear all;
close all;
warning off
plot_counter=1;

num=[1 -.8];
den=conv([1 -.85],[1 -.2]);
T_s=0.1;
d0=4;

A=den;
B=num;