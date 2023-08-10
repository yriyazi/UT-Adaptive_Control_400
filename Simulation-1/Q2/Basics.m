last2digit=03;
m_1=0.1;
m_2=m_1;
b_2=(1+last2digit)/15; k_1=b_2;
b_1=(10+last2digit)/100; k_2=b_1;
k_3=0.5*k_1;
A=[0                1           0               0
  -(k_3+k_1)/m_1    -b_1/m_1    k_3/m_1         b_1/m_1
  0                 0           0               1
  k_3/m_2           b_1/m_2     -(k_3+k_2)/m_2  -b_2/m_2];
B=[0;1;0;0];
C=[0 0 1 0];
D=0;
[b,a]= ss2tf(A,B,C,D);
sys=tf(b,a)
fb = bandwidth(sys)
T_s=0.05*2*pi/fb;
sysd = c2d(sys,T_s,'zoh')
[c,d]=tfdata(sysd,'v')