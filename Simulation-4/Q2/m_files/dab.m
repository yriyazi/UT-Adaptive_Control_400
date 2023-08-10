function [x,y] = dab(a,b,c)
% DAB   Solves the Diophantine-Aryabhatta-Bezout identity
%
%       [X,Y] = DAB(A,B,C)
%
%       AX + BY = C, where A, B, C, X and Y are polynomials
%       and deg Y = deg A - 1.

% Mats Lilja      LastEditDate : Wed Mar  7 16:26:14 1990
% Copyright (c) 1990 by Mats Lilja and Department of Automatic Control,
% Lund Institute of Technology, Lund, SWEDEN
global SIL
na = length(a);
nb = length(b);
nc = length(c);
ny = na - 1;
if ny<1,
  x = c/a;
  y = 0;
  return;
end;
nx = nc - ny;
c = [zeros(1,nb-nx-1) c];
nc = length(c);
nx = nc - ny;
if nx<1,
  x = 0;
  y = c/b;
  return;
end;
b = [zeros(1,nx-nb+1) b];
za = zeros(1,nx-1);
zb = zeros(1,ny-1);
ma = toeplitz([a za],[a(1) za]);
mb = toeplitz([b zb],[b(1) zb]);
m = [ma mb];
SIL = [SIL; det(m)] ;
if rank(m)<min(size(m)),
  disp('Singular problem due to common factors in A and B');
end;
xy = c/m';
x = xy(1:nx);
y = xy(nx+1:nc);
