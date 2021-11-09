function [kmax,lmax,X_rel,jj,kk,M,I] = cube(r,total_mass)

%input parameters:
% r = diagonal length from center of mass to corner
% M = total mass

%outputs:
% X_rel(k,:) = coordinates of node k relative to the center of mass
% jj(l),kk(l) = indices of points connected by link l
% M = matrix of point masses
% I = moment of inertia for the cube

kmax = 8;       %8 nodes
lmax = 12;      %12 links
s = -sqrt(2) + 2*sqrt(0.5 + r^2);   %side length

%create position matrix
sign = [1,1,1; 1,-1,1; -1,-1,1; -1,1,1;
    1,1,-1; 1,-1,-1; -1,-1,-1; -1,1,-1];
X_rel = sign.*s/2;

%fill jj and kk for each link
jj=zeros(lmax,1); kk=zeros(lmax,1);
jj(1) = 1; kk(1) = 2;
jj(2) = 2; kk(2) = 3;
jj(3) = 3; kk(3) = 4;
jj(4) = 4; kk(4) = 1;
jj(5) = 5; kk(5) = 6;
jj(6) = 6; kk(6) = 7;
jj(7) = 7; kk(7) = 8;
jj(8) = 8; kk(8) = 5;
jj(9) = 1; kk(9) = 5;
jj(10) = 2; kk(10) = 6;
jj(11) = 3; kk(11) = 7;
jj(12) = 4; kk(12) = 8;

M = [(total_mass/8)*ones(kmax,1)];  %mass of each point
%moment of inertia tensor
I = diag([(1/6)*total_mass*s^2 (1/6)*total_mass*s^2 (1/6)*total_mass*s^2]);