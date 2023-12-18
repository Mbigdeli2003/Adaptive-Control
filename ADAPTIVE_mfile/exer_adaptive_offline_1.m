%% exer_adaptive_offline

clc
clear all
close all
format long
%% define tranfer function
s=tf('s');
G=(s+0.5)/((3*s+1)*(s^2+0.7*s+1));
H=(s+0.5)/((3*s-1)*(s^2+0.7*s+1));
%% discrisized the system with 20 times bandwidth
fg=bandwidth(G);
fh=bandwidth(H);
Tg=1/(20*fg);
Th=1/(20*fh);
Gd=c2d(G,Tg);
Hd=c2d(H,Th);
%% define idpoly system
%G num & den
Ag=Gd.den{1,1};
Bg=Gd.num{1,1};
% idpoly G system
m0_G=idpoly(Ag,Bg)
% H num & den
Ah=Hd.den{1,1};
Bh=Hd.num{1,1};
% idpoly H system
m0_H=idpoly(Ah,Bh)
%% define inputs
%random input
n=10
u_1=rand(n,1);
% impulse input
p=10;%pulse value
j=5;k=4;
u_2=[zeros(1,j) p zeros(1,k)]';
%  u_2=[ ones(1,10)]';
%% sim system
%% sim system with random input
%G 
y_G_1=sim(m0_G,u_1);
%H
y_H_1=sim(m0_H,u_1);
%% sim system with impulse input
%G 
y_G_2=sim(m0_G,u_2);
%H
y_H_2=sim(m0_H,u_2);
%% system offline identification with RLS

%% system offline identification with random inputs
x0=[0 0 0]';
  y22=-1*[x0(1,1);y_G_1(1:size(u_1)-1,1)];
  y33=-1*[x0(1:2,1);y_G_1(1:size(u_1)-2,1)];
  y44=-1*[x0(1:3,1);y_G_1(1:size(u_1)-3,1)];
  u11=[0;u_1(1:size(u_1)-1,1)];
  u22=[zeros(2,1);u_1(1:size(u_1)-2,1)];
  u33=[zeros(3,1);u_1(1:size(u_1)-3,1)];
  phi_1=[y22 y33 y44 u11 u22 u33];
  %teta estimation
  teta_1_G=(inv(phi_1'*phi_1))*phi_1'*y_G_1
  %% system offline identification H with random inputs
x0=[0 0 0]';
  y22=-1*[x0(1,1);y_H_1(1:size(u_1)-1,1)];
  y33=-1*[x0(1:2,1);y_H_1(1:size(u_1)-2,1)];
  y44=-1*[x0(1:3,1);y_H_1(1:size(u_1)-3,1)];
  u11=[0;u_1(1:size(u_1)-1,1)];
  u22=[zeros(2,1);u_1(1:size(u_1)-2,1)];
  u33=[zeros(3,1);u_1(1:size(u_1)-3,1)];
  phi_1=[y22 y33 y44 u11 u22 u33];
  %teta estimation
  teta_1_H=(inv(phi_1'*phi_1))*phi_1'*y_H_1
%% system offline identification G with impulse inputs
x0=[0 0 0]';
  y22=-1*[x0(1,1);y_G_2(1:size(u_2)-1,1)];
  y33=-1*[x0(1:2,1);y_G_2(1:size(u_2)-2,1)];
  y44=-1*[x0(1:3,1);y_G_2(1:size(u_2)-3,1)];
  u11=[0;u_2(1:size(u_2)-1,1)];
  u22=[zeros(2,1);u_2(1:size(u_2)-2,1)];
  u33=[zeros(3,1);u_2(1:size(u_2)-3,1)];
  phi_2=[y22 y33 y44 u11 u22 u33];
  %teta estimation
  teta_2_G=((inv(phi_2'*phi_2))*(phi_2'))*y_G_2
%% system offline identification H with impulse inputs
x0=[0 0 0]';
  y22=-1*[x0(1,1);y_H_2(1:size(u_2)-1,1)];
  y33=-1*[x0(1:2,1);y_H_2(1:size(u_2)-2,1)];
  y44=-1*[x0(1:3,1);y_H_2(1:size(u_2)-3,1)];
  u11=[0;u_2(1:size(u_2)-1,1)];
  u22=[zeros(2,1);u_2(1:size(u_2)-2,1)];
  u33=[zeros(3,1);u_2(1:size(u_2)-3,1)];
  phi_2=[y22 y33 y44 u11 u22 u33];
  
  %teta estimation
  teta_2_H=((inv(phi_2'*phi_2))*(phi_2'))*y_H_2
  























