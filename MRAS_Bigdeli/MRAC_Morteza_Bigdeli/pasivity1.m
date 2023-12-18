%% passivity MRAS
clc
clear all
close all
N = 15;
gama = 70;
r = 1;
%% Inputs
T = 2*pi/10;
% Square Wave Input
%T = 20;
%% System
alpha = 2;
s = tf('s');
G = (s+3)/(s^2+2*s+4);
bm1 = 2; am1 = 8; am2 = 2;
figure(1);
nyquist(G);%grid on
% state space
S = ss(G);
A = S.a;
B = S.b;
C = S.c;
D = S.d;

Q = eye(2);
P =  lyap(A',Q); %%Solving continuous-time Lyapunov equation
p1 = P(1);
p2 = P(2);
%% Simulation
sim('passivity1.mdl',N);
figure(3)
subplot(413)
plot(alpha*ones(1,length(teta)),'r','linewidth',3)
hold on
plot(teta,'linewidth',2)
grid on
legend('teta')

%figure(4)
subplot(411)
hold on
plot(t,Ym,'r','linewidth',3)
plot(t,Y,'linewidth',2)
legend('y_m ','y')
grid on

grid on
subplot(412)
plot(t,e,'linewidth',2)
title(' error ')
grid on

subplot(414),plot(t,U,'linewidth',2)
title(' Comtrol signal (u) ')
grid on
