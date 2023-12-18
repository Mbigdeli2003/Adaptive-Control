%% STR-dynamic feedback-pole placemnet-nozero pole cancelled
clc
clear all
close all
format long
syms q
%% deifine Uc
%define T 3Ts=3*8.93
tic
t=0:150;
% Uc=zeros(9*length(T));
u1=2*ones(1,length(t));
u2=-1*ones(1,length(t));
u3=ones(1,length(t));
Uc=[zeros(1,3) u1 u2 u3];
%% define Tf
s=tf('s');
%unstable system
H=(s+0.5)/((3*s-1)*(s^2+1*s+1.2));
%desired system
G=(s+0.5)/((3*s+1)*(s^2+1*s+1.2));
%% dicrtisize sys with 20*bandwidth
fh=bandwidth(H);
fg=bandwidth(G);
Th=0.2;
Hd=c2d(H,Th);
Tg=Th;
Gd=c2d(G,Tg);
% H num & den
Ah=Hd.den{1,1};
Bh=Hd.num{1,1};
B=Bh(2:end);
Ag=Gd.den{1,1};
Bg=Gd.num{1,1};
B_g=Bg(2:end);
%system polynaminals defining with q
A_q=poly2sym(Ah,q);
B_q=poly2sym(B,q);
%% define desired system no pole zero cancellation
Am=Ag;
%finding Beta Beta=Am(1)/B(1)
Beta=(sum(Am))/(sum(B));
Bm=Beta*B;
%system polynaminals defining with q
Am_q=poly2sym(Am,q);
Bm_q=poly2sym(Bm,q);
%% control with diophantine
B_plus=1;
B_minus=Bm_q;
%defining S R A0
a0=0.02;a1=0.01;%defining A0 coeffs
syms r0 r1 r2 s0 s1 s2
R_q=(r0*q^2)+(r1*q)+r2;
S_q=(s0*q^2)+(s1*q)+s2;
A0_q=(q+a0)*(q+a1);
%% diphontine Eq solve
Dioph=(A_q*R_q)+(B_q*S_q)-(A0_q*Am_q);
coef=coeffs(Dioph,q);
% sort equation with degree repectively q0...q^5
eq1=coef(1,1);eq2=coef(1,2);eq3=coef(1,3);eq4=coef(1,4);
eq5=coef(1,5);eq6=coef(1,6);
% diphontine eq, finding parameters
Dioph_ans=solve(eq1,eq2,eq3,eq4,eq5,eq6);
%% finding R & S & T value
param_Cont=[Dioph_ans.r0 Dioph_ans.r1 Dioph_ans.r2 Dioph_ans.s0 Dioph_ans.s1 Dioph_ans.s2];
RS=sym2poly(param_Cont);
R=RS(1,1:3);S=RS(1,4:6);
T=sym2poly(A0_q);

%% testing controller obtaiend parameters Ru=TUc-Sy
yh=zeros(1,length(Uc));
yg=zeros(1,length(Uc));
Y=zeros(1,length(Uc));
U=zeros(1,length(Uc));
% for i=4:length(Uc)
%     % main system
%   yh(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*Uc(i-1))+(Bh(1,3)*Uc(i-2))+(Bh(1,4)*Uc(i-3));
% end
  for i=4:length(Uc)
      for j=3:length(Uc)
    U(j)=-(R(1,2)*U(j-1)+R(1,3)*U(j-2))+(T(1,1)*Uc(j)+T(1,2)*Uc(j-1)+T(1,3)*Uc(j-2))-(S(1,1)*Y(j)+S(1,2)*Y(j-1)+S(1,3)*Y(j-2));
     
   Y(i)=-(Ah(1,2)*Y(i-1))-Ah(1,3)*Y(i-2)-(Ah(1,4)*Y(i-3))+(Bh(1,2)*U(i-1))+(Bh(1,3)*U(i-2))+(Bh(1,4)*U(i-3));
      end
  end
 for i=4:length(Uc) 
  %  desired system
yg(i)=-(Ag(1,2)*yg(i-1))-Ag(1,3)*yg(i-2)-(Ag(1,4)*yg(i-3))+(Bg(1,2)*Uc(i-1))+(Bg(1,3)*Uc(i-2))+(Bg(1,4)*Uc(i-3));
%system with controller
% Y(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*U(i-1))+(Bh(1,3)*U(i-2))+(Bh(1,4)*U(i-3));
 end
 %% SSE
 SSE=((sum(2.3*Y-Uc))^2)*10
 %% plotting
%% plotting
 figure;
 subplot(2,1,1)
 plot(Uc,'b','linewidth',3);axis([0 500 -4 6]);hold on
 plot(2.3*Y,'g-','linewidth',3);axis([0 500 -4 6]);hold on
 plot(2.5*yg,'r--','linewidth',2.5);legend('Uc','Y','Yd');axis([0 500 -4 6]);grid on
subplot(2,1,2);
plot(U,'b','linewidth',1);legend('control signal')
grid on

toc

































