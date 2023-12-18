%% STR Stochastic MA_nonadaptive
clc
clear all
close all
%% define Tf
s=tf('s');
%unstable system
H=(s+0.5)/((3*s-1)*(s^2+1*s+1.2));
%% dicrtisize sys with 1/(10*Ts)
Th=0.2;
Hd=c2d(H,Th);
Tg=Th;
% H num & den
Ah=Hd.den{1,1};
Bh=Hd.num{1,1};
B=Bh(2:end);
B_plus=B/B(1,1);
B_minus=deconv(B,B_plus);
%% Define e
 N=1000;
varN=0.001;% noise variance
e=random('normal',0,varN,N,1);
e=0.01*randn(N,1);
%% controller
Deg_A=length(Ah);
Deg_B_plus=length(B_plus);
d=Deg_A-Deg_B_plus;
c=[0.1 -0.2 -0.3];
C=poly(c);
q=poly(zeros(1,d-1));
qB = conv(q,B_plus);
qBC = conv(qB,C);
%% silvester matrix and diophontine equation
E=[1          0         0       0       0       0;
   Ah(1,2)    1         0       B(1,2)  0       0;
   Ah(1,3)   Ah(1,2)   1        B(1,2)  B(1,1)  0;
   Ah(1,4)   Ah(1,3)   Ah(1,2)  B(1,3)  B(1,2)  B(1,1);
   0         Ah(1,4)   Ah(1,3)  0       B(1,3)  B(1,2);
   0         0         Ah(1,4)  0       0       B(1,3)];
%% R and S Parameter
er=inv(E)*qBC';
R=[er(1) er(2) er(3)];
S=[er(4) er(5) er(6)];
R_1=deconv(R,B_plus);
%% output result
U=zeros(1,length(e));
Y=zeros(1,length(e));
A_L=zeros(1,length(e));

for i=4:length(e)
  Y(i)=-(Ah(1,2)*Y(i-1)+Ah(1,3)*Y(i-2)+Ah(1,4)*Y(i-3))+B(1)*U(i-1)+B(2)*U(i-2)+B(3)*U(i-3)+C(1)*e(i)+C(2)*e(i-1)+C(3)*e(i-2)+C(4)*e(i-3);
    U(i)=(-R(1,2)*U(i-1)-R(1,3)*U(i-2)-S(1,1)*Y(i)-S(1,2)*Y(i-1)-S(1,3)*Y(i-2))/R(1,1);      

  A_L(i)=A_L(i-1)+Y(i).^2;
end
%% Plotting
% plot accumulated loss
figure;
plot(A_L,'b--','linewidth',2);legend('Acumulated Loss');grid on
% plot Y and U
figure;
subplot(2,1,1);
plot(Y,'linewidth',2);legend('Y output');grid on
subplot(2,1,2);
plot(U,'linewidth',2);legend('U Control signal');grid on

figure;
plot(Y-e','linewidth',2);legend('Y(t)-e(t)');axis([0 1000 -10 10]);grid on


%% Variance

Varianvce_Y=var(Y)
Variance_e=var(e)















