%% STR Stochastic MV Nonadaptive
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
%% Define e
N=1000;
varN=0.001;% noise variance
e=random('normal',0,varN,N,1);
% e=0.01*randn(N,1)
%% controller
Deg_A=length(Ah);
Deg_B=length(B);
d0=Deg_A-Deg_B;
c=[0.1 -0.2 -0.3];
C=poly(c);
q=poly(zeros(1,d0-1));
B2=conv(q,C);
F1=deconv(C,Ah);
F=F1;
G1=C-conv(F,Ah);
G=[G1(2) G1(3) G1(4)];
BF = conv(B,F);
Y=zeros(1,N);
U=(zeros(1,N));
A_L=zeros(1,N);

for i=4:N
    Y(i)=-(Ah(1,2)*Y(i-1)+Ah(1,3)*Y(i-2)+Ah(1,4)*Y(i-3))+B(1)*U(i-1)+B(2)*U(i-2)+B(3)*U(i-3)+C(1)*e(i)+C(2)*e(i-1)+C(3)*e(i-2)+C(4)*e(i-3);
    U(i)=-(BF(1,2)*U(i-1)+BF(1,3)*U(i-2)+G(1,1)*Y(i)+G(1,2)*Y(i-1)+G(1,3)*Y(i-2))/(BF(1));
    
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
plot(Y-e','linewidth',2);legend('Y(t)-e(t)');grid on


%% Variance

Varianvce_Y=var(Y)
Variance_e=var(e)





















