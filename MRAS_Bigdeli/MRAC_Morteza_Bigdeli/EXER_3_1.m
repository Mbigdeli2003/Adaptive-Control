%% MIT Rule
clc
clear all
close all
warning off
UC=1;
gama1=5;
gama2=5;
sim ('exer_3_1',500);
subplot(3,1,1);
plot(time,Y,'r--','linewidth',4.5);hold on
plot(time,Yd,'g--','linewidth',4);hold on;
plot(time,Uc,'b','linewidth',2);grid on
legend('Y','Yd','Uc')
subplot(3,1,2);
plot(time,U,'linewidth',2),legend('U control signal');grid on
subplot(3,1,3);
plot(time,ERROR,'r','linewidth',2.5),legend('ERROR');grid on
figure(2)
c=ones(length(time),1);
subplot(211),plot(time,to,'b-*','linewidth',2),legend('t0');grid on
subplot(212),plot(time,so,'r-*','linewidth',2),legend('s0');grid on

figure(3)
plot(to,so,'b--','linewidth',2.5);xlabel('t0');ylabel('s0');grid on