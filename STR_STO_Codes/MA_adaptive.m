%% STR Stochastic MA_adaptive
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
%% controller and identification
Deg_A=length(Ah);
Deg_B_plus=length(B_plus);
d=Deg_A-Deg_B_plus;
c=[0.1 -0.2 -0.3];
C=poly(c);
Y=zeros(1,length(e));
U=zeros(1,length(e));
beta=1e-2;
phi_H= beta*ones(1,9);

k =beta* ones(length(phi_H(1,:)));
p =beta* ones(length(phi_H(1,:)));
alpha =1e3;
p(:,:,3) = alpha*eye(length(phi_H(1,:)));
teta_H = [beta*ones(9,1),beta*ones(9,1),beta*ones(9,1)];
phi_H= ones(1,9);

K = zeros(length(phi_H(1,:)));
P = zeros(length(phi_H(1,:)));
alpha =1e6;

P(:,:,3) = alpha*eye(length(phi_H(1,:)));
teta_H_1 = [beta*ones(9,1),beta*ones(9,1),beta*ones(9,1)];
A_L=zeros(1,length(e));










for i=4:length(e)
     
%      if  teta_H(1:3,i-1)'-Ah(1:3)==0 & teta_H(4:6,i-1)'-Bh(1:3)==0
 Y(i)=-(Ah(1,2)*Y(i-1))-Ah(1,3)*Y(i-2)-(Ah(1,4)*Y(i-3))+(Bh(1,2)*U(i-1))+(Bh(1,3)*U(i-2))+(Bh(1,4)*U(i-3))+C(1,1)*e(i)+C(1,2)*e(i-1)+C(1,3)*e(i-2)+C(1,4)*e(i-3);   
 A_L(i) = A_L(i-1)+Y(i).^2;
 phi_H(i,:)=[-Y(i-1) -Y(i-2) -Y(i-3) U(i-1) U(i-2) U(i-3) e(i-1) e(i-2) e(i-3)];
  k(:,i)=p(:,:,i-1)*phi_H(i,:)'/(1+phi_H(i,:)*p(:,:,i-1)*phi_H(i,:)');
    p(:,:,i)=(eye(length(phi_H(i,:)))-k(:,i)*phi_H(i,:))*p(:,:,i-1);
    teta_H(:,i)=teta_H(:,i-1)+k(:,i)*(Y(i)-phi_H(i,:)*teta_H(:,i-1));
A_hat(i,1:3)=teta_H(1:3,i);
A_hat=[1 A_hat(i,1:3)];
B_hat(i,1:3)=teta_H(4:6,i);
B_hat=B_hat(i,1:3);
C_hat(i,1:3)=teta_H(7:9,i);
C_hat=[1 C_hat(i,1:3)];


    
    B_hat_plus = B_hat/B_hat(1);
     d = length(A_hat)-length(B_hat_plus); % d=deg A- deg B_pluse=1
     q = poly(zeros(1,d-1));
    qB = conv(q,B_hat_plus);
    qBC = conv(qB,C_hat);
    %% Silvester Matrix and dioph equation
    E11=[1           0           0          0           0           0      ;
         A_hat(2)    1           0          B_hat(1)    0           0      ;
         A_hat(3)    A_hat(2)    1          B_hat(2)    B_hat(1)    0      ;
         A_hat(4)    A_hat(3)    A_hat(2)   B_hat(3)    B_hat(2)    B_hat(1)  ;
           0         A_hat(4)    A_hat(3)     0         B_hat(3)    B_hat(2)  ;
           0          0          A_hat(4)     0        0            B_hat(3)] ;
 e1=inv(E11)*qBC';
 R=[e1(1,1) e1(2,1) e1(3,1)];
 S=[e1(4,1) e1(5,1) e1(6,1)];
 R1 = deconv(R,B_hat_plus);
 U(i)=(-R(1,2)*U(i-1)-R(1,3)*U(i-2)-S(1,1)*Y(i)-S(1,2)*Y(i-1)-S(1,3)*Y(i-2))/R(1); 
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


Variance_Y=var(Y)
Varinace_e=var(e)
figure;
subplot(2,3,1)
plot(4:i,Ah(1,2),'r*');
hold on
plot(teta_H(1,:),'b--');legend('a1')
subplot(2,3,2)
plot(4:i,Ah(1,3),'r*')
hold on
plot(teta_H(2,:),'b--');legend('a2')
subplot(2,3,3)
plot(4:i,Ah(1,4),'r*')
hold on
plot(teta_H(3,:),'b--');legend('a3')
subplot(2,3,4)
plot(4:i,Bh(1,2),'r*')
hold on
plot(teta_H(4,:),'b--');legend('b1')
subplot(2,3,5)
plot(4:i,Bh(1,3),'r*')
hold on
plot(teta_H(5,:),'b--');legend('b2')
subplot(2,3,6)
plot(4:i,Bh(1,4),'r*')
hold on
plot(teta_H(6,:),'b--');legend('b3')

















