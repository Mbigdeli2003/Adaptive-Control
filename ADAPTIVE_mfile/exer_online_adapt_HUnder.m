%% Exer_2_adaptive_online
clc
clear all
close all
format long
%% define tranfer function
s=tf('s');
G=(s+0.5)/((3*s+1)*(s^2+0.7*s+1));
H=(s+0.5)/((3*s-1)*(s^2+0.7*s+1));
%% discrete system with 20 times bandwidth sapmilg frequency
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
AG=[Ag(1,2:4) Bg(1,2:4)]';
% idpoly G system
m0_G=idpoly(Ag,Bg)
% H num & den
Ah=Hd.den{1,1};
Bh=Hd.num{1,1};
AH=[Ah(1,2:4) Bh(1,2:4)]';
% idpoly H system
m0_H=idpoly(Ah,Bh)

%% system H with random inputs
n=1000;
%  u_1=[ zeros(1,3) rand(1,n)]';
yh=zeros(n,1);
% u = [0;0;u_1];
%random input
%    u_1=[ones(1,n)];
% impulse input
%  u_1=[0 ones(1,n)]'
t=0:0.1:100;
w=10;
% u_1=sin(w*t);
p=1;%pulse value
j=500;k=500;
%   u_1=[zeros(1,j) p zeros(1,k)];
%   u_2=ones(10,1);

%sinus inputs
t=0:0.1:1;
w=10;
u_3=sin(w*t)';
u_1=rand(1,n);
u_1=[ zeros(1,3) u_1]';


phi_H= zeros(1,6);

k = zeros(length(phi_H(1,:)));
p = zeros(length(phi_H(1,:)));
alpha =1e10;
p(:,:,3) = alpha*eye(length(phi_H(1,:)));
teta_H = [zeros(6,1),zeros(6,1),zeros(6,1)];

phi_H_1= zeros(1,4);

K = zeros(length(phi_H_1(1,:)));
P = zeros(length(phi_H_1(1,:)));
alpha =1e10;
P(:,:,3) = alpha*eye(length(phi_H_1(1,:)));
teta_H_1 = [zeros(4,1),zeros(4,1),zeros(4,1)];


%% adaptive identification system H
for i=4:1:max(size(u_1))
 yh(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*u_1(i-1))+(Bh(1,3)*u_1(i-2))+(Bh(1,4)*u_1(i-3));   
  phi_H(i,:)=[-yh(i-1) -yh(i-2) -yh(i-3) u_1(i-1) u_1(i-2) u_1(i-3)];
    k(:,i)=p(:,:,i-1)*phi_H(i,:)'/(1+phi_H(i,:)*p(:,:,i-1)*phi_H(i,:)');
    p(:,:,i)=(eye(length(phi_H(i,:)))-k(:,i)*phi_H(i,:))*p(:,:,i-1);
    teta_H(:,i)=teta_H(:,i-1)+k(:,i)*(yh(i)-phi_H(i,:)*teta_H(:,i-1));
%%%%%%
 yh(i)=-(Ah(1,2)*yh(i-1))-Ah(1,3)*yh(i-2)-(Ah(1,4)*yh(i-3))+(Bh(1,2)*u_1(i-1))+(Bh(1,3)*u_1(i-2))+(Bh(1,4)*u_1(i-3));   
  phi_H_1(i,:)=[-yh(i-1) -yh(i-2)  u_1(i-1) u_1(i-2)];
    K(:,i)=P(:,:,i-1)*phi_H_1(i,:)'/(1+phi_H_1(i,:)*P(:,:,i-1)*phi_H_1(i,:)');
    P(:,:,i)=(eye(length(phi_H_1(i,:)))-K(:,i)*phi_H_1(i,:))*P(:,:,i-1);
    teta_H_1(:,i)=teta_H_1(:,i-1)+K(:,i)*(yh(i)-phi_H_1(i,:)*teta_H_1(:,i-1));



end
init ='z';
yu = sim (m0_H,u_1);
%% plotting system parameters
% figure;
% 
% subplot(2,3,1)
% plot(4:i,Ah(1,2),'r*');
% hold on
% plot(teta_H(1,:),'b--');legend('a1')
% subplot(2,3,2)
% plot(4:i,Ah(1,3),'r*')
% hold on
% plot(teta_H(2,:),'b--');legend('a2')
% subplot(2,3,3)
% plot(4:i,Ah(1,4),'r*')
% hold on
% plot(teta_H(3,:),'b--');legend('a3')
% subplot(2,3,4)
% plot(4:i,Bh(1,2),'r*')
% hold on
% plot(teta_H(4,:),'b--');legend('b1')
% subplot(2,3,5)
% plot(4:i,Bh(1,3),'r*')
% hold on
% plot(teta_H(5,:),'b--');legend('b2')
% subplot(2,3,6)
% plot(4:i,Bh(1,4),'r*')
% hold on
% plot(teta_H(6,:),'b--');legend('b3')
% 
%  figure
%  plot(u_1,'r');legend('step input')
% figure
% subplot(2,1,1)
% plot(yu,'b');legend('real output')
% subplot(2,1,2)
% plot(phi_H*teta_H,'black--')
% ;legend('phi*teta')
figure;

subplot(2,2,1)
plot(4:i,Ah(1,2),'r*');
hold on
plot(teta_H_1(1,:),'b--');legend('a1')
subplot(2,2,2)
plot(4:i,Ah(1,3),'r*')
hold on
plot(teta_H_1(2,:),'b--');legend('a2')
% subplot(2,3,3)
% plot(4:i,Ag(1,4),'r*')
% hold on
% plot(teta_G(3,:),'b--');legend('a3')
subplot(2,2,3)
plot(4:i,Bh(1,2),'r*')
hold on
plot(teta_H_1(3,:),'b--');legend('b1')
subplot(2,2,4)
plot(4:i,Bh(1,3),'r*')
hold on
plot(teta_H_1(4,:),'b--');legend('b2')
% subplot(2,3,6)
% plot(4:i,Bg(1,4),'r*')
% hold on
% plot(teta_G(6,:),'b--');legend('b3')

 figure
 plot(u_1,'r');legend('random input')
figure
subplot(2,1,1)
plot(yu,'b');legend('real output')
subplot(2,1,2)
plot(phi_H_1*teta_H_1,'black--')
;legend('phi*teta')




